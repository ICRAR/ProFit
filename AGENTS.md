# AGENTS.md — ProFit

Guidance for AI/LLM agents working in this repository. For human-facing docs see
`README.md` and the vignettes.

## What this package is

`ProFit` fits projected 2D profiles (Sersic, Moffat, Ferrer, CoreSersic, King,
BrokenExponential, point sources, sky) to galaxy images. It is a thin R layer
over a vendored C++ library, `libprofit`, which does all of the image
generation. The R side handles parameter bookkeeping, priors, and likelihoods;
the C++ side handles profile integration and convolution.

Projected usage is overwhelmingly through MCMC/optimisers, which call
`profitLikeModel()` tens of thousands of times for a single fit. **Almost all
performance questions in this package are questions about the cost of one
`profitLikeModel()` call.**

## Repository layout

| Path | Contents |
|------|----------|
| `R/profitMakeModel.R` | R entry point; assembles the model list passed to C++. Hot-path-adjacent. |
| `R/profitLikeModel.R` | The likelihood. Called once per sampler iteration. |
| `R/profitSetupData.R` | One-off setup: PSF parsing, `calcregion`, `convopt`, buffers. |
| `R/profitRemakeModelList.R` | Maps the free-parameter vector onto the modellist. |
| `R/profitBenchmark*.R` | Benchmarking helpers, including convolver/integrator selection. |
| `src/r_profit.cpp` | `.Call` boundary. Reads R lists into `Model`/`Profile` objects. |
| `src/libprofit/src/profit/` | Vendored libprofit. `radial.cpp`, `sersic.cpp`, `convolve.cpp`, `image.*`, `model.cpp` are the files that matter for speed. |
| `configure` | Detects OpenCL / FFTW / OpenMP and generates `src/Makevars` from `src/Makevars.in`. |

`src/Makevars` and `src/libprofit/src/profit/config.h` are **generated**; edit
`src/Makevars.in` and `configure` instead. The same applies to the OpenCL
kernel headers in `src/libprofit/src/profit/cl/`.

## Building and testing

```sh
R CMD INSTALL .            # respects ./configure; needed after C++ changes
```

`tests/` holds a small set of base-R regression scripts (no testthat
dependency), run automatically by `R CMD check`. They are deliberately narrow:
they cover the PSF padding/point-source geometry, the pixel-scale handling,
the evaluation plumbing (thread propagation, buffer sizing and reuse) and the
convolver choice, because those are the things that have actually broken. Everything else is verified by example and by
vignette, so when changing numerical code, compare old vs new output
explicitly (see "Verifying numerical changes" below) rather than assuming.

`libprofit` is vendored rather than linked, so a deliberate upstream sync means
replacing files wholesale — check `src/libprofit/VERSION` and the git log
messages that mention `libprofit` versions.

## The performance model (read this before optimising anything)

Rough shares of wall time for a typical fit (`profitLikeModel`, converged
sampler, `200x200` to `1000x1000`, `25x25` PSF):

- **Convolution** — usually the largest single item. Brute-force is the
  default, FFTW is typically 1.5–2.3× faster at these sizes, and much more for
  larger PSFs.
- **Radial profile integration** — second. Dominated by sub-pixel integration
  near the centre (`subsample_pixel`), which is recursive.
- **Likelihood arithmetic** — can be 10–20% at megapixel sizes. Its cost is
  mostly data movement and `dnorm()`, not the model.
- **R bookkeeping** (`profitRemakeModellist`, `ParmOff`, list assembly) — small
  per call but non-zero; measure before touching it.

Do not optimise by reading the code alone. The cost distribution depends
strongly on image size, PSF size and profile type. Measure first:

```r
bench = function(f, n = 20) { f(); unname(system.time(for (i in 1:n) f())["elapsed"]) / n }
```

Sampling the running process is the fastest way to attribute cost; the
built-in `stats` plumbing (below) is the alternative.

### Where the time goes, concretely

- `AssociativeBruteForceConvolver::convolve_impl` (`convolve.cpp`) — the SIMD
  dot-product convolver. Default when no convolver is supplied.
- `RadialProfile::subsample_pixel` (`radial.cpp`) — recursive sub-pixel
  integration. **This function is called per pixel, allocates on the heap per
  call, and recurses**; it is the highest-leverage target in the integrator.
- `SersicProfile::eval_function` (`sersic.cpp`) — templated specialisations
  avoid `pow()` for common `nser * (box + 2)` values. Preserve them.
- `CoreSersicProfile::evaluate_at` — recomputes several position-independent
  `pow()` terms at every sample unless cached in `initial_calculations()`.

### Cheap wins found by measurement (re-check before assuming they still hold)

1. **Convolver choice — fixed in 2.8.2.** `profitSetupData()` used to default to
   `"brute"` unconditionally, via `benchconvmethods[1]`. It now uses
   `.profitAutoConvolver()`: FFTW when the PSF area is at least 15x15 pixels,
   brute force below that. Worth ~1.3x on its own at 200-1500 px with a 25x25
   PSF, up to ~20x for a 51x51 one, and bit-identical output. The cutoff is on
   PSF *area*, not side length, because convolution work per output pixel
   scales with kernel area either way and the image dimensions cancel out.
   FFTW's own cost is dominated by how well the extended FFT dimensions
   (`src + psf - 1` per axis) factorise — `Model::evaluate` derives its
   crop=false bookkeeping from those, so they deliberately are *not* rounded
   up to nice sizes, and a size with a large prime factor can be several times
   slower than brute force. That is why the rule is deliberately conservative
   rather than trying to predict the crossover per size. `benchconvmethods`
   now defaults to `NULL`; set it to force a method, or `nbenchconv > 0` to
   benchmark properly.
2. **Thread propagation — fixed in 2.8.2.** `profitMakeConvolver()` defaults
   `omp_threads = 1L`, and `profitSetupData()` used not to pass its
   `omp_threads` argument through to the convolver, so a user could set
   `omp_threads = 8` and still get single-threaded convolution. The convolver
   is now created with `omp_threads` from `profitSetupData()`. Worth 2.6-2.9x
   on `profitLikeModel` at 200-1500 px, bit-identical output. Note the
   convolver is built *once* at setup time: setting `Data$omp_threads`
   afterwards threads the profiles but not the convolver.
3. **`subsample_pixel` overheads — fixed in 2.8.2.** It declared a
   `std::vector<unsigned int>` sized `resolution * resolution` that was never
   used (a heap allocation per subsampled pixel per recursion level), and
   recomputed a loop-invariant `delta_y_prof` inside the innermost loop. Both
   are gone, and the point buffer is now reserved rather than grown. Combined
   with invariant hoisting in `coresersic.cpp` this gave 10-22% on CoreSersic
   (largest for compact/bright cores) and a few percent on Sersic. Output
   identical to ~1e-16 relative, which is the floor for this kind of
   reordering — anything larger is a bug.
4. **Redundant image copies in the crop path — fixed in 2.8.2.**
   `Model::evaluate` crops in place (`Image::crop_in_place`) instead of
   `image = image.crop(...)`, which is also what allows a caller-supplied
   buffer to survive cropping. `Convolver::mask_and_crop` masked through
   `operator&`, which made a second copy of the already-cropped result; it now
   masks in place, halving the allocation traffic there. It deliberately still
   copies *out* of `img` rather than cropping `img` in place: for the FFT
   convolver `img` is the `ext_src` member, reused across every evaluation, and
   cropping it in place would force it to be reallocated each time. Note
   `crop_in_place` moves rows within existing storage, so a cropped `Image`
   keeps its original capacity; don't assume `size()` reflects capacity.
   Caveat: `mask_and_crop` is only reached by the FFT and OpenCL convolvers —
   brute force applies the mask inline — so with the default auto-selected
   convolver at small PSF sizes none of this is on the hot path.
5. **R likelihood stage — fixed in 2.8.2.** At megapixel sizes this was
   20-33% of `profitLikeModel`, and *most of that was `sum(dnorm(...))`*: R's
   `dnorm` was spending its time on per-element argument checks rather than
   arithmetic (105 ms of a 158 ms stage at 2000x2000, versus 26 ms for the
   closed form). It is now the closed form of
   `sum(dnorm(x, 0, scat_scale, log=TRUE))`, which agrees with `dnorm` to
   machine precision. The `scat_scale != 1` case must keep the *full* form —
   `cutsig` is already `(data - model) / sigma`, so the scatter enters both the
   squared term and the normalisation; dropping it from the squared term is a
   silent wrong-answer bug (see the note below). Worth 25-33% of the call down
   to 12-20%, i.e. 4-15% of `profitLikeModel` overall depending on size.

   Two things that **do not** work here, so don't retry them:

   * *Caching `Data$image[region_which]` and `Data$sigma[region_which]` on the
     `Data` object.* `profitLikeModel(parm, Data)` receives `Data` by value, so
     writing `Data$like_cutim = ...` modifies a local copy that is discarded on
     return. The cache is never populated, let alone reused. A cache that does
     work has to live outside the call (an environment, or a cache attached to
     `Data` at setup time), and attaching one at setup risks going stale when
     `Data$sigma` is edited in place — which is a supported pattern.
   * *`crossprod(x)` instead of `sum(x*x)`.* Once the `(cutim - cutmod)/cutsig`
     division is in the same expression the two are indistinguishable; all the
     time is in the divide and the memory traffic, not the reduction.

### Things that look like wins but are not

- Reallocating the `Image` per evaluation is *not* usually the bottleneck; the
  allocation itself is cheap relative to integration and convolution.
- The R-side likelihood arithmetic is a rounding error below ~200×200 (a few
  percent at 500×500). It only becomes worth touching at ~1000×1000 and above.
- `omp_2d_for` already falls back to a serial loop when `threads <= 1`; adding
  `#pragma omp` in more places does not help without also passing threads down.

## Known correctness hazards

All three were live bugs and were fixed in 2.8.2, each with a regression test
in `tests/`. Worth re-reading anyway: all three were invisible on the default
configuration (square PSF, square pixels, default buffer) and had survived for
a long time.

- **Image buffer reuse — fixed in 2.8.2.** `profitLikeModel` used to read
  `Data$model_buff_image` while `profitSetupData` wrote `Data$model_image_buff`,
  so the preallocated buffer was never passed down. `profitSetupData` also
  sized it with a dummy `null`-profile evaluation that has no convolution, so
  it came out at `dim(image)` rather than the drawing dimensions. The buffer is
  now `dim(calcregion)`, which is exactly what
  `Model::analyze_expansion_requirements` computes, and it is wired up. This
  only became reusable once `Model::evaluate` cropped in place;
  `_R_profit_make_model` also has to check dimensions before returning the
  caller's matrix, because a cropped-in-place buffer keeps its storage but
  shrinks its dimensions.
- **`PsfProfile::evaluate` (`psf.cpp`) axis swap — fixed in 2.8.2.** It mixed
  up `offset.x`/`offset.y` when computing `end_x` and `origin_y`, and used
  `psf_height`/`psf_width` as the wrong per-axis bound. With a non-square PSF
  and non-zero padding (i.e. any fit with an extended profile) this misplaced
  point sources by several pixels and lost ~17% of their flux. Covered by
  `tests/test-psf-profile.R`. Square PSFs never showed it, which is why it
  survived.
- **`RadialProfile::evaluate_cpu` scale swap — fixed in 2.8.2.** It computed
  the pixel's `y` from `scale.x` and its `x` from `scale.y`. Invisible at the
  default square scale (`scale_x == scale_y`), which is why it survived, but on
  a rectangular grid it transposes the two axes and the image is wrong by up to
  ~2.4x relative error, flux included. `scale_x`/`scale_y` are now reachable
  from R (`profitMakeModel()` and `profitSetupData()`, defaulting to 1) so this
  is testable; `tests/test-pixel-scale.R` compares the rendering against a
  direct analytic Sersic evaluation. Note that `xcen`/`ycen` are in *model*
  coordinates, so with `scale_y = 2` a centre of 50 lands at pixel row ~25 —
  that is correct, not a bug.

### The `scat_scale` trap in the likelihood

`cutsig` is already divided by `sigma`, so `dnorm(cutsig, sd = scat_scale)`
squares it *again* by `scat_scale`. Any closed-form rewrite of the `"norm"`
likelihood has to scale `cutsig` by `scat_scale` in the squared term as well as
adding the `log(scat_scale)` normalisation. Getting this wrong is invisible
whenever `scat_scale == 1`, which is the default for every fit that does not
set `log_scat_scale`, so it will pass any test that does not exercise that
option — `tests/test-likelihood.R` does.

## Code conventions

- R: base pipe `|>`, `=` for assignment as used throughout existing files,
  two-space indent, snake_case for arguments actually matches existing code
  (`finesample`, `magzero`); follow the surrounding file rather than
  introducing a new style.
- C++: C++17 (`-std=gnu++17`), tabs for indentation in `src/libprofit`,
  `std::` qualification, `override` on virtuals.
- Prefer `std::` over `using namespace` in new code; `r_profit.cpp` is the one
  file that already does `using namespace std`.
- Never change public R API signatures without a version bump in `DESCRIPTION`.

## Verifying numerical changes

Any change to integration or convolution must be checked numerically, not just
for speed. The pattern that works:

```r
a <- old_build()   # same model, previous code
b <- new_build()   # same model, changed code
max(abs(a - b)) / max(abs(a))   # expect ~1e-16 for algebraic refactors, 0 for none
```

For refactors that only reorder floating-point operations, expect differences
around `1e-16` relative. Anything larger in a refactor is a bug, not noise.
For changes to integration accuracy or the sub-sampling heuristics, compare
enclosed flux to a high-accuracy reference instead of comparing images
directly.

## Memory

`Data$image`, `Data$sigma`, `Data$mask`, `Data$calcregion` and the buffers are
all full-image arrays held for the whole fit. `finesample > 1` multiplies
drawing dimensions by the factor in each axis, so a `1000x1000` image at
`finesample = 2` draws `2048x2048` doubles (about 32 MB) per profile buffer.
For large mosaics this, not the profile maths, is what limits how many fits
can run concurrently. Consider `Rprofmem` before adding allocations to the
per-evaluation path.

## Remote setup note

`DESCRIPTION` `Remotes:` lists several `asgr/*` packages. `Rfits`, `magicaxis`,
`celestial` and `ParmOff` are hard dependencies; without them the package will
not load, and `.Call`-level benchmarks cannot run.