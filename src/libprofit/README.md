# libprofit

[![Travis Build Status](https://travis-ci.org/SKA-ScienceDataProcessor/dfms.svg?branch=master)](https://travis-ci.org/ICRAR/libprofit)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/mr42dlpq9vaeqnln?svg=true)](https://ci.appveyor.com/project/rtobar/libprofit)
[![Coverage Status](https://coveralls.io/repos/github/ICRAR/libprofit/badge.svg?branch=master)](https://coveralls.io/github/ICRAR/libprofit?branch=master)
[![Docs](https://readthedocs.org/projects/libprofit/badge/?version=latest)](https://libprofit.readthedocs.io/en/latest/)
[![SonarQube Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=libprofit&metric=alert_status)](https://sonarcloud.io/dashboard?id=libprofit)

libprofit is a low-level C++ library that produces images based on different luminosity profiles.

libprofit currently supports the following profiles:
**sersic**, **coresersic**, **moffat**, **ferrer**, **king**, **brokenexp**, **sky** and **psf**.
With time we intend to add more profiles to the library. Users can also provide their own profiles.

For more information read [libprofit's documentation](https://libprofit.readthedocs.io/).
