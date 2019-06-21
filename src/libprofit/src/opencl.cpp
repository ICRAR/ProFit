/**
 * OpenCL utility methods for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <map>
#include <streambuf>
#include <sstream>
#include <string>
#include <mutex>
#include <vector>

#include "profit/crc.h"
#include "profit/exceptions.h"
#include "profit/opencl_impl.h"
#include "profit/utils.h"

// The individual kernel sources, generated from the coresponding .cl files
#ifdef PROFIT_OPENCL
#include "profit/cl/common-float.h"
#include "profit/cl/common-double.h"
#include "profit/cl/sersic-float.h"
#include "profit/cl/sersic-double.h"
#include "profit/cl/moffat-float.h"
#include "profit/cl/moffat-double.h"
#include "profit/cl/ferrer-float.h"
#include "profit/cl/ferrer-double.h"
#include "profit/cl/king-float.h"
#include "profit/cl/king-double.h"
#include "profit/cl/brokenexponential-float.h"
#include "profit/cl/brokenexponential-double.h"
#include "profit/cl/coresersic-float.h"
#include "profit/cl/coresersic-double.h"
#include "profit/cl/convolve-float.h"
#include "profit/cl/convolve-double.h"
#endif // PROFIT_OPENCL

namespace profit {

OpenCL_command_times &OpenCL_command_times::operator+=(const OpenCL_command_times &other) {
	submit += other.submit;
	wait += other.wait;
	exec += other.exec;
	return *this;
}

const OpenCL_command_times OpenCL_command_times::operator+(const OpenCL_command_times &other) const {
	OpenCL_command_times t1 = *this;
	t1 += other;
	return t1;
}

// Simple implementation of public methods for non-OpenCL builds
#ifndef PROFIT_OPENCL
std::map<int, OpenCL_plat_info> get_opencl_info() {
	return std::map<int, OpenCL_plat_info>();
}

OpenCLEnvPtr get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double, bool enable_profiling)
{
	UNUSED(platform_idx);
	UNUSED(device_idx);
	UNUSED(use_double);
	UNUSED(enable_profiling);
	return nullptr;
}

#else

// Functions to read the duration of OpenCL events (queue->submit and start->end)
template <cl_int S, cl_int E>
static inline
nsecs_t _cl_duration(const cl::Event &evt) {
	auto start = evt.getProfilingInfo<S>();
	auto end = evt.getProfilingInfo<E>();
	if( start > end ) {
		return 0;
	}
	return end - start;
}

static
nsecs_t _cl_submit_time(const cl::Event &evt) {
	return _cl_duration<CL_PROFILING_COMMAND_QUEUED, CL_PROFILING_COMMAND_SUBMIT>(evt);
}

static
nsecs_t _cl_wait_time(const cl::Event &evt) {
	return _cl_duration<CL_PROFILING_COMMAND_SUBMIT, CL_PROFILING_COMMAND_START>(evt);
}

static
nsecs_t _cl_exec_time(const cl::Event &evt) {
	return _cl_duration<CL_PROFILING_COMMAND_START, CL_PROFILING_COMMAND_END>(evt);
}

OpenCL_command_times cl_cmd_times(const cl::Event &evt) {
	return {_cl_submit_time(evt), _cl_wait_time(evt), _cl_exec_time(evt)};
}

static cl_ver_t get_opencl_version(const std::string &version)
{
	// Version string should be of type "OpenCL<space><major_version.minor_version><space><platform-specific information>"

	if( version.find("OpenCL ") != 0) {
		throw opencl_error(std::string("OpenCL version string doesn't start with 'OpenCL ': ") + version);
	}

	auto next_space = version.find(" ", 7);
	auto opencl_version = version.substr(7, next_space);
	auto dot_idx = opencl_version.find(".");
	if( dot_idx == opencl_version.npos ) {
		throw opencl_error("OpenCL version doesn't contain a dot: " + opencl_version);
	}

	auto major = stoui(opencl_version.substr(0, dot_idx));
	auto minor = stoui(opencl_version.substr(dot_idx+1, opencl_version.npos));
	return major * 100U + minor * 10U;
}

static cl_ver_t get_opencl_version(const cl::Platform &platform) {
	return get_opencl_version(platform.getInfo<CL_PLATFORM_VERSION>());
}

static cl_ver_t get_opencl_version(const cl::Device &device) {
	return get_opencl_version(device.getInfo<CL_DEVICE_VERSION>());
}

static
bool supports_double(const cl::Device &device)
{
	if (get_opencl_version(device) < 120) {
		std::string extensions = device.getInfo<CL_DEVICE_EXTENSIONS>();
		return extensions.find("cl_khr_fp64") != std::string::npos;
	}

	return device.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG>() != 0;
}

static
std::map<int, OpenCL_plat_info> _get_opencl_info() {

	std::vector<cl::Platform> all_platforms;
	try {
		cl::Platform::get(&all_platforms);
	} catch (const cl::Error &e) {
		// No platform found by ICD loader, we tolerate that
		if (e.err() != CL_PLATFORM_NOT_FOUND_KHR) {
			throw;
		}
	}

	std::map<int, OpenCL_plat_info> pinfo;
	unsigned int pidx = 0;
	for(auto platform: all_platforms) {
		std::vector<cl::Device> devices;

		try {
			platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
		} catch (const cl::Error &e) {
			// If no devices are found we simply return an empty devices list.
			if (e.err() != CL_DEVICE_NOT_FOUND) {
				throw;
			}
		}

		std::map<int, OpenCL_dev_info> dinfo;
		unsigned int didx = 0;
		for(auto device: devices) {
			dinfo[didx++] = OpenCL_dev_info{
				device.getInfo<CL_DEVICE_NAME>(),
				get_opencl_version(device),
				supports_double(device)
			};
		}

		std::string name = platform.getInfo<CL_PLATFORM_NAME>();
		pinfo[pidx++] = OpenCL_plat_info{name, get_opencl_version(platform), dinfo};
	}

	return pinfo;
}

std::map<int, OpenCL_plat_info> get_opencl_info() {

	// Wrap cl::Error exceptions
	try {
		return _get_opencl_info();
	} catch(const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

class invalid_cache_entry : public std::exception {
};

class KernelCache {

	typedef std::vector<unsigned char> Binary;
	typedef std::vector<Binary> Binaries;
	typedef std::pair<std::string, uint32_t> SourceInformation;

public:
	cl::Program get_program(const cl::Context &context, const cl::Device &device);

private:
	std::string get_entry_name_for(const cl::Device &device);
	cl::Program build(const cl::Context &context, const cl::Device &device, const SourceInformation &source_info);
	cl::Program from_cache(const cl::Context &context, const std::string &cache_entry_name, const SourceInformation &source_info);
	void to_cache(const std::string &cache_entry_name, const SourceInformation &source_info, const cl::Program &program);

	static void init_sources();

	// Lazy-initialized via std::call_once
	static void _init_sources();
	static SourceInformation float_only_sources;
	static SourceInformation all_sources;
	static std::once_flag init_sources_flag;
};

KernelCache::SourceInformation KernelCache::float_only_sources;
KernelCache::SourceInformation KernelCache::all_sources;
std::once_flag KernelCache::init_sources_flag;
void KernelCache::init_sources() {
	std::call_once(init_sources_flag, _init_sources);
}

void KernelCache::_init_sources() {

	auto sources = common_float + sersic_float + moffat_float + ferrer_float +
	               king_float + brokenexponential_float + coresersic_float + convolve_float;
	float_only_sources = std::make_pair(sources, crc32(sources));

	sources += common_double + sersic_double + moffat_double + ferrer_double +
	           king_double + brokenexponential_double + coresersic_double + convolve_double;
	all_sources = std::make_pair(sources, crc32(sources));
}

static
std::string &valid_fname(std::string &&name)
{
	std::string::size_type pos = 0;
	while ((pos = name.find_first_of("/;: ", pos)) != name.npos) {
		name.replace((pos++), 1, "_");
	}
	return name;
}

std::string KernelCache::get_entry_name_for(const cl::Device &device)
{
	cl::Platform plat(device.getInfo<CL_DEVICE_PLATFORM>());
	auto plat_part = valid_fname(plat.getInfo<CL_PLATFORM_NAME>()) + "_" + std::to_string(get_opencl_version(plat));
	auto dev_part = valid_fname(device.getInfo<CL_DEVICE_NAME>()) + "_" + std::to_string(get_opencl_version(device));
	auto the_dir = create_dirs(get_profit_home(), {std::string("opencl_cache"), plat_part});
	return the_dir + "/" + dev_part;
}

cl::Program KernelCache::build(const cl::Context &context, const cl::Device &device, const SourceInformation &source_info)
{

	// Create a program with all the relevant kernels
	// The source of the kernels is kept in a different file,
	// but gets #included here, and thus gets embedded in the resulting
	// shared library (instead of, for instance, loading the sources from a
	// particular location on disk at runtime)

	cl::Program program(context, {std::get<0>(source_info)});
	try {
		program.build({device});
	} catch (const cl::Error &e) {
		throw opencl_error("Error building program: " + program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
	}

	return program;
}

cl::Program KernelCache::from_cache(const cl::Context &context, const std::string &cache_entry_name, const SourceInformation &source_info)
{

	if (!file_exists(cache_entry_name)) {
		throw invalid_cache_entry();
	}

	std::ifstream input(cache_entry_name, std::ios::binary);

	// Source length. If different from expected, cache entry is invalid
	Binary::size_type ssize;
	input.read(reinterpret_cast<char *>(&ssize), sizeof(ssize));
	if (ssize != std::get<0>(source_info).size()) {
		throw invalid_cache_entry();
	}

	// Source checksum. If different from expected, cache entry is invalid
	uint32_t bchecksum;
	input.read(reinterpret_cast<char *>(&bchecksum), sizeof(bchecksum));
	if (bchecksum != std::get<1>(source_info)) {
		throw invalid_cache_entry();
	}

	// Cache entry is valid, return its contents (size + data for each of the binaries)
	Binaries::size_type nbinaries;
	input.read(reinterpret_cast<char *>(&nbinaries), sizeof(nbinaries));
	Binaries binaries;
	for(unsigned int i = 0; i < nbinaries; i++) {
		Binary::size_type bsize;
		input.read(reinterpret_cast<char *>(&bsize), sizeof(bsize));
		Binary binary(bsize);
		input.read(reinterpret_cast<char *>(binary.data()), bsize);
		binaries.push_back(std::move(binary));
	}

	std::vector<cl_int> binaries_status;
	auto program = cl::Program(context, context.getInfo<CL_CONTEXT_DEVICES>(), binaries, &binaries_status, nullptr);
	for(auto binary_status: binaries_status) {
		if (binary_status != CL_SUCCESS) {
			std::ostringstream os;
			os << "Error when loading binary from kernel cache: " << binary_status;
			throw opencl_error(os.str());
		}
	}

	// build it and good-bye
	program.build();
	return program;
}

void KernelCache::to_cache(const std::string &cache_entry_name, const SourceInformation &source_info, const cl::Program &program)
{
	std::ofstream output(cache_entry_name, std::ios::binary);

	// Source length
	auto ssize = std::get<0>(source_info).size();
	output.write(reinterpret_cast<char *>(&ssize), sizeof(ssize));

	// Source checksum
	auto bchecksum = std::get<1>(source_info);
	output.write(reinterpret_cast<char *>(&bchecksum), sizeof(bchecksum));

	// Binaries (#, then size + data for each)
	auto binaries = program.getInfo<CL_PROGRAM_BINARIES>();
	auto nbinaries = binaries.size();
	output.write(reinterpret_cast<char *>(&nbinaries), sizeof(nbinaries));
	for(auto binary: binaries) {
		auto bsize = binary.size();
		output.write(reinterpret_cast<char *>(&bsize), sizeof(bsize));
		output.write(reinterpret_cast<const char *>(binary.data()), bsize);
	}
}


cl::Program KernelCache::get_program(const cl::Context &context, const cl::Device &device)
{

	// Make sure we know all sources and their checksums
	init_sources();

	bool device_supports_double = supports_double(device);
	SourceInformation sources_for_device = float_only_sources;
	if (device_supports_double) {
		sources_for_device = all_sources;
	}

	// Act as a cache! Load compiled program from file if it exists,
	// otherwise build it from source
	auto cache_entry = get_entry_name_for(device);
	try {
		return from_cache(context, cache_entry, sources_for_device);
	} catch (const invalid_cache_entry &e) {
		auto program = build(context, device, sources_for_device);
		to_cache(cache_entry, sources_for_device, program);
		return program;
	}

}


// Lazy initialization of singleton via static local variable
KernelCache get_cache() {
	static KernelCache cache;
	return cache;
}

static
cl::Platform get_platform(unsigned int platform_idx)
{
	std::vector<cl::Platform> all_platforms;
	if( cl::Platform::get(&all_platforms) != CL_SUCCESS ) {
		throw opencl_error("Error while getting OpenCL platforms");
	}
	if( all_platforms.size() == 0 ){
		throw opencl_error("No platforms found. Check OpenCL installation");
	}

	if( platform_idx >= all_platforms.size() ) {
		std::ostringstream ss;
		ss << "OpenCL platform index " << platform_idx << " must be < " << all_platforms.size();
		throw invalid_parameter(ss.str());
	}

	return all_platforms[platform_idx];
}

static
cl::Device get_device(const cl::Platform &platform, unsigned int device_idx, bool use_double)
{
	std::vector<cl::Device> all_devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if( all_devices.size() == 0 ){
		throw opencl_error("No devices found. Check OpenCL installation");
	}
	if( device_idx >= all_devices.size() ) {
		std::ostringstream ss;
		ss << "OpenCL device index " << device_idx << " must be < " << all_devices.size();
		throw invalid_parameter(ss.str());
	}

	auto device = all_devices[device_idx];
	if( use_double && !supports_double(device)) {
		throw opencl_error("Double precision requested but not supported by device");
	}
	return device;
}

static
OpenCLEnvPtr _get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double, bool enable_profiling) {

	auto platform = get_platform(platform_idx);
	auto device = get_device(platform, device_idx, use_double);

	cl::Context context(device);

	// Check if there is an entry in the cache for this platform/device
	// This considers the version information of each, so if we update the
	// platform we need to regenerate the cache entry
	KernelCache cache = get_cache();
	auto program = cache.get_program(context, device);

	cl::CommandQueue queue(context, device, enable_profiling ? CL_QUEUE_PROFILING_ENABLE : 0);

	return std::make_shared<OpenCLEnvImpl>(device, get_opencl_version(device), context, queue, program, use_double, enable_profiling);
}

OpenCLEnvPtr get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double, bool enable_profiling) {

	// Wrap cl::Error exceptions
	try {
		return _get_opencl_environment(platform_idx, device_idx, use_double, enable_profiling);
	} catch(const cl::Error &e) {
		std::ostringstream os;
		os << "OpenCL error: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

unsigned long OpenCLEnvImpl::max_local_memory() {
	return device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
}

unsigned int OpenCLEnvImpl::compute_units() {
	return device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
}

cl::Event OpenCLEnvImpl::queue_write(const cl::Buffer &buffer, const void *data, const std::vector<cl::Event>* wait_evts) {
	cl::Event wevt;
	queue.enqueueWriteBuffer(buffer, CL_FALSE, 0, buffer.getInfo<CL_MEM_SIZE>(), data, wait_evts, &wevt);
	return wevt;
}

cl::Event OpenCLEnvImpl::queue_kernel(const cl::Kernel &kernel, const cl::NDRange global, const std::vector<cl::Event>* wait_evts, const cl::NDRange &local) {
	cl::Event kevt;
	queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, wait_evts, &kevt);
	return kevt;
}

cl::Event OpenCLEnvImpl::queue_read(const cl::Buffer &buffer, void *data, const std::vector<cl::Event>* wait_evts) {
	cl::Event revt;
	queue.enqueueReadBuffer(buffer, CL_FALSE, 0, buffer.getInfo<CL_MEM_SIZE>(), data, wait_evts, &revt);
	return revt;
}

cl::Kernel OpenCLEnvImpl::get_kernel(const std::string &name) {
	return cl::Kernel(program, name.c_str());
}

#endif /* PROFIT_OPENCL */

} // namespace profit