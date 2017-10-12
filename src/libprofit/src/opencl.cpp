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
#include <streambuf>
#include <sstream>
#include <string>
#include <vector>
#include <sys/time.h>

#include "profit/exceptions.h"
#include "profit/opencl.h"

#ifdef PROFIT_OPENCL

namespace profit {

OpenCL_command_times::OpenCL_command_times() :
	submit(0), exec(0)
{
	// no-op
}

OpenCL_command_times &OpenCL_command_times::operator+=(const OpenCL_command_times &other) {
	submit += other.submit;
	exec += other.exec;
	return *this;
}

const OpenCL_command_times OpenCL_command_times::operator+(const OpenCL_command_times &other) const {
	OpenCL_command_times t1;
	t1 += other;
	return t1;
}

OpenCL_times::OpenCL_times() :
	kernel_prep(0), nwork_items(0),
	writing_times(), reading_times(), filling_times(), kernel_times(),
	total(0)
{
	// no-op
}

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
nsecs_t _cl_exec_time(const cl::Event &evt) {
	return _cl_duration<CL_PROFILING_COMMAND_START, CL_PROFILING_COMMAND_END>(evt);
}

OpenCL_command_times cl_cmd_times(const cl::Event &evt) {
	OpenCL_command_times times;
	times.submit = _cl_submit_time(evt);
	times.exec = _cl_exec_time(evt);
	return times;
}

static cl_ver_t get_opencl_version(const cl::Platform &platform) {

	std::string version = platform.getInfo<CL_PLATFORM_VERSION>();

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

	unsigned long major = stoul(opencl_version.substr(0, dot_idx));
	unsigned long minor = stoul(opencl_version.substr(dot_idx+1, opencl_version.npos));
	return major*100u + minor*10u;
}

static
std::map<int, OpenCL_plat_info> _get_opencl_info() {

	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);

	std::map<int, OpenCL_plat_info> pinfo;
	unsigned int pidx = 0;
	for(auto platform: all_platforms) {
		std::vector<cl::Device> devices;
		platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

		std::map<int, OpenCL_dev_info> dinfo;
		unsigned int didx = 0;
		for(auto device: devices) {
			dinfo[didx] = OpenCL_dev_info{
				device.getInfo<CL_DEVICE_NAME>(),
				device.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG>() != 0
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

static
OpenCLEnvPtr _get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double, bool enable_profiling) {

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

	cl::Platform platform = all_platforms[platform_idx];

	//get default device of the default platform
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

	cl::Device device = all_devices[device_idx];

	if( use_double ) {
		auto config = device.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG>();
		if( config == 0 ) {
			throw opencl_error("Double precision requested but not supported by device");
		}
	}

	// Create a program with all the relevant kernels
	// The source of the kernels is kept in a different file,
	// but gets #included here, and thus gets embedded in the resulting
	// shared library (instead of, for instance, loading the sources from a
	// particular location on disk at runtime)
	const char *common_float =
#include "cl/common-float.cl"
	;
	const char *common_double =
#include "cl/common-double.cl"
	;
	const char *sersic_float =
#include "cl/sersic-float.cl"
	;
	const char *sersic_double =
#include "cl/sersic-double.cl"
	;
	const char *moffat_float =
#include "cl/moffat-float.cl"
	;
	const char *moffat_double =
#include "cl/moffat-double.cl"
	;
	const char *ferrer_float =
#include "cl/ferrer-float.cl"
	;
	const char *ferrer_double =
#include "cl/ferrer-double.cl"
	;
	const char *king_float =
#include "cl/king-float.cl"
	;
	const char *king_double =
#include "cl/king-double.cl"
	;
	const char *brokenexp_float =
#include "cl/brokenexponential-float.cl"
	;
	const char *brokenexp_double =
#include "cl/brokenexponential-double.cl"
	;
	const char *coresersic_float =
#include "cl/coresersic-float.cl"
	;
	const char *coresersic_double =
#include "cl/coresersic-double.cl"
	;
	const char *convolve_float =
#include "cl/convolve-float.cl"
	;
	const char *convolve_double =
#include "cl/convolve-double.cl"
	;

	cl::Program::Sources sources;
	sources.push_back(common_float);
	sources.push_back(sersic_float);
	sources.push_back(moffat_float);
	sources.push_back(ferrer_float);
	sources.push_back(king_float);
	sources.push_back(brokenexp_float);
	sources.push_back(coresersic_float);
	sources.push_back(convolve_float);
	if( use_double ) {
		sources.push_back(common_double);
		sources.push_back(sersic_double);
		sources.push_back(moffat_double);
		sources.push_back(ferrer_double);
		sources.push_back(king_double);
		sources.push_back(brokenexp_double);
		sources.push_back(coresersic_double);
		sources.push_back(convolve_double);
	}

	cl::Context context(device);
	cl::Program program(context, sources);
	try {
		program.build({device});
	} catch (const cl::Error &e) {
		throw opencl_error("Error building program: " + program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
	}

	cl::CommandQueue queue(context, device, enable_profiling ? CL_QUEUE_PROFILING_ENABLE : 0);

	return std::make_shared<OpenCL_env>(device, get_opencl_version(platform), context, queue, program, use_double, enable_profiling);
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

unsigned long OpenCL_env::max_local_memory() {
	return device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
}

unsigned int OpenCL_env::compute_units() {
	return device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
}

cl::Event OpenCL_env::queue_write(const cl::Buffer &buffer, const void *data, const std::vector<cl::Event>* wait_evts) {
	cl::Event wevt;
	queue.enqueueWriteBuffer(buffer, CL_FALSE, 0, buffer.getInfo<CL_MEM_SIZE>(), data, wait_evts, &wevt);
	return wevt;
}

cl::Event OpenCL_env::queue_kernel(const cl::Kernel &kernel, const cl::NDRange global, const std::vector<cl::Event>* wait_evts, const cl::NDRange &local) {
	cl::Event kevt;
	queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local, wait_evts, &kevt);
	return kevt;
}

cl::Event OpenCL_env::queue_read(const cl::Buffer &buffer, void *data, const std::vector<cl::Event>* wait_evts) {
	cl::Event revt;
	queue.enqueueReadBuffer(buffer, CL_FALSE, 0, buffer.getInfo<CL_MEM_SIZE>(), data, wait_evts, &revt);
	return revt;
}

cl::Kernel OpenCL_env::get_kernel(const std::string &name) {
	return cl::Kernel(program, name.c_str());
}

}

#endif /* PROFIT_OPENCL */
