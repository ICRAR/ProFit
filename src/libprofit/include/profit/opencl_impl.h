/**
 * Internal header file for OpenCL functionality
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

#ifndef PROFIT_OPENCL_IMPL_H
#define PROFIT_OPENCL_IMPL_H

#include "profit/opencl.h"

#ifdef PROFIT_OPENCL
/* Quickly fail for OpenCL < 1.1 */
# if !defined(PROFIT_OPENCL_MAJOR) || !defined(PROFIT_OPENCL_MINOR)
#  error "No OpenCL version specified"
# elif PROFIT_OPENCL_MAJOR < 1 || (PROFIT_OPENCL_MAJOR == 1 && PROFIT_OPENCL_MINOR < 1 )
#  error "libprofit requires at minimum OpenCL >= 1.1"
# endif

/* MacOS 10.14 (Mojave) started deprecating OpenCL entirely */
#if defined(__APPLE__) || defined(__MACOSX)
# define CL_SILENCE_DEPRECATION
#endif

/* We use exceptions in our code */
# define CL_HPP_ENABLE_EXCEPTIONS

/* MacOS, and maybe other platforms might not define this */
#if !defined(CL_PLATFORM_NOT_FOUND_KHR)
# define CL_PLATFORM_NOT_FOUND_KHR -1001
#endif

/* Define the target OpenCL version based on the given major/minor version */
# define PASTE(x,y) x ## y ## 0
# define MAKE_VERSION(x,y) PASTE(x,y)
# define CL_HPP_TARGET_OPENCL_VERSION  MAKE_VERSION(PROFIT_OPENCL_MAJOR, PROFIT_OPENCL_MINOR)
# define CL_TARGET_OPENCL_VERSION  MAKE_VERSION(PROFIT_OPENCL_MAJOR, PROFIT_OPENCL_MINOR)
# define CL_HPP_MINIMUM_OPENCL_VERSION PROFIT_OPENCL_TARGET_VERSION

/*
 * GCC 6 gives lots of "ignoring attributes on template arguments" warnings
 * because the built-in types have alignment attributes, and when these types
 * are used as template parameters (e.g., when creating a vector<cl_int>) these
 * attribute information is lost. This isn't an issue in reality though, and
 * therefore we can safely disable this warning.
 *
 * Clang also complains about missing brances during the initialization of
 * cl_image_desc structures because some systems define an anonymous union for
 * the last member (which according to the spec is cl_mem buffer). Again,
 * there's nothing we can do about this structure definition, and therefore we
 * simply shut up the warnings.
 */
#if defined __GNUC__ && __GNUC__>=6
# pragma GCC diagnostic ignored "-Wignored-attributes"
#elif defined __clang__
# pragma clang diagnostic ignored "-Wmissing-braces"
#endif
#include "profit/cl/cl2.hpp"

namespace profit
{

/**
 * Returns the time spent submitting and executing the command associated to the
 * given event as an OpenCL_comand_times structure.
 * @param evt An event associated to a command
 * @return A structure holding the times spent submitting and executing the
 *         command associated to the given event
 *         (i.e., CL_PROFILING_COMMAND_SUBMIT - CL_PROFILING_COMMAND_QUEUED, and
 *          CL_PROFILING_COMMAND_END - CL_PROFILING_COMMAND_START).
 */
OpenCL_command_times cl_cmd_times(const cl::Event &evt);

class OpenCLEnvImpl;
typedef std::shared_ptr<OpenCLEnvImpl> OpenCLEnvImplPtr;

/**
 * An OpenCL environment
 *
 * This structure holds all the required information to make libprofit work
 * against a given device in a particular platform.
 */
class OpenCLEnvImpl : public OpenCLEnv {

public:

	OpenCLEnvImpl(cl::Device device, cl_ver_t version, cl::Context context,
	              cl::CommandQueue queue, cl::Program program,
	              bool use_double, bool use_profiling) :
		use_double(use_double), use_profiling(use_profiling),
		device(device), version(version), context(context), queue(queue), program(program)
	{ }

	static OpenCLEnvImplPtr fromOpenCLEnvPtr(const OpenCLEnvPtr &ptr) {
		return std::static_pointer_cast<OpenCLEnvImpl>(ptr);
	}

	/**
	 * Whether double floating-point precision has been requested on this device
	 * or not.
	 */
	bool use_double;

	/**
	 * Whether profiling information should be gathered or not
	 */
	bool use_profiling;

	// Implementing OpenCL_env's interface
	cl_ver_t get_version() override {
		return version;
	}

	// Implementing OpenCL_env's interface
	std::string get_platform_name() override {
		return cl::Platform(device.getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>();
	}

	// Implementing OpenCL_env's interface
	std::string get_device_name() override {
		return device.getInfo<CL_DEVICE_NAME>();
	}

	/**
	 * Returns the amount of memory, in bytes, that each OpenCL Compute Unit
	 * has.
	 */
	unsigned long max_local_memory();

	/**
	 * Returns the number of Computer Units available in the device wrapped
	 * by this OpenCL environment.
	 */
	unsigned int compute_units();

	/**
	 * Returns a buffer that can hold `n_elements` elements of type `T`.
	 * The buffer is created with the given `flags`.
	 */
	template <typename T>
	cl::Buffer get_buffer(int flags, cl::size_type n_elements) {
		return cl::Buffer(context, flags, sizeof(T) * n_elements);
	}

	/**
	 * Queues a write of `data` into `buffer` and returns the generated event.
	 */
	cl::Event queue_write(const cl::Buffer &buffer, const void *data, const std::vector<cl::Event>* wait_evts = NULL);

	/**
	 * Queues the execution of `kernel` in the `global` NDRange. The execution
	 * waits on `wait_evts` before commencing. It also uses the `local` NDRange
	 * to control the work group sizes.
	 */
	cl::Event queue_kernel(const cl::Kernel &kernel, const cl::NDRange global,
	                       const std::vector<cl::Event>* wait_evts = NULL,
	                       const cl::NDRange &local = cl::NullRange);

	/**
	 * Queues a read of `buffer` into `data` and returns the generated event.
	 */
	cl::Event queue_read(const cl::Buffer &buffer, void *data, const std::vector<cl::Event>* wait_evts = NULL);

#if CL_HPP_TARGET_OPENCL_VERSION >= 120
	template <typename PatternType>
	cl::Event queue_fill(const cl::Buffer &buffer, PatternType pattern, const std::vector<cl::Event>* wait_evts = NULL) {
		cl::Event fevt;
		queue.enqueueFillBuffer(buffer, pattern, 0, buffer.getInfo<CL_MEM_SIZE>(), wait_evts, &fevt);
		return fevt;
	}
#endif

	/**
	 * Get a reference to the named kernel
	 */
	cl::Kernel get_kernel(const std::string &name);

private:

	/** The device to be used throughout OpenCL operations */
	cl::Device device;

	/** The OpenCL supported by the platform this device belongs to */
	cl_ver_t version;

	/** The OpenCL context used throughout the OpenCL operations */
	cl::Context context;

	/** The queue set up against this device to be used by libprofit */
	cl::CommandQueue queue;

	/**
	 * The set of kernels and routines compiled against this device and
	 * required by libprofit
	 */
	cl::Program program;

};

} /* namespace profit */

#endif /* PROFIT_OPENCL */

#endif /* PROFIT_OPENCL_IMPL_H */