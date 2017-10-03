/**
 * Header file for OpenCL functionality
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

#ifdef PROFIT_OPENCL

#ifndef PROFIT_OPENCL_H
#define PROFIT_OPENCL_H

#include <map>
#include <memory>
#include <string>

#include "profit/config.h"
#include "profit/common.h"

/*
 * We only give specific inputs to cl2.hpp when building the library,
 * but not when including this file from an external program
 */
#ifdef PROFIT_BUILD

/* Quickly fail for OpenCL < 1.1 */
# if !defined(PROFIT_OPENCL_MAJOR) || !defined(PROFIT_OPENCL_MINOR)
#  error "No OpenCL version specified"
# elif PROFIT_OPENCL_MAJOR < 1 || (PROFIT_OPENCL_MAJOR == 1 && PROFIT_OPENCL_MINOR < 1 )
#  error "libprofit requires at minimum OpenCL >= 1.1"
# endif

/* We use exceptions in our code */
# define CL_HPP_ENABLE_EXCEPTIONS

#endif /* PROFIT_BUILD */

/* Define the target OpenCL version based on the given major/minor version */
# define PASTE(x,y) x ## y ## 0
# define MAKE_VERSION(x,y) PASTE(x,y)
# define CL_HPP_TARGET_OPENCL_VERSION  MAKE_VERSION(PROFIT_OPENCL_MAJOR, PROFIT_OPENCL_MINOR)
# define CL_HPP_MINIMUM_OPENCL_VERSION 110

/*
 * GCC 6 gives lots of "ignoring attributes on template arguments" warnings
 * Until the cl2.hpp header files doesn't get a proper fix (it's taken from
 * the official Kronos github) we simply turn the warnings off.
 */
#if defined __GNUC__ && __GNUC__>=6
# pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include "profit/cl/cl2.hpp"

namespace profit
{

/**
 * A datatype for storing an OpenCL version.
 * It should have the form major*100 + minor*10 (e.g., 120 for OpenCL 1.2)
 */
typedef unsigned int cl_ver_t;

/**
 * A structure holding two times associated with OpenCL commands:
 * submission and execution
 */
struct OpenCL_command_times {
	OpenCL_command_times();
	nsecs_t submit;
	nsecs_t exec;
	OpenCL_command_times &operator+=(const OpenCL_command_times &other);
	const OpenCL_command_times operator+(const OpenCL_command_times &other) const;
};

/**
 * A structure holding a number of OpenCL command times (filling, writing,
 * kernel and reading) plus other OpenCL-related times.
 */
struct OpenCL_times {
	OpenCL_times();
	nsecs_t kernel_prep;
	unsigned int nwork_items;
	OpenCL_command_times writing_times;
	OpenCL_command_times reading_times;
	OpenCL_command_times filling_times;
	OpenCL_command_times kernel_times;
	nsecs_t total;
};

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

/**
 * An OpenCL environment
 *
 * This structure holds all the required information to make libprofit work
 * against a given device in a particular platform.
 */
class OpenCL_env {

public:

	OpenCL_env(cl::Device device, cl_ver_t version, cl::Context context,
	           cl::CommandQueue queue, cl::Program program,
	           bool use_double, bool use_profiling) :
		use_double(use_double), use_profiling(use_profiling),
		device(device), version(version), context(context), queue(queue), program(program)
	{ }

	/**
	 * Whether double floating-point precision has been requested on this device
	 * or not.
	 */
	bool use_double;

	/**
	 * Whether profiling information should be gathered or not
	 */
	bool use_profiling;

	/**
	 * Returns the maximum OpenCL version supported by the underlying device.
	 */
	cl_ver_t get_version() {;
		return version;
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

/**
 * A structure holding information about a specific OpenCL device
 */
typedef struct _OpenCL_dev_info {

	/** The name of the device */
	std::string name;

	/** Whether or not this device supports double floating-point precision */
	bool double_support;

} OpenCL_dev_info;

/**
 * An structure holding information about a specific OpenCL platform.
 */
typedef struct _OpenCL_plat_info {

	/** The name of the platform */
	std::string name;

	/** The supported OpenCL version */
	cl_ver_t supported_opencl_version;

	/** A map containing information about all devices on this platform */

	std::map<int, OpenCL_dev_info> dev_info;
} OpenCL_plat_info;

/**
 * Queries the system about the OpenCL supported platforms and devices and returns
 * the information the caller.
 *
 * @return A map keyed by index, containing the information of each of the
 *         OpenCL platforms found on this system.
 */
std::map<int, OpenCL_plat_info> get_opencl_info();

/**
 * Prepares an OpenCL working space for using with libprofit.
 *
 * This method will get the requested device on the requested platform, compile
 * the libprofit OpenCL kernel sources to be used against it, and set up a queue
 * on the device.
 *
 * @param platform_idx The index of the platform to use
 * @param device_idx The index of device to use in the platform
 * @param use_double Whether double floating-point support should be used in
 *        the device or not.
 * @param enable_profiling Whether OpenCL profiling capabilities should be
 *        turned on in the OpenCL Queue created within this envinronment.
 * @return A pointer to a OpenCL_env structure, which contains the whole set of
 *         elements required to work with the requested device.
 */
std::shared_ptr<OpenCL_env> get_opencl_environment(
	unsigned int platform_idx,
	unsigned int device_idx,
	bool use_double,
	bool enable_profiling);

} /* namespace profit */

#endif /* PROFIT_OPENCL_H */

#endif /* PROFIT_OPENCL */
