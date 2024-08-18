/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__RESULTS_REAL_TIME
#define INQ__INTERFACE__RESULTS_REAL_TIME

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <interface/actions.hpp>
#include <input/environment.hpp>
#include <real_time/results.hpp>

#include <utility>

namespace inq {
namespace interface {

struct {		

	constexpr auto name() const {
		return "results real-time";
	}

	constexpr auto one_line() const {
		return "Get information about the results obtained from a real-time calculation";
	}

	constexpr auto help() const {
		
		return R""""(

The 'results real-time' command
==================

This command queries the results obtained from a real-time
calculation. Without arguments, it prints the values calculated.

The options allows you to query a specific value. In this case only
the value will be printed without any other text, so it is suitable
for easy parsing in scripting. The values are returned in atomic
units.

These are the available subcommands:

- Shell:  `results real-time`

  When no arguments are given, print the values calculated.

  Example: `inq results real-time`.


- Shell:  `results real-time total-steps`

  Returns the total number of real-time simulation steps done.

  Example: `inq results real-time total-steps`.


- Shell:  `results real-time total-time`

  Returns the total simulated time (in atomic units).

  Example: `inq results real-time total-time`.

- Shell:  `results real-time time [step]`

  Returns the time values. If not additional arguments are passed, inq
  prints the whole series for each time step. Alternatively, you can
  pass a step index to get the energy value.

  Note that for the moment inq uses a uniform time integration, so the
  time is just the step index times the time-step.

  Examples: `inq results real-time time`.
            `inq results real-time time 99`.


- Shell:  `results real-time total-energy [step]`

  Returns the values of the total energy during the propagation. If
  not additional arguments are passed, inq prints the whole series for
  each time step in a two column format, the first column is the time
  and the second one is the energy. This output is suitable to view
  on a plotting program like gnuplot.

  Alternatively, you can pass a step index to get the energy value for
  that step.

  Examples: `inq results real-time total-energy`
            `inq results real-time total-energy 43`.


)"""";
	}

private:

	static auto load() {
		try { return real_time::results::load(".inq/default_results_real_time"); }
		catch(...){
			actions::error(input::environment::global().comm(), "Cannot find real-time results, run a real-time simulation first");
			exit(1);
		}
	}
	
public:
	
	static void status() {
		auto res = load();
		if(input::environment::global().comm().root()) std::cout << res;
	}
	
	void operator()() const {
		status();
	}

	static auto total_steps() {
		return load().total_steps;
	}
	
	static auto total_time() {
		return load().total_time;
	}

	static auto time() {
		return load().time;
	}
	
	static auto total_energy() {
		return load().total_energy;
	}

	static auto dipole() {
		auto && res = load();
		if(res.dipole.size() == 0)	actions::error(input::environment::global().comm(), "The dipole was not calculated during the real-time simulation");
		return res.dipole;
	}
	
	static auto current() {
		auto && res = load();
		if(res.current.size() == 0)	actions::error(input::environment::global().comm(), "The current was not calculated during the real-time simulation.");
		return res.current;
	}
	
private:

	template <typename ArgsType, typename ArrayType> 
	void array_output_scalar(ArgsType args, ArrayType const & array, std::string const & label, std::string const & command_name) const {
		
		if(args.size() == 1) {
			
			auto time_array = time();
			if(input::environment::global().comm().root()) {
				printf("%-30s\t%-30s\n", "#time [atu]", label.c_str());
				for(auto ii = 0ul; ii < time_array.size(); ii++) printf("%-30.20e\t%-30.20e\n", time_array[ii], array[ii]);
			}
			actions::normal_exit();
			
		} else if (args.size() == 2) {
			auto index = utils::str_to<long>(args[1]);
			if(index < 0 or index >= (long) array.size()) actions::error(input::environment::global().comm(), "Invalid index ", index, " in the '" + command_name + "' command");
			if(input::environment::global().comm().root()) printf("%-30.20e\n", array[index]);
			actions::normal_exit();
		} else {
			actions::error(input::environment::global().comm(), "Invalid syntax in the '" + command_name + "' command");
		}
	}
	
	template <typename ArgsType, typename ArrayType> 
	void array_output_vector(ArgsType args, ArrayType const & array, std::string const & label, std::string const & command_name) const {
		
		if(args.size() == 1) {
			
			auto time_array = time();
			if(input::environment::global().comm().root()) {
				printf("%-30s\t%-30s\t%-30s\t%-30s\t\n", "#time [atu]", ("x-" + label).c_str(), ("y-" + label).c_str(), ("z-" + label).c_str());
				for(auto ii = 0ul; ii < time_array.size(); ii++) printf("%-30.20e\t%-30.20e\t%-30.20e\t%-30.20e\n", time_array[ii], array[ii][0], array[ii][1], array[ii][2]);
			}
			actions::normal_exit();
			
		} else if (args.size() == 2 or args.size() == 3) {
			auto index = utils::str_to<long>(args[1]);
			if(index < 0 or index >= (long) array.size()) actions::error(input::environment::global().comm(), "Invalid index ", index, " in the '" + command_name + "' command");
			
			if(args.size() == 2) {
				if(input::environment::global().comm().root()) printf("%.20e\t%.20e\t%.20e\n", array[index][0], array[index][1], array[index][2]);
				actions::normal_exit();
			}
			
			auto idir = utils::str_to_index(args[2]);
			if(idir == -1) actions::error(input::environment::global().comm(), "Invalid coordinate index in the '" + command_name + "' command");
			if(input::environment::global().comm().root()) printf("%.20e\n", array[index][idir]);
			actions::normal_exit();
			
		} else {
			actions::error(input::environment::global().comm(), "Invalid syntax in the '" + command_name + "' command");
		}
	}
	
public:
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0){
			operator()();
			actions::normal_exit();
		}

		if(args.size() == 1 and args[0] == "total-steps"){
			if(input::environment::global().comm().root()) printf("%ld\n", total_steps());
			actions::normal_exit();
		}
		
		if(args.size() == 1 and args[0] == "total-time"){
			if(input::environment::global().comm().root()) printf("%.20e\n", total_time());
			actions::normal_exit();
		}

		if(args[0] == "time"){
			auto time_array = time();
			if(args.size() == 1) {
				if(input::environment::global().comm().root()) {
					for(auto & ti : time_array) printf("%.20e\n", ti);
				}
			} else if (args.size() == 2) {
				if(input::environment::global().comm().root()) printf("%.20e\n", time_array[utils::str_to<long>(args[1])]);
			} else {
				actions::error(input::environment::global().comm(), "Invalid syntax in the 'results real-time time' command");
			}
			actions::normal_exit();
		}

		if(args[0] == "total-energy")   array_output_scalar(args, total_energy(),    "total-energy [Ha]",    "result real-time energy");
		if(args[0] == "dipole")         array_output_vector(args, dipole(),          "dipole [au]",          "result real-time dipole");
		if(args[0] == "current")        array_output_vector(args, current(),         "current [au]",         "result real-time current");
		
		actions::error(input::environment::global().comm(), "Invalid syntax in the 'results real-time' command");
	}
		
#ifdef INQ_PYTHON_INTERFACE
	template <class PythonModule>
	void python_interface(PythonModule & module) const {
		namespace py = pybind11;
		using namespace pybind11::literals;
 
		auto sub = module.def_submodule("real_time", help());

		sub.def("status",       &status);
		sub.def("total_steps",  &total_steps);
		sub.def("total_time",   &total_time);
		sub.def("time",         &time);
		sub.def("total_energy", &total_energy);

		sub.def("dipole", []() {
			
			auto dipole_multi = dipole();
			
			py::array_t<double, py::array::c_style> dipole_array({(long) dipole_multi.size(), 3l});
			
			auto arr = dipole_array.mutable_unchecked();
			
			for (py::ssize_t iter = 0; iter < arr.shape(0); iter++) {
				for (py::ssize_t idir = 0; idir < arr.shape(1); idir++) {
					arr(iter, idir) = dipole_multi[iter][idir];
				}
			}
		
			return dipole_array;
		});

		sub.def("current", []() {
			
			auto current_multi = current();
			
			py::array_t<double, py::array::c_style> current_array({(long) current_multi.size(), 3l});
			
			auto arr = current_array.mutable_unchecked();
			
			for (py::ssize_t iter = 0; iter < arr.shape(0); iter++) {
				for (py::ssize_t idir = 0; idir < arr.shape(1); idir++) {
					arr(iter, idir) = current_multi[iter][idir];
				}
			}
		
			return current_array;
		});

	}
#endif
	
} const results_real_time;

}
}
#endif

#ifdef INQ_INTERFACE_RESULTS_REAL_TIME_UNIT_TEST
#undef INQ_INTERFACE_RESULTS_REAL_TIME_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
