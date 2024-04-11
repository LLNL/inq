/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__UTIL
#define INQ__INTERFACE__UTIL

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace interface {

struct {

	constexpr auto name() const {
		return "util";
	}

	constexpr auto one_line() const {
		return "Miscelaneous utility commands";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'util' command
==================

This command provided some simple utilities for use in scripts or for
the testing of inq.

These are the available subcommands:

- `util match calc <math-expression>`

  Calculates the result of the given mathematical expression. This is
  useful if you need to calculate a value in a script or double check
  an expression you are given to inq in an argument.

  Example: `inq util calc "1/sqrt(2.0) + exp(-1/2)`.


- `util match <value1> <value2> ... <tolerance>`

  Checks if two set of values match within a certain tolerance. If they match,
  the command will run successfully, with a 0 return code. If the match
  fails, then the command will fail with a return value of 1.

  When more than two values are passed, it is assumed these are two
  sets given in sequential order. For example, if 6 values are passed
  inq assumes than the first 3 are one set and the second 3 another
  one. Then comparisons are made between the first value with the
  fourth one, the second with the fifth, and the third with the sixth
  one.

  Examples: `inq util match 1.0 2.0 1e-5`
            `inq util match  1.0 0.0 0.0  1.001 1e-12 1e-23  1e-2`

- `util test-data`

  Returns the path where inq install the data files used for
  tests. This is not really useful for users, just developers.

  Example: `inq util test-data`


)"""";
	}

	auto calc(std::string const & expr) const {
		return utils::str_to<double>(expr);		
	}

	bool match(double value, double reference, double tolerance) const {
    
    auto diff = fabs(reference - value);
    
    if(diff > tolerance){
			if(input::environment::global().comm().root()) {
				tfm::format(std::cout, "\nMatch: FAILED\n");
				tfm::format(std::cout, "  calculated value = %.12f\n",  value);
				tfm::format(std::cout, "  reference value  = %.12f\n",  reference);
				tfm::format(std::cout, "  difference       = %.1e\n",   diff);
				tfm::format(std::cout, "  tolerance        = %.1e\n\n", tolerance);
			}
      return false;
    } else {
      if(input::environment::global().comm().root()) tfm::format(std::cout, "Match: SUCCESS (value = %.12f , diff = %.1e)\n", value, diff);
      return true;
    }
  }

	auto test_data() const {
		return config::path::unit_tests_data();
	}

	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

		using utils::str_to;

		if(args.size() == 2 and args[0] == "calc"){
			if(input::environment::global().comm().root()) printf("%.20e\n", calc(args[1]));
			actions::normal_exit();
		}

		if(args[0] == "match"){
			if(args.size()%2 == 1) actions::error(input::environment::global().comm(), "Invalid number of arguments in the 'util match' command");

			auto tol = str_to<double>(args[args.size() - 1]);
			auto ncomp = args.size()/2 - 1;

			auto no_error = true;
			for(auto icomp = 0ul; icomp < ncomp; icomp++){
				auto val = str_to<double>(args[icomp + 1]);
				auto ref = str_to<double>(args[icomp + ncomp + 1]);
				no_error = match(val, ref, tol) and no_error;
			}
			
			if(no_error) {
				actions::normal_exit();
			} else {
				actions::error_exit();
			}
		}
			
		if(args.size() == 1 and args[0] == "test-data"){
			if(input::environment::global().comm().root()) std::cout << test_data() << std::endl;
			actions::normal_exit();
		}

		actions::error(input::environment::global().comm(), "Invalid syntax in the 'util' command");
	}
	
} const util ;

}
}
#endif

#ifdef INQ_INTERFACE_UTIL_UNIT_TEST
#undef INQ_INTERFACE_UTIL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
