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

	std::string name() const {
		return "util";
	}

	std::string one_line() const {
		return "Miscelaneous utility commands";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'util' command
==================

This command provided some simple utilities used for the testing of
inq. They are not very useful for users.

These are the available subcommands:

- `util match <value1> <value2> <tolerance>`

  Checks if two values match within a certain tolerance. If they match,
  the command will run successfully, with a 0 return code. If the match
  fails, then the command will fail with a return value of 1.

  Example: `inq util match 1.0 2.0 1e-5`.


- `util test-data`

  Returns the path where inq install the data files used for
  tests. This is not really useful for users, just developers.

  Example: `inq util test-data`


)"""";
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
		
		if(args.size() == 4 and args[0] == "match"){
      auto val = atof(args[1].c_str());
      auto ref = atof(args[2].c_str());
      auto tol = atof(args[3].c_str());

      if(match(val, ref, tol)){
        exit(0);
      } else {
        exit(1);
      }
		}

		if(args.size() == 1 and args[0] == "test-data"){
			std::cout << test_data() << std::endl;
			exit(0);
		}
		
		if(input::environment::global().comm().root()) std::cerr << "Error: Invalid syntax in the 'util' command" << std::endl;
		exit(1);
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
