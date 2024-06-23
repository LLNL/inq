/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__HISTORY
#define INQ__INTERFACE__HISTORY

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
		return "history";
	}

	constexpr auto one_line() const {
		return "Manage inq command history";
	}
	
	constexpr auto help() const {
		return R""""(

The 'history' command
==================

This command provides access to the history of the inq commands
executed in the current directory. This command is only available for
the shell interface.

These are the available subcommands:

- `history`

  Prints a full list of the inq commands that have been executed in the
  current directory.

  Shell example: `inq history`


)"""";
	}

	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

    if(args.size() == 0) {
      std::ifstream hist_file(".inq_history");
      if(hist_file.is_open()) std::cout << hist_file.rdbuf();
    }                   

		actions::error(input::environment::global().comm(), "Invalid syntax in the 'history' command");
	}
	
} const history ;

}
}
#endif

#ifdef INQ_INTERFACE_HISTORY_UNIT_TEST
#undef INQ_INTERFACE_HISTORY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
