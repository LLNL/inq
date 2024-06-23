/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__HISTORY
#define INQ__INTERFACE__HISTORY

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

- `history clear`

  Clears the list of saved commands.

  Shell example: `inq history clear`

- `history script`

  Prints a bash script containing the list of commands that reproduce
  the calculation in the current directory. The recommended script
  headers are included. Any command issues before `inq clear` will not
  be included.

  Shell example: `inq history script > my_calculation.sh`


)"""";
	}

	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

    if(args.size() == 0) {
      std::ifstream hist_file(".inq_history");
      if(hist_file.is_open()) std::cout << hist_file.rdbuf();
      actions::normal_exit();
    }                   

		if(args.size() == 1 and args[0] == "clear") {
      std::ofstream hist_file(".inq_history");
			actions::normal_exit();
		}

    if(args.size() == 1 and args[0] == "script") {

      std::cout << "#!/bin/bash\n\n"
                << "set -e #make the script fail if a command fails\n"
                << "set -x #output commands to the terminal\n\n"
                << "inq clear\n";

      std::ifstream hist_file(".inq_history");

      std::list<std::string> file_lines;
      while(hist_file) {
        std::string line;
        std::getline(hist_file, line);

        if(utils::lowercase(line).rfind("inq history", 0) == 0) continue;
				
        file_lines.push_back(line);
      }

      for(auto & line : file_lines) std::cout << line << '\n';
      
      actions::normal_exit();
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
