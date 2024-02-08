/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>
#include <interface/aggregate.hpp>

using namespace inq;

int main(int argc, char* argv[]) {
	using namespace std::string_literals;
	using interface::operator+;
 
	std::map<std::string, std::string> dictionary = {
    { "ground_state"s,     "ground-state"s     },
    {	"groundstate"s,      "ground-state"s     },
    {	"hartree_fock"s,     "hartree-fock"s     },
		{	"hartreefock"s,      "hartree-fock"s     },
		{	"kpoint"s,           "kpoints"s          },
		{	"kpoints"s,          "kpoints"s          },
		{	"k-point"s,          "kpoints"s          },
		{	"k-points"s,         "kpoints"s          },	
		{ "exact_exchange"s,   "exact-exchange"s   },
		{ "exactexchange"s,    "exact-exchange"s   },
		{ "extra_electrons"s,  "extra-electrons"s  },
		{ "extraelectrons"s,   "extra-electrons"s  },
		{ "extra_states"s,     "extra-states"s     },
		{ "extrastates"s,      "extra-states"s     },
		{ "max_steps"s,        "max-steps"s        },
		{ "maxsteps"s,         "max-steps"s        },
		{ "mix"s,              "mixing"s           },
		{ "non_collinear"s,    "non-collinear"s    },
		{ "noncollinear"s,     "non-collinear"s    },
		{ "non_interacting"s,  "non-interacting"s  },
		{ "noninteracting"s,   "non-interacting"s  },
		{ "non_local"s,        "non-local"s        },
		{ "nonlocal"s,         "non-local"s        },
		{ "gridshifted"s,      "shifted-grid"s     },
		{ "grid-shifted"s,     "shifted-grid"s     },
		{ "shiftedgrid"s,      "shifted-grid"s     },
		{ "tol"s         ,     "tolerance"s        }
	};
	
	auto comm = input::environment::global().comm(); //Initialize MPI 

	auto all_commands =
		interface::item(interface::cell)
		+ interface::item(interface::clear)
		+ interface::item(interface::electrons)
		+ interface::item(interface::energy)
		+ interface::item(interface::ground_state)
		+ interface::item(interface::ions)
		+ interface::item(interface::kpoints)
		+ interface::item(interface::run)
		+ interface::item(interface::theory)
		+ interface::item(interface::util);

	auto all_helpers =
		interface::item(interface::units);
		
	if(argc == 1){
		if(comm.root()) {
			std::cout << "\n";
			std::cout << "Usage: inq <command> [arguments]\n\n";
			std::cout << "The following commands are available:\n";
			std::cout << "  " << "help"                         << "\t\t" << "Prints detailed information about other commands\n";
			std::cout << all_commands.list();
			std::cout << "\n";
			std::cout << "And the following options:\n";
			std::cout << "  -q,--quiet    Run silently, do not print information unless explicitly asked to.\n";
			std::cout << std::endl;
		}
		exit(0);
	}

	auto quiet = false;

	std::vector<std::string> args;
	for(int iarg = 1; iarg < argc; iarg++) {
		auto arg = std::string(argv[iarg]);

		if(arg == "-q" or arg == "--quiet") {
			quiet = true;
			continue;
		}

		auto lower = true;
		if(args.size() > 0 and args.back() == "file") lower = false; //do not convert filenames to lowercase
		if(lower) arg = utils::lowercase(arg);

		//convert spelling
		auto search = dictionary.find(arg);
		if(search != dictionary.end()) arg = search->second;

		//convert spelling for words with a space
		if(iarg + 1 < argc){
			auto fusion = utils::lowercase(arg + argv[iarg + 1]);
			auto search = dictionary.find(fusion);
			if(search != dictionary.end()) {
				arg = search->second;
				iarg++;
			}
		}
		
		args.emplace_back(arg);
	}
	
	auto command = args[0];
	args.erase(args.begin());

	all_commands.execute(command, args, quiet);
	
	if(command == "help") {
		if(args.size() == 0){
			if(comm.root()) {
				std::cout << "\n";
				std::cout << "Usage: inq help <command>\n\n";
				std::cout << "The 'help' command prints detailed information about other inq commands:\n\n";
				std::cout << all_commands.list();
				std::cout << "\nThere is also some additional help topics you can read:\n\n";
				std::cout << all_helpers.list();
				std::cout << std::endl;
			}
			exit(0);
		}

		command = args[0];
		args.erase(args.begin());

		if(comm.root()) {
			all_commands.help(command);
			all_helpers.help(command);
			
			std::cerr << "inq error: unknown help item '" << command << "'." << std::endl;
			exit(1);
		} else {
			exit(0);
		}
	}
	
	if(comm.root()) std::cerr << "inq error: unknown command '" << command << "'." << std::endl;
	exit(1);
	
	fftw_cleanup();
  return 0;
}
