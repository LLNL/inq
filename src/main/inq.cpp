/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>
#include <utils/lowercase.hpp>

int main(int argc, char* argv[]) {
	using namespace std::string_literals;
 
	std::map<std::string, std::string> dictionary = {
    { "ground_state"s,     "ground-state"s     },
    {	"groundstate"s,      "ground-state"s     },
    {	"hartree_fock"s,     "hartree-fock"s     },
		{	"hartreefock"s,      "hartree-fock"s     },
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
		{ "tol"s         ,     "tolerance"s        }
	};
	
	using namespace inq;
	
	input::environment::global(); //Initialize MPI 

	if(argc == 1){
		std::cout << "Usage: inq <command> [arguments]\n\n";
		std::cout << "The following commands are available:\n";
		std::cout << "  " << interface::clear.name()        << "\t\t" << interface::clear       .one_line() << '\n';
		std::cout << "  " << interface::cell.name()         << "\t\t" << interface::cell        .one_line() << '\n';
		std::cout << "  " << interface::ions.name()         << "\t\t" << interface::ions        .one_line() << '\n';
		std::cout << "  " << interface::theory.name()       << "\t\t" << interface::theory      .one_line() << '\n';
		std::cout << "  " << interface::electrons.name()    << "\t\t" << interface::electrons   .one_line() << '\n';
		std::cout << "  " << interface::ground_state.name() << "\t\t" << interface::ground_state.one_line() << '\n';
		std::cout << "  " << interface::run.name()          << "\t\t" << interface::run         .one_line() << '\n';
		std::cout << "  " << interface::energy.name()       << "\t\t" << interface::energy      .one_line() << '\n';
		std::cout << "  " << interface::util.name()         << "\t\t" << interface::util        .one_line() << '\n';
		std::cout << "\n";
		std::cout << "And the following options:\n";
		std::cout << "  -q,--quiet    Run silently, do not print information unless explicitly asked to.";
		std::cout << std::endl;
		exit(1);
	}

	auto quiet = false;

	std::vector<std::string> args;
	for(int iarg = 1; iarg < argc; iarg++) {
		auto arg = std::string(argv[iarg]);

		if(arg == "-q" or arg == "--quiet") {
			quiet = true;
			continue;
		}

		arg = utils::lowercase(arg);

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
	
	if(command == interface::clear       .name()) interface::clear       .command(args, quiet);
	if(command == interface::cell        .name()) interface::cell        .command(args, quiet);
	if(command == interface::ions        .name()) interface::ions        .command(args, quiet);
	if(command == interface::theory      .name()) interface::theory      .command(args, quiet);
	if(command == interface::electrons   .name()) interface::electrons   .command(args, quiet);
	if(command == interface::ground_state.name()) interface::ground_state.command(args, quiet);
	if(command == interface::run         .name()) interface::run         .command(args, quiet);
	if(command == interface::energy      .name()) interface::energy      .command(args, quiet);
	if(command == interface::util        .name()) interface::util        .command(args, quiet);	
	
	std::cerr << "inq error: unknown command '" << command << "'." << std::endl;
	exit(1);
	
	fftw_cleanup();
  return 0;
}
