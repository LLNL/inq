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
 
	std::unordered_map<std::string, std::string> aliases = {
		{ "calculator"s,               "calc"s                     },
		{ "functionals"s,              "functional"s               },
		{	"groundstate"s,              "ground-state"s             },
		{	"hartreefock"s,              "hartree-fock"s             },
		{	"kpoint"s,                   "kpoints"s                  },
		{	"kpoints"s,                  "kpoints"s                  },
		{	"k-point"s,                  "kpoints"s                  },
		{	"k-points"s,                 "kpoints"s                  },
		{ "exactexchange"s,            "exact-exchange"s           },
		{ "extraelectrons"s,           "extra-electrons"s          },
		{ "extrastates"s,              "extra-states"s             },
		{ "filename"s,                 "file"s                     },
		{ "freq"s,                     "frequency"s                },
		{ "maxsteps"s,                 "max-steps"s                },
		{ "mix"s,                      "mixing"s                   },
		{ "noncollinear"s,             "non-collinear"s            },
		{ "noninteracting"s,           "non-interacting"s          },
		{ "nonlocal"s,                 "non-local"s                },
		{ "perturbation"s,             "perturbations"s            },
		{ "realtime"s,                 "real-time"s                },
		{ "result"s,                   "results"s                  },
		{ "results-ground-state"s,     "results ground-state"s     },
		{ "resultsgroundstate"s,       "results ground-state"s     },
		{ "resultsground-state"s,      "results ground-state"s     },
		{ "results-groundstate"s,      "results ground-state"s     },
		{ "resultsgs"s,                "results ground-state"s     },
		{ "results-gs"s,               "results ground-state"s     },
		{ "results-real-time"s,        "results real-time"s        },
		{ "resultsrealtime"s,          "results real-time"s        },
		{ "resultsreal-time"s,         "results real-time"s        },
		{ "results-realtime"s,         "results real-time"s        },
		{ "resultsrt"s,                "results real-time"s        },
		{ "results-rt"s,               "results real-time"s        },
		{ "grid-shifted"s,             "shifted-grid"s             },
		{ "shiftedgrid"s,              "shifted-grid"s             },
		{ "timestep"s,                 "time-step"s                },
		{ "totaltime"s,                "total-time"s               },
		{ "totalsteps"s,               "total-steps"s              },
		{ "tol"s,                      "tolerance"s                },
		{ "utils"s,                    "util"s                     }
	};
	
	auto comm = input::environment::global().comm(); //Initialize MPI 

	auto all_commands =
		interface::item(interface::cell)
		+ interface::item(interface::clear)
		+ interface::item(interface::electrons)
		+ interface::item(interface::ground_state)
		+ interface::item(interface::ions)
		+ interface::item(interface::kpoints)
		+ interface::item(interface::perturbations)		
		+ interface::item(interface::results_ground_state)
		+ interface::item(interface::results_real_time)
		+ interface::item(interface::real_time)
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
			std::cout << interface::list_item("help", "Prints detailed information about other commands");
			std::cout << all_commands.list();
			std::cout << "\n";
			std::cout << "And the following options:\n";
			std::cout << interface::list_item("-q,--quiet", "Run silently, do not print information unless explicitly asked to");
			std::cout << interface::list_item("-d,--debug", "Print debug information (useful for inq developers)");
			std::cout << std::endl;
		}
		exit(0);
	}

	auto quiet = false;
	auto debug = false;
	
	auto uniformize = [](auto arg){
		arg = utils::lowercase(arg);
		std::replace(arg.begin(), arg.end(), '_', '-'); //replace underscores with dashes
		return arg;
	};
	
	std::vector<std::string> args;
	for(int iarg = 1; iarg < argc; iarg++) {
		auto arg = std::string(argv[iarg]);

		if(arg == "-q" or arg == "--quiet") {
			quiet = true;
			continue;
		}

		if(arg == "-d" or arg == "--debug") {
			debug = true;
			continue;
		}

		//if it's a filename, don't do anything to it
		if(args.size() > 0 and args.back() == "file") {
			args.emplace_back(arg);
			continue;
		}

		arg = uniformize(arg);

		//process aliases
		auto search = aliases.find(arg);
		if(search != aliases.end()) arg = search->second;

		//process aliases for words with one or two spaces
		if(iarg + 1 < argc){
			auto fusion = uniformize(arg + argv[iarg + 1]);
			auto search = aliases.find(fusion);
			if(search != aliases.end()) {
				arg = search->second;
				iarg++;
			}
		}
		
		if(iarg + 2 < argc){
			auto fusion = uniformize(arg + argv[iarg + 1] + argv[iarg + 2]);
			auto search = aliases.find(fusion);
			if(search != aliases.end()) {
				arg = search->second;
				iarg += 2;
			}
		}
		
		args.emplace_back(arg);
	}

	if(debug) {
		std::cout << "Processed arguments: ";
		for(auto const & arg : args){
			std::cout << "|" << arg;
		}
		std::cout << "|" << std::endl;
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
