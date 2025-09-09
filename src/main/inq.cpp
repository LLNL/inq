/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>
#include <interface/aggregate.hpp>
#include <interface/aliases.hpp>
#include <interface/runtime_options.hpp>

using namespace inq;

int main(int argc, char* argv[]) {
	using interface::operator+;
 
	auto comm = input::environment::global().comm(); //Initialize MPI 

	auto all_commands =
		interface::item(interface::cell)
		+ interface::item(interface::clear)
		+ interface::item(interface::electrons)
		+ interface::item(interface::ground_state)
		+ interface::item(interface::history)
		+ interface::item(interface::ions)
		+ interface::item(interface::kpoints)
		+ interface::item(interface::perturbations)
		+ interface::item(interface::results_ground_state)
		+ interface::item(interface::results_real_time)
		+ interface::item(interface::real_time)
		+ interface::item(interface::run)
		+ interface::item(interface::species)
		+ interface::item(interface::spectrum)
		+ interface::item(interface::status)
		+ interface::item(interface::theory)
		+ interface::item(interface::util);

	auto all_helpers =
		interface::item(interface::units)
		+ interface::item(interface::results);
		
	interface::history_file.add_entry(argc, argv);

	interface::runtime_options run_opts;
	
	auto uniformize = [](auto arg){
		arg = utils::lowercase(arg);
		std::replace(arg.begin(), arg.end(), '_', '-'); //replace underscores with dashes
		return arg;
	};

	std::vector<std::string> args;
	for(int iarg = 1; iarg < argc; iarg++) {
		auto arg = std::string(argv[iarg]);

		if(arg == "-h" or arg == "--help") {
			args = {"help"};
			break;
		}
		
		if(arg == "-q" or arg == "--quiet") {
			run_opts.quiet = true;
			continue;
		}

		if(arg == "-d" or arg == "--debug") {
			run_opts.debug = true;
			continue;
		}

		//if it's a filename, don't do anything to it
		if(args.size() > 0 and args.back() == "file") {
			args.emplace_back(arg);
			continue;
		}

		arg = uniformize(arg);

		auto aliases = interface::aliases();
		
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

	if(run_opts.debug) {
		std::cout << "Processed arguments: ";
		for(auto const & arg : args){
			std::cout << "|" << arg;
		}
		std::cout << "|" << std::endl;
	}

	if(args.size() == 0){
		args = {"status"};
	}

	auto command = args[0];
	args.erase(args.begin());

	all_commands.execute(command, args, run_opts);
	
	if(command == "help") {
		if(args.size() == 0){
			if(comm.root()) {
				std::cout << "\n";
				std::cout << "Usage: inq <command> [arguments]\n\n";
				std::cout << "The following commands are available:\n";
				std::cout << interface::list_item("help", "Prints detailed information about other commands");
				std::cout << all_commands.list();
				std::cout << "\n";
				std::cout << "And the following options:\n";
				std::cout << interface::list_item("-h,--help",  "Prints this help dialog");
				std::cout << interface::list_item("-q,--quiet", "Run silently, do not print information unless explicitly asked to");
				std::cout << interface::list_item("-d,--debug", "Print debug information (useful for inq developers)");
				std::cout << "\nTo get more information about any command use: inq help <command>\n";
				std::cout << "\nBesides commands, there is also some additional help topics you can read with 'help':\n\n";
				std::cout << all_helpers.list();
				std::cout << std::endl;
			}
			interface::actions::normal_exit();
		}

		command = args[0];
		args.erase(args.begin());

		if(comm.root()) {
			all_commands.help(command);
			all_helpers.help(command);

			interface::actions::error(comm, "Unknown help item '" + command + "'.");
		} else {
			interface::actions::normal_exit();
		}
	}
	
	interface::actions::error(comm, "Unknown command '" + command + "'.");
	
	fftw_cleanup();
  return 0;
}
