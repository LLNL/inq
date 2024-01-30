/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char* argv[]) {

	using namespace inq;
	
	input::environment::global(); //Initialize MPI 

	if(argc == 1){
		std::cout << "Usage: inq <command> [arguments]\n\n";
		std::cout << "The following commands are available:\n";
		std::cout << "  " << interface::clear.name() << "\t\t" << interface::clear.one_line() << '\n';
		std::cout << "  " << interface::cell.name()  << "\t\t" << interface::cell.one_line()  << '\n';		
		std::cout << "  " << interface::ions.name()  << "\t\t" << interface::ions.one_line()  << '\n';
		std::cout << std::endl;
		exit(1);
	}

	auto quiet = false;

	auto command = std::string(argv[1]);

	std::vector<std::string> args;
	for(int iarg = 2; iarg < argc; iarg++) args.emplace_back(argv[iarg]);
	
	if(command == interface::clear.name()) interface::clear.command(args, quiet);
	if(command == interface::cell.name())  interface::cell.command (args, quiet);
	if(command == interface::ions.name()) {

		if(args.size() == 0) {
			interface::ions();
			exit(0);
		}

		if(args[0] == "clear"){

			if(args.size() != 1) {
				std::cerr << "The 'ions clear' command doesn't take arguments." << std::endl;
				exit(1);
			}
			interface::ions.clear();
			if(not quiet) interface::ions();
			exit(0);
		}
 
		if(args[0] == "add"){

			if(args.size() != 6) {
				std::cerr << "Wrong arguments for ions add.\nUse: inq ions add <symbol> <pos_x> <pos_y> <pos_z> <units>" << std::endl;
				exit(1);
			}

			auto symbol = args[1];
			auto units = magnitude::length::parse(args[5]);
			auto xx = atof(args[2].c_str())*units;
			auto yy = atof(args[3].c_str())*units;
			auto zz = atof(args[4].c_str())*units;
			
			interface::ions.add(symbol, {xx, yy, zz});
			if(not quiet) interface::ions();
			exit(0);
		}
 
		std::cerr << "Invalid syntax in the ions command" << std::endl;
		exit(1);
	}
	
	std::cerr << "inq error: unknown command '" << command << "'." << std::endl;
	exit(1);
	
	fftw_cleanup();
  return 0;
}
