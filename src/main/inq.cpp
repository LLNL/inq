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
		std::cout << "  " << interface::clear.name()        << "\t\t" << interface::clear       .one_line() << '\n';
		std::cout << "  " << interface::cell.name()         << "\t\t" << interface::cell        .one_line() << '\n';
		std::cout << "  " << interface::ions.name()         << "\t\t" << interface::ions        .one_line() << '\n';
		std::cout << "  " << interface::theory.name()       << "\t\t" << interface::theory      .one_line() << '\n';
		std::cout << "  " << interface::electrons.name()    << "\t\t" << interface::electrons   .one_line() << '\n';
		std::cout << "  " << interface::ground_state.name() << "\t\t" << interface::ground_state.one_line() << '\n';
		std::cout << "  " << interface::run.name()          << "\t\t" << interface::run         .one_line() << '\n';
		std::cout << std::endl;
		exit(1);
	}

	auto quiet = false;

	std::vector<std::string> args;
	for(int iarg = 1; iarg < argc; iarg++) {
		auto arg = std::string(argv[iarg]);
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
	
	std::cerr << "inq error: unknown command '" << command << "'." << std::endl;
	exit(1);
	
	fftw_cleanup();
  return 0;
}
