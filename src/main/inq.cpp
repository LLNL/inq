/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/interface.hpp>
#include <inq/inq.hpp>

static auto parse_periodicity(std::string per_string){
	
	std::transform(per_string.begin(), per_string.end(), per_string.begin(), ::tolower);
	
	if(per_string == "finite" or per_string == "0" or per_string == "0d") {
		return 0;
	} else if (per_string == "wire" or per_string == "1" or per_string == "1d"){
		return 1;
	} else if (per_string == "slab" or per_string == "2" or per_string == "2d"){
		return 2;
	} else if (per_string == "periodic" or per_string == "3" or per_string == "3d"){
		return 3;
	} else {
		throw std::runtime_error("inq error: unknown periodicity '" + per_string + "'.");
	}
}

int main(int argc, char* argv[]) {

	using namespace inq;
	
	input::environment::global(); //Initialize MPI 

	if(argc == 1){
		std::cout << "Usage: inq <command> [arguments]\n\n";
		std::cout << "The following commands are available:\n";
		std::cout << "  clear      Removes any inq information from the current directory.\n";
		std::cout << "  cell       Defines the simulation cell.\n";
		std::cout << "  ions       Defines the ions in the simulation.\n";
		std::cout << std::endl;
		exit(1);
	}

	auto quiet = false;

	auto command = std::string(argv[1]);
	
	if(command == "clear") {
		interface::clear();
		exit(0);
	}

	if(command == "cell") {

		if(argc == 2) {
			interface::cell();
			exit(0);
		}

		if(argv[2] == std::string("cubic")){
			if(argc != 5 and argc != 6) {
				std::cerr << "Wrong arguments for a cubic cell definition.\nUse: inq cell cubic <lattice_parameter> <units> [periodicity]" << std::endl;
				exit(1);
			}

			auto aa = atof(argv[3])*magnitude::length::parse(argv[4]);

			int per = 3;
			if(argc == 6) per = parse_periodicity(argv[5]);

			interface::cell.cubic(aa, per);
			if(not quiet) interface::cell();
			exit(0);
		}

		std::cerr << "Invalid syntax in the cell command" << std::endl;
		exit(1);
	}

	if(command == "ions") {

		if(argc == 2) {
			interface::ions();
			exit(0);
		}

		if(argv[2] == std::string("clear")){

			if(argc != 3) {
				std::cerr << "The 'ions clear' command doesn't take arguments." << std::endl;
				exit(1);
			}
			interface::ions.clear();
			if(not quiet) interface::ions();
			exit(0);
		}
 
		if(argv[2] == std::string("add")){

			if(argc != 8) {
				std::cerr << "Wrong arguments for ions add.\nUse: inq ions add <symbol> <pos_x> <pos_y> <pos_z> <units>" << std::endl;
				exit(1);
			}

			auto symbol = std::string(argv[3]);
			auto units = magnitude::length::parse(argv[7]);
			auto xx = atof(argv[4])*units;
			auto yy = atof(argv[5])*units;
			auto zz = atof(argv[6])*units;
			
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
