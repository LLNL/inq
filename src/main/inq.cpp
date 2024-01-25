/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/interface.hpp>
#include <inq/inq.hpp>

int main(int argc, char* argv[]) {

	using namespace inq;
	
	input::environment::global(); //Initialize MPI 

	if(argc == 1){
		std::cout << "Usage: inq [command] [arguments]\n\n";
		std::cout << "The following commands are available:\n";
		std::cout << "  clear      Removes any inq information from the current directory.\n";
		std::cout << std::endl;
		exit(1);
	}
	
	if(argv[1] == std::string("clear")) interface::clear();
	
	fftw_cleanup();
  return 0;
}
