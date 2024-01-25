/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/interface.hpp>
#include <inq/inq.hpp>

int main(int argc, char* argv[]) {

	inq::input::environment::global(); //Initialize MPI 

	if(argc == 1){
		std::cout << "usage: " << argv[0] << " [command] [arguments]" << std::endl;
		exit(1);
	}
	
	fftw_cleanup();
  return 0;
}
