/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <input/environment.hpp>
#include <fftw3.h>

int main( int argc, char* argv[] ) {

	inq::input::environment::global(); //Initialize MPI 
	
  int result = Catch::Session().run( argc, argv );

	fftw_cleanup();
	
  return result;
}
