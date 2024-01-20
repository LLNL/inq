/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__INTERFACE
#define INPUT__INTERFACE

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace interface {

void cell_cubic(quantity<magnitude::length> const aa, int periodicity = 3){
  systems::ions ions(systems::cell::cubic(aa).periodicity(periodicity));
  ions.save(input::environment::global().comm(), ".default_ions");
}

void cell() {
  auto cell = systems::ions::load(".default_ions").cell();
  if(input::environment::global().comm().root()) std::cout << cell;
}

void ions_add(input::species const & sp, vector3<quantity<magnitude::length>> const & pos){
  auto ions = systems::ions::load(".default_ions");
  ions.insert(sp, pos);
  ions.save(input::environment::global().comm(), ".default_ions");
}

void ions(){
  auto ions = systems::ions::load(".default_ions");		
  if(input::environment::global().comm().root()) std::cout << ions;
}


}
}
#endif

#ifdef INQ_INPUT_INTERFACE_UNIT_TEST
#undef INQ_INPUT_INTERFACE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
