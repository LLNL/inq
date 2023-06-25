/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__PARSE_XYZ
#define INQ__INPUT__PARSE_XYZ

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <vector>
#include <cmath>

#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <systems/ions.hpp>
#include <input/species.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace input {

auto parse_xyz(const std::string & xyz_file_name, quantity<magnitude::length> unit = magnitude::operator""_angstrom(1.0)){

	using namespace inq;
 
	std::vector<systems::ions::atom> geo;

	std::ifstream xyz_file(xyz_file_name.c_str());
	
	assert(xyz_file.is_open());

	int natoms;
	std::string comment_line;
  
	xyz_file >> natoms;
  
	std::getline(xyz_file, comment_line);
	std::getline(xyz_file, comment_line);
  
	std::string atom_name;
	vector3<double> atom_position;
  
	for(int iatom = 0; iatom < natoms; iatom++){
		xyz_file >> atom_name >> atom_position;
		geo.push_back(systems::ions::atom{input::species{atom_name}, atom_position*unit.in_atomic_units()});
	}
  
	xyz_file.close();
  
	assert(unsigned(natoms) == geo.size());

	return geo;
}

}
}
#endif

#ifdef INQ_INPUT_PARSE_XYZ_UNIT_TEST
#undef INQ_INPUT_PARSE_XYZ_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <config/path.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  auto geo = input::parse_xyz(config::path::unit_tests_data() + "benzene.xyz");
  
  CHECK(geo.size() == 12);

  CHECK(geo[2].species() == "C");
  CHECK(geo[2].position()[0] == 2.2846788549_a);
  CHECK(geo[2].position()[1] == -1.3190288178_a);
  CHECK(geo[2].position()[2] == 0.0_a);

  CHECK(geo[11].species() == "H");
  CHECK(geo[11].position()[0] == -4.0572419367_a);
  CHECK(geo[11].position()[1] == 2.343260364_a);
  CHECK(geo[11].position()[2] == 0.0_a);

}
#endif
