/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__PARSE_XYZ
#define INQ__INPUT__PARSE_XYZ

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <vector>
#include <cmath>

#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <input/atom.hpp>
#include <input/species.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace input {

auto parse_xyz(const std::string & xyz_file_name, quantity<magnitude::length> unit = magnitude::operator""_angstrom(1.0)){

	using namespace inq;
 
	std::vector<input::atom> geo;

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
		geo.push_back(atom_name | atom_position*unit.in_atomic_units());
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

  CHECK(geo[2].species() == pseudo::element("C"));
  CHECK(geo[2].position()[0] == 2.2846788549_a);
  CHECK(geo[2].position()[1] == -1.3190288178_a);
  CHECK(geo[2].position()[2] == 0.0_a);

  CHECK(geo[11].species() == pseudo::element("H"));
  CHECK(geo[11].position()[0] == -4.0572419367_a);
  CHECK(geo[11].position()[1] == 2.343260364_a);
  CHECK(geo[11].position()[2] == 0.0_a);

  geo.push_back("Cl" | vector3<double>(-3.0, 4.0, 5.0));

  CHECK(geo.size() == 13);
  CHECK(geo[12].species() == pseudo::element("Cl"));
  CHECK(geo[12].position()[0] == -3.0_a);
  CHECK(geo[12].position()[1] == 4.0_a);
  CHECK(geo[12].position()[2] == 5.0_a);

}
#endif
