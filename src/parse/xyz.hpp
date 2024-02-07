/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARSE__XYZ
#define INQ__PARSE__XYZ

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <vector>
#include <cmath>

#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace parse {

class xyz {

	std::vector<input::species> atoms_;
	std::vector<vector3<double>> positions_;

	public:
	
	xyz(const std::string & xyz_file_name, quantity<magnitude::length> unit = magnitude::operator""_angstrom(1.0)){
		
		std::ifstream xyz_file(xyz_file_name);
		
		if(not xyz_file.is_open()) throw std::runtime_error("error: failed to open XYZ file '" + xyz_file_name + "'.");
		
		int natoms;
		std::string comment_line;
		
		xyz_file >> natoms;
		
		std::getline(xyz_file, comment_line);
		std::getline(xyz_file, comment_line);
		
		std::string atom_name;
		vector3<double> atom_position;
		
		for(int iatom = 0; iatom < natoms; iatom++){
			xyz_file >> atom_name >> atom_position;
			atoms_.emplace_back(atom_name);
			positions_.push_back(atom_position*unit.in_atomic_units());
		}
		
		xyz_file.close();
		
		assert(unsigned(natoms) == atoms_.size());
		assert(unsigned(natoms) == positions_.size());
	}

	auto size() const {
		return long(atoms_.size());
	}
	
	auto & atoms() const {
		return atoms_;
	}

	auto & positions() const {
		return positions_;
	}
	
};

}
}
#endif

#ifdef INQ_PARSE_XYZ_UNIT_TEST
#undef INQ_PARSE_XYZ_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <config/path.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  auto xyz_file = parse::xyz(config::path::unit_tests_data() + "benzene.xyz");
  
  CHECK(xyz_file.size() == 12);

  CHECK(xyz_file.atoms()[2] == "C");
  CHECK(xyz_file.positions()[2][0] == 2.2846788549_a);
  CHECK(xyz_file.positions()[2][1] == -1.3190288178_a);
  CHECK(xyz_file.positions()[2][2] == 0.0_a);

  CHECK(xyz_file.atoms()[11] == "H");
  CHECK(xyz_file.positions()[11][0] == -4.0572419367_a);
  CHECK(xyz_file.positions()[11][1] == 2.343260364_a);
  CHECK(xyz_file.positions()[11][2] == 0.0_a);

}
#endif
