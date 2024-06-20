/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PARSE__POSCAR
#define INQ__PARSE__POSCAR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>


#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace parse {

class poscar {

	std::vector<vector3<double>> lattice_vectors_;
	std::vector<input::species> atoms_;
	std::vector<vector3<double>> positions_;
	
public:
	
	poscar(const std::string & poscar_file_name):
		lattice_vectors_(3)
	{

		//This follows the specification from https://www.vasp.at/wiki/index.php/POSCAR
		
		using namespace inq::magnitude;

		std::ifstream poscar_file(poscar_file_name);

		if(not poscar_file.is_open()) throw std::runtime_error("error: failed to open POSCAR file '" + poscar_file_name + "'.");

		{
			std::string comment_line;
			std::getline(poscar_file, comment_line);
		}
		
		double scaling_factor;
		
		poscar_file >> scaling_factor;
		
		for(int idir = 0; idir < 3; idir++)	{
			poscar_file >> lattice_vectors_[idir];
			lattice_vectors_[idir] *= (scaling_factor*1.0_A).in_atomic_units();
		}

		std::vector<std::string> species;
		
		std::string species_line;
		std::getline(poscar_file, species_line);
		std::getline(poscar_file, species_line);

		std::istringstream iss(species_line + " ");

		std::string species_name;
		while(iss) {
			iss >> species_name;
			if(iss.eof()) break;
			species.push_back(species_name);
			input::species sp(species.back());
			if(not sp.valid()) throw std::runtime_error("Cannot read the species from POSCAR file \'" + poscar_file_name +
																									"\'. Make sure your file contains the optional \'Species names\' line (see\n https://www.vasp.at/wiki/index.php/POSCAR for details).");
		}
		
		std::vector<int> species_num(species.size());

		for(unsigned ispecies = 0; ispecies < species.size(); ispecies++) poscar_file >> species_num[ispecies];
		
		std::string tail;
		std::getline(poscar_file, tail);
		
		std::string coordinate_type;
		std::getline(poscar_file, coordinate_type);

		if(coordinate_type[0] == 's' or coordinate_type[0] == 'S') std::getline(poscar_file, coordinate_type); //skip the Optional selective dynamics line

		if(coordinate_type[0] == 'C' or coordinate_type[0] == 'c' or coordinate_type[0] == 'K' or coordinate_type[0] == 'k') {
			// Cartesian
			for(unsigned ispecies = 0; ispecies < species_num.size(); ispecies++){
				for(int iatom = 0; iatom < species_num[ispecies]; iatom++){
					vector3<double> pos;
					poscar_file >> pos;
					atoms_.emplace_back(input::species(species[ispecies]));
					positions_.emplace_back(scaling_factor*in_atomic_units(1.0_A*pos));
					std::getline(poscar_file, tail);
				}
			}

		} else {
			// Direct 
			auto cell = systems::cell(lattice_vectors_[0], lattice_vectors_[1], lattice_vectors_[2]);
			for(unsigned ispecies = 0; ispecies < species_num.size(); ispecies++){
				for(int iatom = 0; iatom < species_num[ispecies]; iatom++){
					vector3<double, contravariant> pos;
					poscar_file >> pos;
					atoms_.emplace_back(input::species(species[ispecies]));
					positions_.emplace_back(cell.metric().to_cartesian(pos));
					std::getline(poscar_file, tail);					
				}
			}
		}
		
	}
	
	auto size() const {
		return long(atoms_.size());
	}	
	
	auto & lattice() const {
		return lattice_vectors_;
	}
	
	auto & atoms() const {
		return atoms_;
	}

	auto & positions() const {
		return positions_;
	}

	auto cell() const {
		return systems::cell(lattice_vectors_[0], lattice_vectors_[1], lattice_vectors_[2]);
	}
	
};


}
}
#endif

#ifdef INQ_PARSE_POSCAR_UNIT_TEST
#undef INQ_PARSE_POSCAR_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <config/path.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("BN"){
	
		parse::poscar vasp_file(config::path::unit_tests_data() + "bn.poscar");

		CHECK(vasp_file.lattice()[0][0] == 0.0_a);
		CHECK(vasp_file.lattice()[0][1] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[0][2] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[1][0] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[1][1] == 0.0_a);
		CHECK(vasp_file.lattice()[1][2] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[2][0] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[2][1] == 3.3731611325_a);
		CHECK(vasp_file.lattice()[2][2] == 0.0_a);		
		
		CHECK(vasp_file.size() == 2);

		CHECK(vasp_file.atoms()[0] == "B");
		CHECK(vasp_file.atoms()[1] == "N");

		CHECK(vasp_file.positions()[0][0] == 0.0_a);
		CHECK(vasp_file.positions()[0][1] == 0.0_a);
		CHECK(vasp_file.positions()[0][2] == 0.0_a);
		CHECK(vasp_file.positions()[1][0] == 1.6865805662_a);
		CHECK(vasp_file.positions()[1][1] == 1.6865805662_a);
		CHECK(vasp_file.positions()[1][2] == 1.6865805662_a);		
	}

	SECTION("Al"){
		
		parse::poscar vasp_file(config::path::unit_tests_data() + "al.poscar");

		CHECK(vasp_file.lattice()[0][0] == 7.6458319003_a);
		CHECK(vasp_file.lattice()[0][1] == 0.0_a);
		CHECK(vasp_file.lattice()[0][2] == 0.0_a);
		CHECK(vasp_file.lattice()[1][0] == 0.0_a);
		CHECK(vasp_file.lattice()[1][1] == 7.6458319003_a);
		CHECK(vasp_file.lattice()[1][2] == 0.0_a);
		CHECK(vasp_file.lattice()[2][0] == 0.0_a);
		CHECK(vasp_file.lattice()[2][1] == 0.0_a);
		CHECK(vasp_file.lattice()[2][2] == 7.6458319003_a);		
		
		CHECK(vasp_file.size() == 4);

		CHECK(vasp_file.atoms()[0] == "Al");
		CHECK(vasp_file.atoms()[1] == "Al");
		CHECK(vasp_file.atoms()[2] == "Al");		
		CHECK(vasp_file.atoms()[3] == "Al");

		CHECK(vasp_file.positions()[0][0] == 0.0_a);
		CHECK(vasp_file.positions()[0][1] == 0.0_a);
		CHECK(vasp_file.positions()[0][2] == 0.0_a);
		CHECK(vasp_file.positions()[1][0] == 3.8229159501_a);
		CHECK(vasp_file.positions()[1][1] == 3.8229159501_a);
		CHECK(vasp_file.positions()[1][2] == 0.0_a);
		CHECK(vasp_file.positions()[2][0] == 0.0_a);
		CHECK(vasp_file.positions()[2][1] == 3.8229159501_a);
		CHECK(vasp_file.positions()[2][2] == 3.8229159501_a);
		CHECK(vasp_file.positions()[3][0] == 3.8229159501_a);
		CHECK(vasp_file.positions()[3][1] == 0.0_a);
		CHECK(vasp_file.positions()[3][2] == 3.8229159501_a);
	}

	SECTION("Ni"){
		
		parse::poscar vasp_file(config::path::unit_tests_data() + "POSCAR");

		CHECK(vasp_file.lattice()[0][0] == 3.33536661_a);
		CHECK(vasp_file.lattice()[0][1] == 3.33536661_a);
		CHECK(vasp_file.lattice()[0][2] == 0.0_a);
		CHECK(vasp_file.lattice()[1][0] == -3.33536661_a);
		CHECK(vasp_file.lattice()[1][1] == 3.33536661_a);
		CHECK(vasp_file.lattice()[1][2] == 0.0_a);
		CHECK(vasp_file.lattice()[2][0] == 0.0_a);
		CHECK(vasp_file.lattice()[2][1] == 0.0_a);
		CHECK(vasp_file.lattice()[2][2] == 33.3536660997_a);		
		
		CHECK(vasp_file.size() == 5);

		CHECK(vasp_file.atoms()[0] == "Ni");
		CHECK(vasp_file.atoms()[1] == "Ni");
		CHECK(vasp_file.atoms()[2] == "Ni");
		CHECK(vasp_file.atoms()[3] == "Ni");		
		CHECK(vasp_file.atoms()[4] == "Ni");

		CHECK(vasp_file.positions()[0][0] == 0.0_a);
		CHECK(vasp_file.positions()[0][1] == 0.0_a);
		CHECK(vasp_file.positions()[0][2] == 0.0_a);

		CHECK(vasp_file.positions()[1][0] == 0.0_a);
		CHECK(vasp_file.positions()[1][1] == 3.33536661_a);
		CHECK(vasp_file.positions()[1][2] == 3.33536661_a);

		CHECK(vasp_file.positions()[2][0] == 0.0_a);
		CHECK(vasp_file.positions()[2][1] == 0.0_a);
		CHECK(vasp_file.positions()[2][2] == 6.6707332199_a);

		CHECK(vasp_file.positions()[3][0] == 0.0_a);
		CHECK(vasp_file.positions()[3][1] == 3.33536661_a);
		CHECK(vasp_file.positions()[3][2] == 10.0060998299_a);

		CHECK(vasp_file.positions()[4][0] == 0.0_a);
		CHECK(vasp_file.positions()[4][1] == 0.0_a);
		CHECK(vasp_file.positions()[4][2] == 13.3414664399_a);
		
	}

	SECTION("POSCAR with spaces in species"){
		
		parse::poscar vasp_file(config::path::unit_tests_data() + "co.POSCAR");
		CHECK(vasp_file.size() == 128);
		
	}

}
#endif
