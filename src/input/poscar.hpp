/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__POSCAR
#define INQ__INPUT__POSCAR

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
#include <string>
#include <iostream>
#include <sstream>


#include <math/vector3.hpp>

#include <pseudopod/element.hpp>
#include <input/atom.hpp>
#include <input/species.hpp>
#include <ions/unit_cell.hpp>
#include <magnitude/length.hpp>

namespace inq {
namespace input {

class poscar {

	std::vector<math::vector3<double>> lattice_vectors_;
	std::vector<input::atom> geo_;
	
public:
	
	poscar(const std::string & poscar_file_name):
		lattice_vectors_(3)
	{

		using namespace inq::magnitude;

		std::ifstream poscar_file(poscar_file_name.c_str());
		
		assert(poscar_file.is_open());

		{
			std::string comment_line;
			std::getline(poscar_file, comment_line);
		}
		
		double scaling_factor;
		
		poscar_file >> scaling_factor;

		std::cout << scaling_factor << std::endl;
		
		for(int idir = 0; idir < 3; idir++)	{
			poscar_file >> lattice_vectors_[idir];
			lattice_vectors_[idir] *= (scaling_factor*1.0_A).in_atomic_units();
		}

		std::vector<std::string> species;

		std::string species_line;
		std::getline(poscar_file, species_line);
		std::getline(poscar_file, species_line);

		std::istringstream iss(species_line);

		std::string species_name;
		while (std::getline(iss, species_name, ' ')){
			species.push_back(species_name);
		}

		std::vector<int> species_num(species.size());

		auto num_atoms = 0;
		for(unsigned ispecies = 0; ispecies < species.size(); ispecies++){
			poscar_file >> species_num[ispecies];
			num_atoms += species_num[ispecies];
		}

		std::string coordinate_type;

		poscar_file >> coordinate_type;

		if(coordinate_type == "Direct") {

			ions::unit_cell cell(lattice_vectors_);
			
			for(unsigned ispecies = 0; ispecies < species_num.size(); ispecies++){
				for(int iatom = 0; iatom < species_num[ispecies]; iatom++){
					math::vector3<double, math::contravariant> pos;
					poscar_file >> pos;
					geo_.emplace_back(input::species(species[ispecies]), cell.metric().to_cartesian(pos));
				}
			}
			
		} else if(coordinate_type == "Cartesian") {
			for(unsigned ispecies = 0; ispecies < species_num.size(); ispecies++){
				for(int iatom = 0; iatom < species_num[ispecies]; iatom++){
					math::vector3<double> pos;
					poscar_file >> pos;
					geo_.emplace_back(input::species(species[ispecies]), pos);
				}
			}
		} else {
			throw std::runtime_error("Unsupported POSCAR file");
		}

	}

	auto num_atoms() const {
		return long(geo_.size());
	}	
	
	auto & lattice_vectors() const {
		return lattice_vectors_;
	}
	
	auto & geo() const {
		return geo_;
	}
	
};


}
}

#ifdef INQ_INPUT_POSCAR_UNIT_TEST
#undef INQ_INPUT_POSCAR_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <config/path.hpp>

TEST_CASE("function ions::poscar", "[inq::input::poscar]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("BN"){
	
		input::poscar vasp_file(config::path::unit_tests_data() + "bn.poscar");

		CHECK(vasp_file.lattice_vectors()[0][0] == 0.0_a);
		CHECK(vasp_file.lattice_vectors()[0][1] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[0][2] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[1][0] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[1][1] == 0.0_a);
		CHECK(vasp_file.lattice_vectors()[1][2] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[2][0] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[2][1] == 3.3731611325_a);
		CHECK(vasp_file.lattice_vectors()[2][2] == 0.0_a);		
		
		CHECK(vasp_file.num_atoms() == 2);

		CHECK(vasp_file.geo()[0].species() == "B");
		CHECK(vasp_file.geo()[1].species() == "N");

		CHECK(vasp_file.geo()[0].position()[0] == 0.0_a);
		CHECK(vasp_file.geo()[0].position()[1] == 0.0_a);
		CHECK(vasp_file.geo()[0].position()[2] == 0.0_a);
		CHECK(vasp_file.geo()[1].position()[0] == 1.6865805662_a);
		CHECK(vasp_file.geo()[1].position()[1] == 1.6865805662_a);
		CHECK(vasp_file.geo()[1].position()[2] == 1.6865805662_a);		
	}

}


#endif

#endif
