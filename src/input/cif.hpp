/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__CIF
#define INQ__INPUT__CIF

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <gemmi/cif.hpp>
#include <gemmi/symmetry.hpp>
#include <vector>

#include <input/atom.hpp>
#include <magnitude/length.hpp>
#include <math/vector3.hpp>

namespace inq {
namespace input {

class cif {

	std::vector<vector3<double>> lattice_vectors_;
	std::vector<input::atom> geo_;

  static auto to_double(std::string const * strptr, std::string const & cif_file_name){
    if(strptr == NULL) throw std::runtime_error("Error: cannot read cell parameters from CIF file '" + cif_file_name + "'.");
    return stod(*strptr);
  }
  
public:
	
	cif(const std::string & cif_file_name):
		lattice_vectors_(3)
	{
    namespace cif = gemmi::cif;
    cif::Document doc = cif::read_file(cif_file_name);
		
		//THE CELL
		auto length_a = to_double(doc.sole_block().find_value("_cell_length_a"), cif_file_name);
    auto length_b = to_double(doc.sole_block().find_value("_cell_length_b"), cif_file_name);
    auto length_c = to_double(doc.sole_block().find_value("_cell_length_c"), cif_file_name);
    auto angle_alpha = to_double(doc.sole_block().find_value("_cell_angle_alpha"), cif_file_name);
    auto angle_beta = to_double(doc.sole_block().find_value("_cell_angle_beta"), cif_file_name);
    auto angle_gamma = to_double(doc.sole_block().find_value("_cell_angle_gamma"), cif_file_name);

    angle_alpha *= M_PI/180.0;
    angle_beta  *= M_PI/180.0;
    angle_gamma *= M_PI/180.0;

    using namespace inq::magnitude;
    auto angstrom = in_atomic_units(1.0_A);
    
    lattice_vectors_[0] = angstrom*length_a*vector3<double>({1.0, 0.0, 0.0});
    lattice_vectors_[1] = angstrom*length_b*vector3<double>({cos(angle_gamma), sin(angle_gamma), 0.0});

    auto cosbeta = cos(angle_beta);
    auto c1 = (cos(angle_alpha) - cos(angle_beta)*cos(angle_gamma))/sin(angle_gamma);
    auto c2 = sqrt(1.0 - cosbeta*cosbeta - c1*c1);
    lattice_vectors_[2] = angstrom*length_c*vector3<double>({cosbeta, c1, c2});

		ions::unit_cell cell(lattice_vectors_);
		
		//SYMMETRIES
		auto symmetries = doc.sole_block().find_loop("_space_group_symop_operation_xyz");
		if(symmetries.length() == 0) symmetries = doc.sole_block().find_loop("_symmetry_equiv_pos_as_xyz");

		
		//THE ATOMS
		auto atom_symbols = doc.sole_block().find_loop("_atom_site_label");
		if(atom_symbols.length() == 0) throw std::runtime_error("Error: cannot find the atoms in CIF file '" + cif_file_name + "'.");

		auto natoms = atom_symbols.length();
		
		auto atom_x = doc.sole_block().find_loop("_atom_site_fract_x");
		auto atom_y = doc.sole_block().find_loop("_atom_site_fract_y");
		auto atom_z = doc.sole_block().find_loop("_atom_site_fract_z");

		if(natoms != atom_x.length() or natoms != atom_y.length() or natoms != atom_z.length()){
			throw std::runtime_error("Error: read the atomic coordinates in CIF file '" + cif_file_name + "'.");
		}

		std::vector<std::string> symbols;
		std::vector<vector3<double, contravariant>> positions;
		
		for(int iatom = 0; iatom < natoms; iatom++){
			auto symbol = atom_symbols[iatom];
			symbol.erase(std::remove_if(symbol.begin(), symbol.end(), ::isdigit), symbol.end()); //remove the numbers
			auto pos = vector3<double, contravariant>{stod(atom_x[iatom]), stod(atom_y[iatom]), stod(atom_z[iatom])};
			pos = cell.position_in_cell(pos);
			
			symbols.emplace_back(symbol);
			positions.emplace_back(pos);

			for(auto & symm: symmetries){
				symm.erase(std::remove(symm.begin(), symm.end(), '\''), symm.end());
				auto op = gemmi::parse_triplet(symm);
				auto symm_pos = vector3<double, contravariant>{op.apply_to_xyz(std::array<double, 3>{pos[0], pos[1], pos[2]})};
				symm_pos = cell.position_in_cell(symm_pos);
				
				auto duplicated = false;
				
				for(int jatom = 0; jatom < (int) symbols.size(); jatom++){
					if(cell.metric().norm(symm_pos - positions[jatom]) > 1e-8) continue;
					
					if(symbol != symbols[jatom]) throw std::runtime_error("Error: the file '" + cif_file_name + "' contains two different atoms in the same position.");
					duplicated = true;
					break;
				}

				if(duplicated) continue;

				symbols.emplace_back(symbol);
				positions.emplace_back(symm_pos);
			}
		}

		assert(symbols.size() == positions.size());
		
		for(auto iatom = 0; iatom < (int) symbols.size(); iatom++){
			geo_.emplace_back(input::species(symbols[iatom]), cell.metric().to_cartesian(positions[iatom]));
		}
		
	}
	
	auto num_atoms() const {
		return long(geo_.size());
	}	
	
	auto & lattice() const {
		return lattice_vectors_;
	}
	
	auto & atoms() const {
		return geo_;
	}

	auto cell() const {
		return ions::unit_cell(lattice());
	}
	
};


}
}
#endif

#ifdef INQ_INPUT_CIF_UNIT_TEST
#undef INQ_INPUT_CIF_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <config/path.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("Al"){
		
		input::cif cif_file(config::path::unit_tests_data() + "Al.cif");

		CHECK(cif_file.lattice()[0][0] == 7.634890386_a);
		CHECK(cif_file.lattice()[0][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[0][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][1] == 7.634890386_a);
		CHECK(cif_file.lattice()[1][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][2] == 7.634890386_a);		
		
		CHECK(cif_file.num_atoms() == 4);

		CHECK(cif_file.atoms()[0].species() == "Al");
		CHECK(cif_file.atoms()[1].species() == "Al");
		CHECK(cif_file.atoms()[2].species() == "Al");		
		CHECK(cif_file.atoms()[3].species() == "Al");

		CHECK(cif_file.atoms()[0].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[2] == Approx(0.0).margin(1e-12));
		
		CHECK(cif_file.atoms()[1].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[1] == -3.817445193_a);
		CHECK(cif_file.atoms()[1].position()[2] == -3.817445193_a);

		CHECK(cif_file.atoms()[2].position()[0] == -3.817445193_a);
		CHECK(cif_file.atoms()[2].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[2].position()[2] == -3.817445193_a);

		CHECK(cif_file.atoms()[3].position()[0] == -3.817445193_a);
		CHECK(cif_file.atoms()[3].position()[1] == -3.817445193_a);
		CHECK(cif_file.atoms()[3].position()[2] == Approx(0.0).margin(1e-12));
	}
	
	SECTION("SrTiO3"){
		
		input::cif cif_file(config::path::unit_tests_data() + "9002806.cif");

		CHECK(cif_file.lattice()[0][0] == 10.4316661532_a);
		CHECK(cif_file.lattice()[0][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[0][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][1] == 10.4316661532_a);
		CHECK(cif_file.lattice()[1][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][2] == 14.7525249371_a);		
		
		CHECK(cif_file.num_atoms() == 20);

		//Conversions for units and cell convention to compare with obabel
		//generated coordinates
		//
		//  1.38005 ->  2.6079165383
		//  1.95168 ->  3.6881312343
		//  2.76010 -> -5.2158330766
		//  3.90335 -> -7.3762624686
		//  4.14015 -> -2.6079165383
		//  5.85503 -> -3.6881312343
		
		//Sr         0.00000        0.00000        1.95168
		CHECK(cif_file.atoms()[0].species() == "Sr");
		CHECK(cif_file.atoms()[0].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[2] == 3.6881312343_a);

		//Sr         0.00000        0.00000        5.85503
		CHECK(cif_file.atoms()[1].species() == "Sr");
		CHECK(cif_file.atoms()[1].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[2] == -3.6881312343_a);

		//Sr         2.76010        2.76010        5.85503
		CHECK(cif_file.atoms()[2].species() == "Sr");
		CHECK(cif_file.atoms()[2].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[2].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[2].position()[2] == -3.6881312343_a);

		//Sr         2.76010        2.76010        1.95167
		CHECK(cif_file.atoms()[3].species() == "Sr");
		CHECK(cif_file.atoms()[3].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[3].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[3].position()[2] ==  3.6881312343_a);

		//Ti         2.76010        0.00000        0.00000
		CHECK(cif_file.atoms()[4].species() == "Ti");
		CHECK(cif_file.atoms()[4].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[4].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[4].position()[2] == Approx(0.0).margin(1e-12));

		//Ti         2.76010        0.00000        3.90335
		CHECK(cif_file.atoms()[5].species() == "Ti");
		CHECK(cif_file.atoms()[5].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[5].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[5].position()[2] == -7.3762624686_a);

		//Ti         0.00000        2.76010        3.90335
		CHECK(cif_file.atoms()[6].species() == "Ti");
		CHECK(cif_file.atoms()[6].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[6].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[6].position()[2] == -7.3762624686_a);

		//Ti         0.00000        2.76010        0.00000
		CHECK(cif_file.atoms()[7].species() == "Ti");
		CHECK(cif_file.atoms()[7].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[7].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[7].position()[2] == Approx(0.0).margin(1e-12));

		//O          0.00000        2.76010        1.95168
		CHECK(cif_file.atoms()[8].species() == "O");
		CHECK(cif_file.atoms()[8].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[8].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[8].position()[2] == 3.6881312343_a);

		//O          0.00000        2.76010        5.85503
		CHECK(cif_file.atoms()[9].species() == "O");
		CHECK(cif_file.atoms()[9].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[9].position()[1] == -5.2158330766_a);
		CHECK(cif_file.atoms()[9].position()[2] == -3.6881312343_a);

		//O          2.76010        0.00000        5.85503
		CHECK(cif_file.atoms()[10].species() == "O");
		CHECK(cif_file.atoms()[10].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[10].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[10].position()[2] == -3.6881312343_a);

		//O          2.76010        0.00000        1.95167
		CHECK(cif_file.atoms()[11].species() == "O");
		CHECK(cif_file.atoms()[11].position()[0] == -5.2158330766_a);
		CHECK(cif_file.atoms()[11].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[11].position()[2] ==  3.6881312343_a);

		//O          4.14015        1.38005        0.00000
		CHECK(cif_file.atoms()[12].species() == "O");
		CHECK(cif_file.atoms()[12].position()[0] == -2.6079165383_a);
		CHECK(cif_file.atoms()[12].position()[1] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[12].position()[2] == Approx(0.0).margin(1e-12));

		//4.14015        1.38005        3.90335
		CHECK(cif_file.atoms()[13].species() == "O");
		CHECK(cif_file.atoms()[13].position()[0] == -2.6079165383_a);
		CHECK(cif_file.atoms()[13].position()[1] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[13].position()[2] == -7.3762624686_a);

		//O          1.38005        4.14015        3.90335
		CHECK(cif_file.atoms()[14].species() == "O");
		CHECK(cif_file.atoms()[14].position()[0] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[14].position()[1] == -2.6079165383_a);
		CHECK(cif_file.atoms()[14].position()[2] == -7.3762624686_a);

		//O          1.38005        1.38005        3.90335
		CHECK(cif_file.atoms()[15].species() == "O");
		CHECK(cif_file.atoms()[15].position()[0] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[15].position()[1] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[15].position()[2] == -7.3762624686_a);

		//O          4.14015        4.14015        3.90335
		CHECK(cif_file.atoms()[16].species() == "O");
		CHECK(cif_file.atoms()[16].position()[0] == -2.6079165383_a);
		CHECK(cif_file.atoms()[16].position()[1] == -2.6079165383_a);
		CHECK(cif_file.atoms()[16].position()[2] == -7.3762624686_a);

		//O          4.14015        4.14015        0.00000
		CHECK(cif_file.atoms()[17].species() == "O");
		CHECK(cif_file.atoms()[17].position()[0] == -2.6079165383_a);
		CHECK(cif_file.atoms()[17].position()[1] == -2.6079165383_a);
		CHECK(cif_file.atoms()[17].position()[2] == Approx(0.0).margin(1e-12));

		//O          1.38005        1.38005        0.00000
		CHECK(cif_file.atoms()[18].species() == "O");
		CHECK(cif_file.atoms()[18].position()[0] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[18].position()[1] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[18].position()[2] == Approx(0.0).margin(1e-12));

		//O          1.38005        4.14015        0.00000
		CHECK(cif_file.atoms()[19].species() == "O");
		CHECK(cif_file.atoms()[19].position()[0] ==  2.6079165383_a);
		CHECK(cif_file.atoms()[19].position()[1] == -2.6079165383_a);
		CHECK(cif_file.atoms()[19].position()[2] == Approx(0.0).margin(1e-12));
		
	}

	SECTION("Ca2PI symmetrized"){
		
		input::cif cif_file(config::path::unit_tests_data() + "Ca2PI_symm.cif");

		CHECK(cif_file.lattice()[0][0] == 8.1469916149_a);
		CHECK(cif_file.lattice()[0][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[0][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][0] == -4.0734958074_a);
		CHECK(cif_file.lattice()[1][1] ==  7.0555017029_a);
		CHECK(cif_file.lattice()[1][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][2] == 42.0773092856_a);		
		
		CHECK(cif_file.num_atoms() == 12);
		
		CHECK(cif_file.atoms()[0].species() == "Ca");
		CHECK(cif_file.atoms()[0].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[2] == 9.6727962369_a);
		
		CHECK(cif_file.atoms()[1].species() == "Ca");
		CHECK(cif_file.atoms()[1].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[2] == -9.6727962369_a);

		CHECK(cif_file.atoms()[2].species() == "Ca");
		CHECK(cif_file.atoms()[2].position()[0] ==  -4.0734958074_a);
		CHECK(cif_file.atoms()[2].position()[1] ==   2.3518339010_a);
		CHECK(cif_file.atoms()[2].position()[2] == -18.3787432869_a);
		
		CHECK(cif_file.atoms()[3].species() == "Ca");
		CHECK(cif_file.atoms()[3].position()[0] ==  -4.0734958074_a);
		CHECK(cif_file.atoms()[3].position()[1] ==   2.3518339010_a);
		CHECK(cif_file.atoms()[3].position()[2] ==   4.3529735250_a);
		
		CHECK(cif_file.atoms()[4].species() == "Ca");
		CHECK(cif_file.atoms()[4].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[4].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[4].position()[2] ==  -4.3529735250_a);
		
		CHECK(cif_file.atoms()[5].species() == "Ca");
		CHECK(cif_file.atoms()[5].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[5].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[5].position()[2] ==  18.3787432869_a);

		CHECK(cif_file.atoms()[6].species() == "P");
		CHECK(cif_file.atoms()[6].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[6].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[6].position()[2] == -21.0386546428_a);

		CHECK(cif_file.atoms()[7].species() == "P");
		CHECK(cif_file.atoms()[7].position()[0] == -4.0734958074_a);
		CHECK(cif_file.atoms()[7].position()[1] ==  2.3518339010_a);
		CHECK(cif_file.atoms()[7].position()[2] == -7.0128848809_a);

		CHECK(cif_file.atoms()[8].species() == "P");
		CHECK(cif_file.atoms()[8].position()[0] ==  4.0734958074_a);
		CHECK(cif_file.atoms()[8].position()[1] == -2.3518339010_a);
		CHECK(cif_file.atoms()[8].position()[2] ==  7.0128848809_a);

		CHECK(cif_file.atoms()[9].species() == "I");
		CHECK(cif_file.atoms()[9].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[9].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[9].position()[2] == Approx(0.0).margin(1e-12));

		CHECK(cif_file.atoms()[10].species() == "I");
		CHECK(cif_file.atoms()[10].position()[0] == -4.0734958074_a);
		CHECK(cif_file.atoms()[10].position()[1] ==  2.3518339010_a);
		CHECK(cif_file.atoms()[10].position()[2] == 14.0257697619_a);

		CHECK(cif_file.atoms()[11].species() == "I");
		CHECK(cif_file.atoms()[11].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[11].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[11].position()[2] == -14.0257697619_a);
		
	}
	
	SECTION("Ca2PI not symmetrized"){
		
		input::cif cif_file(config::path::unit_tests_data() + "Ca2PI.cif");

		CHECK(cif_file.lattice()[0][0] == 8.1469916149_a);
		CHECK(cif_file.lattice()[0][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[0][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][0] == -4.0734958074_a);
		CHECK(cif_file.lattice()[1][1] ==  7.0555017029_a);
		CHECK(cif_file.lattice()[1][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][2] == 42.0773092856_a);		
		
		CHECK(cif_file.num_atoms() == 12);

		//the order her is to match the symmetrized test above
		CHECK(cif_file.atoms()[1].species() == "Ca");
		CHECK(cif_file.atoms()[1].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[2] == 9.6727962369_a);

		CHECK(cif_file.atoms()[4].species() == "Ca");
		CHECK(cif_file.atoms()[4].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[4].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[4].position()[2] == -9.6727962369_a);

		CHECK(cif_file.atoms()[3].species() == "Ca");
		CHECK(cif_file.atoms()[3].position()[0] ==  -4.0734958074_a);
		CHECK(cif_file.atoms()[3].position()[1] ==   2.3518339010_a);
		CHECK(cif_file.atoms()[3].position()[2] == -18.3787432869_a);

		CHECK(cif_file.atoms()[0].species() == "Ca");
		CHECK(cif_file.atoms()[0].position()[0] ==  -4.0734958074_a);
		CHECK(cif_file.atoms()[0].position()[1] ==   2.3518339010_a);
		CHECK(cif_file.atoms()[0].position()[2] ==   4.3529735250_a);

		CHECK(cif_file.atoms()[5].species() == "Ca");
		CHECK(cif_file.atoms()[5].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[5].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[5].position()[2] ==  -4.3529735250_a);

		CHECK(cif_file.atoms()[2].species() == "Ca");
		CHECK(cif_file.atoms()[2].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[2].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[2].position()[2] ==  18.3787432869_a);

		CHECK(cif_file.atoms()[7].species() == "P");
		CHECK(cif_file.atoms()[7].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[7].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[7].position()[2] == -21.0386546428_a);

		CHECK(cif_file.atoms()[8].species() == "P");
		CHECK(cif_file.atoms()[8].position()[0] == -4.0734958074_a);
		CHECK(cif_file.atoms()[8].position()[1] ==  2.3518339010_a);
		CHECK(cif_file.atoms()[8].position()[2] == -7.0128848809_a);
		
		CHECK(cif_file.atoms()[6].species() == "P");
		CHECK(cif_file.atoms()[6].position()[0] ==  4.0734958074_a);
		CHECK(cif_file.atoms()[6].position()[1] == -2.3518339010_a);
		CHECK(cif_file.atoms()[6].position()[2] ==  7.0128848809_a);

		CHECK(cif_file.atoms()[9].species() == "I");
		CHECK(cif_file.atoms()[9].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[9].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[9].position()[2] == Approx(0.0).margin(1e-12));

		CHECK(cif_file.atoms()[10].species() == "I");
		CHECK(cif_file.atoms()[10].position()[0] == -4.0734958074_a);
		CHECK(cif_file.atoms()[10].position()[1] ==  2.3518339010_a);
		CHECK(cif_file.atoms()[10].position()[2] == 14.0257697619_a);

		CHECK(cif_file.atoms()[11].species() == "I");
		CHECK(cif_file.atoms()[11].position()[0] ==   4.0734958074_a);
		CHECK(cif_file.atoms()[11].position()[1] ==  -2.3518339010_a);
		CHECK(cif_file.atoms()[11].position()[2] == -14.0257697619_a);
		
	}
	
	SECTION("Na"){
		
		input::cif cif_file(config::path::unit_tests_data() + "Na.cif");

		//These lattice vectors match openbabel
		CHECK(cif_file.lattice()[0][0] == 17.7976863062_a);
		CHECK(cif_file.lattice()[0][1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[0][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[1][0] == 16.3923048366_a);
		CHECK(cif_file.lattice()[1][1] == 6.9318092872_a);
		CHECK(cif_file.lattice()[1][2] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.lattice()[2][0] == 16.3923048366_a);
		CHECK(cif_file.lattice()[2][1] ==  3.3234384423_a);
		CHECK(cif_file.lattice()[2][2] ==  6.0831518898_a);		
		
		CHECK(cif_file.num_atoms() == 3);
		
		CHECK(cif_file.atoms()[0].species() == "Na");
		CHECK(cif_file.atoms()[0].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[2] == Approx(0.0).margin(1e-12));

		CHECK(cif_file.atoms()[1].species() == "Na");
		CHECK(cif_file.atoms()[1].position()[0] == 11.2403978126_a);
		CHECK(cif_file.atoms()[1].position()[1] ==  2.2789211505_a);
		CHECK(cif_file.atoms()[1].position()[2] ==  1.3517980130_a);

		CHECK(cif_file.atoms()[2].species() == "Na");
		CHECK(cif_file.atoms()[2].position()[0] == -11.2403978126_a);
		CHECK(cif_file.atoms()[2].position()[1] ==  -2.2789211505_a);
		CHECK(cif_file.atoms()[2].position()[2] ==  -1.3517980130_a);
	}
	
}
#endif
