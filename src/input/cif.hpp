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
#include <vector>

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
    auto length_a = to_double(doc.sole_block().find_value("_cell_length_a"), cif_file_name);
    auto length_b = to_double(doc.sole_block().find_value("_cell_length_b"), cif_file_name);
    auto length_c = to_double(doc.sole_block().find_value("_cell_length_c"), cif_file_name);
    auto angle_alpha = to_double(doc.sole_block().find_value("_cell_angle_alpha"), cif_file_name);
    auto angle_beta = to_double(doc.sole_block().find_value("_cell_angle_beta"), cif_file_name);
    auto angle_gamma = to_double(doc.sole_block().find_value("_cell_angle_gamma"), cif_file_name);

		//    std::cout << length_a << '\t' << length_b << '\t' << length_c << std::endl;    
		//    std::cout << angle_alpha << '\t' << angle_beta << '\t' << angle_gamma << std::endl;

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

		auto atom_symbols = doc.sole_block().find_loop("_atom_site_type_symbol");
		if(atom_symbols.length() == 0) atom_symbols = doc.sole_block().find_loop("_atom_site_label");
		if(atom_symbols.length() == 0) throw std::runtime_error("Error: cannot find the atoms in CIF file '" + cif_file_name + "'.");

		auto natoms = atom_symbols.length();
		
		auto atom_x = doc.sole_block().find_loop("_atom_site_fract_x");
		auto atom_y = doc.sole_block().find_loop("_atom_site_fract_y");
		auto atom_z = doc.sole_block().find_loop("_atom_site_fract_z");

		if(natoms != atom_x.length() or natoms != atom_y.length() or natoms != atom_z.length()){
			throw std::runtime_error("Error: read the atomic coordinates in CIF file '" + cif_file_name + "'.");
		}
		
		ions::unit_cell cell(lattice_vectors_);

		for(int iatom = 0; iatom < natoms; iatom++){
			auto symbol = atom_symbols[iatom];
			symbol.erase(std::remove_if(symbol.begin(), symbol.end(), ::isdigit), symbol.end()); //remove the numbers
			auto pos = vector3<double, contravariant>{stod(atom_x[iatom]), stod(atom_y[iatom]), stod(atom_z[iatom])};
			geo_.emplace_back(input::species(symbol), cell.metric().to_cartesian(pos));
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
		CHECK(cif_file.atoms()[1].position()[1] == 3.817445193_a);
		CHECK(cif_file.atoms()[1].position()[2] == 3.817445193_a);

		CHECK(cif_file.atoms()[2].position()[0] == 3.817445193_a);
		CHECK(cif_file.atoms()[2].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[2].position()[2] == 3.817445193_a);

		CHECK(cif_file.atoms()[3].position()[0] == 3.817445193_a);
		CHECK(cif_file.atoms()[3].position()[1] == 3.817445193_a);
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
		/*
		CHECK(cif_file.num_atoms() == 4);

		CHECK(cif_file.atoms()[0].species() == "Al");
		CHECK(cif_file.atoms()[1].species() == "Al");
		CHECK(cif_file.atoms()[2].species() == "Al");
		CHECK(cif_file.atoms()[3].species() == "Al");

		CHECK(cif_file.atoms()[0].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[0].position()[2] == Approx(0.0).margin(1e-12));
		
		CHECK(cif_file.atoms()[1].position()[0] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[1].position()[1] == 3.817445193_a);
		CHECK(cif_file.atoms()[1].position()[2] == 3.817445193_a);

		CHECK(cif_file.atoms()[2].position()[0] == 3.817445193_a);
		CHECK(cif_file.atoms()[2].position()[1] == Approx(0.0).margin(1e-12));
		CHECK(cif_file.atoms()[2].position()[2] == 3.817445193_a);

		CHECK(cif_file.atoms()[3].position()[0] == 3.817445193_a);
		CHECK(cif_file.atoms()[3].position()[1] == 3.817445193_a);
		CHECK(cif_file.atoms()[3].position()[2] == Approx(0.0).margin(1e-12));
		*/
	}
	
}
#endif
