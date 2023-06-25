/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__SYSTEMS__IONS
#define INQ__SYSTEMS__IONS

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <spglib.h>

#include <input/species.hpp>
#include <ions/unit_cell.hpp>
#include <gpu/array.hpp>


namespace inq {
namespace systems {

class ions {

	inq::ions::unit_cell cell_;
	std::vector<input::species> atoms_;
	std::vector<vector3<double>> coordinates_;
	std::vector<vector3<double>> velocities_;	

	template <typename PositionType>
	void add_atom(input::species const & element, PositionType const & position){
		atoms_.push_back(element);
		coordinates_.push_back(in_atomic_units(position));
		velocities_.push_back(vector3<double>(0.0, 0.0, 0.0));					
	}

public:

	ions(inq::ions::unit_cell arg_cell_input):
		cell_(std::move(arg_cell_input)){
	}

	auto & atoms() const {
		return atoms_;
	}
	
	auto & coordinates() const {
		return coordinates_;
	}
	
	auto & coordinates() {
		return coordinates_;
	}

	auto & velocities() const {
		return velocities_;
	}
    
	auto & velocities() {
		return velocities_;
	}
	
	auto symmetry_string() const{

		assert(size() > 0);
		
		char symbol[11];
		
		std::vector<int> types(size());
		std::vector<double> positions(3*size());
		
		for(int iatom = 0; iatom < size(); iatom++){
			types[iatom] = atoms()[iatom].atomic_number();
			auto pos = cell_.metric().to_contravariant(cell_.position_in_cell(coordinates()[iatom]));
			positions[3*iatom + 0] = pos[0];
			positions[3*iatom + 1] = pos[1];
			positions[3*iatom + 2] = pos[2];
		}

		double amat[9];
		amat[0] = cell_.lattice(0)[0];
		amat[1] = cell_.lattice(0)[1];
		amat[2] = cell_.lattice(0)[2];
		amat[3] = cell_.lattice(1)[0];
		amat[4] = cell_.lattice(1)[1];
		amat[5] = cell_.lattice(1)[2];
		amat[6] = cell_.lattice(2)[0];
		amat[7] = cell_.lattice(2)[1];
		amat[8] = cell_.lattice(2)[2];
		
		auto symnum = spg_get_international(symbol, reinterpret_cast<double (*)[3]>(amat), reinterpret_cast<double (*)[3]>(positions.data()), types.data(), size(), 1e-4);
		return symbol + std::string(" (number ") + std::to_string(symnum) + std::string(")");
	}
	
	auto & cell() const {
		return cell_;
	}

	auto insert(std::string const & symbol, vector3<quantity<magnitude::length>> const & pos){
		add_atom(input::species(pseudo::element(symbol)), pos);
	}
	
	auto insert(input::species const & sp, vector3<quantity<magnitude::length>> const & pos){
		add_atom(sp, pos);
	}

	template <class ContainerType>
	auto insert(ContainerType const & container){
		for(auto atom : container) add_atom(atom.species(), atom.position());
	}

	auto insert_fractional(input::species const & sp, vector3<double, contravariant> const & pos){
		add_atom(sp, cell_.metric().to_cartesian(pos));
	}

	auto insert_fractional(std::string const & symbol, vector3<double, contravariant> const & pos){
		insert_fractional(input::species(pseudo::element(symbol)), pos);
	}

	int size() const {
		return (long) coordinates_.size();
	}
	
	template <class output_stream>
	void info(output_stream & out) const {
		out << "GEOMETRY:" << std::endl;
		out << "  Number of atoms = " << size() << std::endl;
		out << std::endl;
	}
	
	template<class OStream>
	friend OStream& operator<<(OStream& os, ions const& self){
		self.info(os);
		return os;
	}

    
};

}
}
#endif

#ifdef INQ_SYSTEMS_IONS_UNIT_TEST
#undef INQ_SYSTEMS_IONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("Create empty and add an atom"){
		
		auto dcc = 1.42_A;
		auto aa = sqrt(3)*dcc;
		auto lz = 10.0_b;
		systems::ions ions(ions::unit_cell::lattice(aa*vector3{1.0, 0.0, 0.0}, aa*vector3{-1.0/2.0, sqrt(3.0)/2.0, 0.0}, {0.0_b, 0.0_b, lz}).periodicity(2));
		
		CHECK(ions.cell().periodicity() == 2);
		
    CHECK(ions.size() == 0);

    ions.insert(pseudo::element("Xe"), {1000.0_b, -200.0_b, 6.0_b});

    CHECK(ions.size() == 1);
    CHECK(ions.atoms()[0].atomic_number() == 54);
    CHECK(ions.atoms()[0] == pseudo::element(54));
    CHECK(ions.atoms()[0].charge() == -54.0_a);
    CHECK(ions.atoms()[0].mass() == 239333.5935636_a);
    CHECK(ions.coordinates()[0][0] == 1000.0_a);
    CHECK(ions.coordinates()[0][1] == -200.0_a);
    CHECK(ions.coordinates()[0][2] == 6.0_a);
		CHECK(ions.velocities()[0][0] == 0.0_a);
    CHECK(ions.velocities()[0][1] == 0.0_a);
    CHECK(ions.velocities()[0][2] == 0.0_a);
		
    ions.coordinates()[0][0] += 8;  
    
    CHECK(ions.coordinates()[0][0] == 1008.0_a);
		CHECK(ions.velocities().size() == ions.coordinates().size());
		
  }
	
	SECTION("Read an xyz file"){
		
		systems::ions ions(ions::unit_cell::cubic(66.6_A).finite());

    ions.insert(input::parse_xyz(config::path::unit_tests_data() + "benzene.xyz"));
		
    CHECK(ions.size() == 12);
    
    CHECK(ions.atoms()[2] == pseudo::element("C"));
    CHECK(ions.atoms()[2].charge() == -6.0_a);
    CHECK(ions.atoms()[2].mass() == 21892.1617296_a);
    CHECK(ions.coordinates()[2][0] == 2.2846788549_a);
    CHECK(ions.coordinates()[2][1] == -1.3190288178_a);
    CHECK(ions.coordinates()[2][2] == 0.0_a);

    CHECK(ions.atoms()[11] == pseudo::element("H"));
    CHECK(ions.atoms()[11].charge() == -1.0_a);
    CHECK(ions.atoms()[11].mass() == 1837.17994584_a);
    CHECK(ions.coordinates()[11][0] == -4.0572419367_a);
    CHECK(ions.coordinates()[11][1] == 2.343260364_a);
    CHECK(ions.coordinates()[11][2] == 0.0_a);
		CHECK(ions.velocities()[11][0] == 0.0_a);
    CHECK(ions.velocities()[11][1] == 0.0_a);
    CHECK(ions.velocities()[11][2] == 0.0_a);

		CHECK(ions.velocities().size() == ions.coordinates().size());
		
    ions.insert(pseudo::element("Cl"), {-3.0_b, 4.0_b, 5.0_b});

    CHECK(ions.size() == 13);
    CHECK(ions.atoms()[12].atomic_number() == 17);
    CHECK(ions.atoms()[12] == pseudo::element(17));
    CHECK(ions.atoms()[12].charge() == -17.0_a);
    CHECK(ions.atoms()[12].mass() == 64614.105771_a);
    CHECK(ions.coordinates()[12][0] == -3.0_a);
    CHECK(ions.coordinates()[12][1] == 4.0_a);
    CHECK(ions.coordinates()[12][2] == 5.0_a);
		CHECK(ions.velocities()[12][0] == 0.0_a);
    CHECK(ions.velocities()[12][1] == 0.0_a);
    CHECK(ions.velocities()[12][2] == 0.0_a);

		CHECK(ions.velocities().size() == ions.coordinates().size());
		
  }

}
#endif
