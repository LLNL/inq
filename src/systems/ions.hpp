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

#include <ions/geometry.hpp>
#include <ions/unit_cell.hpp>
#include <gpu/array.hpp>
#include <systems/box.hpp>

namespace inq {
namespace systems {

class ions {

public:

	ions(inq::ions::unit_cell arg_cell_input):
		cell_(std::move(arg_cell_input)){
	}
	
	auto symmetry_string() const{

		assert(geo_.num_atoms() > 0);
		
		char symbol[11];
		
		std::vector<int> types(geo_.num_atoms());
		std::vector<double> positions(3*geo_.num_atoms());
		
		for(int iatom = 0; iatom < geo_.num_atoms(); iatom++){
			types[iatom] = geo_.atoms()[iatom].atomic_number();
			auto pos = cell_.metric().to_contravariant(cell_.position_in_cell(geo_.coordinates()[iatom]));
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
		
		auto symnum = spg_get_international(symbol, reinterpret_cast<double (*)[3]>(amat), reinterpret_cast<double (*)[3]>(positions.data()), types.data(), geo_.num_atoms(), 1e-4);
		return symbol + std::string(" (number ") + std::to_string(symnum) + std::string(")");
	}

	auto & geo() const {
		return geo_;
	}
	auto & geo() {
		return geo_;
	}

	auto & cell() const {
		return cell_;
	}

	auto insert(std::string const & symbol, vector3<quantity<magnitude::length>> const & pos){
		geo_.add_atom(input::species(pseudo::element(symbol)), pos);
	}
	
	auto insert(input::species const & sp, vector3<quantity<magnitude::length>> const & pos){
		geo_.add_atom(sp, pos);
	}

	template <class ContainerType>
	auto insert(ContainerType const & container){
		for(auto atom : container) geo_.add_atom(atom.species(), atom.position());
	}

	auto insert_fractional(input::species const & sp, vector3<double, contravariant> const & pos){
		geo_.add_atom(sp, cell_.metric().to_cartesian(pos));
	}

	auto insert_fractional(std::string const & symbol, vector3<double, contravariant> const & pos){
		insert_fractional(input::species(pseudo::element(symbol)), pos);
	}
	
	inq::ions::unit_cell cell_;
	inq::ions::geometry geo_;
	
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
}
#endif
