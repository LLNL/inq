/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__GEOMETRY
#define INQ__IONS__GEOMETRY

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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

#include <math/vector3.hpp>
#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <input/atom.hpp>
#include <config/path.hpp>

#include <vector>
#include <cassert>

namespace inq {
namespace ions {

class geometry {

public:
    
	enum class error {
										FILE_NOT_FOUND
	};
      
	geometry(){
	}

	template <class container_type>
	geometry(const container_type & atom_container){

		atoms_.reserve(atom_container.size());
		coordinates_.reserve(atom_container.size());
			
		for(auto it = atom_container.begin(); it != atom_container.end(); it++){
			atoms_.push_back(it->species());
			coordinates_.push_back(it->position());
			velocities_.push_back(math::vector3<double>(0.0, 0.0, 0.0));			
		}
			
	}
    
	int num_atoms() const {
		return (long) coordinates_.size();
	}

	void add_atom(const input::species & element, const math::vector3<double> & position){
		atoms_.push_back(element);
		coordinates_.push_back(position);
		velocities_.push_back(math::vector3<double>(0.0, 0.0, 0.0));					
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
	
	template <class output_stream>
	void info(output_stream & out) const {
		out << "GEOMETRY:" << std::endl;
		out << "  Number of atoms = " << num_atoms() << std::endl;
		out << std::endl;
	}
	
	template<class OStream>
	friend OStream& operator<<(OStream& os, geometry const& self){
		self.info(os);
		return os;
	}

private:

	std::vector<input::species> atoms_;
	std::vector<math::vector3<double>> coordinates_;
	std::vector<math::vector3<double>> velocities_;	
    
};
}
}
  
#ifdef INQ_IONS_GEOMETRY_UNIT_TEST
#undef INQ_IONS_GEOMETRY_UNIT_TEST

#include <catch2/catch.hpp>

#include <input/parse_xyz.hpp>

TEST_CASE("Class ions::geometry", "[geometry]") {

	using namespace inq;
	using namespace Catch::literals;

  SECTION("Create empty and add an atom"){
    ions::geometry geo;

    CHECK(geo.num_atoms() == 0);

    geo.add_atom(pseudo::element("Xe"), math::vector3<double>(1000.0, -200.0, 6.0));

    CHECK(geo.num_atoms() == 1);
    CHECK(geo.atoms()[0].atomic_number() == 54);
    CHECK(geo.atoms()[0] == pseudo::element(54));
    CHECK(geo.atoms()[0].charge() == -54.0_a);
    CHECK(geo.atoms()[0].mass() == 239333.5935636_a);
    CHECK(geo.coordinates()[0][0] == 1000.0_a);
    CHECK(geo.coordinates()[0][1] == -200.0_a);
    CHECK(geo.coordinates()[0][2] == 6.0_a);
		CHECK(geo.velocities()[0][0] == 0.0_a);
    CHECK(geo.velocities()[0][1] == 0.0_a);
    CHECK(geo.velocities()[0][2] == 0.0_a);
		
    geo.coordinates()[0][0] += 8;  
    
    CHECK(geo.coordinates()[0][0] == 1008.0_a);

		assert(geo.velocities().size() == geo.coordinates().size());
				
  }    
 
  SECTION("Read an xyz file"){

    ions::geometry geo(input::parse_xyz(config::path::unit_tests_data() + "benzene.xyz"));

    CHECK(geo.num_atoms() == 12);
    
    CHECK(geo.atoms()[2] == pseudo::element("C"));
    CHECK(geo.atoms()[2].charge() == -6.0_a);
    CHECK(geo.atoms()[2].mass() == 21892.1617296_a);
    CHECK(geo.coordinates()[2][0] == 2.2846788549_a);
    CHECK(geo.coordinates()[2][1] == -1.3190288178_a);
    CHECK(geo.coordinates()[2][2] == 0.0_a);

    CHECK(geo.atoms()[11] == pseudo::element("H"));
    CHECK(geo.atoms()[11].charge() == -1.0_a);
    CHECK(geo.atoms()[11].mass() == 1837.17994584_a);
    CHECK(geo.coordinates()[11][0] == -4.0572419367_a);
    CHECK(geo.coordinates()[11][1] == 2.343260364_a);
    CHECK(geo.coordinates()[11][2] == 0.0_a);
		CHECK(geo.velocities()[11][0] == 0.0_a);
    CHECK(geo.velocities()[11][1] == 0.0_a);
    CHECK(geo.velocities()[11][2] == 0.0_a);

		assert(geo.velocities().size() == geo.coordinates().size());
		
    geo.add_atom(pseudo::element("Cl"), math::vector3<double>(-3.0, 4.0, 5.0));

    CHECK(geo.num_atoms() == 13);
    CHECK(geo.atoms()[12].atomic_number() == 17);
    CHECK(geo.atoms()[12] == pseudo::element(17));
    CHECK(geo.atoms()[12].charge() == -17.0_a);
    CHECK(geo.atoms()[12].mass() == 64614.105771_a);
    CHECK(geo.coordinates()[12][0] == -3.0_a);
    CHECK(geo.coordinates()[12][1] == 4.0_a);
    CHECK(geo.coordinates()[12][2] == 5.0_a);
		CHECK(geo.velocities()[12][0] == 0.0_a);
    CHECK(geo.velocities()[12][1] == 0.0_a);
    CHECK(geo.velocities()[12][2] == 0.0_a);

		assert(geo.velocities().size() == geo.coordinates().size());
		
  }

}
#endif

#endif
