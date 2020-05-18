/* -*- indent-tabs-mode: t -*- */

#ifndef IONS_GEOMETRY
#define IONS_GEOMETRY

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

#include <math/vec3d.hpp>
#include <pseudopod/element.hpp>
#include <input/species.hpp>
#include <input/atom.hpp>

#include <vector>
#include <cassert>

namespace ions {

  class geometry {

  public:
    
    enum class error {
      FILE_NOT_FOUND
    };
      
    geometry(){
    }

		geometry(const char * xyz_file_name){
			geometry(std::string(xyz_file_name));
		}
		
    // Generates a geometry from an xyz file
    geometry(const std::string & xyz_file_name){
      std::ifstream xyz_file(xyz_file_name.c_str());

      if(!xyz_file.is_open()) throw error::FILE_NOT_FOUND;

      int natoms;
      std::string comment_line;
      
      xyz_file >> natoms;
      
      std::getline(xyz_file, comment_line);
      std::getline(xyz_file, comment_line);
      
      std::string atom_name;
      math::vec3d atom_position;
      
      for(int iatom = 0; iatom < natoms; iatom++){
      xyz_file >> atom_name >> atom_position;
      add_atom(pseudo::element(atom_name), atom_position*1.8897261);
      }
      
      xyz_file.close();
      
      assert(natoms == num_atoms());
      
    }

		template <class container_type>
		geometry(const container_type & atom_container){

			atoms_.reserve(atom_container.size());
			coordinates_.reserve(atom_container.size());
			
			for(auto it = atom_container.begin(); it != atom_container.end(); it++){
				atoms_.push_back(it->species());
				coordinates_.push_back(it->position());
			}
			
    }
    
    int num_atoms() const { return coordinates_.size(); }

    void add_atom(const input::species & element, const math::vec3d & position){
      atoms_.push_back(element);
      coordinates_.push_back(position);
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

    template <class output_stream>
    void info(output_stream & out) const {
      out << "GEOMETRY:" << std::endl;
      out << "  Number of atoms = " << num_atoms() << std::endl;
      out << std::endl;
    }
    
  private:

    std::vector<input::species> atoms_;
    std::vector<math::vec3d> coordinates_;
    
  };
  
#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class ions::geometry", "[geometry]") {

  using namespace Catch::literals;

  SECTION("Create empty and add an atom"){
    ions::geometry geo;

    CHECK(geo.num_atoms() == 0);

    geo.add_atom(pseudo::element("Xe"), math::vec3d(1000.0, -200.0, 6.0));

    CHECK(geo.num_atoms() == 1);
    CHECK(geo.atoms()[0].atomic_number() == 54);
    CHECK(geo.atoms()[0] == pseudo::element(54));
    CHECK(geo.atoms()[0].charge() == -54.0_a);
    CHECK(geo.atoms()[0].mass() == 131.2936_a);
    CHECK(geo.coordinates()[0][0] == 1000.0_a);
    CHECK(geo.coordinates()[0][1] == -200.0_a);
    CHECK(geo.coordinates()[0][2] == 6.0_a);

    geo.coordinates()[0][0] += 8;  
    
    CHECK(geo.coordinates()[0][0] == 1008.0_a);
  }    
 
  SECTION("Read an xyz file"){
    ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");

    CHECK(geo.num_atoms() == 12);
    
    CHECK(geo.atoms()[2] == pseudo::element("C"));
    CHECK(geo.atoms()[2].charge() == -6.0_a);
    CHECK(geo.atoms()[2].mass() == 12.0096_a);
    CHECK(geo.coordinates()[2][0] == 2.2846788549_a);
    CHECK(geo.coordinates()[2][1] == -1.3190288178_a);
    CHECK(geo.coordinates()[2][2] == 0.0_a);

    CHECK(geo.atoms()[11] == pseudo::element("H"));
    CHECK(geo.atoms()[11].charge() == -1.0_a);
    CHECK(geo.atoms()[11].mass() == 1.00784_a);
    CHECK(geo.coordinates()[11][0] == -4.0572419367_a);
    CHECK(geo.coordinates()[11][1] == 2.343260364_a);
    CHECK(geo.coordinates()[11][2] == 0.0_a);

    geo.add_atom(pseudo::element("Cl"), math::vec3d(-3.0, 4.0, 5.0));

    CHECK(geo.num_atoms() == 13);
    CHECK(geo.atoms()[12].atomic_number() == 17);
    CHECK(geo.atoms()[12] == pseudo::element(17));
    CHECK(geo.atoms()[12].charge() == -17.0_a);
    CHECK(geo.atoms()[12].mass() == 35.446_a);
    CHECK(geo.coordinates()[12][0] == -3.0_a);
    CHECK(geo.coordinates()[12][1] == 4.0_a);
    CHECK(geo.coordinates()[12][2] == 5.0_a);
    
  }

  SECTION("Try to read a non-existent file"){
    CHECK_THROWS(ions::geometry("/this_file_should_not_exist,_i_hope_it_doesnt"));
  }
  
}
#endif

}

#endif
