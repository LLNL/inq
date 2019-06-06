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

#include <config/path.hpp>
#include <math/d3vector.hpp>
#include <pseudo/element.hpp>

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

    // Generates a geometry from an xyz file
    geometry(const std::string & xyz_file_name){
      std::ifstream xyz_file(xyz_file_name.c_str());

      if(!xyz_file.is_open()) throw error::FILE_NOT_FOUND;

      int num_atoms;
      std::string comment_line;
      
      xyz_file >> num_atoms;
      
      std::getline(xyz_file, comment_line);
      std::getline(xyz_file, comment_line);
      
      std::string atom_name;
      math::d3vector atom_position;
      
      for(int iatom = 0; iatom < num_atoms; iatom++){
      xyz_file >> atom_name >> atom_position;
      add_atom(pseudo::element(atom_name), atom_position*1.8897261);
      }
      
      xyz_file.close();
      
      assert(num_atoms == number_of_atoms());
      
    }
    
    int number_of_atoms() const { return coordinates_.size(); }

    void add_atom(const pseudo::element & element, const math::d3vector & position){
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
  private:

    std::vector<pseudo::element> atoms_;
    std::vector<math::d3vector> coordinates_;
    
  };
  
#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class ions::geometry", "[geometry]") {

  using namespace Catch::literals;

  SECTION("Create empty and add an atom"){
    ions::geometry geo;

    REQUIRE(geo.number_of_atoms() == 0);

    geo.add_atom(pseudo::element("Xe"), math::d3vector(1000.0, -200.0, 6.0));

    REQUIRE(geo.number_of_atoms() == 1);
    REQUIRE(geo.atoms()[0].atomic_number() == 54);
    REQUIRE(geo.atoms()[0] == pseudo::element(54));
    REQUIRE(geo.atoms()[0].charge() == -54.0_a);
    REQUIRE(geo.atoms()[0].mass() == 131.2936_a);
    REQUIRE(geo.coordinates()[0][0] == 1000.0_a);
    REQUIRE(geo.coordinates()[0][1] == -200.0_a);
    REQUIRE(geo.coordinates()[0][2] == 6.0_a);

    geo.coordinates()[0][0] += 8;  
    
    REQUIRE(geo.coordinates()[0][0] == 1008.0_a);
  }    
 
  SECTION("Read an xyz file"){
    ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");

    REQUIRE(geo.number_of_atoms() == 12);
    
    REQUIRE(geo.atoms()[2] == pseudo::element("C"));
    REQUIRE(geo.atoms()[2].charge() == -6.0_a);
    REQUIRE(geo.atoms()[2].mass() == 12.0096_a);
    REQUIRE(geo.coordinates()[2][0] == 2.2846788549_a);
    REQUIRE(geo.coordinates()[2][1] == -1.3190288178_a);
    REQUIRE(geo.coordinates()[2][2] == 0.0_a);

    REQUIRE(geo.atoms()[11] == pseudo::element("H"));
    REQUIRE(geo.atoms()[11].charge() == -1.0_a);
    REQUIRE(geo.atoms()[11].mass() == 1.00784_a);
    REQUIRE(geo.coordinates()[11][0] == -4.0572419367_a);
    REQUIRE(geo.coordinates()[11][1] == 2.343260364_a);
    REQUIRE(geo.coordinates()[11][2] == 0.0_a);

    geo.add_atom(pseudo::element("Cl"), math::d3vector(-3.0, 4.0, 5.0));

    REQUIRE(geo.number_of_atoms() == 13);
    REQUIRE(geo.atoms()[12].atomic_number() == 17);
    REQUIRE(geo.atoms()[12] == pseudo::element(17));
    REQUIRE(geo.atoms()[12].charge() == -17.0_a);
    REQUIRE(geo.atoms()[12].mass() == 35.446_a);
    REQUIRE(geo.coordinates()[12][0] == -3.0_a);
    REQUIRE(geo.coordinates()[12][1] == 4.0_a);
    REQUIRE(geo.coordinates()[12][2] == 5.0_a);
    
  }

  SECTION("Try to read a non-existent file"){
    REQUIRE_THROWS(ions::geometry("/this_file_should_not_exist,_i_hope_it_doesnt"));
  }
  

  
}
#endif

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
