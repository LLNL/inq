#ifndef IONS_GEOMETRY
#define IONS_GEOMETRY

#include "position.hpp"
#include <pseudo/element.hpp>

#include <map>
#include <vector>
#include <cassert>

namespace ions {

  class geometry {

  public:
    
    enum class error {
      FILE_NOT_FOUND,
      INDEX_OUT_OF_RANGE
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
      Position atom_position;
      
      for(int iatom = 0; iatom < num_atoms; iatom++){
      xyz_file >> atom_name >> atom_position;
      add_atom(pseudo::element(atom_name), atom_position*1.8897261);
      }
      
      xyz_file.close();
      
      assert(num_atoms == number_of_atoms());
      
    }
    
    int number_of_atoms() const { return coordinates.size(); }

    void add_atom(const pseudo::element & element, const Position & position){

      auto element_location = element_list.insert(std::pair<int, pseudo::element>(element.atomic_number(), element));
      elements.push_back(element_location.first);
      coordinates.push_back(position);

    }

    class Atom {

    public:

      const pseudo::element & element() const {
	return geo->elements[index]->second;
      }
      
      double charge() const {
	return geo->elements[index]->second.charge();
      }

      double mass() const {
	return geo->elements[index]->second.mass();
      }

      const Position & position() const {
	return geo->coordinates[index];
      }
      
    private:
     
      Atom(const int atom_index, const geometry * atom_geo)
	:index(atom_index), geo(atom_geo){
      }
 
      const int index;
      const geometry * geo;
      
      friend class geometry;
    };
    

    Atom atom(const int atom_index) const {
      if(atom_index >= number_of_atoms()) throw error::INDEX_OUT_OF_RANGE;
      return Atom(atom_index, this);
    }
    
  private:

    typedef std::map<int, pseudo::element> ElementsContainer;
    
    ElementsContainer element_list;
    std::vector<ElementsContainer::const_iterator> elements;
    std::vector<Position> coordinates;
    
  };
  
#ifdef UNIT_TEST
#include <catch.hpp>

TEST_CASE("Class ions::geometry", "[geometry]") {

  using namespace Catch::literals;

  SECTION("Read an xyz file"){
    ions::geometry geo(SHARE_DIR + std::string("unit_tests_data/benzene.xyz"));

    REQUIRE(geo.number_of_atoms() == 12);
    REQUIRE_THROWS(geo.atom(12));
    REQUIRE_THROWS(geo.atom(425));
    
    REQUIRE(geo.atom(2).element() == pseudo::element("C"));
    REQUIRE(geo.atom(2).charge() == -6.0_a);
    REQUIRE(geo.atom(2).mass() == 12.0096_a);
    REQUIRE(geo.atom(2).position().x() == 2.2846788549_a);
    REQUIRE(geo.atom(2).position().y() == -1.3190288178_a);
    REQUIRE(geo.atom(2).position().z() == 0.0_a);

    REQUIRE(geo.atom(11).element() == pseudo::element("H"));
    REQUIRE(geo.atom(11).charge() == -1.0_a);
    REQUIRE(geo.atom(11).mass() == 1.00784_a);
    REQUIRE(geo.atom(11).position().x() == -4.0572419367_a);
    REQUIRE(geo.atom(11).position().y() == 2.343260364_a);
    REQUIRE(geo.atom(11).position().z() == 0.0_a);

    geo.add_atom(pseudo::element("Cl"), ions::Position(-3.0, 4.0, 5.0));

    REQUIRE(geo.number_of_atoms() == 13);
    REQUIRE(geo.atom(12).element().atomic_number() == 17);
    REQUIRE(geo.atom(12).element() == pseudo::element(17));
    REQUIRE(geo.atom(12).charge() == -17.0_a);
    REQUIRE(geo.atom(12).mass() == 35.446_a);
    REQUIRE(geo.atom(12).position().x() == -3.0_a);
    REQUIRE(geo.atom(12).position().y() == 4.0_a);
    REQUIRE(geo.atom(12).position().z() == 5.0_a);
    
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
