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
      FILE_NOT_FOUND
    };
      
    geometry(){
    }

    // Generates a geometry from an xyz file
    geometry(const std::string & xyz_file_name){
      std::ifstream xyz_file(xyz_file_name.c_str());

      if(!xyz_file.is_open()){
	throw error::FILE_NOT_FOUND;
      }
      
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
      return Atom(atom_index, this);
    }
    
  private:

    typedef std::map<int, pseudo::element> ElementsContainer;
    
    ElementsContainer element_list;
    std::vector<ElementsContainer::const_iterator> elements;
    std::vector<Position> coordinates;
    
  };
  
#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class ions::geometry", "[geometry]") {

  using Catch::Matchers::WithinULP;

  SECTION("Try to read a non-existent file"){
    REQUIRE_THROWS(ions::geometry("/this_file_should_not_exist,_i_hope_it_doesnt"));
  }
  
  SECTION("Read an xyz file"){
    ions::geometry geo(SHARE_DIR + std::string("unit_test_data/benzene.xyz"));
  }
  
}
#endif

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
