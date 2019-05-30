#ifndef ST_GEOMETRY
#define ST_GEOMETRY

#include "position.hpp"
#include <pseudo/element.hpp>

#include <map>
#include <vector>
#include <cassert>

namespace ions {

  class geometry {

  public:

    geometry(){
    }

    // Generates a geometry from an xyz file
    geometry(const std::string & xyz_file_name){
      std::ifstream xyz_file(xyz_file_name.c_str());
      
      int num_atoms;
      std::string comment_line;
      
      xyz_file >> num_atoms;
      
      std::getline(xyz_file, comment_line);
      std::getline(xyz_file, comment_line);
      
      std::string atom_name;
      Position atom_position;
      
      for(int iatom = 0; iatom < num_atoms; iatom++){
      xyz_file >> atom_name >> atom_position;
      add_atom(pseudopotential::element(atom_name), atom_position*1.8897261);
      }
      
      xyz_file.close();
      
      assert(num_atoms == number_of_atoms());
      
    }
    
    int number_of_atoms() const { return coordinates.size(); }

    void add_atom(const pseudopotential::element & element, const Position & position){

      auto element_location = element_list.insert(std::pair<int, pseudopotential::element>(element.atomic_number(), element));
      elements.push_back(element_location.first);
      coordinates.push_back(position);

    }

    class Atom {

    public:

      const pseudopotential::element & element() const {
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

    typedef std::map<int, pseudopotential::element> ElementsContainer;
    
    ElementsContainer element_list;
    std::vector<ElementsContainer::const_iterator> elements;
    std::vector<Position> coordinates;
    
  };
  
}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
