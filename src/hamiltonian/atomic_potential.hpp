#ifndef ATOMIC_POTENTIAL_HPP
#define ATOMIC_POTENTIAL_HPP

#include <pseudo/set.hpp>
#include <pseudo/pseudopotential.hpp>

#include <unordered_map>

namespace hamiltonian {

  class atomic_potential {

  public:

    enum class error {
      PSEUDOPOTENTIAL_NOT_FOUND
    };

    template <class atom_array>
    atomic_potential(const int natoms, const atom_array & atom_list):
      pseudo_set_(config::path::share() + "pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"){

      nelectrons_ = 0.0;
      for(int iatom = 0; iatom < natoms; iatom++){
	if(!pseudo_set_.has(atom_list[iatom])) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 

	auto insert = pseudopotential_list_.emplace(atom_list[iatom].symbol(), pseudo::pseudopotential(pseudo_set_.file_path(atom_list[iatom])));

	auto & pseudo = insert.first->second;
	
	nelectrons_ += pseudo.valence_charge();
	
      }
      
    }

    int number_of_species() const {
      return pseudopotential_list_.size();
    }

    const double & number_of_electrons() const {
      return nelectrons_;
    }
    
  private:

    double nelectrons_;
    pseudo::set pseudo_set_;
    std::unordered_map<std::string, pseudo::pseudopotential> pseudopotential_list_;
        
  };

}

#ifdef UNIT_TEST
#include <catch.hpp>

TEST_CASE("Class hamiltonian::atomic_potential", "[atomic_potential]") {

  using namespace Catch::literals;
  using pseudo::element;

  SECTION("Non-existing element"){
    std::vector<element> el_list({element("P"), element("X")});

    REQUIRE_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list));
  }
  
  SECTION("Duplicated element"){
    std::vector<element> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.size(), el_list.begin());

    REQUIRE(pot.number_of_species() == 1);
    REQUIRE(pot.number_of_electrons() == 10.0_a);
    
  }

  SECTION("Empty list"){
    std::vector<element> el_list;
    
    hamiltonian::atomic_potential pot(el_list.size(), el_list);

    REQUIRE(pot.number_of_species() == 0);
    REQUIRE(pot.number_of_electrons() == 0.0_a);
  }

  SECTION("CNOH"){
    element el_list[] = {element("C"), element("N"), element("O"), element("H")};

    hamiltonian::atomic_potential pot(4, el_list);

    REQUIRE(pot.number_of_species() == 4);
    REQUIRE(pot.number_of_electrons() == 16.0_a);
  }
  
}

#endif
  
#endif
