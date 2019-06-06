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

    template <class species_list_it>
    atomic_potential(const species_list_it species_begin, const species_list_it species_end):
      pseudo_set_(config::path::share() + "pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"){

      for(auto it = species_begin; it != species_end; it++){
	if(!pseudo_set_.has(*it)) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 

	pseudopotential_list_.emplace(it->symbol(), pseudo::pseudopotential(pseudo_set_.file_path(*it)));
      }
      
    }

    int number_of_species() const {
      return pseudopotential_list_.size();
    }
    
  private:

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

    REQUIRE_THROWS(hamiltonian::atomic_potential(el_list.begin(), el_list.end()));
  }
  
  SECTION("Duplicated element"){
    std::vector<element> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.begin(), el_list.end());

    REQUIRE(pot.number_of_species() == 1);
    
  }

  SECTION("Empty list"){
    std::vector<element> el_list;
    
    hamiltonian::atomic_potential pot(el_list.begin(), el_list.end());

    REQUIRE(pot.number_of_species() == 0);
  }

  SECTION("CNOH"){
    std::vector<element> el_list({element("C"), element("N"), element("O"), element("H")});

    hamiltonian::atomic_potential pot(el_list.begin(), el_list.end());

    REQUIRE(pot.number_of_species() == 4);
  }
  
}

#endif
  
#endif
