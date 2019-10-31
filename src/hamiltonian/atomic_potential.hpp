/* -*- indent-tabs-mode: t -*- */

#ifndef ATOMIC_POTENTIAL_HPP
#define ATOMIC_POTENTIAL_HPP

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

#include <pseudopod/set.hpp>
#include <pseudopod/pseudopotential.hpp>
#include <basis/spherical_grid.hpp>
#include <math/array.hpp>
#include <solvers/poisson.hpp>

#include <unordered_map>

namespace hamiltonian {

  class atomic_potential {

  public:

    enum class error {
      PSEUDOPOTENTIAL_NOT_FOUND
    };

    template <class atom_array>
    atomic_potential(const int natoms, const atom_array & atom_list):
			sep_(0.5), //this is the default from qball, but it can be optimized for the grid. Check AtomsSet.cc:1102
      pseudo_set_("pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"){

      nelectrons_ = 0.0;
      for(int iatom = 0; iatom < natoms; iatom++){
				if(!pseudo_set_.has(atom_list[iatom])) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 

				auto file_path = pseudo_set_.file_path(atom_list[iatom]);
				if(atom_list[iatom].has_file()) file_path = atom_list[iatom].file_path();

				auto insert = pseudopotential_list_.emplace(atom_list[iatom].symbol(), pseudo::pseudopotential(file_path, sep_));
				
				auto & pseudo = insert.first->second;
				
				nelectrons_ += pseudo.valence_charge();
				
      }
      
    }
		
    int num_species() const {
      return pseudopotential_list_.size();
    }

    const double & num_electrons() const {
      return nelectrons_;
    }
		
		template <class element_type>
    const pseudo::pseudopotential & pseudo_for_element(const element_type & el) const {
      return pseudopotential_list_.at(el.symbol());
    }

    template <class basis_type, class cell_type, class geo_type>
    auto local_potential(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

      basis::field<basis_type, double> potential(basis);

			for(long ii = 0; ii < potential.basis().size(); ii++) potential[ii] = 0.0;
			
      for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, ps.short_range_potential_radius());
				
				//DATAOPERATIONS LOOP 1D (random access output)
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					auto rr = length(basis.rvector(sphere.points()[ipoint]) - atom_position);
					auto sr_potential = ps.short_range_potential().value(rr);
					potential.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]] += sr_potential;
				}
				
      }

			return potential;			
    }

		
    template <class basis_type, class cell_type, class geo_type>
    auto ionic_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

      basis::field<basis_type, double> density(basis);
			
			density = 0.0;
			
			for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, sep_.long_range_density_radius());

				//DATAOPERATIONS LOOP 1D (random access output)
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					double rr = length(basis.rvector(sphere.points()[ipoint]) - atom_position);
					density.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]]
						+= ps.valence_charge()*sep_.long_range_density(rr);
				}
      }

			return density;			
    }
    
    template <class output_stream>
    void info(output_stream & out) const {
      out << "ATOMIC POTENTIAL:" << std::endl;
      out << "  Number of species   = " << num_species() << std::endl;
      out << "  Number of electrons = " << num_electrons() << std::endl;
      out << std::endl;
    }

		auto & range_separation() const {
			return sep_;
		}
    
  private:

		const math::erf_range_separation sep_;
    double nelectrons_;
    pseudo::set pseudo_set_;
    std::unordered_map<std::string, pseudo::pseudopotential> pseudopotential_list_;
        
  };

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/geometry.hpp>

TEST_CASE("Class hamiltonian::atomic_potential", "[atomic_potential]") {

  using namespace Catch::literals;
	using pseudo::element;
  using input::species;

  SECTION("Non-existing element"){
    std::vector<species> el_list({element("P"), element("X")});

    REQUIRE_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list));
  }
  
  SECTION("Duplicated element"){
    std::vector<species> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.size(), el_list.begin());

    REQUIRE(pot.num_species() == 1);
    REQUIRE(pot.num_electrons() == 10.0_a);
    
  }

  SECTION("Empty list"){
    std::vector<species> el_list;
    
    hamiltonian::atomic_potential pot(el_list.size(), el_list);

    REQUIRE(pot.num_species() == 0);
    REQUIRE(pot.num_electrons() == 0.0_a);
  }

  SECTION("CNOH"){
    species el_list[] = {element("C"), element("N"), element("O"), element("H")};

    hamiltonian::atomic_potential pot(4, el_list);

    REQUIRE(pot.num_species() == 4);
    REQUIRE(pot.num_electrons() == 16.0_a);
  }

  SECTION("Construct from a geometry"){
    
    ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");

    hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());

    REQUIRE(pot.num_species() == 2);
    REQUIRE(pot.num_electrons() == 30.0_a);
  }
  
}

#endif
  
#endif
