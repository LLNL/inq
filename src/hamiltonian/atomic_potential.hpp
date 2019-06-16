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

#include <pseudo/set.hpp>
#include <pseudo/pseudopotential.hpp>
#include <basis/spherical_grid.hpp>
#include <multi/array.hpp>
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
      pseudo_set_(config::path::share() + "pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"){

      nelectrons_ = 0.0;
      for(int iatom = 0; iatom < natoms; iatom++){
	if(!pseudo_set_.has(atom_list[iatom])) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 

	auto insert = pseudopotential_list_.emplace(atom_list[iatom].symbol(), pseudo::pseudopotential(pseudo_set_.file_path(atom_list[iatom])));

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

    const pseudo::pseudopotential & pseudo_for_element(const pseudo::element & el) const {
      return pseudopotential_list_.at(el.symbol());
    }

    template <class basis_type, class cell_type, class geo_type, class array_dim3>
    void local_potential(const basis_type & basis, const cell_type & cell, const geo_type & geo, array_dim3 & potential) const {

      array_dim3 density(extensions(potential), 0.0);

      for(int iatom = 0; iatom < geo.num_atoms(); iatom++){

	auto atom_position = geo.coordinates()[iatom];
	
	auto & ps = pseudo_for_element(geo.atoms()[iatom]);
	basis::spherical_grid sphere(basis, cell, atom_position, ps.long_range_density_radius());

	for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
	  double rr = length(basis.rvector(sphere.points()[ipoint]) - atom_position);
	  density[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]] += ps.long_range_density(rr);
	}

      }

      solvers::poisson psolver;

      psolver.solve(basis, density, potential);
      
    }
    
    template <class output_stream>
    void info(output_stream & out) const {
      out << "ATOMIC POTENTIAL:" << std::endl;
      out << "  Number of species   = " << num_species() << std::endl;
      out << "  Number of electrons = " << num_electrons() << std::endl;
      out << std::endl;
    }    
    
  private:

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

  SECTION("Non-existing element"){
    std::vector<element> el_list({element("P"), element("X")});

    REQUIRE_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list));
  }
  
  SECTION("Duplicated element"){
    std::vector<element> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.size(), el_list.begin());

    REQUIRE(pot.num_species() == 1);
    REQUIRE(pot.num_electrons() == 10.0_a);
    
  }

  SECTION("Empty list"){
    std::vector<element> el_list;
    
    hamiltonian::atomic_potential pot(el_list.size(), el_list);

    REQUIRE(pot.num_species() == 0);
    REQUIRE(pot.num_electrons() == 0.0_a);
  }

  SECTION("CNOH"){
    element el_list[] = {element("C"), element("N"), element("O"), element("H")};

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
