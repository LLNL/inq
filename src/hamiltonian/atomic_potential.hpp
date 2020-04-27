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
#include <utils/partition.hpp>

#include <unordered_map>

#include <mpi3/environment.hpp>

namespace hamiltonian {

  class atomic_potential {

  public:

    enum class error {
      PSEUDOPOTENTIAL_NOT_FOUND
    };

    template <class atom_array>
    atomic_potential(const int natoms, const atom_array & atom_list, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			sep_(0.625), //this is the default from octopus
      pseudo_set_("pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"),
			comm_(comm),
			part_(natoms, comm_)
		{

			has_nlcc_ = false;
      nelectrons_ = 0.0;
      for(int iatom = 0; iatom < natoms; iatom++){
				if(!pseudo_set_.has(atom_list[iatom])) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 

				auto file_path = pseudo_set_.file_path(atom_list[iatom]);
				if(atom_list[iatom].has_file()) file_path = atom_list[iatom].file_path();

				auto insert = pseudopotential_list_.emplace(atom_list[iatom].symbol(), pseudo::pseudopotential(file_path, sep_));
				
				auto & pseudo = insert.first->second;
				
				nelectrons_ += pseudo.valence_charge();
				has_nlcc_ = has_nlcc_ or pseudo.has_nlcc_density();
				
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

			for(long ii = 0; ii < potential.basis().size(); ii++) potential.linear()[ii] = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, ps.short_range_potential_radius());

				//DATAOPERATIONS LOOP + GPU::RUN 1D (random access output)
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					auto rr = sphere.distance()[ipoint];
					auto sr_potential = ps.short_range_potential().value(rr);
					potential.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]] += sr_potential;
				}
				
      }

			if(part_.parallel()){
				comm_.all_reduce_in_place_n(static_cast<double *>(potential.linear().data()), potential.linear().size(), std::plus<>{});
			}

			return potential;			
    }
		
    template <class basis_type, class cell_type, class geo_type>
    basis::field<basis_type, double> ionic_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

      basis::field<basis_type, double> density(basis);
			
			density = 0.0;

			for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, sep_.long_range_density_radius());

				//DATAOPERATIONS LOOP + GPU::RUN 1D (random access output)
#ifdef HAVE_CUDA
				//OPTIMIZATION: this should be done in parallel for atoms too
				gpu::run(sphere.size(),
								 [dns = begin(density.cubic()), pts = begin(sphere.points()),
									chrg = ps.valence_charge(),
									sp = sep_, distance = begin(sphere.distance())] __device__
								 (auto ipoint){
									 double rr = distance[ipoint];
									 dns[pts[ipoint][0]][pts[ipoint][1]][pts[ipoint][2]] += chrg*sp.long_range_density(rr);
								 });
#else
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					auto rr = sphere.distance()[ipoint];
					density.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]]
						+= ps.valence_charge()*sep_.long_range_density(rr);
				}
#endif
      }

			if(part_.parallel()){
				comm_.all_reduce_in_place_n(static_cast<double *>(density.linear().data()), density.linear().size(), std::plus<>{});
			}
			
			return density;			
    }
    
    template <class basis_type, class cell_type, class geo_type>
    auto atomic_electronic_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

      basis::field<basis_type, double> density(basis);

			for(long ii = 0; ii < density.basis().size(); ii++) density.linear()[ii] = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);

				//TODO: implement the case when the pseudo does not have the density
				assert(ps.has_electronic_density());

				basis::spherical_grid sphere(basis, cell, atom_position, ps.electronic_density_radius());

				//DATAOPERATIONS LOOP + GPU::RUN 1D (random access output)
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					auto rr = sphere.distance()[ipoint];
					auto density_val = ps.electronic_density().value(rr);
					density.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]] += density_val;
				}
				
      }

			if(part_.parallel()){
				comm_.all_reduce_in_place_n(static_cast<double *>(density.linear().data()), density.linear().size(), std::plus<>{});
			}

			return density;			
    }

		auto has_nlcc() const {
			return has_nlcc_;
		}

		template <class basis_type, class cell_type, class geo_type>
    auto nlcc_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

      basis::field<basis_type, double> density(basis);

			for(long ii = 0; ii < density.basis().size(); ii++) density.linear()[ii] = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);

				if(not ps.has_nlcc_density()) continue;

				assert(has_nlcc());
				
				basis::spherical_grid sphere(basis, cell, atom_position, ps.nlcc_density_radius());

				//DATAOPERATIONS LOOP + GPU::RUN 1D (random access output)
				for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
					auto rr = sphere.distance()[ipoint];
					auto density_val = ps.nlcc_density().value(rr);
					density.cubic()[sphere.points()[ipoint][0]][sphere.points()[ipoint][1]][sphere.points()[ipoint][2]] += density_val;
				}
				
      }

			if(part_.parallel()){
				comm_.all_reduce_in_place_n(static_cast<double *>(density.linear().data()), density.linear().size(), std::plus<>{});
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
		mutable boost::mpi3::communicator comm_;
		utils::partition part_;
		bool has_nlcc_;
        
  };

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/geometry.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::atomic_potential", "[hamiltonian::atomic_potential]") {

  using namespace Catch::literals;
	using pseudo::element;
  using input::species;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	SECTION("Non-existing element"){
    std::vector<species> el_list({element("P"), element("X")});

    CHECK_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list));
  }
  
  SECTION("Duplicated element"){
    std::vector<species> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.size(), el_list.begin());

    CHECK(pot.num_species() == 1);
    CHECK(pot.num_electrons() == 10.0_a);
    
  }

  SECTION("Empty list"){
    std::vector<species> el_list;
    
    hamiltonian::atomic_potential pot(el_list.size(), el_list);

    CHECK(pot.num_species() == 0);
    CHECK(pot.num_electrons() == 0.0_a);

		CHECK(not pot.has_nlcc());
		
  }

  SECTION("CNOH"){
    species el_list[] = {element("C"), element("N"), element("O"), element("H")};

    hamiltonian::atomic_potential pot(4, el_list);

    CHECK(pot.num_species() == 4);
    CHECK(pot.num_electrons() == 16.0_a);
  }

  SECTION("Construct from a geometry"){

    ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");

    hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms(), comm);

    CHECK(pot.num_species() == 2);
    CHECK(pot.num_electrons() == 30.0_a);

		double ll = 20.0;
		auto cell = input::cell::cubic(ll, ll, ll);
		basis::real_space rs(cell, input::basis::cutoff_energy(20.0));

		rs.info(std::cout);
		
		auto vv = pot.local_potential(rs, cell, geo);

		CHECK(operations::integral(vv) == -45.5544154295_a);

		CHECK(vv.cubic()[5][3][0] == -1.574376555_a);
		CHECK(vv.cubic()[3][1][0] == -0.258229883_a);
							 
		auto id = pot.ionic_density(rs, cell, geo);

		CHECK(operations::integral(id) == -30.0000000746_a);
		CHECK(id.cubic()[5][3][0] == -0.9448936487_a);
		CHECK(id.cubic()[3][1][0] == -0.2074502252_a);

		auto nn = pot.atomic_electronic_density(rs, cell, geo);
		
		CHECK(operations::integral(nn) == 29.9562520003_a);
		CHECK(nn.cubic()[5][3][0] == 0.1330589609_a);
		CHECK(nn.cubic()[3][1][0] == 0.1846004508_a);

		CHECK(pot.has_nlcc());
		
		auto nlcc = pot.nlcc_density(rs, cell, geo);
		
		CHECK(operations::integral(nlcc) == 3.0083012065_a);
		CHECK(nlcc.cubic()[5][3][0] == 0.6248217151_a);
		CHECK(nlcc.cubic()[3][1][0] == 0.0007040027_a);
		
  }
  
}

#endif
  
#endif
