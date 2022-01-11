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

#include <utils/raw_pointer_cast.hpp>
#include <config/path.hpp>

#include <pseudopod/set.hpp>
#include <pseudopod/pseudopotential.hpp>


#include <basis/spherical_grid.hpp>
#include <basis/double_grid.hpp>
#include <math/array.hpp>
#include <operations/integral.hpp>
#include <solvers/poisson.hpp>
#include <utils/partition.hpp>

#include <unordered_map>

#include <mpi3/environment.hpp>

#include <utils/profiling.hpp>


namespace inq {
namespace hamiltonian {

  class atomic_potential {

  public:

		using pseudopotential_type = pseudo::pseudopotential<math::array<double, 1>>;
		
    enum class error {
      PSEUDOPOTENTIAL_NOT_FOUND
    };

    template <class atom_array>
    atomic_potential(const int natoms, const atom_array & atom_list, double gcutoff, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			sep_(0.625), //this is the default from octopus
      pseudo_set_("pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"),
			comm_(comm),
			part_(natoms, comm_)
		{

			CALI_CXX_MARK_FUNCTION;
			
			has_nlcc_ = false;
      nelectrons_ = 0.0;

      for(int iatom = 0; iatom < natoms; iatom++){
				if(!pseudo_set_.has(atom_list[iatom])) throw error::PSEUDOPOTENTIAL_NOT_FOUND; 
				
				auto map_ref = pseudopotential_list_.find(atom_list[iatom].symbol());
				
				if(map_ref == pseudopotential_list_.end()){
					
					auto file_path = pseudo_set_.file_path(atom_list[iatom]);
					if(atom_list[iatom].has_file()) file_path = atom_list[iatom].file_path();

					//sorry for this, emplace has a super ugly syntax
					auto insert = pseudopotential_list_.emplace(std::piecewise_construct, std::make_tuple(atom_list[iatom].symbol()), std::make_tuple(file_path, sep_, gcutoff, atom_list[iatom].filter_pseudo()));
					map_ref = insert.first;
					
				}
				
				auto & pseudo = map_ref->second;
				
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
    const pseudopotential_type & pseudo_for_element(const element_type & el) const {
      return pseudopotential_list_.at(el.symbol());
    }

    template <class basis_type, class cell_type, class geo_type>
    basis::field<basis_type, double> local_potential(const basis_type & basis, const cell_type & cell, const geo_type & geo, int single_atom = -1) const {

			CALI_CXX_MARK_SCOPE("atomic_potential::local_potential");
			
      basis::field<basis_type, double> potential(basis);

			potential = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){

				if(single_atom >= 0 and single_atom != iatom) continue;
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, ps.short_range_potential_radius());

				if(not basis.double_grid().enabled()){

					gpu::run(sphere.size(),
									 [pot = begin(potential.cubic()),
										sph = sphere.ref(),
										spline = ps.short_range_potential().cbegin()] GPU_LAMBDA (auto ipoint){
										 auto rr = sph.distance(ipoint);
										 auto potential_val = spline.value(rr);
										 gpu::atomic::add(&pot[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]], potential_val);
									 });

				} else {

					CALI_CXX_MARK_SCOPE("atomic_potential::double_grid");
					
					gpu::run(sphere.size(),
									 [pot = begin(potential.cubic()),
										sph = sphere.ref(),
										spline = ps.short_range_potential().cbegin(),
										dg = basis.double_grid().ref(),
										spac = basis.rspacing()] GPU_LAMBDA (auto ipoint){
										 gpu::atomic::add(&pot[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]],
																			dg.value([spline] GPU_LAMBDA (auto pos) { return spline.value(length(pos)); }, spac, sph.point_pos(ipoint)));
									 });
				
				}
			}
			
			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("atomic_potential::local_potential::reduce");
				comm_.all_reduce_in_place_n(raw_pointer_cast(potential.linear().data_elements()), potential.linear().size(), std::plus<>{});
			}
			
			return potential;			
		}
		
		template <class basis_type, class cell_type, class geo_type>
    basis::field<basis_type, double> ionic_density(const basis_type & basis, const cell_type & cell, const geo_type & geo, int single_atom = -1) const {

			CALI_CXX_MARK_FUNCTION;
	
      basis::field<basis_type, double> density(basis);
			
			density = 0.0;

			for(auto iatom = part_.start(); iatom < part_.end(); iatom++){

				if(single_atom >= 0 and single_atom != iatom) continue;
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);
				basis::spherical_grid sphere(basis, cell, atom_position, sep_.long_range_density_radius());

				//OPTIMIZATION: this should be done in parallel for atoms too
				gpu::run(sphere.size(),
								 [dns = begin(density.cubic()),
									sph = sphere.ref(),
									chrg = ps.valence_charge(),
									sp = sep_] GPU_LAMBDA (auto ipoint){
									 double rr = sph.distance(ipoint);
									 gpu::atomic::add(&dns[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]], chrg*sp.long_range_density(rr));
								 });
      }
			
			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("ionic_density::reduce");
				comm_.all_reduce_in_place_n(raw_pointer_cast(density.linear().data_elements()), density.linear().size(), std::plus<>{});
			}
			
			return density;			
    }
    
    template <class basis_type, class cell_type, class geo_type>
    basis::field<basis_type, double> atomic_electronic_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

			CALI_CXX_MARK_FUNCTION;

      basis::field<basis_type, double> density(basis);

			density = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);

				if(ps.has_electronic_density()){

					basis::spherical_grid sphere(basis, cell, atom_position, ps.electronic_density_radius());
					
					gpu::run(sphere.size(),
									 [dens = begin(density.cubic()),
										sph = sphere.ref(),
										spline = ps.electronic_density().cbegin()] GPU_LAMBDA (auto ipoint){
										 auto rr = sph.distance(ipoint);
										 auto density_val = spline.value(rr);
										 gpu::atomic::add(&dens[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]], density_val);
									 });

				} else {

					//just some crude guess for now
					basis::spherical_grid sphere(basis, cell, atom_position, 3.0);
					
					gpu::run(sphere.size(),
									 [dens = begin(density.cubic()),
										sph = sphere.ref(),
										zval = ps.valence_charge()] GPU_LAMBDA (auto ipoint){
										 auto rr = sph.distance(ipoint);
										 gpu::atomic::add(&dens[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]], zval/(M_PI)*exp(-2.0*rr));
									 });
					
				}
      }

			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("atomic_electronic_density::reduce");
				comm_.all_reduce_in_place_n(raw_pointer_cast(density.linear().data_elements()), density.linear().size(), std::plus<>{});
			}

			return density;			
    }

		auto has_nlcc() const {
			return has_nlcc_;
		}

		template <class basis_type, class cell_type, class geo_type>
    basis::field<basis_type, double> nlcc_density(const basis_type & basis, const cell_type & cell, const geo_type & geo) const {

			CALI_CXX_MARK_FUNCTION;
	
      basis::field<basis_type, double> density(basis);

			density = 0.0;
			
      for(auto iatom = part_.start(); iatom < part_.end(); iatom++){
				
				auto atom_position = geo.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(geo.atoms()[iatom]);

				if(not ps.has_nlcc_density()) continue;

				assert(has_nlcc());
				
				basis::spherical_grid sphere(basis, cell, atom_position, ps.nlcc_density_radius());

				gpu::run(sphere.size(),
								 [dens = begin(density.cubic()),
									sph = sphere.ref(),
									spline = ps.nlcc_density().cbegin()] GPU_LAMBDA (auto ipoint){
									 auto rr = sph.distance(ipoint);
									 auto density_val = spline.value(rr);
									 gpu::atomic::add(&dens[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]], density_val);
								 });
				
      }

			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("nlcc_density::reduce");
				comm_.all_reduce_in_place_n(raw_pointer_cast(density.linear().data_elements()), density.linear().size(), std::plus<>{});
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

		pseudo::math::erf_range_separation const sep_;
    double nelectrons_;
    pseudo::set pseudo_set_;
    std::unordered_map<std::string, pseudopotential_type> pseudopotential_list_;
		mutable boost::mpi3::communicator comm_;
		inq::utils::partition part_;
		bool has_nlcc_;

  };

}
}

#ifdef INQ_HAMILTONIAN_ATOMIC_POTENTIAL_UNIT_TEST
#undef INQ_HAMILTONIAN_ATOMIC_POTENTIAL_UNIT_TEST

#include <input/parse_xyz.hpp>
#include <catch2/catch_all.hpp>
#include <ions/geometry.hpp>
#include <basis/real_space.hpp>
#include <systems/box.hpp>

TEST_CASE("Class hamiltonian::atomic_potential", "[hamiltonian::atomic_potential]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using pseudo::element;
  using input::species;

	double const gcut = 0.785;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	SECTION("Non-existing element"){
    std::vector<species> el_list({element("P"), element("X")});

    CHECK_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list, gcut));
  }
  
  SECTION("Duplicated element"){
    std::vector<species> el_list({element("N"), element("N")});

    hamiltonian::atomic_potential pot(el_list.size(), el_list.begin(), gcut);

    CHECK(pot.num_species() == 1);
    CHECK(pot.num_electrons() == 10.0_a);
    
  }

  SECTION("Empty list"){
    std::vector<species> el_list;
    
    hamiltonian::atomic_potential pot(el_list.size(), el_list, gcut);

    CHECK(pot.num_species() == 0);
    CHECK(pot.num_electrons() == 0.0_a);

		CHECK(not pot.has_nlcc());
		
  }

  SECTION("CNOH"){
    species el_list[] = {element("C"), element("N"), element("O"), element("H")};

    hamiltonian::atomic_potential pot(4, el_list, gcut);

    CHECK(pot.num_species() == 4);
    CHECK(pot.num_electrons() == 16.0_a);
  }

  SECTION("Construct from a geometry"){

	ions::geometry geo(input::parse_xyz(config::path::unit_tests_data() + "benzene.xyz"));

	auto cell = systems::box::cubic(20.0_b).cutoff_energy(20.0_Ha);
	basis::real_space rs(cell);

	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms(), rs.gcutoff(), comm);

    CHECK(pot.num_species() == 2);
    CHECK(pot.num_electrons() == 30.0_a);

		rs.info(std::cout);
		
		auto vv = pot.local_potential(rs, cell, geo);

		CHECK(operations::integral(vv) == -45.5744357466_a);

		CHECK(vv.cubic()[5][3][0] == -1.6226427555_a);
		CHECK(vv.cubic()[3][1][0] == -0.2739253316_a);
							 
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
