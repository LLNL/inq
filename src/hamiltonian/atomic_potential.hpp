/* -*- indent-tabs-mode: t -*- */

#ifndef ATOMIC_POTENTIAL_HPP
#define ATOMIC_POTENTIAL_HPP

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <stdexcept>

#include <utils/raw_pointer_cast.hpp>
#include <config/path.hpp>

#include <pseudopod/set.hpp>
#include <pseudopod/pseudopotential.hpp>

#include <basis/spherical_grid.hpp>
#include <basis/double_grid.hpp>
#include <gpu/array.hpp>
#include <operations/integral.hpp>
#include <parallel/partition.hpp>
#include <solvers/poisson.hpp>
#include <states/ks_states.hpp>

#include <unordered_map>

#include <mpi3/environment.hpp>

#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

	class atomic_potential {

	public:

		using pseudopotential_type = pseudo::pseudopotential<gpu::array<double, 1>>;
		
		template <class atom_array>
		atomic_potential(const int natoms, const atom_array & atom_list, double gcutoff, bool double_grid):
			sep_(0.625), //this is the default from octopus
			natoms_(natoms),
			pseudo_set_("pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/"),
			double_grid_(double_grid)
		{

			CALI_CXX_MARK_FUNCTION;

			gcutoff *= double_grid_.spacing_factor(); 
			
			has_nlcc_ = false;
			nelectrons_ = 0.0;

			for(int iatom = 0; iatom < natoms; iatom++){
				if(!pseudo_set_.has(atom_list[iatom])) throw std::runtime_error("inq error: pseudopotential for element " + atom_list[iatom].symbol() + " not found.");
				
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

		template <class CommType, class basis_type, class ions_type>
		basis::field<basis_type, double> local_potential(CommType & comm, const basis_type & basis, const ions_type & ions, int single_atom = -1) const {

			CALI_CXX_MARK_SCOPE("atomic_potential::local_potential");
			
			basis::field<basis_type, double> potential(basis);

			parallel::partition part(natoms_, comm);
			
			potential.fill(0.0);
			
			for(auto iatom = part.start(); iatom < part.end(); iatom++){

				if(single_atom >= 0 and single_atom != iatom) continue;
				
				auto atom_position = ions.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(ions.atoms()[iatom]);
				basis::spherical_grid sphere(basis, atom_position, ps.short_range_potential_radius());

				if(not double_grid_.enabled()){

					gpu::run(sphere.size(),
									 [pot = begin(potential.cubic()),
										sph = sphere.ref(),
										spline = ps.short_range_potential().function()] GPU_LAMBDA (auto ipoint){
										 auto rr = sph.distance(ipoint);
										 auto potential_val = spline(rr);
										 gpu::atomic::add(&pot[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]], potential_val);
									 });

				} else {

					CALI_CXX_MARK_SCOPE("atomic_potential::double_grid");
					
					gpu::run(sphere.size(),
									 [pot = begin(potential.cubic()),
										sph = sphere.ref(),
										spline = ps.short_range_potential().function(),
										dg = double_grid_.ref(),
										spac = basis.rspacing(), metric = basis.cell().metric()] GPU_LAMBDA (auto ipoint){
										 gpu::atomic::add(&pot[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]],
																			dg.value([spline] GPU_LAMBDA (auto pos) { return spline(pos.length()); }, spac, metric.to_cartesian(sph.point_pos(ipoint))));
									 });
				}
			}
			
			potential.all_reduce(comm);
			return potential;			
		}

		////////////////////////////////////////////////////////////////////////////////////
		
		template <class CommType, class basis_type, class ions_type>
		basis::field<basis_type, double> ionic_density(CommType & comm, const basis_type & basis, const ions_type & ions, int single_atom = -1) const {

			CALI_CXX_MARK_FUNCTION;

			parallel::partition part(natoms_, comm);
			
			basis::field<basis_type, double> density(basis);
			
			density.fill(0.0);

			for(auto iatom = part.start(); iatom < part.end(); iatom++){

				if(single_atom >= 0 and single_atom != iatom) continue;
				
				auto atom_position = ions.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(ions.atoms()[iatom]);
				basis::spherical_grid sphere(basis, atom_position, sep_.long_range_density_radius());

				//OPTIMIZATION: this should be done in parallel for atoms too
				gpu::run(sphere.size(),
								 [dns = begin(density.cubic()),
									sph = sphere.ref(),
									chrg = ps.valence_charge(),
									sp = sep_] GPU_LAMBDA (auto ipoint){
									 double rr = sph.distance(ipoint);
									 gpu::atomic::add(&dns[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]], chrg*sp.long_range_density(rr));
								 });
			}

			density.all_reduce(comm);
			return density;			
		}

		////////////////////////////////////////////////////////////////////////////////////
		
		template <class CommType, class basis_type, class ions_type>
		basis::field_set<basis_type, double> atomic_electronic_density(CommType & comm, const basis_type & basis, const ions_type & ions, states::ks_states const & states) const {

			CALI_CXX_MARK_FUNCTION;

			parallel::partition part(natoms_, comm);

			auto nspin = states.num_density_components();
			basis::field_set<basis_type, double> density(basis, nspin);

			density.fill(0.0);

			double polarization = 1.0;
			if(nspin == 2) polarization = 0.6;
			
			for(auto iatom = part.start(); iatom < part.end(); iatom++){
				
				auto atom_position = ions.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(ions.atoms()[iatom]);

				if(ps.has_electronic_density()){

					basis::spherical_grid sphere(basis, atom_position, ps.electronic_density_radius());
					
					gpu::run(nspin, sphere.size(),
									 [dens = begin(density.hypercubic()), sph = sphere.ref(), spline = ps.electronic_density().function(), polarization] GPU_LAMBDA (auto ispin, auto ipoint){
										 auto rr = sph.distance(ipoint);
										 auto density_val = spline(rr);
										 auto pol = polarization;
										 if(ispin == 1) pol = 1.0 - pol;
										 gpu::atomic::add(&dens[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]][ispin], pol*density_val);
									 });

				} else {

					//just some crude guess for now
					basis::spherical_grid sphere(basis, atom_position, 3.0);
					
					gpu::run(nspin, sphere.size(),
									 [dens = begin(density.hypercubic()), sph = sphere.ref(), zval = ps.valence_charge(), polarization] GPU_LAMBDA (auto ispin, auto ipoint){
										 auto rr = sph.distance(ipoint);
										 auto pol = polarization;
										 if(ispin == 1) pol = 1.0 - pol;
										 gpu::atomic::add(&dens[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]][ispin], pol*zval/(M_PI)*exp(-2.0*rr));
									 });
					
				}
			}

			density.all_reduce(comm);
			return density;			
		}

		////////////////////////////////////////////////////////////////////////////////////

		auto has_nlcc() const {
			return has_nlcc_;
		}

		////////////////////////////////////////////////////////////////////////////////////

		template <class CommType, class basis_type, class ions_type>
		basis::field<basis_type, double> nlcc_density(CommType & comm, const basis_type & basis, const ions_type & ions) const {

			CALI_CXX_MARK_FUNCTION;

			parallel::partition part(natoms_, comm);
			
			basis::field<basis_type, double> density(basis);

			density.fill(0.0);
			
			for(auto iatom = part.start(); iatom < part.end(); iatom++){
				
				auto atom_position = ions.coordinates()[iatom];
				
				auto & ps = pseudo_for_element(ions.atoms()[iatom]);

				if(not ps.has_nlcc_density()) continue;

				assert(has_nlcc());
				
				basis::spherical_grid sphere(basis, atom_position, ps.nlcc_density_radius());

				gpu::run(sphere.size(),
								 [dens = begin(density.cubic()),
									sph = sphere.ref(),
									spline = ps.nlcc_density().function()] GPU_LAMBDA (auto ipoint){
									 auto rr = sph.distance(ipoint);
									 auto density_val = spline(rr);
									 gpu::atomic::add(&dens[sph.grid_point(ipoint)[0]][sph.grid_point(ipoint)[1]][sph.grid_point(ipoint)[2]], density_val);
								 });
				
			}

			density.all_reduce(comm);
			return density;			
		}
		
		template <class output_stream>
		void info(output_stream & out) const {
			out << "ATOMIC POTENTIAL:" << std::endl;
			out << "	Number of species		= " << num_species() << std::endl;
			out << "	Number of electrons = " << num_electrons() << std::endl;
			out << std::endl;
		}

		auto & range_separation() const {
			return sep_;
		}

		auto & double_grid() const {
			return double_grid_;
		}
		
	private:

		pseudo::math::erf_range_separation sep_;
		int natoms_;
		double nelectrons_;
		pseudo::set pseudo_set_;
		std::unordered_map<std::string, pseudopotential_type> pseudopotential_list_;
		bool has_nlcc_;
		basis::double_grid double_grid_;

	};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ATOMIC_POTENTIAL_UNIT_TEST
#undef INQ_HAMILTONIAN_ATOMIC_POTENTIAL_UNIT_TEST

#include <input/parse_xyz.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;
	using pseudo::element;
	using input::species;

	double const gcut = 0.785;

	parallel::communicator comm{boost::mpi3::environment::get_self_instance()};
	
	SECTION("Non-existing element"){
		std::vector<species> el_list({element("P"), element("X")});
	
		CHECK_THROWS(hamiltonian::atomic_potential(el_list.size(), el_list, gcut, /*double_grid = */ false));
	}
	
	SECTION("Duplicated element"){
		std::vector<species> el_list({element("N"), element("N")});

		hamiltonian::atomic_potential pot(el_list.size(), el_list.begin(), gcut, /*double_grid = */ false);

		CHECK(pot.num_species() == 1);
		CHECK(pot.num_electrons() == 10.0_a);
		
	}

	SECTION("Empty list"){
		std::vector<species> el_list;
		
		hamiltonian::atomic_potential pot(el_list.size(), el_list, gcut, false);

		CHECK(pot.num_species() == 0);
		CHECK(pot.num_electrons() == 0.0_a);

		CHECK(not pot.has_nlcc());
		
	}

	SECTION("CNOH"){
		species el_list[] = {element("C"), element("N"), element("O"), element("H")};

		hamiltonian::atomic_potential pot(4, el_list, gcut, /*double_grid = */ false);

		CHECK(pot.num_species() == 4);
		CHECK(pot.num_electrons() == 16.0_a);
	}

	SECTION("Construct from a geometry"){

		auto ions = systems::ions(ions::unit_cell::cubic(20.0_b));

		ions.insert(input::parse_xyz(config::path::unit_tests_data() + "benzene.xyz"));
		
		basis::real_space rs(ions.cell(), /*spacing = */ 0.49672941, comm);
		
		hamiltonian::atomic_potential pot(ions.num_atoms(), ions.atoms(), rs.gcutoff(), /*double_grid = */ false);
		
		CHECK(pot.num_species() == 2);
		CHECK(pot.num_electrons() == 30.0_a);
		
		rs.info(std::cout);
		
		auto vv = pot.local_potential(comm, rs, ions);
		
		CHECK(operations::integral(vv) == -45.5744357466_a);
		
		CHECK(vv.cubic()[5][3][0] == -1.6226427555_a);
		CHECK(vv.cubic()[3][1][0] == -0.2739253316_a);
		
		auto id = pot.ionic_density(comm, rs, ions);
		
		CHECK(operations::integral(id) == -30.0000000746_a);
		CHECK(id.cubic()[5][3][0] == -0.9448936487_a);
		CHECK(id.cubic()[3][1][0] == -0.2074502252_a);

		states::ks_states unp(states::ks_states::spin_config::UNPOLARIZED, 11.0);
		
		auto nn_unp = pot.atomic_electronic_density(comm, rs, ions, unp);

		CHECK(nn_unp.set_size() == 1);		
		CHECK(operations::integral_sum(nn_unp) == 29.9562520003_a);
		CHECK(nn_unp.hypercubic()[5][3][0][0] == 0.1330589609_a);
		CHECK(nn_unp.hypercubic()[3][1][0][0] == 0.1846004508_a);

		states::ks_states pol(states::ks_states::spin_config::POLARIZED, 11.0);
		
		auto nn_pol = pot.atomic_electronic_density(comm, rs, ions, pol);

		CHECK(nn_pol.set_size() == 2);
		CHECK(operations::integral_sum(nn_pol) == 29.9562519176_a);
		CHECK(nn_pol.hypercubic()[5][3][0][0] == Approx(1.2*0.066529473));
		CHECK(nn_pol.hypercubic()[3][1][0][0] == Approx(1.2*0.0923002217));
		CHECK(nn_pol.hypercubic()[5][3][0][1] == Approx(0.8*0.066529473));
		CHECK(nn_pol.hypercubic()[3][1][0][1] == Approx(0.8*0.0923002217));
		
		CHECK(pot.has_nlcc());
		
		auto nlcc = pot.nlcc_density(comm, rs, ions);
		
		CHECK(operations::integral(nlcc) == 3.0083012065_a);
		CHECK(nlcc.cubic()[5][3][0] == 0.6248217151_a);
		CHECK(nlcc.cubic()[3][1][0] == 0.0007040027_a);
	
	}
	
}
#endif
