/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__KS_HAMILTONIAN
#define INQ__HAMILTONIAN__KS_HAMILTONIAN

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <basis/field.hpp>
#include <multi/adaptors/fftw.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <hamiltonian/exchange_operator.hpp>
#include <hamiltonian/projector.hpp>
#include <hamiltonian/projector_all.hpp>
#include <hamiltonian/projector_fourier.hpp>
#include <hamiltonian/scalar_potential.hpp>
#include <input/environment.hpp>
#include <ions/geometry.hpp>
#include <operations/space.hpp>
#include <operations/laplacian.hpp>
#include <operations/gradient.hpp>
#include <states/ks_states.hpp>
#include <states/orbital_set.hpp>

#include <utils/profiling.hpp>

#include <list>
#include <unordered_map>

namespace inq {
namespace hamiltonian {

template <typename PotentialType>
class ks_hamiltonian {
		
public:

	using potential_type = PotentialType;
	
	void update_projectors(const basis::real_space & basis, const atomic_potential & pot, const ions::geometry & geo){
			
		CALI_CXX_MARK_FUNCTION;

		std::list<projector> projectors;
			
		projectors_fourier_map_.clear();			
			
		for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
			if(non_local_in_fourier_){
				auto insert = projectors_fourier_map_.emplace(geo.atoms()[iatom].symbol(), projector_fourier(basis, pot.pseudo_for_element(geo.atoms()[iatom])));
				insert.first->second.add_coord(basis.cell().metric().to_contravariant(geo.coordinates()[iatom]));
			} else {
				projectors.emplace_back(basis, pot.pseudo_for_element(geo.atoms()[iatom]), geo.coordinates()[iatom], iatom);
				if(projectors.back().empty()) projectors.pop_back(); 
			}
		}

		projectors_all_ = projector_all(projectors);
	}
		
	exchange_operator exchange;

	////////////////////////////////////////////////////////////////////////////////////////////
		
	ks_hamiltonian(const basis::real_space & basis, ions::brillouin const & bzone, states::ks_states const & states, const atomic_potential & pot, bool fourier_pseudo, const ions::geometry & geo,
								 const int num_hf_orbitals, const double exchange_coefficient, bool use_ace = false):
		exchange(basis.cell(), bzone, exchange_coefficient, use_ace),
		scalar_potential_(states.num_density_components(), basis),
		uniform_vector_potential_({0.0, 0.0, 0.0}),
		non_local_in_fourier_(fourier_pseudo),
		states_(states)
	{
		for(auto & pot : scalar_potential_) pot.fill(0.0);
		update_projectors(basis, pot, geo);
	}

	////////////////////////////////////////////////////////////////////////////////////////////
		
	void non_local(const states::orbital_set<basis::fourier_space, complex> & phi, states::orbital_set<basis::fourier_space, complex> & vnlphi) const {

		if(not non_local_in_fourier_) return;
			
		for(auto it = projectors_fourier_map_.cbegin(); it != projectors_fourier_map_.cend(); ++it){
			it->second(phi, vnlphi);
		}
	}
		
	////////////////////////////////////////////////////////////////////////////////////////////
		
	auto non_local(const states::orbital_set<basis::real_space, complex> & phi) const {

		CALI_CXX_MARK_FUNCTION;
 
		if(non_local_in_fourier_) {

			auto phi_fs = operations::space::to_fourier(phi);
			states::orbital_set<basis::fourier_space, complex> vnlphi_fs(phi_fs.skeleton());

			vnlphi_fs.fill(0.0);
			non_local(phi_fs, vnlphi_fs);
			return operations::space::to_real(vnlphi_fs);
					
		} else {
				
			auto proj = projectors_all_.project(phi, phi.kpoint());
				
			states::orbital_set<basis::real_space, complex> vnlphi(phi.skeleton());
			vnlphi.fill(0.0);

			projectors_all_.apply(proj, vnlphi, phi.kpoint());
			
			return vnlphi;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	auto operator()(const states::orbital_set<basis::real_space, complex> & phi) const {
			
		CALI_CXX_MARK_SCOPE("hamiltonian_real");

		auto proj = projectors_all_.project(phi, phi.kpoint() + uniform_vector_potential_);
			
		auto phi_fs = operations::space::to_fourier(phi);
		
		auto hphi_fs = operations::laplacian(phi_fs, -0.5, -2.0*phi.basis().cell().metric().to_contravariant(phi.kpoint() + uniform_vector_potential_));

		non_local(phi_fs, hphi_fs);
			
		auto hphi = operations::space::to_real(hphi_fs);

		hamiltonian::scalar_potential_add(scalar_potential_[phi.spin_index()], 0.5*phi.basis().cell().metric().norm(phi.kpoint() + uniform_vector_potential_), phi, hphi);
		exchange(phi, hphi);

		projectors_all_.apply(proj, hphi, phi.kpoint() + uniform_vector_potential_);

		return hphi;
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	auto operator()(const states::orbital_set<basis::fourier_space, complex> & phi) const{
			
		CALI_CXX_MARK_SCOPE("hamiltonian_fourier");

		auto phi_rs = operations::space::to_real(phi);

		auto proj = projectors_all_.project(phi_rs, phi.kpoint() + uniform_vector_potential_);
			
		auto hphi_rs = hamiltonian::scalar_potential(scalar_potential_[phi.spin_index()], 0.5*phi.basis().cell().metric().norm(phi.kpoint() + uniform_vector_potential_), phi_rs);
		
		exchange(phi_rs, hphi_rs);
 
		projectors_all_.apply(proj, hphi_rs, phi.kpoint() + uniform_vector_potential_);
			
		auto hphi = operations::space::to_fourier(hphi_rs);

		operations::laplacian_add(phi, hphi, -0.5, -2.0*phi.basis().cell().metric().to_contravariant(phi.kpoint() + uniform_vector_potential_));
		non_local(phi, hphi);

		return hphi;
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	auto momentum(const states::orbital_set<basis::real_space, complex> & phi) const{
		CALI_CXX_MARK_FUNCTION;

		return operations::gradient(phi, /* factor = */ 1.0, /*shift = */ phi.kpoint() + uniform_vector_potential_);
	}
		
	////////////////////////////////////////////////////////////////////////////////////////////
		
	auto & projectors_all() const {
		return projectors_all_;
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	template <class output_stream>
	void info(output_stream & out) const {
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	auto & scalar_potential() {
		return scalar_potential_;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////

		auto & uniform_vector_potential() const {
			return uniform_vector_potential_;
		}
		auto & uniform_vector_potential() {
			return uniform_vector_potential_;
		}

		vector3<double, covariant> vp_accel={0,0,0};
		vector3<double, covariant> vp_velocity={0,0,0};
		vector3<double, covariant> vp_induced={0,0,0};

	////////////////////////////////////////////////////////////////////////////////////////////
		
private:
		
	std::vector<basis::field<basis::real_space, PotentialType>> scalar_potential_;
	vector3<double, covariant> uniform_vector_potential_;
	projector_all projectors_all_;		
	bool non_local_in_fourier_;
	std::unordered_map<std::string, projector_fourier> projectors_fourier_map_;
	std::vector<std::unordered_map<std::string, projector_fourier>::iterator> projectors_fourier_;
	states::ks_states states_;
		
	template <typename Perturbation>
	friend class self_consistency;
	
};
}
}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
#endif

#ifdef INQ_HAMILTONIAN_KS_HAMILTONIAN_UNIT_TEST
#undef INQ_HAMILTONIAN_KS_HAMILTONIAN_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
  
	parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto set_comm = basis::set_subcomm(cart_comm);
	auto basis_comm = basis::basis_subcomm(cart_comm);	

	systems::box box = systems::box::cubic(10.0_b).cutoff_energy(20.0_Ha);
	auto ions = systems::ions(box);

  basis::real_space rs(box, basis_comm);

	SECTION("Basis"){
		
		CHECK(rs.size() == 8000);
		CHECK(rs.rspacing()[0] == 0.5_a);
		CHECK(rs.rspacing()[1] == 0.5_a);	
		CHECK(rs.rspacing()[2] == 0.5_a);
		CHECK(rs.volume_element() == 0.125_a);
	}
	
	hamiltonian::atomic_potential pot(ions.geo().num_atoms(), ions.geo().atoms(), rs.gcutoff());
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

  states::orbital_set<basis::real_space, complex> phi(rs, st.num_states(), 1, vector3<double, covariant>{0.0, 0.0, 0.0}, 0, cart_comm);

	auto bzone = ions::brillouin(ions, input::kpoints::gamma());
	
	hamiltonian::ks_hamiltonian<double> ham(rs, bzone, st, pot, false, ions.geo(), st.num_states(), 0.0);

	SECTION("Constant function"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					ham.scalar_potential()[0].cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.hypercubic()[ix][iy][iz][ist] = 1.0;
					}
				}
			}
		}
		
		auto hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						diff += fabs(hphi.hypercubic()[ix][iy][iz][ist] - 0.0);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();
		
		CHECK(diff < 1e-14);
		
	}
	
	SECTION("Plane wave"){
		
		double kk = 2.0*M_PI/rs.rlength()[0];
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					ham.scalar_potential()[0].cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto ixg = rs.cubic_part(0).local_to_global(ix);
						auto iyg = rs.cubic_part(1).local_to_global(iy);
						auto izg = rs.cubic_part(2).local_to_global(iz);	
						auto istg = phi.set_part().local_to_global(ist);
						
						double xx = rs.point_op().rvector_cartesian(ixg, iyg, izg)[0];
						phi.hypercubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
					}
				}
			}
		}

		auto hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						auto istg = phi.set_part().local_to_global(ist);
						diff += fabs(hphi.hypercubic()[ix][iy][iz][ist] - 0.5*istg.value()*kk*istg.value()*kk*phi.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff < 1e-14);
		
	}

	SECTION("Harmonic oscillator"){

		double ww = 2.0;

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					auto ixg = rs.cubic_part(0).local_to_global(ix);
					auto iyg = rs.cubic_part(1).local_to_global(iy);
					auto izg = rs.cubic_part(2).local_to_global(iz);	
					
					double r2 = rs.point_op().r2(ixg, iyg, izg);
					ham.scalar_potential()[0].cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.hypercubic()[ix][iy][iz][ist] = exp(-ww*r2);
					}
					
				}
			}
		}

		auto hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						diff += fabs(hphi.hypercubic()[ix][iy][iz][ist] - 1.5*ww*phi.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff == 0.0051420503_a);
		
	}


	SECTION("Plane wave - fourier"){
		
		double kk = 2.0*M_PI/rs.rlength()[0];
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					ham.scalar_potential()[0].cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto ixg = rs.cubic_part(0).local_to_global(ix);
						auto iyg = rs.cubic_part(1).local_to_global(iy);
						auto izg = rs.cubic_part(2).local_to_global(iz);	
						auto istg = phi.set_part().local_to_global(ist);
						
						double xx = rs.point_op().rvector_cartesian(ixg, iyg, izg)[0];
						phi.hypercubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
					}
				}
			}
		}

		auto hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto istg = phi.set_part().local_to_global(ist);

						diff += fabs(hphi.hypercubic()[ix][iy][iz][ist] - 0.5*istg.value()*kk*istg.value()*kk*phi.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff < 1e-14);
		
	}

	SECTION("Harmonic oscillator - fourier"){

		double ww = 2.0;

		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					auto ixg = rs.cubic_part(0).local_to_global(ix);
					auto iyg = rs.cubic_part(1).local_to_global(iy);
					auto izg = rs.cubic_part(2).local_to_global(iz);	
					
					double r2 = rs.point_op().r2(ixg, iyg, izg);
					ham.scalar_potential()[0].cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.hypercubic()[ix][iy][iz][ist] = exp(-ww*r2);
					}
					
				}
			}
		}

		auto hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						diff += fabs(hphi.hypercubic()[ix][iy][iz][ist] - 1.5*ww*phi.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff == 0.0051420503_a);
		
	}
	
}
#endif
