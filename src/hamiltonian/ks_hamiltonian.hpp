/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__KS_HAMILTONIAN
#define INQ__HAMILTONIAN__KS_HAMILTONIAN

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
#include <states/ks_states.hpp>
#include <states/orbital_set.hpp>

#include <utils/profiling.hpp>

#include <list>
#include <unordered_map>

namespace inq {
namespace hamiltonian {
  template <class basis_type>
  class ks_hamiltonian {
		
  public:

		void update_projectors(const basis_type & basis, const ions::UnitCell & cell, const atomic_potential & pot, const ions::geometry & geo){
			
			CALI_CXX_MARK_FUNCTION;
			
			projectors_.clear();
			projectors_fourier_map_.clear();			
			
			for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
				if(non_local_in_fourier_){
					auto insert = projectors_fourier_map_.emplace(geo.atoms()[iatom].symbol(),
																												projector_fourier(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom])));
					insert.first->second.add_coord(geo.coordinates()[iatom]);
				} else {
					projectors_.emplace_back(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom]), geo.coordinates()[iatom], iatom);
					if(projectors_.back().empty()) projectors_.pop_back(); 
				}
			}

			projectors_all_ = projector_all(projectors_);
			
		}
		
		basis::field<basis::real_space, double> scalar_potential;
		exchange_operator exchange;
		
    ks_hamiltonian(const basis_type & basis, const ions::UnitCell & cell, const atomic_potential & pot, bool fourier_pseudo, const ions::geometry & geo,
									 const int num_hf_orbitals, const double exchange_coefficient, boost::mpi3::cartesian_communicator<2> comm):
			scalar_potential(basis),
			exchange(basis, num_hf_orbitals, exchange_coefficient, std::move(comm)),
			non_local_in_fourier_(fourier_pseudo)
		{
			
			scalar_potential = pot.local_potential(basis, cell, geo);
			update_projectors(basis, cell, pot, geo);

    }

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void non_local(const basis::field_set<basis::fourier_space, complex> & phi, basis::field_set<basis::fourier_space, complex> & vnlphi) const {

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

				vnlphi_fs.fields() = 0.0;
				non_local(phi_fs.fields(), vnlphi_fs.fields());
				return operations::space::to_real(vnlphi_fs);
					
			} else {

				auto proj = projectors_all_.project(phi.fields(), phi.kpoint());
				
				states::orbital_set<basis::real_space, complex> vnlphi(phi.skeleton());
				vnlphi.fields() = 0.0;

				projectors_all_.apply(proj, vnlphi.fields(), phi.kpoint());
			
				return vnlphi;
							
			}
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    auto operator()(const states::orbital_set<basis::real_space, complex> & phi) const {
			
			CALI_CXX_MARK_SCOPE("hamiltonian_real");

			auto proj = projectors_all_.project(phi.fields(), phi.kpoint());
			
			auto phi_fs = operations::space::to_fourier(phi.fields());
		
			auto hphi_fs = operations::laplacian(phi_fs, -0.5, -2.0*phi.kpoint());

			non_local(phi_fs, hphi_fs);
			
			auto hphi = states::orbital_set<basis::real_space, complex>{operations::space::to_real(hphi_fs), phi.kpoint()};

			hamiltonian::scalar_potential_add(scalar_potential, 0.5*norm(phi.kpoint()), phi, hphi);
			exchange(phi.fields(), hphi.fields());

			projectors_all_.apply(proj, hphi.fields(), phi.kpoint());

			return hphi;
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    auto operator()(const states::orbital_set<basis::fourier_space, complex> & phi) const{
			
			CALI_CXX_MARK_SCOPE("hamiltonian_fourier");

			auto phi_rs = operations::space::to_real(phi);

			auto proj = projectors_all_.project(phi_rs.fields(), phi.kpoint());
			
			auto hphi_rs = hamiltonian::scalar_potential(scalar_potential, 0.5*norm(phi.kpoint()), phi_rs);
		
			exchange(phi_rs.fields(), hphi_rs.fields());
 
			projectors_all_.apply(proj, hphi_rs.fields(), phi.kpoint());
			
			auto hphi = operations::space::to_fourier(hphi_rs.fields());

			operations::laplacian_add(phi.fields(), hphi, -0.5, -2.0*phi.kpoint());
			non_local(phi.fields(), hphi);

			return states::orbital_set<basis::fourier_space, complex>{std::move(hphi), phi.kpoint()};
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////
		
		int num_projectors() const {
			int nn = 0;
			for(auto it = projectors_.cbegin(); it != projectors_.cend(); ++it){
				nn += it->num_projectors();
			}
			return nn;			
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		auto & projectors() const {
			return projectors_;			
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    template <class output_stream>
    void info(output_stream & out) const {
      out << "HAMILTONIAN:" << std::endl;
      out << "  Total number of projectors = " << num_projectors() << std::endl;
      out << std::endl;
    }	
		
  private:

		std::list<projector> projectors_;
		projector_all projectors_all_;		
		bool non_local_in_fourier_;
		std::unordered_map<std::string, projector_fourier> projectors_fourier_map_;
		std::vector<std::unordered_map<std::string, projector_fourier>::iterator> projectors_fourier_;
		
  };

}
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_HAMILTONIAN_KS_HAMILTONIAN_UNIT_TEST
#undef INQ_HAMILTONIAN_KS_HAMILTONIAN_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::ks_hamiltonian", "[hamiltonian::ks_hamiltonian]"){

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
  using math::vector3;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

	ions::geometry geo;
	systems::box box = systems::box::cubic(10.0_b).cutoff_energy(20.0_Ha);
  basis::real_space rs(box, basis_comm);

	SECTION("Basis"){
		
		CHECK(rs.size() == 8000);
		CHECK(rs.rspacing()[0] == 0.5_a);
		CHECK(rs.rspacing()[1] == 0.5_a);	
		CHECK(rs.rspacing()[2] == 0.5_a);
		CHECK(rs.volume_element() == 0.125_a);
	}
	
	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms(), rs.gcutoff(), set_comm);
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

  states::orbital_set<basis::real_space, complex> phi(rs, st.num_states(), cart_comm);

	hamiltonian::ks_hamiltonian<basis::real_space> ham(rs, box, pot, false, geo, st.num_states(), 0.0, cart_comm);

	SECTION("Constant function"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					ham.scalar_potential.cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = 1.0;
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
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.0);
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

					ham.scalar_potential.cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);	
						auto istg = phi.set_part().local_to_global(ist);
						
						double xx = rs.point_op().rvector(ixg, iyg, izg)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
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
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*istg.value()*kk*istg.value()*kk*phi.cubic()[ix][iy][iz][ist]);
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

					auto ixg = rs.cubic_dist(0).local_to_global(ix);
					auto iyg = rs.cubic_dist(1).local_to_global(iy);
					auto izg = rs.cubic_dist(2).local_to_global(iz);	
					
					double r2 = rs.point_op().r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-ww*r2);
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
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 1.5*ww*phi.cubic()[ix][iy][iz][ist]);
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

					ham.scalar_potential.cubic()[ix][iy][iz] = 0.0;
					
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);	
						auto istg = phi.set_part().local_to_global(ist);
						
						double xx = rs.point_op().rvector(ixg, iyg, izg)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
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

						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*istg.value()*kk*istg.value()*kk*phi.cubic()[ix][iy][iz][ist]);
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

					auto ixg = rs.cubic_dist(0).local_to_global(ix);
					auto iyg = rs.cubic_dist(1).local_to_global(iy);
					auto izg = rs.cubic_dist(2).local_to_global(iz);	
					
					double r2 = rs.point_op().r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-ww*r2);
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
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 1.5*ww*phi.cubic()[ix][iy][iz][ist]);
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

#endif
