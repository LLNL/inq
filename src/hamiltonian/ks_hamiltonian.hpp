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
#include <states/ks_states.hpp>
#include <multi/adaptors/fftw.hpp>
#include <hamiltonian/projector.hpp>
#include <hamiltonian/projector_fourier.hpp>
#include <hamiltonian/exchange_operator.hpp>
#include <operations/space.hpp>
#include <operations/laplacian.hpp>

namespace inq {
namespace hamiltonian {
  template <class basis_type>
  class ks_hamiltonian {
		
  public:

		basis::field<basis::real_space, double> scalar_potential;
		exchange_operator exchange;
		
    ks_hamiltonian(const basis_type & basis, const ions::UnitCell & cell, const atomic_potential & pot, bool fourier_pseudo, const ions::geometry & geo,
									 const int num_hf_orbitals, const double exchange_coefficient, boost::mpi3::cartesian_communicator<2> const & comm = {boost::mpi3::environment::get_world_instance(), {}}):
			scalar_potential(basis),
			exchange(basis, num_hf_orbitals, exchange_coefficient, comm),
			non_local_in_fourier_(fourier_pseudo)
		{
			
			scalar_potential = pot.local_potential(basis, cell, geo);
			
			for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
				if(non_local_in_fourier_){
					auto insert = projectors_fourier_map_.emplace(geo.atoms()[iatom].symbol(),
																												projector_fourier(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom])));
					insert.first->second.add_coord(geo.coordinates()[iatom]);
				} else {
					projectors_.push_back(projector(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom]), geo.coordinates()[iatom]));
				}
			}

    }

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void non_local(const basis::field_set<basis::real_space, complex> & phi, basis::field_set<basis::real_space, complex> & vnlphi) const {
			for(unsigned iproj = 0; iproj < projectors_.size(); iproj++) projectors_[iproj](phi, vnlphi);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void non_local(const basis::field_set<basis::fourier_space, complex> & phi, basis::field_set<basis::fourier_space, complex> & vnlphi) const {
			for(auto it = projectors_fourier_map_.cbegin(); it != projectors_fourier_map_.cend(); ++it){
				it->second(phi, vnlphi);
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////
		
		auto non_local(const basis::field_set<basis::real_space, complex> & phi) const {

			if(non_local_in_fourier_) {

				auto phi_fs = operations::space::to_fourier(phi);
				basis::field_set<basis::fourier_space, complex> vnlphi_fs(phi_fs.skeleton());

				vnlphi_fs = 0.0;
				non_local(phi_fs, vnlphi_fs);
				return operations::space::to_real(vnlphi_fs);
					
			} else {

				basis::field_set<basis::real_space, complex> vnlphi(phi.skeleton());
				vnlphi = 0.0;
				non_local(phi, vnlphi);
				return vnlphi;
							
			}
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void fourier_space_terms(const basis::field_set<basis::fourier_space, complex> & phi, basis::field_set<basis::fourier_space, complex> & hphi) const {
			operations::laplacian_add(phi, hphi);
			if(non_local_in_fourier_) non_local(phi, hphi);
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void real_space_terms(const basis::field_set<basis::real_space, complex> & phi, basis::field_set<basis::real_space, complex> & hphi) const {
			//the non local potential in real space
			if(not non_local_in_fourier_) non_local(phi, hphi);

			assert(scalar_potential.linear().num_elements() == phi.basis().local_size());
			
			//the scalar local potential in real space
			//DATAOPERATIONS LOOP + GPU:RUN 2D
#ifdef ENABLE_CUDA
			gpu::run(phi.local_set_size(), phi.basis().local_size(),
							 [pot = begin(scalar_potential.linear()), it_hphi = begin(hphi), it_phi = begin(hphi)] __device__
							 (auto ist, auto ip){
								 it_hphi[ip][ist] += pot[ip]*it_phi[ip][ist];
							 });
							 
#else
			for(long ip = 0; ip < phi.basis().local_size(); ip++){
				double vv  = scalar_potential.linear()[ip];
				for(int ist = 0; ist < phi.local_set_size(); ist++) hphi.matrix()[ip][ist] += vv*phi.matrix()[ip][ist];
			}
#endif
			
			exchange(phi, hphi);
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    auto operator()(const basis::field_set<basis::real_space, complex> & phi) const{
      
			auto phi_fs = operations::space::to_fourier(phi);
			
			basis::field_set<basis::fourier_space, complex> hphi_fs(phi_fs.skeleton());
			
			hphi_fs = 0.0;
			
			fourier_space_terms(phi_fs, hphi_fs);
	
			auto hphi = operations::space::to_real(hphi_fs);

			real_space_terms(phi, hphi);

			return hphi;
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
    auto operator()(const basis::field_set<basis::fourier_space, complex> & phi) const{

			auto phi_rs = operations::space::to_real(phi);

			basis::field_set<basis::real_space, complex> hphi_rs(phi_rs.skeleton());

			hphi_rs = 0.0;
						
			real_space_terms(phi_rs, hphi_rs);
		
			auto hphi = operations::space::to_fourier(hphi_rs);

			fourier_space_terms(phi, hphi);

			return hphi;
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		int num_projectors() const {
			int nn = 0;
			for(unsigned iproj = 0; iproj < projectors_.size(); iproj++){
				nn += projectors_[iproj].num_projectors();
			}
			return nn;			
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    template <class output_stream>
    void info(output_stream & out) const {
      out << "HAMILTONIAN:" << std::endl;
      out << "  Total number of projectors = " << num_projectors() << std::endl;
      out << std::endl;
    }	
		
  private:

		std::vector<projector> projectors_;
		bool non_local_in_fourier_;
		std::unordered_map<std::string, projector_fourier> projectors_fourier_map_;
		std::vector<std::unordered_map<std::string, projector_fourier>::iterator> projectors_fourier_;
		
  };

}
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::ks_hamiltonian", "[hamiltonian::ks_hamiltonian]"){

	using namespace inq;
	using namespace Catch::literals;
  using math::vec3d;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);

	SECTION("Basis"){
		
		CHECK(rs.size() == 8000);
		CHECK(rs.rspacing()[0] == 0.5_a);
		CHECK(rs.rspacing()[1] == 0.5_a);	
		CHECK(rs.rspacing()[2] == 0.5_a);
		CHECK(rs.volume_element() == 0.125_a);
	}
	
	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms(), rs.gcutoff(), set_comm);
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

  basis::field_set<basis::real_space, complex> phi(rs, st.num_states(), cart_comm);
	basis::field_set<basis::real_space, complex> hphi(rs, st.num_states(), cart_comm);

	hamiltonian::ks_hamiltonian<basis::real_space> ham(rs, cell, pot, false, geo, st.num_states(), 0.0, cart_comm);

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
		
		hphi = ham(phi);
		
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
						
						double xx = rs.rvector(ixg, iyg, izg)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg*kk*xx), sin(istg*kk*xx));
					}
				}
			}
		}

		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						auto istg = phi.set_part().local_to_global(ist);
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*istg*kk*istg*kk*phi.cubic()[ix][iy][iz][ist]);
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
					
					double r2 = rs.r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-0.5*ww*r2);
					}
					
				}
			}
		}

		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - ww*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff == 0.0055687279_a);
		
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
						
						double xx = rs.rvector(ixg, iyg, izg)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg*kk*xx), sin(istg*kk*xx));
					}
				}
			}
		}

		hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){

						auto istg = phi.set_part().local_to_global(ist);

						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*istg*kk*istg*kk*phi.cubic()[ix][iy][iz][ist]);
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
					
					double r2 = rs.r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-0.5*ww*r2);
					}
					
				}
			}
		}

		hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.local_set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - ww*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		diff /= hphi.set_size()*hphi.basis().size();

		CHECK(diff == 0.0055687279_a);
		
	}
	
}

#endif

#endif
