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
#include <hamiltonian/atomic_potential.hpp>
#include <hamiltonian/exchange_operator.hpp>
#include <hamiltonian/projector.hpp>
#include <hamiltonian/projector_fourier.hpp>
#include <hamiltonian/scalar_potential.hpp>
#include <input/environment.hpp>
#include <ions/geometry.hpp>
#include <operations/space.hpp>
#include <operations/laplacian.hpp>

#include <utils/profiling.hpp>

#include <future>
#include <list>
#include <unordered_map>

namespace inq {
namespace hamiltonian {
  template <class basis_type>
  class ks_hamiltonian {
		
  public:

		basis::field<basis::real_space, double> scalar_potential;
		exchange_operator exchange;
		
    ks_hamiltonian(const basis_type & basis, const ions::UnitCell & cell, const atomic_potential & pot, bool fourier_pseudo, const ions::geometry & geo,
									 const int num_hf_orbitals, const double exchange_coefficient, boost::mpi3::cartesian_communicator<2> const & comm):
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
					projectors_.emplace_back(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom]), geo.coordinates()[iatom], iatom);
					if(projectors_.back().empty()) projectors_.pop_back(); 
				}
			}

    }

		////////////////////////////////////////////////////////////////////////////////////////////

    ~ks_hamiltonian(){
			
			//we have to empty the map by hand, to preserve the order of the destruction
			while(not projectors_.empty()){
				projectors_.erase(projectors_.begin());
			}

    }

		////////////////////////////////////////////////////////////////////////////////////////////
		
		auto non_local_projection(const basis::field_set<basis::real_space, complex> & phi) const {

			CALI_CXX_MARK_FUNCTION;
			
			using proj_type = decltype(projectors_.cbegin()->project(phi));
			
			auto policy = std::launch::deferred;
			if(input::environment::threaded()) policy = std::launch::async;
				
			std::vector<std::future<proj_type>> projections;

			if(not non_local_in_fourier_) {
				for(auto it = projectors_.cbegin(); it != projectors_.cend(); ++it){
					projections.emplace_back(std::async(policy, [it, &phi]{ return it->project(phi);}));
				}
			}

			return projections;
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////		

		template <typename ProjType>
		void non_local_apply(ProjType && projections, basis::field_set<basis::real_space, complex> & vnlphi) const {

			if(non_local_in_fourier_) return;
				
			CALI_CXX_MARK_FUNCTION;
			
			auto projit = projections.begin();
			for(auto it = projectors_.cbegin(); it != projectors_.cend(); ++it){
				it->apply(projit->get(), vnlphi);
				++projit;
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
		void non_local(const basis::field_set<basis::fourier_space, complex> & phi, basis::field_set<basis::fourier_space, complex> & vnlphi) const {

			if(not non_local_in_fourier_) return;
			
			for(auto it = projectors_fourier_map_.cbegin(); it != projectors_fourier_map_.cend(); ++it){
				it->second(phi, vnlphi);
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////
		
		auto non_local(const basis::field_set<basis::real_space, complex> & phi) const {

			CALI_CXX_MARK_FUNCTION;
 
			if(non_local_in_fourier_) {

				auto phi_fs = operations::space::to_fourier(phi);
				basis::field_set<basis::fourier_space, complex> vnlphi_fs(phi_fs.skeleton());

				vnlphi_fs = 0.0;
				non_local(phi_fs, vnlphi_fs);
				return operations::space::to_real(vnlphi_fs);
					
			} else {

				auto && proj = non_local_projection(phi);
				
				basis::field_set<basis::real_space, complex> vnlphi(phi.skeleton());
				vnlphi = 0.0;

				non_local_apply(std::move(proj), vnlphi);
			
				return vnlphi;
							
			}
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////

    auto operator()(const basis::field_set<basis::real_space, complex> & phi) const{

			CALI_CXX_MARK_SCOPE("hamiltonian_real");

			auto && proj = non_local_projection(phi);
			
			auto phi_fs = operations::space::to_fourier(phi);
		
			auto hphi_fs = operations::laplacian(phi_fs);

			non_local(phi_fs, hphi_fs);
			
			auto hphi = operations::space::to_real(hphi_fs);

			hamiltonian::scalar_potential_add(scalar_potential, phi, hphi);
			exchange(phi, hphi);

			non_local_apply(std::move(proj), hphi);

			return hphi;
			
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		
    auto operator()(const basis::field_set<basis::fourier_space, complex> & phi) const{

			CALI_CXX_MARK_SCOPE("hamiltonian_fourier");
			
			auto phi_rs = operations::space::to_real(phi);

			auto && proj = non_local_projection(phi_rs);
			
			auto hphi_rs = hamiltonian::scalar_potential(scalar_potential, phi_rs);
		
			exchange(phi_rs, hphi_rs);

			non_local_apply(std::move(proj), hphi_rs);
			
			auto hphi = operations::space::to_fourier(hphi_rs);

			operations::laplacian_add(phi, hphi);
			non_local(phi, hphi);
			
			return hphi;
			
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
	using namespace Catch::literals;
  using math::vector3;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);	

  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
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
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
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
					
					double r2 = rs.r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-ww*r2);
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
						
						double xx = rs.rvector(ixg, iyg, izg)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(istg.value()*kk*xx), sin(istg.value()*kk*xx));
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
					
					double r2 = rs.r2(ixg, iyg, izg);
					ham.scalar_potential.cubic()[ix][iy][iz] = 0.5*ww*ww*r2;

					for(int ist = 0; ist < phi.local_set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-ww*r2);
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
