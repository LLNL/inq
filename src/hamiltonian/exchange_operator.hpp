/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__EXCHANGE_OPERATOR
#define INQ__HAMILTONIAN__EXCHANGE_OPERATOR

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

#include <states/orbital_set.hpp>
#include <basis/real_space.hpp>
#include <solvers/poisson.hpp>

namespace inq {
namespace hamiltonian {
  class exchange_operator {
		
  public:

		exchange_operator(const basis::real_space & basis, const int num_hf_orbitals, const double exchange_coefficient, boost::mpi3::cartesian_communicator<2> comm):
			hf_occupations(num_hf_orbitals),
			hf_orbitals(basis, num_hf_orbitals, std::move(comm)),
			exchange_coefficient_(exchange_coefficient){

		}

		void operator()(const states::orbital_set<basis::real_space, complex> & phi, states::orbital_set<basis::real_space, complex> & exxphi) const {

			if(exchange_coefficient_ != 0.0){

				CALI_CXX_MARK_SCOPE("hartree_fock_exchange");

				basis::field_set<basis::real_space, complex> rhoij(phi.fields().basis(), phi.fields().set_size());
				
				for(int jj = 0; jj < hf_orbitals.set_size(); jj++){

					{ CALI_CXX_MARK_SCOPE("hartree_fock_exchange_gen_dens");
						gpu::run(phi.fields().set_size(), phi.fields().basis().size(),
										 [rho = begin(rhoij.matrix()), hfo = begin(hf_orbitals.matrix()), ph = begin(phi.fields().matrix()), jj] GPU_LAMBDA (auto ist, auto ipoint){ 
											 rho[ipoint][ist] = conj(hfo[ipoint][jj])*ph[ipoint][ist];
										 });
					}
					
					poisson_solver_.in_place(rhoij);

					{ CALI_CXX_MARK_SCOPE("hartree_fock_exchange_mul_pot");
						gpu::run(phi.fields().set_size(), phi.fields().basis().size(),
										 [pot = begin(rhoij.matrix()), hfo = begin(hf_orbitals.matrix()), exph = begin(exxphi.fields().matrix()), occ = begin(hf_occupations), jj, coeff = exchange_coefficient_]
										 GPU_LAMBDA (auto ist, auto ipoint){
											 exph[ipoint][ist] -= 0.5*coeff*occ[jj]*hfo[ipoint][jj]*pot[ipoint][ist];
										 });
					}
					
				}
				
			}
		}

		auto operator()(const states::orbital_set<basis::real_space, complex> & phi) const {
			states::orbital_set<basis::real_space, complex> exxphi(phi.skeleton());
			exxphi.fields() = 0.0;
			operator()(phi, exxphi);
			return exxphi;
		}

		auto enabled() const {
			return exchange_coefficient_ != 0.0;
		}
		
		math::array<double, 1> hf_occupations;
		states::orbital_set<basis::real_space, complex> hf_orbitals;

	private:
		solvers::poisson poisson_solver_;
		double exchange_coefficient_;

  };

}
}

#ifdef INQ_HAMILTONIAN_EXCHANGE_OPERATOR_UNIT_TEST
#undef INQ_HAMILTONIAN_EXCHANGE_OPERATOR_UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::exchange", "[hamiltonian::exchange]"){

	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;
  /*
  auto ecut = 20.0_Ha;
  double ll = 10.0;
	*/
	/*
	ions::geometry geo;
  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

  states::orbital_set<basis::real_space, complex> phi(rs, st.num_states());
	states::orbital_set<basis::real_space, complex> hphi(rs, st.num_states());
	
	hamiltonian::exchange<basis::real_space> ham(rs, cell, pot, geo, st.num_states(), 0.0);
	*/
}

#endif

#endif
