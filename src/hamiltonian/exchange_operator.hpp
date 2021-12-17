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

				// Hartree-Fock exchange
				for(int ii = 0; ii < phi.fields().set_size(); ii++){
					for(int jj = 0; jj < phi.fields().set_size(); jj++){
						
						basis::field<basis::real_space, complex> rhoij(phi.fields().basis());

						//DATAOPERATIONS LOOP 1D
						for(long ipoint = 0; ipoint < phi.fields().basis().size(); ipoint++) rhoij.linear()[ipoint] = conj(hf_orbitals.matrix()[ipoint][jj])*phi.fields().matrix()[ipoint][ii];
						
						//OPTIMIZATION: this could be done in place
						auto potij = poisson_solver_(rhoij);
						
						//DATAOPERATIONS LOOP 1D
						for(long ipoint = 0; ipoint < phi.fields().basis().size(); ipoint++) {
							exxphi.fields().matrix()[ipoint][ii] -= 0.5*exchange_coefficient_*hf_occupations[jj]*hf_orbitals.matrix()[ipoint][jj]*potij.linear()[ipoint];
						}
						
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
#include <catch2/catch.hpp>
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
