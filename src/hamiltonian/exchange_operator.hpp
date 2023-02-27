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

#include <basis/real_space.hpp>
#include <operations/overlap.hpp>
#include <operations/overlap_diagonal.hpp>
#include <operations/rotate.hpp>
#include <parallel/arbitrary_partition.hpp>
#include <parallel/array_iterator.hpp>
#include <solvers/cholesky.hpp>
#include <solvers/poisson.hpp>
#include <states/orbital_set.hpp>

#include <optional>

namespace inq {
namespace hamiltonian {
  class exchange_operator {
		
  public:

		exchange_operator(double const exchange_coefficient, bool const use_ace):
			exchange_coefficient_(exchange_coefficient),
			use_ace_(use_ace){
		}

		//////////////////////////////////////////////////////////////////////////////////

		template <class ElectronsType>
		double update(ElectronsType const & el){
			if(not enabled()) return 0.0;

			CALI_CXX_MARK_SCOPE("exchage_operator::update");

			auto part = parallel::arbitrary_partition(el.max_local_set_size()*el.lot_size(), el.states_comm_);

			occupations_.reextent(part.local_size());
			kpoints_.reextent(part.local_size());
			
			if(not orbitals_.has_value()) orbitals_.emplace(el.states_basis_, part, el.states_basis_comm_);
			
			auto iphi = 0;
			auto ist = 0;
			for(auto & phi : el.lot()){

				occupations_({ist, ist + phi.local_set_size()}) = el.occupations()[iphi];
				kpoints_({ist, ist + phi.local_set_size()}).fill(phi.kpoint());
				orbitals_->matrix()({0, phi.basis().local_size()}, {ist, ist + phi.local_set_size()}) = phi.matrix();

				iphi++;
				ist += phi.local_set_size();			
			}

			if(not ace_orbitals_.has_value()) ace_orbitals_.emplace(el.states_basis_, part, el.states_basis_comm_);
			
			iphi = 0;
			ist = 0;
			for(auto & phi : el.lot()){

				auto exxphi = direct(phi, -1.0);
				ace_orbitals_->matrix()({0, phi.basis().local_size()}, {ist, ist + phi.local_set_size()}) = exxphi.matrix();
				
				iphi++;
				ist += phi.local_set_size();
			}

			auto exx_matrix = operations::overlap(*ace_orbitals_, *orbitals_);
			double energy = -0.5*real(operations::sum_product(occupations_, exx_matrix.diagonal()));
			el.lot_states_comm_.all_reduce_n(&energy, 1);
				
			solvers::cholesky(exx_matrix.array());
			operations::rotate_trs(exx_matrix, *ace_orbitals_);

			return energy;
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		auto direct(const states::orbital_set<basis::real_space, complex> & phi, double scale = 1.0) const {
			states::orbital_set<basis::real_space, complex> exxphi(phi.skeleton());
			exxphi.fill(0.0);
			direct(phi, exxphi, scale);
			return exxphi;
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		template <class HFType, class HFOccType, class KptType, class PhiType, class ExxphiType>
		void block_exchange(double factor, HFType const & hf, HFOccType const & hfocc, KptType const & kpt, PhiType const & phi, ExxphiType & exxphi) const {

			auto nst = phi.local_set_size();
			auto nhf = (~hf).size();
			basis::field_set<basis::real_space, complex> rhoij(phi.basis(), nst);
			
			for(int jj = 0; jj < nhf; jj++){
				
				{ CALI_CXX_MARK_SCOPE("hartree_fock_exchange_gen_dens");
					gpu::run(nst, phi.basis().local_size(),
									 [rho = begin(rhoij.matrix()), hfo = begin(hf), ph = begin(phi.matrix()), jj] GPU_LAMBDA (auto ist, auto ipoint){ 
										 rho[ipoint][ist] = conj(hfo[ipoint][jj])*ph[ipoint][ist];
									 });
				}

				poisson_solver_.in_place(rhoij, -phi.kpoint() + kpt[jj]);
				
				{ CALI_CXX_MARK_SCOPE("hartree_fock_exchange_mul_pot");
					gpu::run(nst, exxphi.basis().local_size(),
									 [pot = begin(rhoij.matrix()), hfo = begin(hf), exph = begin(exxphi.matrix()), occ = begin(hfocc), jj, factor]
									 GPU_LAMBDA (auto ist, auto ipoint){
										 exph[ipoint][ist] += factor*occ[jj]*hfo[ipoint][jj]*pot[ipoint][ist];
									 });
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		void direct(const states::orbital_set<basis::real_space, complex> & phi, states::orbital_set<basis::real_space, complex> & exxphi, double scale = 1.0) const {
			if(not enabled()) return;
			
			CALI_CXX_MARK_SCOPE("hartree_fock_exchange");
			
			double factor = -0.5*scale*exchange_coefficient_;

			if(not orbitals_->set_part().parallel()){
				block_exchange(factor, orbitals_->matrix(), occupations_, kpoints_, phi, exxphi);
			} else {

				auto occ_it = parallel::array_iterator(orbitals_->set_part(), orbitals_->set_comm(), occupations_);
				auto kpt_it = parallel::array_iterator(orbitals_->set_part(), orbitals_->set_comm(), kpoints_);
				for(auto hfo_it = orbitals_->par_set_begin(); hfo_it != orbitals_->par_set_end(); ++hfo_it){
					block_exchange(factor, hfo_it.matrix(), *occ_it, *kpt_it, phi, exxphi);
					++occ_it;
					++kpt_it;
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		auto ace(const states::orbital_set<basis::real_space, complex> & phi) const {
			states::orbital_set<basis::real_space, complex> exxphi(phi.skeleton());
			exxphi.fill(0.0);
			ace(phi, exxphi);
			return exxphi;
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		auto operator()(const states::orbital_set<basis::real_space, complex> & phi) const {
			states::orbital_set<basis::real_space, complex> exxphi(phi.skeleton());
			exxphi.fill(0.0);
			operator()(phi, exxphi);
			return exxphi;
		}

		//////////////////////////////////////////////////////////////////////////////////

		void operator()(const states::orbital_set<basis::real_space, complex> & phi, states::orbital_set<basis::real_space, complex> & exxphi) const {
			if(not enabled()) return;

			if(use_ace_) ace(phi, exxphi);
			else direct(phi, exxphi);
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		void ace(const states::orbital_set<basis::real_space, complex> & phi, states::orbital_set<basis::real_space, complex> & exxphi) const {			
			if(not enabled()) return;
			namespace blas = boost::multi::blas;

			auto olap = operations::overlap(*ace_orbitals_, phi);
			operations::rotate(olap, *ace_orbitals_, exxphi, -1.0, 1.0);
		}

		//////////////////////////////////////////////////////////////////////////////////
		
		bool enabled() const {
			return fabs(exchange_coefficient_) > 1.0e-14;
		}

		//////////////////////////////////////////////////////////////////////////////////

	private:
		math::array<double, 1> occupations_;
		math::array<vector3<double, covariant>, 1> kpoints_;		
		std::optional<basis::field_set<basis::real_space, complex, parallel::arbitrary_partition>> orbitals_;
		std::optional<basis::field_set<basis::real_space, complex, parallel::arbitrary_partition>> ace_orbitals_;
		solvers::poisson poisson_solver_;
		double exchange_coefficient_;
		bool use_ace_;
		
  };

}
}
#endif

#ifdef INQ_HAMILTONIAN_EXCHANGE_OPERATOR_UNIT_TEST
#undef INQ_HAMILTONIAN_EXCHANGE_OPERATOR_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif
