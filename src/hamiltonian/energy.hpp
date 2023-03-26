/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__ENERGY
#define INQ__HAMILTONIAN__ENERGY

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

#include <operations/shift.hpp>

#include <tinyformat/tinyformat.h>

namespace inq {
namespace hamiltonian {

	class energy {

		double ion_ = 0.0;
		double eigenvalues_ = 0.0;
		double external_ = 0.0;
		double nonlocal_ = 0.0;
		double hartree_ = 0.0;
		double xc_ = 0.0;
		double nvxc_ = 0.0;
		double hf_exchange_ = 0.0;

	public:
		
		energy() = default;

		template <typename HamType, typename ElType>
		auto calculate(HamType const & ham, ElType & el) {

			CALI_CXX_MARK_SCOPE("energy::calculate");

			auto normres = math::array<complex, 2>({el.kpin().size(), el.max_local_set_size()});
			
			eigenvalues_ = 0.0;
			nonlocal_ = 0.0;
			hf_exchange_ = 0.0;
			
			int iphi = 0;
			for(auto & phi : el.kpin()){
				
				auto residual = ham(phi);
				el.eigenvalues()[iphi] = operations::overlap_diagonal_normalized(residual, phi, operations::real_part{});
				operations::shift(-1.0, el.eigenvalues()[iphi], phi, residual);
				
				normres[iphi] = operations::overlap_diagonal(residual);
				auto nl_me = operations::overlap_diagonal_normalized(ham.non_local(phi), phi);
				auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(phi), phi);
				
				auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
				
				eigenvalues_ += operations::sum(el.occupations()[iphi], el.eigenvalues()[iphi], energy_term);
				nonlocal_ += operations::sum(el.occupations()[iphi], nl_me, energy_term);
				hf_exchange_ += 0.5*operations::sum(el.occupations()[iphi], exchange_me, energy_term);
				
				iphi++;
			}

			el.kpin_states_comm().all_reduce_n(&eigenvalues_, 1);
			el.kpin_states_comm().all_reduce_n(&nonlocal_, 1);
			el.kpin_states_comm().all_reduce_n(&hf_exchange_, 1);

			return normres;
		}
	
		auto kinetic() const {
			return eigenvalues_ - 2.0*hartree_ - nvxc_ - 2.0*hf_exchange_ - external_ - nonlocal_;
		}
		
		auto total() const {
			return kinetic() + hartree_ + external_ + nonlocal_ + xc_ + hf_exchange_ + ion_;
		}

		auto & eigenvalues() const {
			return eigenvalues_;
		}
		
		void eigenvalues(double const & val) {
			eigenvalues_ = val;
		}

		auto & hartree() const {
			return hartree_;
		}

		void hartree(double const & val) {
			hartree_ = val;
		}

		auto & external() const {
			return external_;
		}

		void external(double const & val) {
			external_ = val;
		}

		auto & nonlocal() const {
			return nonlocal_;
		}

		void nonlocal(double const & val) {
			nonlocal_ = val;
		}

		auto & xc() const {
			return xc_;
		}

		void xc(double const & val) {
			xc_ = val;
		}

		auto & nvxc() const {
			return nvxc_;
		}

		void nvxc(double const & val) {
			nvxc_ = val;
		}

		auto & hf_exchange() const {
			return hf_exchange_;
		}

		void hf_exchange(double const & val) {
			hf_exchange_ = val;
		}

		auto & ion() const {
			return ion_;
		}

		void ion(double const & val) {
			ion_ = val;
		}
		
		template <class out_type>
		void print(out_type & out) const {

			tfm::format(out, "\n");
			tfm::format(out, "  total          = %20.12f\n", total());			
			tfm::format(out, "  kinetic        = %20.12f\n", kinetic());
			tfm::format(out, "  eigenvalues    = %20.12f\n", eigenvalues_);
			tfm::format(out, "  hartree        = %20.12f\n", hartree_);
			tfm::format(out, "  external       = %20.12f\n", external_);
			tfm::format(out, "  nonlocal       = %20.12f\n", nonlocal_);
			tfm::format(out, "  xc             = %20.12f\n", xc_);
			tfm::format(out, "  intnvxc        = %20.12f\n", nvxc_);
			tfm::format(out, "  HF exchange    = %20.12f\n", hf_exchange_);
			tfm::format(out, "  ion            = %20.12f\n", ion_);
			tfm::format(out, "\n");

		}
		
		template<class OStream>
		friend OStream& operator<<(OStream& os, energy const& self){
			self.print(os);
			return os;
		}
		
	};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ENERGY_UNIT_TEST
#undef INQ_HAMILTONIAN_ENERGY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

