/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__ENERGY
#define INQ__HAMILTONIAN__ENERGY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <operations/shift.hpp>

#include <tinyformat/tinyformat.h>

namespace inq {
namespace hamiltonian {

	class energy {

		double ion_ = 0.0;
		double ion_kinetic_ = 0.0;		
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

			auto normres = gpu::array<complex, 2>({static_cast<gpu::array<complex, 2>::size_type>(el.kpin().size()), el.max_local_set_size()});
			
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
				
				auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
				
				eigenvalues_ += operations::sum(el.occupations()[iphi], el.eigenvalues()[iphi], energy_term);
				nonlocal_ += operations::sum(el.occupations()[iphi], nl_me, energy_term);

				if(ham.exchange.enabled()){
					CALI_CXX_MARK_SCOPE("energy::calculate::exchange");

					auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(phi), phi);
					hf_exchange_ += 0.5*operations::sum(el.occupations()[iphi], exchange_me, energy_term);
				}

				iphi++;
			}

			if(el.kpin_states_comm().size() > 1){	
				CALI_CXX_MARK_SCOPE("energy::calculate::reduce");

				el.kpin_states_comm().all_reduce_n(&eigenvalues_, 1);
				el.kpin_states_comm().all_reduce_n(&nonlocal_, 1);
				el.kpin_states_comm().all_reduce_n(&hf_exchange_, 1);
			}

			return normres;
		}
	
		auto kinetic() const {
			return eigenvalues_ - 2.0*hartree_ - nvxc_ - 2.0*hf_exchange_ - external_ - nonlocal_;
		}
		
		auto total() const {
			return kinetic() + hartree_ + external_ + nonlocal_ + xc_ + hf_exchange_ + ion_ + ion_kinetic_;
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

		auto & ion_kinetic() const {
			return ion_kinetic_;
		}

		void ion_kinetic(double const & val) {
			ion_kinetic_ = val;
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

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

