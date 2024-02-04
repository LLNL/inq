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

			auto energy_term = [](auto occ, auto ev){ return occ*real(ev); };
			
			eigenvalues_ = 0.0;
			nonlocal_ = 0.0;
			hf_exchange_ = 0.0;
			
			int iphi = 0;
			for(auto & phi : el.kpin()){

				{
					CALI_CXX_MARK_SCOPE("energy::calculate::eigenvalues");
					auto residual = ham(phi);
					el.eigenvalues()[iphi] = operations::overlap_diagonal_normalized(residual, phi, operations::real_part{});
					operations::shift(-1.0, el.eigenvalues()[iphi], phi, residual);
					normres[iphi] = operations::overlap_diagonal(residual);
					eigenvalues_ += operations::sum(el.occupations()[iphi], el.eigenvalues()[iphi], energy_term);
				}

				{
					CALI_CXX_MARK_SCOPE("energy::calculate::nonlocal");
					auto nl_me = operations::overlap_diagonal_normalized(ham.non_local(phi), phi);
					nonlocal_ += operations::sum(el.occupations()[iphi], nl_me, energy_term);
				}

				if(ham.exchange.enabled()){
					CALI_CXX_MARK_SCOPE("energy::calculate::exchange");
					auto exchange_me = operations::overlap_diagonal_normalized(ham.exchange(phi), phi);
					hf_exchange_ += 0.5*operations::sum(el.occupations()[iphi], exchange_me, energy_term);
				}

				iphi++;
			}

			if(el.kpin_states_comm().size() > 1){	
				CALI_CXX_MARK_SCOPE("energy::calculate::reduce");

				double red[3] = {eigenvalues_, nonlocal_, hf_exchange_};
				el.kpin_states_comm().all_reduce_n(red, 3);
				eigenvalues_ = red[0];
				nonlocal_    = red[1];
				hf_exchange_ = red[2];
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

		void save(parallel::communicator & comm, std::string const & dirname) const {
			auto error_message = "INQ error: Cannot save the energy to directory '" + dirname + "'.";
			
			comm.barrier();

			auto exception_happened = true;
			if(comm.root()) {
			
				try { std::filesystem::create_directories(dirname); }
				catch(...) {
					comm.broadcast_value(exception_happened);
					throw std::runtime_error(error_message);
				}

				utils::save_value(comm, dirname + "/ion",          ion_,         error_message);
				utils::save_value(comm, dirname + "/ion_kinetic",  ion_kinetic_, error_message);
				utils::save_value(comm, dirname + "/eigenvalues",  eigenvalues_, error_message);
				utils::save_value(comm, dirname + "/external",     external_,    error_message);
				utils::save_value(comm, dirname + "/nonlocal",     nonlocal_,    error_message);
				utils::save_value(comm, dirname + "/hartree",      hartree_,     error_message);
				utils::save_value(comm, dirname + "/xc",           xc_,          error_message);
				utils::save_value(comm, dirname + "/nvxc",         nvxc_,        error_message);
				utils::save_value(comm, dirname + "/hf_exchange_", hf_exchange_, error_message);
				
				exception_happened = false;
				comm.broadcast_value(exception_happened);
			
			} else {
				comm.broadcast_value(exception_happened);
				if(exception_happened) throw std::runtime_error(error_message);
			}
		
			comm.barrier();
			
		}
		
		template<class OStream>
		friend OStream & operator<<(OStream & out, energy const & self){

			tfm::format(out, "\n");
			tfm::format(out, "  total          = %20.12f\n", self.total());
			tfm::format(out, "  kinetic        = %20.12f\n", self.kinetic());
			tfm::format(out, "  eigenvalues    = %20.12f\n", self.eigenvalues_);
			tfm::format(out, "  hartree        = %20.12f\n", self.hartree());
			tfm::format(out, "  external       = %20.12f\n", self.external());
			tfm::format(out, "  nonlocal       = %20.12f\n", self.nonlocal());
			tfm::format(out, "  xc             = %20.12f\n", self.xc());
			tfm::format(out, "  intnvxc        = %20.12f\n", self.nvxc());
			tfm::format(out, "  HF exchange    = %20.12f\n", self.hf_exchange());
			tfm::format(out, "  ion            = %20.12f\n", self.ion());
			tfm::format(out, "\n");

			return out;
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

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	hamiltonian::energy en;

	en.ion(1.0);
	en.ion_kinetic(2.0);
	en.eigenvalues(3.0);
	en.external(4.0);
	en.nonlocal(5.0);
	en.hartree(6.0);
	en.xc(7.0);
	en.nvxc(8.0);
	en.hf_exchange(10.0);

	CHECK(en.ion() == 1.0);
	CHECK(en.ion_kinetic() == 2.0);
	CHECK(en.eigenvalues() == 3.0);
	CHECK(en.external() == 4.0);
	CHECK(en.nonlocal() == 5.0);
	CHECK(en.hartree() == 6.0);
	CHECK(en.xc() == 7.0);
	CHECK(en.nvxc() == 8.0);
	CHECK(en.hf_exchange() == 10.0);
	
	en.save(comm, "save_energy");
	
}
#endif

