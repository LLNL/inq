/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__SCF
#define INQ__INPUT__SCF

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <magnitude/energy.hpp>
#include <utils/merge_optional.hpp>

#include <optional>
#include <cassert>

namespace inq {
namespace input {

  class scf {

  public:

    enum class scf_eigensolver { STEEPEST_DESCENT };
    enum class mixing_algo { LINEAR, BROYDEN };

	private:
		
		std::optional<scf_eigensolver> eigensolver_;
    std::optional<double> mixing_;
		std::optional<double> energy_tol_;
		std::optional<mixing_algo> mixing_algo_;
		std::optional<bool> verbose_;
		std::optional<bool> subspace_diag_;
		std::optional<int> scf_steps_;
		std::optional<bool> calc_forces_;

  public:

    auto steepest_descent(){
      scf solver = *this;;
      solver.eigensolver_ = scf_eigensolver::STEEPEST_DESCENT;
      return solver;
    }

    auto eigensolver() const {
      return eigensolver_.value_or(scf_eigensolver::STEEPEST_DESCENT);
    }

    auto mixing(double mixing_factor) {
      scf solver = *this;;
      solver.mixing_ = mixing_factor;
      return solver;
    }

    auto mixing() const {
      return mixing_.value_or(0.3);
    }

		auto energy_tolerance(quantity<magnitude::energy> etol) {
			scf solver = *this;;
      solver.energy_tol_ = etol.in_atomic_units();
      return solver;
    }
				
		auto energy_tolerance() const {
			return energy_tol_.value_or(1e-5);
		}
		
		auto linear_mixing(){
			scf solver = *this;;
      solver.mixing_algo_ = mixing_algo::LINEAR;
      return solver;
		}

		auto broyden_mixing(){
			scf solver = *this;;
      solver.mixing_algo_ = mixing_algo::BROYDEN;
      return solver;
		}
		
		auto mixing_algorithm() const {
			return mixing_algo_.value_or(mixing_algo::BROYDEN);
		}

		auto silent(){
			scf solver = *this;;
      solver.verbose_ = false;
      return solver;
		}
		
		auto verbose_output() const {
			return verbose_.value_or(true);
		}

		auto no_subspace_diag() {
			scf solver = *this;;
      solver.subspace_diag_ = false;
      return solver;
		}
		
		auto subspace_diag() const {
			return subspace_diag_.value_or(true);
		}

		auto scf_steps(int val) {
			scf solver = *this;;
			solver.scf_steps_ = val;
      return solver;
		}

		auto scf_steps() const {
			return scf_steps_.value_or(200);
		}		

		auto calculate_forces() {
			scf solver = *this;;
      solver.calc_forces_ = true;
      return solver;
		}
		
		auto calc_forces() const {
			return calc_forces_.value_or(false);
		}
		
  };
}
}
#endif

#ifdef INQ_INPUT_SCF_UNIT_TEST
#undef INQ_INPUT_SCF_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;

	SECTION("Defaults"){

    input::scf solver;

    CHECK(solver.eigensolver() == input::scf::scf_eigensolver::STEEPEST_DESCENT);
    CHECK(solver.mixing() == 0.3_a);
  }

  SECTION("Composition"){

    auto solver = input::scf{}.calculate_forces().mixing(0.05);

		CHECK(solver.calc_forces());
    CHECK(solver.mixing() == 0.05_a);
  }
}
#endif
