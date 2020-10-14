/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__SCF
#define INQ__INPUT__SCF

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <utils/merge_optional.hpp>

#include <nonstd/optional.hpp>
#include <cassert>

namespace inq {
namespace input {

  class scf {

  public:

    enum class scf_eigensolver { STEEPEST_DESCENT,
                                 CONJUGATE_GRADIENT,
				 DAVIDSON
    };

    enum class mixing_algo { LINEAR,
			     PULAY,
			     BROYDEN
    };

    static auto steepest_descent(){
      scf solver;
      solver.eigensolver_ = scf_eigensolver::STEEPEST_DESCENT;
      return solver;
    }

    static auto conjugate_gradient(){
      scf solver;
      solver.eigensolver_ = scf_eigensolver::CONJUGATE_GRADIENT;
      return solver;
    }

    static auto davidson(){
      scf solver;
      solver.eigensolver_ = scf_eigensolver::DAVIDSON;
      return solver;
    }
		
    auto eigensolver() const {
      return eigensolver_.value_or(scf_eigensolver::STEEPEST_DESCENT);
    }

    static auto mixing(double mixing_factor) {
      scf solver;
      solver.mixing_ = mixing_factor;
      return solver;
    }

    auto mixing() const {
      return mixing_.value_or(0.3);
    }

    enum class mix_field { DENSITY,
                           POTENTIAL
    };
    
    auto static density_mixing() {
      scf solver;
      solver.mix_field_ = mix_field::DENSITY;
      return solver;
    }

    auto static potential_mixing() {
      scf solver;
      solver.mix_field_ = mix_field::POTENTIAL;
      return solver;
    }

    auto mix_field_requested() const {
      return mix_field_.value_or(mix_field::POTENTIAL);
    }
    
    auto mix_density() const {
      return mix_field_requested() == mix_field::DENSITY;
    }

    auto mix_potential() const {
      return mix_field_requested() == mix_field::POTENTIAL;
    }

		auto static energy_tolerance(double etol) {
			scf solver;
      solver.energy_tol_ = etol;
      return solver;
    }
				
		auto energy_tolerance() const {
			return energy_tol_.value_or(1e-7);
		}
		
		auto static linear_mixing(){
			scf solver;
      solver.mixing_algo_ = mixing_algo::LINEAR;
      return solver;
		}

		auto static pulay_mixing(){
			scf solver;
      solver.mixing_algo_ = mixing_algo::PULAY;
      return solver;
		}

		auto static broyden_mixing(){
			scf solver;
      solver.mixing_algo_ = mixing_algo::BROYDEN;
      return solver;
		}
		
		auto mixing_algorithm() const {
			return mixing_algo_.value_or(mixing_algo::BROYDEN);
		}

		auto static silent(){
			scf solver;
      solver.verbose_ = false;
      return solver;
		}
		
		auto verbose_output() const {
			return verbose_.value_or(true);
		}
		
    friend auto operator|(const scf & solver1, const scf & solver2){
			using inq::utils::merge_optional;

			scf rsolver;
			rsolver.eigensolver_	= merge_optional(solver1.eigensolver_, solver2.eigensolver_);
			rsolver.mixing_	= merge_optional(solver1.mixing_, solver2.mixing_);
			rsolver.mix_field_	= merge_optional(solver1.mix_field_, solver2.mix_field_);
			rsolver.energy_tol_	= merge_optional(solver1.energy_tol_, solver2.energy_tol_);
			rsolver.mixing_algo_	= merge_optional(solver1.mixing_algo_, solver2.mixing_algo_);
			rsolver.verbose_	= merge_optional(solver1.verbose_, solver2.verbose_);			
			return rsolver;
		}
    
  private:

    nonstd::optional<scf_eigensolver> eigensolver_;
    nonstd::optional<double> mixing_;
    nonstd::optional<mix_field> mix_field_;
		nonstd::optional<double> energy_tol_;
		nonstd::optional<mixing_algo> mixing_algo_;
		nonstd::optional<bool> verbose_;
		
  };
    
}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_SCF_UNIT_TEST
#undef INQ_INPUT_SCF_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("class input::scf", "[input::scf]") {
  
	using namespace inq;
	using namespace Catch::literals;

	SECTION("Defaults"){

    input::scf solver;

    CHECK(solver.eigensolver() == input::scf::scf_eigensolver::STEEPEST_DESCENT);
    CHECK(solver.mixing() == 0.3_a);
    
  }

  SECTION("Composition"){

    auto solver = input::scf::conjugate_gradient() | input::scf::mixing(0.05);
    
    CHECK(solver.eigensolver() == input::scf::scf_eigensolver::CONJUGATE_GRADIENT);
    CHECK(solver.mixing() == 0.05_a);
    
  }

}

#endif
   
#endif
