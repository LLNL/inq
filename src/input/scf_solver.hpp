/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__SCF_SOLVER
#define INPUT__SCF_SOLVER

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

namespace input {

  class scf_solver {

  public:

    scf_solver(){
		}

    enum class scf_eigensolver { STEEPEST_DESCENT,
                                 CONJUGATE_GRADIENT
    };

    static auto steepest_descent(){
      scf_solver solver;
      solver.eigensolver_ = scf_eigensolver::STEEPEST_DESCENT;
      return solver;
    }

    static auto conjugate_gradient(){
      scf_solver solver;
      solver.eigensolver_ = scf_eigensolver::CONJUGATE_GRADIENT;
      return solver;
    }
		
    auto eigensolver() const {
      return eigensolver_.value_or(scf_eigensolver::STEEPEST_DESCENT);
    }

    static auto mixing(double mixing_factor) {
      scf_solver solver;
      solver.mixing_ = mixing_factor;
      return solver;
    }

    auto mixing() const {
      return mixing_.value_or(0.3);
    }
    
    friend auto operator|(const scf_solver & solver1, const scf_solver & solver2){
			using utils::merge_optional;

			scf_solver rsolver;
			rsolver.eigensolver_	= merge_optional(solver1.eigensolver_, solver2.eigensolver_);
			rsolver.mixing_	= merge_optional(solver1.mixing_, solver2.mixing_);
			return rsolver;
		}
    
  private:

    nonstd::optional<scf_eigensolver> eigensolver_;
    nonstd::optional<double> mixing_;
    
  };
    
}

////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class input::scf_solver", "[input::scf_solver]") {
  
  using namespace Catch::literals;

	SECTION("Defaults"){

    input::scf_solver solver;

    REQUIRE(solver.eigensolver() == input::scf_solver::scf_eigensolver::STEEPEST_DESCENT);
    REQUIRE(solver.mixing() == 0.3_a);
    
  }

  SECTION("Composition"){

    auto solver = input::scf_solver::conjugate_gradient() | input::scf_solver::mixing(0.05);
    
    REQUIRE(solver.eigensolver() == input::scf_solver::scf_eigensolver::CONJUGATE_GRADIENT);
    REQUIRE(solver.mixing() == 0.05_a);
    
  }

}

#endif
   
#endif
