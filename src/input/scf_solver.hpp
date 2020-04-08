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

#include <nonstd/optional.hpp>
#include <cassert>

namespace input {

  class scf_solver {

  public:

    enum class scf_eigensolver { STEEPEST_DESCENT,
                                 CONJUGATE_GRADIENT
    };

    scf_solver(){
		}

    static auto steepest_descent(){
      scf_solver inter;
      inter.eigensolver_ = scf_eigensolver::STEEPEST_DESCENT;
      return inter;
    }

    static auto conjugate_gradient(){
      scf_solver inter;
      inter.eigensolver_ = scf_eigensolver::CONJUGATE_GRADIENT;
      return inter;
    }
		
    auto eigensolver() const {
      return eigensolver_.value_or(scf_eigensolver::STEEPEST_DESCENT);
    }

  private:

    nonstd::optional<scf_eigensolver> eigensolver_;
    
  };
    
}

////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#endif
   
#endif
