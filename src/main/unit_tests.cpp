/* -*- indent-tabs-mode: t -*- */

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

#define UNIT_TEST

#include <catch2/catch.hpp>

#include <pseudopod/spline.hpp>
#include <pseudopod/erf_range_separation.hpp>
#include <pseudopod/element.hpp>
#include <pseudopod/pseudopotential.hpp>
#include <pseudopod/spherical_harmonic.hpp>

#include <utils/distribution.hpp>

#include <gpu/run.hpp>

#include <input/basis.hpp>
#include <input/cell.hpp>
#include <input/species.hpp>

#include <config/path.hpp>

#include <math/d3vector.hpp>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <ions/interaction.hpp>
#include <ions/periodic_replicas.hpp>

#include <basis/grid.hpp>
#include <basis/spherical_grid.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>

#include <states/ks_states.hpp>

#include <hamiltonian/projector.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>

#include <operations/add.hpp>
#include <operations/calculate_density.hpp>
#include <operations/sum.hpp>
#include <operations/overlap.hpp>
#include <operations/orthogonalization.hpp>
#include <operations/diagonalize.hpp>
#include <operations/shift.hpp>
#include <operations/randomize.hpp>

#include <operations/matrix_operator.hpp>

#include <solvers/poisson.hpp>
#include <solvers/linear_mixer.hpp>
#include <solvers/pulay_mixer.hpp>
#include <solvers/linear.hpp>
#include <solvers/steepest_descent.hpp>
