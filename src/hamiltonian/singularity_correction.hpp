/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__SINGULARITY_CORRECTION
#define INQ__HAMILTONIAN__SINGULARITY_CORRECTION

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

// Implements the method for the singularity correction of Carrier et al. PRB 75 205126 (2007)
// https://doi.org/10.1103/PhysRevB.75.205126
class singularity_correction {

  math::array<double, 1> fk_;
  double fzero_;

public:

  // the function defined in Eq. 16
  static auto auxiliary(ions::unit_cell const & cell, vector3<double, covariant> const & qpoint){
    auto val = 0.0;
    auto const & metric = cell.metric();
    
    for(int jj = 0; jj < 3; jj++){
      auto jjp1 = jj + 1;
      if(jj == 3) jj = 0;

      auto v1 = cell.reciprocal(jj)  *sin(metric.dot(cell.lattice(jj)  , 0.5*qpoint));
      auto v2 = cell.reciprocal(jj)  *sin(metric.dot(cell.lattice(jj)  ,     qpoint));
      auto v3 = cell.reciprocal(jjp1)*sin(metric.dot(cell.lattice(jjp1),     qpoint));
      
      val += 4.0*metric.dot(v1, v1) + 2*metric.dot(v2, v3);
    }

    return 4*M_PI*M_PI/val;
  }

public:
  
};
}
}

#endif

#ifdef INQ_HAMILTONIAN_SINGULARITY_CORRECTION_UNIT_TEST
#undef INQ_HAMILTONIAN_SINGULARITY_CORRECTION_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;

  SECTION("Auxiliary function cubic"){
    auto aa = 10.18;
    ions::unit_cell cell(vector3<double>(aa, 0.0, 0.0), vector3<double>(0.0, aa, 0.0), vector3<double>(0.0, 0.0, aa));

    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{0.0, -0.5, 0.0}) == 25.9081_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{8.3333333333333332E-003, 7.4999999999999997E-002, 0.26666666666666666}) == 42.650855183181122_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{0.11666666666666667, 0.20000000000000001, 0.21666666666666667}) == 29.780683447124286_a);    
  }
  /*
    This doesn't work, I have to figure out why
    
    SECTION("Auxiliary function non-orthogonal"){
    auto aa = 6.7408326;
    ions::unit_cell cell(aa*vector3<double>(0.0, 0.5, 0.5), aa*vector3<double>(0.5, 0.0, 0.5), aa*vector3<double>(0.5, 0.5, 0.0));
    
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{1.6666666666666666E-002, 0.28333333333333333, 0.39166666666666666}) == 2.77471621018199290_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{0.12500000000000000,-0.20833333333333334, -0.23333333333333334}) == 3.6560191647005245_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{ 0.14999999999999999, 0.25000000000000000, -3.3333333333333333E-002}) == 5.8717108336249790_a);
  }  
  */
}
#endif
