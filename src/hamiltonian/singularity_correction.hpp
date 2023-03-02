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

  singularity_correction(ions::unit_cell const & cell, ions::brillouin const & bzone):
    fk_(bzone.size())
  {
    
    for(int ik = 0; ik < bzone.size(); ik++){
      
      fk_[ik] = 0.0;
      for(int ik2 = 0; ik2 < bzone.size(); ik2++){
        auto qpoint = bzone.kpoint(ik) - bzone.kpoint(ik2);
        if(cell.metric().norm(qpoint) < 1e-6) continue;
        fk_[ik] += bzone.kpoint_weight(ik2)*auxiliary(cell, qpoint);
      }
      fk_[ik] *= 4.0*M_PI/cell.volume();
    }
    
    auto const nsteps = 7;
    auto const nk = 60;
    
    fzero_ = 0.0;
    auto length = 1.0;
    auto kvol_element = pow(2.0*M_PI/(2.0*nk + 1.0), 3)/cell.volume();
    
    for(int istep = 0; istep < nsteps; istep++){

      for(auto ikx = 0; ikx <= nk; ikx++){
        for(auto iky = -nk; iky <= nk; iky++){
          for(auto ikz = -nk; ikz <= nk; ikz++){          
            if(fabs(ikx) <= nk/3 and fabs(iky) <= nk/3 and fabs(ikz) <= nk/3) continue;
            
            auto qpoint = 2.0*M_PI*vector3<int, covariant>(ikx, iky, ikz)*length/(2.0*nk);
            fzero_ += kvol_element*auxiliary(cell, qpoint);
          }
        }
      }
      if(istep < nsteps - 1){
        length /= 3.0;
        kvol_element /= 27.0;
      }
    }

    fzero_ *= 8.0*M_PI/pow(2.0*M_PI, 3);
    fzero_ += 7.7955541794415*pow(cell.volume(), 2.0/3.0)/M_PI/cell.volume()*length;
  }

  auto fk(int ik) const {
    return fk_[ik];
  }

  auto fzero() const {
    return fzero_;
  }  
  
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
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

  SECTION("Auxiliary function cubic"){
    auto aa = 10.18_b;

		auto box = systems::box::lattice({aa, 0.0_b, 0.0_b}, {0.0_b, aa, 0.0_b}, {0.0_b, 0.0_b, aa}).cutoff_energy(35.0_Ha);
		auto ions = systems::ions(box);
    auto const & cell = ions.cell();    

    auto bzone = ions::brillouin(ions, input::kpoints::grid({2, 2, 2}));
        
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{0.0, -0.5, 0.0}) == 25.9081_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{8.3333333333333332E-003, 7.4999999999999997E-002, 0.26666666666666666}) == 42.650855183181122_a);
    CHECK(hamiltonian::singularity_correction::auxiliary(cell, 2.0*M_PI*vector3<double, covariant>{0.11666666666666667, 0.20000000000000001, 0.21666666666666667}) == 29.780683447124286_a);    
    
    auto sing = hamiltonian::singularity_correction(cell, bzone);

    CHECK(sing.fzero() == 0.30983869660201141_a);

    CHECK(sing.fk(0) == 0.18644848345224296_a);
    CHECK(sing.fk(1) == 0.18644848345224296_a);
    CHECK(sing.fk(2) == 0.18644848345224296_a);
    CHECK(sing.fk(3) == 0.18644848345224296_a);
    CHECK(sing.fk(4) == 0.18644848345224296_a);
    CHECK(sing.fk(5) == 0.18644848345224296_a);
    CHECK(sing.fk(6) == 0.18644848345224296_a);
    CHECK(sing.fk(7) == 0.18644848345224296_a);
    
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
