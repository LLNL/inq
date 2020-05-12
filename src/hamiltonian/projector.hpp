/* -*- indent-tabs-mode: t -*- */

#ifndef HAMILTONIAN_PROJECTOR
#define HAMILTONIAN_PROJECTOR

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo Correa.

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

#include <pseudopod/spherical_harmonic.hpp>
#include <pseudopod/pseudopotential.hpp>

#include <math/array.hpp>
#include <math/vec3d.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#ifdef HAVE_CUDA
#include <multi/adaptors/blas/cuda.hpp> // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>

namespace hamiltonian {

  class projector {

  public:
    projector(const basis::real_space & basis, const ions::UnitCell & cell, pseudo::pseudopotential ps, math::vec3d atom_position):
      sphere_(basis, cell, atom_position, ps.projector_radius()),
      nproj_(ps.num_projectors_lm()),
      matrix_({nproj_, sphere_.size()}),
			kb_coeff_(nproj_) {

			std::vector<double> grid(sphere_.size()), proj(sphere_.size());

			int iproj_lm = 0;
      for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				// interpolate the value of the radial part of the projectors to the sphere points
				ps.projector(iproj_l).value(sphere_.size(), sphere_.distance(), proj);
				
				int l = ps.projector_l(iproj_l);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
						matrix_[iproj_lm][ipoint] = proj[ipoint]*math::spherical_harmonic(l, m, sphere_.point_pos()[ipoint]);
					}
					kb_coeff_[iproj_lm]	= ps.kb_coeff(iproj_l); 
					iproj_lm++;
				}
				
      }

			assert(iproj_lm == ps.num_projectors_lm());
			
    }

    template <class field_set_type>
    void operator()(const field_set_type & phi, field_set_type & vnlphi) const {

			if(nproj_ == 0) return;
			
			using boost::multi::blas::gemm;
			using boost::multi::blas::transposed;
				
			auto sphere_phi = sphere_.gather(phi.cubic());
			
			//DATAOPERATIONS BLAS
			auto projections = gemm(sphere_.volume_element(), matrix_, sphere_phi);
			
			//DATAOPERATIONS GPU::RUN 2D
			gpu::run(phi.set_size(), nproj_,
							 [proj = begin(projections), coeff = begin(kb_coeff_)]
							 GPU_LAMBDA (auto ist, auto iproj){
								 proj[iproj][ist] = proj[iproj][ist]*coeff[iproj];
							 });
			
			//DATAOPERATIONS BLAS
			sphere_phi = gemm(transposed(matrix_), projections);
			
			sphere_.scatter_add(sphere_phi, vnlphi.cubic());
			
    }

    int num_projectors() const {
      return nproj_;
    }
		
    auto kb_coeff(int iproj){
      return kb_coeff_[iproj];
    }
		
  private:

    basis::spherical_grid sphere_;
    int nproj_;
		//OPTIMIZATION: make this matrix real
    math::array<complex, 2> matrix_;
		math::array<double, 1> kb_coeff_;
    
  };
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class hamiltonian::projector", "[hamiltonian::projector]") {
  
  using namespace Catch::literals;
  using math::vec3d;
	
	const math::erf_range_separation sep(0.625);
	

  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	pseudo::pseudopotential ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector proj(rs, cell, ps, vec3d(0.0, 0.0, 0.0));

	CHECK(proj.num_projectors() == 8);
	
	CHECK(proj.kb_coeff(0) ==  7.494508815_a);
	CHECK(proj.kb_coeff(1) ==  0.6363049519_a);
	CHECK(proj.kb_coeff(2) == -4.2939052122_a);
	CHECK(proj.kb_coeff(3) == -4.2939052122_a);
	CHECK(proj.kb_coeff(4) == -4.2939052122_a);
	CHECK(proj.kb_coeff(5) == -1.0069878791_a);
	CHECK(proj.kb_coeff(6) == -1.0069878791_a);
	CHECK(proj.kb_coeff(7) == -1.0069878791_a);
	
}
#endif

#endif

