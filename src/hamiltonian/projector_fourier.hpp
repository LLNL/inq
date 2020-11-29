/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR_FOURIER
#define INQ__HAMILTONIAN__PROJECTOR_FOURIER

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

#include <caliper/cali.h>

#ifdef ENABLE_CUDA
#include <multi/adaptors/blas/cuda.hpp> // must be included before blas.hpp
#endif
#include <multi/adaptors/blas.hpp>

#include <config/path.hpp>

#include <pseudopod/spherical_harmonic.hpp>

#include <math/array.hpp>
#include <math/vector3.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <operations/space.hpp>
#include <hamiltonian/atomic_potential.hpp>

namespace inq {
namespace hamiltonian {

  class projector_fourier {

  public:
    projector_fourier(const basis::real_space & basis, const ions::UnitCell & cell, atomic_potential::pseudopotential_type const & ps):
      nproj_(ps.num_projectors_lm()),
			kb_coeff_(nproj_),
			beta_(basis::fourier_space(basis), nproj_){

			//For the moment we calculate the projectors in real space, centered at zero, and then move them to Fourier space

			basis::field_set<basis::real_space, complex> beta_rs(basis, nproj_);
			
			basis::spherical_grid sphere(beta_rs.basis(), cell, math::vector3<double>(0.0, 0.0, 0.0), 1.5*ps.projector_radius());
			std::vector<double> grid(sphere.size()), proj(sphere.size());

			beta_rs = 0.0;
			
			int iproj_lm = 0;
      for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				// interpolate the value of the radial part of the projectors to the sphere points
				ps.projector(iproj_l).value(sphere.size(), sphere.distance(), proj);
				
				int l = ps.projector_l(iproj_l);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
						auto point = sphere.points()[ipoint];					
						beta_rs.cubic()[point[0]][point[1]][point[2]][iproj_lm] = proj[ipoint]*pseudo::math::spherical_harmonic(l, m, sphere.point_pos()[ipoint]);
					}
					
					kb_coeff_[iproj_lm]	= ps.kb_coeff(iproj_l); 
					iproj_lm++;
				}
				
      }

			assert(iproj_lm == ps.num_projectors_lm());

			beta_ = operations::space::to_fourier(beta_rs);
			
    }
		
		void add_coord(math::vector3<double> const & coord){
			coords_.push_back(coord);
		}

    void operator()(basis::field_set<basis::fourier_space, complex> const & phi, basis::field_set<basis::fourier_space, complex> & vnlphi) const {

			CALI_CXX_MARK_SCOPE("projector_fourier");

			if(nproj_ == 0) return;

			for(unsigned icoord = 0; icoord < coords_.size(); icoord++){
			
				math::array<complex, 2> projections({nproj_, phi.set_part().local_size()}, 0.0);
			
				basis::field<basis::fourier_space, complex> eigr(phi.basis());

				auto point_op = eigr.basis().point_op();
				
				for(int ix = 0; ix < eigr.basis().sizes()[0]; ix++){
					for(int iy = 0; iy < eigr.basis().sizes()[1]; iy++){
						for(int iz = 0; iz < eigr.basis().sizes()[2]; iz++){
							double gr = dot(coords_[icoord], point_op.gvector(ix, iy, iz));
							eigr.cubic()[ix][iy][iz] = exp(complex(0.0, gr));
						}
					}
				}
				
				for(int iproj = 0; iproj < nproj_; iproj++){
					for(long ip = 0; ip < phi.basis().part().local_size(); ip++){
						for(int ist = 0; ist < phi.set_part().local_size(); ist++){
							projections[iproj][ist] += eigr.linear()[ip]*conj(beta_.matrix()[ip][iproj])*phi.matrix()[ip][ist];
						}
					}
				}
				
				//DATAOPERATIONS GPU::RUN 2D
				gpu::run(phi.set_part().local_size(), nproj_,
								 [proj = begin(projections), coeff = begin(kb_coeff_), vol = phi.basis().volume_element()]
								 GPU_LAMBDA (auto ist, auto iproj){
									 proj[iproj][ist] = proj[iproj][ist]*coeff[iproj]*vol;
								 });

				//TODO: reduce projections
				
				for(int iproj = 0; iproj < nproj_; iproj++){
					for(long ip = 0; ip < phi.basis().part().local_size(); ip++){
						for(int ist = 0; ist < phi.set_part().local_size(); ist++){
							vnlphi.matrix()[ip][ist] += conj(eigr.linear()[ip])*beta_.matrix()[ip][iproj]*projections[iproj][ist];
						}
					}
				}
				
			}
			
    }
		
    int num_projectors() const {
      return nproj_;
    }
		
    auto kb_coeff(int iproj){
      return kb_coeff_[iproj];
    }
		
  private:

    int nproj_;
		math::array<double, 1> kb_coeff_;
    basis::field_set<basis::fourier_space, complex> beta_;
		std::vector<math::vector3<double>> coords_;
		
  };
  
}
}

#ifdef INQ_HAMILTONIAN_PROJECTOR_FOURIER_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_FOURIER_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("class hamiltonian::projector_fourier", "[hamiltonian::projector_fourier]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;
	
	pseudo::math::erf_range_separation const sep(0.625);
	
  double ecut = 20.0;
  double ll = 10.0;

  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector_fourier proj(rs, cell, ps);

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

