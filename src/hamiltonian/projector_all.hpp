/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR_ALL
#define INQ__HAMILTONIAN__PROJECTOR_ALL

/*
	Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa.

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

#include <math/array.hpp>
#include <math/vector3.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

class projector_all {

#ifndef ENABLE_CUDA
private:
#else
public:
#endif
	
	template <typename ProjectorsType>
	void constructor(ProjectorsType const & projectors){

		CALI_CXX_MARK_FUNCTION;

		nprojs_ = projectors.size();
		max_sphere_size_ = 0;
		max_nlm_ = 0;
		for(auto it = projectors.cbegin(); it != projectors.cend(); ++it) {
			max_sphere_size_ = std::max(max_sphere_size_, it->sphere_.size());
			max_nlm_ = std::max(max_nlm_, it->nproj_);			
		}

		points_ = decltype(points_)({nprojs_, max_sphere_size_});
    coeff_ = decltype(coeff_)({nprojs_, max_nlm_}, 0.0);
    matrices_ = decltype(matrices_)({nprojs_, max_nlm_, max_sphere_size_});
		comms_ = decltype(comms_)({nprojs_});
		
    auto iproj = 0;
    for(auto it = projectors.cbegin(); it != projectors.cend(); ++it) {
			gpu::run(max_sphere_size_,
							 [poi = begin(points_), sph = it->sphere_.ref(), iproj, npoint = it->sphere_.size()] GPU_LAMBDA (auto ipoint){
								 if(ipoint < unsigned (npoint)){
									 poi[iproj][ipoint] = sph.points(ipoint);
								 } else {
									 poi[iproj][ipoint] = {-1, -1, -1};
								 }
							 });
			
			gpu::run(max_sphere_size_, max_nlm_,
							 [mat = begin(matrices_), itmat = begin(it->matrix_), iproj, np = it->sphere_.size(), nlm = it->nproj_] GPU_LAMBDA (auto ipoint, auto ilm){
								 if(ipoint < (unsigned) np and ilm < (unsigned) nlm) {
									 mat[iproj][ilm][ipoint] = itmat[ilm][ipoint];
								 } else {
									 mat[iproj][ilm][ipoint] = 0.0;								 
								 }
							 });
								 
      coeff_[iproj]({0, it->nproj_}) = it->kb_coeff_;

			comms_[iproj] = it->comm_;
			
      iproj++;
    }
    
	}
	
public:

	projector_all():
		nprojs_(0),
    max_sphere_size_(0),
    max_nlm_(0) {
  }
  

	template <typename ProjectorsType>
	projector_all(ProjectorsType const & projectors){
		constructor(projectors);
	}
  
	math::array<complex, 3> project(basis::field_set<basis::real_space, complex> const & phi) const {
    
		math::array<complex, 3> sphere_phi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		math::array<complex, 3> projections_all({nprojs_, max_nlm_, phi.local_set_size()});


		{ CALI_CXX_MARK_SCOPE("projector::gather");
				
			gpu::run(phi.local_set_size(), max_sphere_size_, nprojs_,
							 [sgr = begin(sphere_phi_all), gr = begin(phi.cubic()), poi = begin(points_)] GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
								 if(poi[iproj][ipoint][0] >= 0){
									 sgr[iproj][ipoint][ist] = gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
								 } else {
									 sgr[iproj][ipoint][ist] = complex(0.0, 0.0);
								 }
							 });
		}

		for(auto iproj = 0; iproj < nprojs_; iproj++){
			CALI_CXX_MARK_SCOPE("projector_gemm_1");

			namespace blas = boost::multi::blas;
			blas::real_doubled(projections_all[iproj]) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sphere_phi_all[iproj])); 
		}
				
    { CALI_CXX_MARK_SCOPE("projector_scal");
				
      gpu::run(phi.local_set_size(), max_nlm_, nprojs_,
               [proj = begin(projections_all), coe = begin(coeff_)]
               GPU_LAMBDA (auto ist, auto ilm, auto iproj){
                 proj[iproj][ilm][ist] = proj[iproj][ilm][ist]*coe[iproj][ilm];
               });
		}

		for(auto iproj = 0; iproj < nprojs_; iproj++){
			CALI_CXX_MARK_SCOPE("projector_mpi_reduce");
			
			if(comms_[iproj].size() == 1) continue;
			comms_[iproj].all_reduce_in_place_n(raw_pointer_cast(&projections_all[iproj][0][0]), max_nlm_*phi.local_set_size(), std::plus<>{});
		}
		
		for(auto iproj = 0; iproj < nprojs_; iproj++) {
			CALI_CXX_MARK_SCOPE("projector_gemm_2");
			namespace blas = boost::multi::blas;
			blas::real_doubled(sphere_phi_all[iproj]) = blas::gemm(1., blas::T(matrices_[iproj]), blas::real_doubled(projections_all[iproj]));
		}

		return sphere_phi_all;
			
	}

	////////////////////////////////////////////////////////////////////////////////////////////		

	template <typename SpherePhiType>
	void apply(SpherePhiType & sphere_vnlphi, basis::field_set<basis::real_space, complex> & vnlphi) const {

		CALI_CXX_MARK_FUNCTION;

		for(auto iproj = 0; iproj < nprojs_; iproj++){
			gpu::run(vnlphi.local_set_size(), max_sphere_size_,
							 [sgr = begin(sphere_vnlphi), gr = begin(vnlphi.cubic()), poi = begin(points_), iproj] GPU_LAMBDA (auto ist, auto ipoint){
								 if(poi[iproj][ipoint][0] >= 0){
									 gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist] += sgr[iproj][ipoint][ist];
								 }
							 });
		}
	}

private:
			
	int nprojs_;
	long max_sphere_size_;
	int max_nlm_;
	math::array<math::vector3<int>, 2> points_;
	math::array<double, 2> coeff_;
	math::array<double, 3> matrices_;
	mutable math::array<boost::mpi3::communicator, 1> comms_;	

  
};
  
}
}

#ifdef INQ_HAMILTONIAN_PROJECTOR_ALL_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_ALL_UNIT_TEST

#include <config/path.hpp>
#include <ions/geometry.hpp>

#include <catch2/catch.hpp>

TEST_CASE("class hamiltonian::projector_all", "[hamiltonian::projector_all]") {
  
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
  using math::vector3;

}

#endif

#endif

