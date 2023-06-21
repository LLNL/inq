/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR_FOURIER
#define INQ__HAMILTONIAN__PROJECTOR_FOURIER

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <utils/profiling.hpp>
#include <config/path.hpp>

#include <pseudopod/math/sharmonic.hpp>

#include <gpu/array.hpp>
#include <math/vector3.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <operations/space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

  class projector_fourier {

  public:
    projector_fourier(const basis::real_space & basis, atomic_potential::pseudopotential_type const & ps):
      nproj_(ps.num_projectors_lm()),
			kb_coeff_(nproj_),
			beta_(basis::fourier_space(basis), nproj_){

			//For the moment we calculate the projectors in real space, centered at zero, and then move them to Fourier space

			basis::field_set<basis::real_space, complex> beta_rs(basis, nproj_);
			
			basis::spherical_grid sphere(beta_rs.basis(), vector3<double>(0.0, 0.0, 0.0), 1.5*ps.projector_radius());

			beta_rs.fill(0.0);
			
			int iproj_lm = 0;
      for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){

				int l = ps.projector_l(iproj_l);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere.size(); ipoint++){
						auto point = sphere.ref().grid_point(ipoint);
						beta_rs.hypercubic()[point[0]][point[1]][point[2]][iproj_lm] = ps.projector(iproj_l).function()(sphere.ref().distance(ipoint))*pseudo::math::sharmonic(l, m, sphere.ref().point_pos(ipoint));
					}
					
					kb_coeff_[iproj_lm]	= ps.kb_coeff(iproj_l); 
					iproj_lm++;
				}
				
      }

			assert(iproj_lm == ps.num_projectors_lm());

			beta_ = operations::space::to_fourier(beta_rs);
			
    }
		
		void add_coord(vector3<double, contravariant> const & coord){
			coords_.push_back(coord);
		}

    void operator()(states::orbital_set<basis::fourier_space, complex> const & phi, states::orbital_set<basis::fourier_space, complex> & vnlphi) const {

			CALI_CXX_MARK_SCOPE("projector_fourier");

			if(nproj_ == 0) return;

			for(unsigned icoord = 0; icoord < coords_.size(); icoord++){
			
				gpu::array<complex, 2> projections({nproj_, phi.set_part().local_size()}, 0.0);
			
				basis::field<basis::fourier_space, complex> eigr(phi.basis());

				auto point_op = eigr.basis().point_op();
				
				for(int ix = 0; ix < eigr.basis().local_sizes()[0]; ix++){
					for(int iy = 0; iy < eigr.basis().local_sizes()[1]; iy++){
						for(int iz = 0; iz < eigr.basis().local_sizes()[2]; iz++){
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
				
				gpu::run(phi.set_part().local_size(), nproj_,
								 [proj = begin(projections), coeff = begin(kb_coeff_), vol = phi.basis().volume_element()]
								 GPU_LAMBDA (auto ist, auto iproj){
									 proj[iproj][ist] = proj[iproj][ist]*coeff[iproj]*vol;
								 });

				phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(projections.data_elements()), projections.num_elements(), std::plus<>{});
				
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
		gpu::array<double, 1> kb_coeff_;
    basis::field_set<basis::fourier_space, complex> beta_;
		std::vector<vector3<double, contravariant>> coords_;
		
  };
  
}
}
#endif

#ifdef INQ_HAMILTONIAN_PROJECTOR_FOURIER_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_FOURIER_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
  
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	pseudo::math::erf_range_separation const sep(0.625);

	systems::box box = systems::box::cubic(10.0_b);
	basis::real_space rs(box, /*spacing = */ 0.49672941, comm);
	
	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector_fourier proj(rs, ps);

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


