/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR
#define INQ__HAMILTONIAN__PROJECTOR

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

namespace inq {
namespace hamiltonian {

  class projector {

  public:
    projector(const basis::real_space & basis, const ions::UnitCell & cell, atomic_potential::pseudopotential_type const & ps, math::vector3<double> atom_position, int iatom):
      sphere_(basis, cell, atom_position, ps.projector_radius()),
      nproj_(ps.num_projectors_lm()),
      matrix_({nproj_, sphere_.size()}),
			kb_coeff_(nproj_),
			comm_(sphere_.create_comm(basis.comm())),
			iatom_(iatom){

			std::vector<double> grid(sphere_.size()), proj(sphere_.size());

			int iproj_lm = 0;
      for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				// interpolate the value of the radial part of the projectors to the sphere points
				ps.projector(iproj_l).value(sphere_.size(), sphere_.distance(), proj);
				
				int l = ps.projector_l(iproj_l);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
						matrix_[iproj_lm][ipoint] = proj[ipoint]*pseudo::math::spherical_harmonic(l, m, sphere_.point_pos()[ipoint]);
					}
					kb_coeff_[iproj_lm]	= ps.kb_coeff(iproj_l); 
					iproj_lm++;
				}
				
      }

			assert(iproj_lm == ps.num_projectors_lm());
			
    }

		auto empty() const {
			return nproj_ == 0 or sphere_.size() == 0;
		}

    template <class field_set_type>
    math::array<typename field_set_type::element_type, 2> project(const field_set_type & phi) const {

			assert(not empty());
			
			CALI_CXX_MARK_SCOPE("projector::project");
				
			auto sphere_phi = sphere_.gather(phi.cubic());

			CALI_MARK_BEGIN("projector_allocation");

			math::array<typename field_set_type::element_type, 2> projections({nproj_, phi.local_set_size()});

			CALI_MARK_END("projector_allocation");
			
			{ CALI_CXX_MARK_SCOPE("projector_gemm_1");
				namespace blas = boost::multi::blas;
				blas::real_doubled(projections) = blas::gemm(sphere_.volume_element(), matrix_, blas::real_doubled(sphere_phi));
			}
			
			{	CALI_CXX_MARK_SCOPE("projector_scal");
				
				//DATAOPERATIONS GPU::RUN 2D
				gpu::run(phi.local_set_size(), nproj_,
								 [proj = begin(projections), coeff = begin(kb_coeff_)]
								 GPU_LAMBDA (auto ist, auto iproj){
									 proj[iproj][ist] = proj[iproj][ist]*coeff[iproj];
								 });
			}
			
			if(comm_.size() > 1){
				CALI_CXX_MARK_SCOPE("projector_mpi_reduce");
				comm_.all_reduce_in_place_n(static_cast<typename field_set_type::element_type *>(projections.data_elements()), projections.num_elements(), std::plus<>{});
			}
			
			{
				CALI_CXX_MARK_SCOPE("projector_gemm_2");
				namespace blas = boost::multi::blas;
				blas::real_doubled(sphere_phi) = blas::gemm(1., blas::T(matrix_), blas::real_doubled(projections));
			}

			
			return sphere_phi;
		}
			
		template <class SpherePhiType, class field_set_type>
    void apply(SpherePhiType const & sphere_vnlphi, field_set_type & vnlphi) const {

			CALI_CXX_MARK_SCOPE("projector::apply");
			
			assert(not empty());
			
			sphere_.scatter_add(sphere_vnlphi, vnlphi.cubic());
		}

		template <typename OcType, typename PhiType, typename GPhiType>
		struct force_term {
			OcType oc;
			PhiType phi;
			GPhiType gphi;
			constexpr auto operator()(int ist, int ip) const {
				return -2.0*oc[ist]*real(phi[ip][ist]*conj(gphi[ip][ist]));
			}
		};
		
    template <class PhiType, typename GPhiType, typename OccsType>
    math::vector3<double> force(PhiType const & phi, GPhiType const & gphi, OccsType const & occs) const {

			math::vector3<double> force{0.0, 0.0, 0.0};

			if(nproj_ == 0) return force;

			CALI_CXX_MARK_SCOPE("projector::force");
				
			using boost::multi::blas::gemm;
			using boost::multi::blas::transposed;
			namespace blas = boost::multi::blas;
				
			auto sphere_phi = sphere_.gather(phi.cubic());
			auto sphere_gphi = sphere_.gather(gphi.cubic());			

			math::array<typename PhiType::element_type, 2> projections({nproj_, phi.local_set_size()});

			if(sphere_.size() > 0) {
				
				blas::real_doubled(projections) = gemm(sphere_.volume_element(), matrix_, blas::real_doubled(sphere_phi));
				
				{
					CALI_CXX_MARK_SCOPE("projector_force_scal"); 
					
					gpu::run(phi.local_set_size(), nproj_,
									 [proj = begin(projections), coeff = begin(kb_coeff_)]
									 GPU_LAMBDA (auto ist, auto iproj){
										 proj[iproj][ist] = proj[iproj][ist]*coeff[iproj];
									 });
				}
				
			} else {
				projections.elements().fill(0.0);
			}

			{	CALI_CXX_MARK_SCOPE("projector::force_mpi_reduce_1");
				comm_.all_reduce_in_place_n(static_cast<typename PhiType::element_type *>(projections.data_elements()), projections.num_elements(), std::plus<>{});
			}

			if(sphere_.size() > 0) {
				namespace blas = boost::multi::blas;
				blas::real_doubled(sphere_phi) = blas::gemm(1., transposed(matrix_), blas::real_doubled(projections));


				{
					CALI_CXX_MARK_SCOPE("projector_force_sum");
					force = gpu::run(gpu::reduce(phi.local_set_size()), gpu::reduce(sphere_.size()),
													 force_term<decltype(begin(occs)), decltype(begin(sphere_phi)), decltype(begin(sphere_gphi))>{begin(occs), begin(sphere_phi), begin(sphere_gphi)});
				}
				
			}
			
			return phi.basis().volume_element()*force;
    }
		
    int num_projectors() const {
      return nproj_;
    }
		
    auto kb_coeff(int iproj){
      return kb_coeff_[iproj];
    }

		auto iatom() const {
			return iatom_;
		}

  private:

    basis::spherical_grid sphere_;
    int nproj_;
    math::array<double, 2> matrix_;
		math::array<double, 1> kb_coeff_;
		mutable boost::mpi3::communicator comm_;
		int iatom_;
    
  };
  
}
}

#ifdef INQ_HAMILTONIAN_PROJECTOR_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_UNIT_TEST

#include <config/path.hpp>
#include <ions/geometry.hpp>

#include <catch2/catch.hpp>

TEST_CASE("class hamiltonian::projector", "[hamiltonian::projector]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;
	
	pseudo::math::erf_range_separation const sep(0.625);

	auto comm = boost::mpi3::environment::get_world_instance();

  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut), comm);

	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector proj(rs, cell, ps, vector3<double>(0.0, 0.0, 0.0));

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

