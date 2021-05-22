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
#include <basis/double_grid.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

class projector {

#ifdef ENABLE_CUDA
public:
#endif
	
	void build(const basis::real_space & basis, atomic_potential::pseudopotential_type const & ps) {

		CALI_CXX_MARK_SCOPE("projector::build");
		
		int iproj_lm = 0;
		for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				
			int l = ps.projector_l(iproj_l);

			if(not basis.double_grid().enabled()) {
				
				// now construct the projector with the spherical harmonics
				gpu::run(sphere_.size(), 2*l + 1,
								 [mat = begin(matrix_), spline = ps.projector(iproj_l).cbegin(), sph = sphere_.ref(), l, iproj_lm, kb_ = begin(kb_coeff_), coe = ps.kb_coeff(iproj_l)] GPU_LAMBDA (auto ipoint, auto m) {
									 
									 if(ipoint == 0) kb_[iproj_lm + m] = coe;
									 mat[iproj_lm + m][ipoint] = spline.value(sph.distance(ipoint))*pseudo::math::spherical_harmonic(l, m - l, sph.point_pos(ipoint));
								 });
				
			} else {

				CALI_CXX_MARK_SCOPE("projector::double_grid");
				
				gpu::run(sphere_.size(), 2*l + 1,
								 [mat = begin(matrix_), spline = ps.projector(iproj_l).cbegin(), sph = sphere_.ref(), l, iproj_lm, kb_ = begin(kb_coeff_), coe = ps.kb_coeff(iproj_l),
									dg = basis.double_grid().ref(), spac = basis.rspacing()] GPU_LAMBDA (auto ipoint, auto m) {
									 
									 if(ipoint == 0) kb_[iproj_lm + m] = coe;
									 mat[iproj_lm + m][ipoint] = dg.value([spline, l, m] GPU_LAMBDA(auto pos) { return spline.value(length(pos))*pseudo::math::spherical_harmonic(l, m - l, pos);}, spac, sph.point_pos(ipoint));
								 });
				
			}

			iproj_lm += 2*l + 1;
			
		}

		assert(iproj_lm == ps.num_projectors_lm());

	}
	
public:
	projector(const basis::real_space & basis, const ions::UnitCell & cell, atomic_potential::pseudopotential_type const & ps, math::vector3<double> atom_position, int iatom):
		sphere_(basis, cell, atom_position, ps.projector_radius()),
		nproj_(ps.num_projectors_lm()),
		matrix_({nproj_, sphere_.size()}),
		kb_coeff_(nproj_),
		comm_(sphere_.create_comm(basis.comm())),
		iatom_(iatom){

		build(basis, ps);

	}

	projector(projector const &) = delete;		

	auto empty() const {
		return nproj_ == 0 or sphere_.size() == 0;
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

		math::array<typename PhiType::element_type, 2> sphere_phi({sphere_.size(), phi.local_set_size()});
		math::array<typename GPhiType::element_type, 2> sphere_gphi({sphere_.size(), phi.local_set_size()});		

		{
			CALI_CXX_MARK_SCOPE("projector_force_gather"); 
			gpu::run(phi.local_set_size(), sphere_.size(),
							 [sphi = begin(sphere_phi), sgphi = begin(sphere_gphi), phic = begin(phi.cubic()), gphic = begin(gphi.cubic()), sph = sphere_.ref()] GPU_LAMBDA (auto ist, auto ipoint){
								 sphi[ipoint][ist] = phic[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]][ist];
								 sgphi[ipoint][ist] = gphic[sph.points(ipoint)[0]][sph.points(ipoint)[1]][sph.points(ipoint)[2]][ist];
							 });
		}
		
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

		if(comm_.size() > 1) {
			CALI_CXX_MARK_SCOPE("projector::force_mpi_reduce_1");
			comm_.all_reduce_in_place_n(raw_pointer_cast(projections.data_elements()), projections.num_elements(), std::plus<>{});
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

	friend class projector_all;
	
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
	using namespace inq::magnitude;
	using namespace Catch::literals;
  using math::vector3;
	
	pseudo::math::erf_range_separation const sep(0.625);

	auto comm = boost::mpi3::environment::get_world_instance();

  auto ecut = 20.0_Ha;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut), comm);

	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector proj(rs, cell, ps, vector3<double>(0.0, 0.0, 0.0), 77);

	CHECK(proj.num_projectors() == 8);

	if(not proj.empty()){
		CHECK(proj.kb_coeff(0) ==  7.494508815_a);
		CHECK(proj.kb_coeff(1) ==  0.6363049519_a);
		CHECK(proj.kb_coeff(2) == -4.2939052122_a);
		CHECK(proj.kb_coeff(3) == -4.2939052122_a);
		CHECK(proj.kb_coeff(4) == -4.2939052122_a);
		CHECK(proj.kb_coeff(5) == -1.0069878791_a);
		CHECK(proj.kb_coeff(6) == -1.0069878791_a);
		CHECK(proj.kb_coeff(7) == -1.0069878791_a);
	}
	
	CHECK(proj.iatom() == 77);
	
}
#endif

#endif

