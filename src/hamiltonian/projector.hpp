/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR
#define INQ__HAMILTONIAN__PROJECTOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <sharmonic.h>

#include <gpu/array.hpp>
#include <math/vector3.hpp>
#include <basis/double_grid.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

class projector {

	basis::spherical_grid sphere_;
	int nproj_;
	gpu::array<double, 2> matrix_;
	gpu::array<double, 1> kb_coeff_;
	int iatom_;

public: // for CUDA
	
	void build(basis::real_space const & basis, basis::double_grid const & double_grid, atomic_potential::pseudopotential_type const & ps) {

		CALI_CXX_MARK_SCOPE("projector::build");
		
		int iproj_lm = 0;
		for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				
			int l = ps.projector_l(iproj_l);

			for(auto mm = 0; mm < 2*l + 1; mm++) kb_coeff_[iproj_lm + mm] = ps.kb_coeff(iproj_l, iproj_l);
			
			if(not double_grid.enabled()) {
				
				// now construct the projector with the spherical harmonics
				gpu::run(sphere_.size(), 2*l + 1,
								 [mat = begin(matrix_),
									spline = ps.projector(iproj_l).function(),
									sph = sphere_.ref(), l, iproj_lm,
									metric = basis.cell().metric()] GPU_LAMBDA (auto ipoint, auto m) {
									 mat[iproj_lm + m][ipoint] = spline(sph.distance(ipoint))*sharmonic::cartesian_real(l, m - l, metric.to_cartesian(sph.point_pos(ipoint)));
								 });
				
			} else {

				CALI_CXX_MARK_SCOPE("projector::double_grid");
				
				gpu::run(sphere_.size(), 2*l + 1,
								 [mat = begin(matrix_), spline = ps.projector(iproj_l).function(), sph = sphere_.ref(), l, iproj_lm,
									dg = double_grid.ref(), spac = basis.rspacing(), metric = basis.cell().metric()] GPU_LAMBDA (auto ipoint, auto m) {
									 mat[iproj_lm + m][ipoint] = dg.value([spline, l, m] GPU_LAMBDA(auto pos) { return spline(pos.length())*sharmonic::cartesian_real(l, m - l, pos);}, spac, metric.to_cartesian(sph.point_pos(ipoint)));
								 });
				
			}

			iproj_lm += 2*l + 1;
			
		}

		assert(iproj_lm == ps.num_projectors_lm());

	}
	
public:
	projector(const basis::real_space & basis, basis::double_grid const & double_grid, atomic_potential::pseudopotential_type const & ps, vector3<double> atom_position, int iatom):
		sphere_(basis, atom_position, ps.projector_radius()),
		nproj_(ps.num_projectors_lm()),
		matrix_({nproj_, sphere_.size()}),
		kb_coeff_(nproj_),
		iatom_(iatom){

		build(basis, double_grid, ps);
	}

	projector(projector const &) = delete;		

	auto empty() const {
		return nproj_ == 0;
	}
	
	auto locally_empty() const {
		return nproj_ == 0 or sphere_.size() == 0;
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
	
	auto & sphere() const {
		return sphere_;
	}

	auto & matrix() const {
		return matrix_;
	}

	friend class projector_all;
	
};
  
}
}
#endif

#ifdef INQ_HAMILTONIAN_PROJECTOR_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_UNIT_TEST

#include <config/path.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	pseudo::math::erf_range_separation const sep(0.625);

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	basis::real_space rs(systems::cell::cubic(10.0_b), /*spacing = */ 0.49672941, comm);
	basis::double_grid dg(false);
	
	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "N.upf", sep, rs.gcutoff());
	
	hamiltonian::projector proj(rs, dg, ps, vector3<double>(0.0, 0.0, 0.0), 77);

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

