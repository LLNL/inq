/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__RELATIVISTIC_PROJECTOR
#define INQ__HAMILTONIAN__RELATIVISTIC_PROJECTOR

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pseudopod/math/sharmonic.hpp>

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

class relativistic_projector {

	basis::spherical_grid sphere_;
	int nproj_;
	gpu::array<double, 3> matrix_;
	gpu::array<double, 1> kb_coeff_;
	int iatom_;

public: // for CUDA
	
	void build(basis::real_space const & basis, basis::double_grid const & double_grid, atomic_potential::pseudopotential_type const & ps) {

		CALI_CXX_MARK_SCOPE("relativistic_projector::build");

		nproj_ = 0.0;
		for(int iproj = 0; iproj < ps.num_projectors_l(); iproj++){
				
			int const jj = std::lround(2.0*ps.projector_j(iproj));

			nproj_ += jj + 1;
		}

		std::cout << "RELATIVISTIC " << ps.full_relativistic() << std::endl;
		
		matrix_.reextent({nproj_, sphere_.size(), 2});
		kb_coeff_.reextent(nproj_);
		
		int iproj_lm = 0;
		for(int iproj = 0; iproj < ps.num_projectors_l(); iproj++){
				
			auto ll = ps.projector_l(iproj);
			int const jj = std::lround(2.0*ps.projector_j(iproj));

			std::cout << "LL = " << ll << " JJ = " << jj/2.0 << std::endl;

			auto sgn = 1.0;
			if(jj == 2*ll - 1.0) sgn = -1.0;

			for(auto mj = -jj; mj <= jj; mj += 2){

				auto den = sqrt(jj - sgn + 1);
				auto cc0 = sgn*sqrt(jj/2.0 - 0.5*sgn + sgn*mj/2.0 + 0.5)/den;
				auto cc1 =     sqrt(jj/2.0 - 0.5*sgn - sgn*mj/2.0 + 0.5)/den;
				auto mm0 = lround(mj/2.0 - 0.5);
				auto mm1 = lround(mj/2.0 + 0.5);

				std::cout << cc0 << '\t' << cc1 << '\t' << mm0 << '\t' << mm1 << std::endl;												
				
				gpu::run(sphere_.size(),
								 [mat = begin(matrix_),
									spline = ps.projector(iproj).function(),
									sph = sphere_.ref(), cc0, cc1, mm0, mm1, ll, iproj_lm,
									metric = basis.cell().metric()] GPU_LAMBDA (auto ipoint) {

									 if(abs(mm0) <= ll) {
										 mat[iproj_lm][ipoint][0] = cc0*spline(sph.distance(ipoint))*pseudo::math::sharmonic(ll, mm0, metric.to_cartesian(sph.point_pos(ipoint)));
									 } else {
										 mat[iproj_lm][ipoint][0] = 0.0;
									 }
									 if(abs(mm1) <= ll) {									 
										 mat[iproj_lm][ipoint][1] = cc1*spline(sph.distance(ipoint))*pseudo::math::sharmonic(ll, mm1, metric.to_cartesian(sph.point_pos(ipoint)));
									 } else {
										 mat[iproj_lm][ipoint][1] = 0.0;
									 }
								 });

				kb_coeff_[iproj_lm] = ps.kb_coeff(iproj);
				
				iproj_lm++;
			}
			
		}
	}
	
public:
	relativistic_projector(const basis::real_space & basis, basis::double_grid const & double_grid, atomic_potential::pseudopotential_type const & ps, vector3<double> atom_position, int iatom):
		sphere_(basis, atom_position, ps.projector_radius()),
		iatom_(iatom){

		build(basis, double_grid, ps);
	}

	relativistic_projector(relativistic_projector const &) = delete;		

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

#ifdef INQ_HAMILTONIAN_RELATIVISTIC_PROJECTOR_UNIT_TEST
#undef INQ_HAMILTONIAN_RELATIVISTIC_PROJECTOR_UNIT_TEST

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
	
	hamiltonian::atomic_potential::pseudopotential_type ps(config::path::unit_tests_data() + "He_fr.upf.gz", sep, rs.gcutoff());
	
	hamiltonian::relativistic_projector proj(rs, dg, ps, vector3<double>(0.0, 0.0, 0.0), 77);

	CHECK(proj.num_projectors() == 10);
	/*
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
	*/
}
#endif

