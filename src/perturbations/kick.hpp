/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__PERTURBATIONS__KICK
#define INQ__PERTURBATIONS__KICK

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <states/orbital_set.hpp>

namespace inq {
namespace perturbations {

class kick {

public:
	
	kick(math::vector3<double> const & arg_kick_field):
		kick_field_(arg_kick_field)
	{
	}	

	template <typename PhiType>
	void zero_step(PhiType & phi) const {

		auto cov_efield = phi.basis().cell().metric().to_covariant(kick_field_);
		
		gpu::run(phi.basis().local_sizes()[2], phi.basis().local_sizes()[1], phi.basis().local_sizes()[0],
						 [pop = phi.basis().point_op(), ph = begin(phi.cubic()), cov_efield, nst = phi.set_part().local_size()] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 
							 auto rr = pop.rvector(ix, iy, iz);
							 auto kick_factor = exp(complex(0.0, dot(cov_efield, rr)));
							 for(int ist = 0; ist < nst; ist++) ph[ix][iy][iz][ist] *= kick_factor;
						 });
	}

	auto has_uniform_electric_field() const {
		return false;
	}	

	auto uniform_electric_field(double time) const {
		return math::vector3<double>{0.0, 0.0, 0.0};
	}
	
private:

	math::vector3<double> kick_field_;
	
};

}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_PERTURBATIONS_KICK_UNIT_TEST
#undef INQ_PERTURBATIONS_KICK_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("perturbations::kick", "[perturbations::kick]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using Catch::Approx;
	using math::vector3;
	
	const int nvec = 12;

	auto ecut = 31.2_Ha;
	double phi_absdif = 0.0;
	double phi_dif = 0.0;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	systems::box box = systems::box::orthorhombic(4.2_b, 3.5_b, 6.4_b).cutoff_energy(ecut);
	
	basis::real_space bas(box, comm);

	basis::field_set<basis::real_space, complex> phi(bas, nvec);

	//Construct a field
	for(int ix = 0; ix < phi.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < phi.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < phi.basis().local_sizes()[2]; iz++){
				for(int ist = 0; ist < phi.set_part().local_size(); ist++){
					phi.cubic()[ix][iy][iz][ist] = complex(cos(ist+(ix+iy+iz)), 1.3*sin(ist+(cos(ix-iy-iz))));
				}
			}
		}
	}

	auto phi_old = phi;

	auto kick = perturbations::kick({0.1, 0.0, 0.0});

	kick.zero_step(phi);

	for(int ix = 0; ix < phi.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < phi.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < phi.basis().local_sizes()[2]; iz++){
				for(int ist = 0; ist < phi.set_part().local_size(); ist++){
					phi_absdif += norm(phi.cubic()[ix][iy][iz][ist]) - norm(phi_old.cubic()[ix][iy][iz][ist]);
					phi_dif += norm(phi.cubic()[ix][iy][iz][ist] - phi_old.cubic()[ix][iy][iz][ist]);
				}
			}
		}
	}

	CHECK(phi_absdif == Approx(0).margin(1.0e-9));
	CHECK(phi_dif > 1.0e-9);

}

#endif
#endif
