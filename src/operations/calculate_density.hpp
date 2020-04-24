/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__CALCULATE_DENSITY
#define OPERATIONS__CALCULATE_DENSITY

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

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <math/complex.hpp>
#include <cstdlib>

namespace operations {

  template<class occupations_array_type, class field_set_type>
  basis::field<typename field_set_type::basis_type, double> calculate_density(const occupations_array_type & occupations, field_set_type & phi){

    basis::field<typename field_set_type::basis_type, double> density(phi.basis(), phi.basis_comm());

    //DATAOPERATIONS LOOP + GPU::RUN 2D
#ifdef HAVE_CUDA

		const auto nst = phi.set_part().local_size();
		auto occupationsp = begin(occupations);
		auto phip = begin(phi.matrix());
		auto densityp = begin(density.linear());
		
		gpu::run(phi.basis().part().local_size(),
						 [=] __device__ (auto ipoint){
							 densityp[ipoint] = 0.0;
							 for(int ist = 0; ist < nst; ist++) densityp[ipoint] += occupationsp[ist]*norm(phip[ipoint][ist]);
						 });
		
#else
		
    for(int ipoint = 0; ipoint < phi.basis().part().local_size(); ipoint++){
			density.linear()[ipoint] = 0.0;
      for(int ist = 0; ist < phi.set_part().local_size(); ist++) density.linear()[ipoint] += occupations[ist]*norm(phi.matrix()[ipoint][ist]);
    }
		
#endif

		if(phi.set_part().parallel()){
			phi.set_comm().all_reduce_in_place_n(static_cast<double *>(density.linear().data()), density.linear().size(), std::plus<>{});
		}
		
    return density;
  }

	template <class FieldType>
	void normalize_density(FieldType & density, const double & total_charge){

		auto qq = operations::integral(density);
		assert(qq > 1e-16);
		for(int i = 0; i < density.basis().part().local_size(); i++) density.linear()[i] *= total_charge/qq;

	}
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::calculate_density", "[operations::calculate_density]") {

	using namespace Catch::literals;

	const int npoint = 100;
	const int nvec = 12;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm);
	
	auto basis_comm = cart_comm.axis(1);
	
	basis::trivial bas(npoint, basis_comm);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);

		math::array<double, 1> occ(aa.set_part().local_size());
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii))*(aa.set_part().local_to_global(jj) + 1);
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj) + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) REQUIRE(dd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii)*nvec*(nvec + 1)));

		operations::normalize_density(dd, 33.3);

		REQUIRE(operations::integral(dd) == 33.3_a);
		
	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);

		math::array<double, 1> occ(nvec);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii))*(aa.set_part().local_to_global(jj) + 1)*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii)));
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj) + 1);

		auto dd = operations::calculate_density(occ, aa);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) REQUIRE(dd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii)*nvec*(nvec + 1)));

		operations::normalize_density(dd, 33.3);

		REQUIRE(operations::integral(dd) == 33.3_a);
		
	}
	
}


#endif

#endif
