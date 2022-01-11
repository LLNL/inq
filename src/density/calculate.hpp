/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__DENSITY__CALCULATE
#define INQ__DENSITY__CALCULATE

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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
#include <operations/transfer.hpp>
#include <systems/electrons.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace density {

template<class occupations_array_type, class field_set_type>
void calculate_add(const occupations_array_type & occupations, field_set_type & phi, basis::field<typename field_set_type::basis_type, double> & density){

	const auto nst = phi.set_part().local_size();
	auto occupationsp = begin(occupations);
	auto phip = begin(phi.matrix());
	auto densityp = begin(density.linear());
		
	gpu::run(phi.basis().part().local_size(),
					 [=] GPU_LAMBDA (auto ipoint){
						 for(int ist = 0; ist < nst; ist++) densityp[ipoint] += occupationsp[ist]*norm(phip[ipoint][ist]);
					 });

}

///////////////////////////////////////////////////////////////

template<class occupations_array_type, class field_set_type, class vector_field_set_type>
void calculate_gradient_add(const occupations_array_type & occupations, field_set_type const & phi, vector_field_set_type const & gphi, basis::field<typename vector_field_set_type::basis_type, math::vector3<double>> & gdensity){

	CALI_CXX_MARK_SCOPE("density::calculate_gradient");
 
	gpu::run(3, phi.basis().part().local_size(),
					 [nst = phi.set_part().local_size(), occs = begin(occupations),
            phip = begin(phi.matrix()), gphip = begin(gphi.matrix()), gdensityp = begin(gdensity.linear())] GPU_LAMBDA (auto idir, auto ip){
						 for(int ist = 0; ist < nst; ist++) gdensityp[ip][idir] += occs[ist]*real(conj(gphip[ip][ist][idir])*phip[ip][ist] + conj(phip[ip][ist])*gphip[ip][ist][idir]);
					 });

}

///////////////////////////////////////////////////////////////

basis::field<basis::real_space, double> calculate(const systems::electrons & elec){
	
	basis::field<basis::real_space, double> density(elec.density_basis_);

	density = 0.0;

	int iphi = 0;
	for(auto & phi : elec.lot()) {
		density::calculate_add(elec.occupations()[iphi], phi, density);
		iphi++;
	}
	
	if(elec.lot_states_comm_.size() > 1) elec.lot_states_comm_.all_reduce_in_place_n(raw_pointer_cast(density.linear().data_elements()), density.linear().size(), std::plus<>{});

	return density;
}

}
}

#ifdef INQ_DENSITY_CALCULATE_UNIT_TEST
#undef INQ_DENSITY_CALCULATE_UNIT_TEST

#include <basis/trivial.hpp>
#include <math/complex.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("function density::calculate", "[density::calculate]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 100;
	const int nvec = 12;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {});
	
	auto basis_comm = cart_comm.axis(1);
	
	basis::trivial bas(npoint, basis_comm);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);

		math::array<double, 1> occ(aa.set_part().local_size());
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii).value())*(aa.set_part().local_to_global(jj).value() + 1);
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj).value() + 1);

		basis::field<basis::trivial, double> dd(bas);
		dd = 0.0;
		
		density::calculate_add(occ, aa, dd);

		aa.set_comm().all_reduce_in_place_n(raw_pointer_cast(dd.linear().data_elements()), dd.linear().size(), std::plus<>{});
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(dd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));

	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, nvec, cart_comm);

		math::array<double, 1> occ(nvec);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii).value())*(aa.set_part().local_to_global(jj).value() + 1)*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii).value()));
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj).value() + 1);

		basis::field<basis::trivial, double> dd(bas);
		dd = 0.0;
		
		density::calculate_add(occ, aa, dd);

		aa.set_comm().all_reduce_in_place_n(raw_pointer_cast(dd.linear().data_elements()), dd.linear().size(), std::plus<>{});
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(dd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));

	}
	
}


#endif

#endif
