/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__DENSITY
#define INQ__OBSERVABLES__DENSITY

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
#include <operations/integral.hpp>
#include <operations/transfer.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace observables {
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

template<class occupations_array_type, class field_set_type, class vector_field_set_type, typename VectorSpace>
void calculate_gradient_add(const occupations_array_type & occupations, field_set_type const & phi, vector_field_set_type const & gphi, basis::field<typename vector_field_set_type::basis_type, math::vector3<double, VectorSpace>> & gdensity){

	CALI_CXX_MARK_SCOPE("density::calculate_gradient");
 
	gpu::run(phi.basis().part().local_size(),
					 [nst = phi.set_part().local_size(), occs = begin(occupations),
            phip = begin(phi.matrix()), gphip = begin(gphi.matrix()), gdensityp = begin(gdensity.linear())] GPU_LAMBDA (auto ip){
						 for(int ist = 0; ist < nst; ist++) gdensityp[ip] += occs[ist]*real(conj(gphip[ip][ist])*phip[ip][ist] + conj(phip[ip][ist])*gphip[ip][ist]);
					 });

}

///////////////////////////////////////////////////////////////

template <typename ElecType>
basis::field<basis::real_space, double> calculate(ElecType & elec){
	
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

///////////////////////////////////////////////////////////////

template <class FieldType>
void normalize(FieldType & density, const double & total_charge){

	CALI_CXX_MARK_FUNCTION;
	
	auto qq = operations::integral(density);
	assert(fabs(qq) > 1e-16);
	for(int i = 0; i < density.basis().part().local_size(); i++) density.linear()[i] *= total_charge/qq;
	
}

///////////////////////////////////////////////////////////////

template <class FieldSetType>
basis::field<basis::real_space, double> total(FieldSetType & spin_density){

	CALI_CXX_MARK_FUNCTION;

	assert(spin_density.set_size() == 1 or spin_density.set_size() == 2 or spin_density.set_size() == 4);
	assert(spin_density.set_size() == spin_density.local_set_size());
	
	basis::field<basis::real_space, double> total_density(spin_density.basis());
	
	gpu::run(spin_density.basis().local_size(),
					 [spi = begin(spin_density.matrix()), tot = begin(total_density.linear()), nspin = spin_density.set_size()] GPU_LAMBDA (auto ip){
						 if(nspin == 1) tot[ip] = spi[ip][0];
						 if(nspin == 2) tot[ip] = spi[ip][0] + spi[ip][1];
						 if(nspin == 4) tot[ip] = spi[ip][0] + spi[ip][1] + spi[ip][2] + spi[ip][3];
					 });
}

}
}
}

#ifdef INQ_OBSERVABLES_DENSITY_UNIT_TEST
#undef INQ_OBSERVABLES_DENSITY_UNIT_TEST

#include <basis/trivial.hpp>
#include <math/complex.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("function observables::density", "[observables::density]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

	const int npoint = 100;
	const int nvec = 12;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	parallel::cartesian_communicator<2> cart_comm(comm, {});
	
	auto basis_comm = basis::basis_subcomm(cart_comm);
	
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
		
		observables::density::calculate_add(occ, aa, dd);

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
		
		observables::density::calculate_add(occ, aa, dd);

		aa.set_comm().all_reduce_in_place_n(raw_pointer_cast(dd.linear().data_elements()), dd.linear().size(), std::plus<>{});
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(dd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));

	}
	
}

TEST_CASE("function observables::density::normalize", "[observables::density::normalize]") {

	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 100;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	basis::trivial bas(npoint, comm);
	
	SECTION("double"){
		
		basis::field<basis::trivial, double> aa(bas);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) aa.linear()[ii] = sqrt(bas.part().local_to_global(ii).value());

		observables::density::normalize(aa, 33.3);

		CHECK(operations::integral(aa) == 33.3_a);
		
	}
	
	SECTION("complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			aa.linear()[ii] = sqrt(bas.part().local_to_global(ii).value())*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii).value()));
		}

		observables::density::normalize(aa, 19.2354);

		CHECK(real(operations::integral(aa)) == 19.2354_a);
		
	}
	
}

#endif

#endif
