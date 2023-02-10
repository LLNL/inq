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
void calculate_add(const occupations_array_type & occupations, field_set_type & phi, basis::field_set<typename field_set_type::basis_type, double> & density){
	
	if(not phi.spinors()){

		assert(phi.spin_index() < density.set_size());
		
		gpu::run(phi.basis().part().local_size(),
						 [nst = phi.set_part().local_size(), occ = begin(occupations), ph = begin(phi.matrix()), den = begin(density.matrix()), ispin = phi.spin_index()] GPU_LAMBDA (auto ipoint){
							 for(int ist = 0; ist < nst; ist++) den[ipoint][ispin] += occ[ist]*norm(ph[ipoint][ist]);
						 });
	} else {
		
		assert(density.set_size() == 4);
		assert(std::get<1>(sizes(phi.spinor_matrix())) == phi.spinor_local_set_size());
		assert(std::get<2>(sizes(phi.spinor_matrix())) == phi.spinor_dim());
		
		gpu::run(phi.basis().part().local_size(),
						 [nst = phi.spinor_local_set_size(), occ = begin(occupations), ph = begin(phi.spinor_matrix()), den = begin(density.matrix())] GPU_LAMBDA (auto ipoint){
							 for(int ist = 0; ist < nst; ist++) {
								 den[ipoint][0] += occ[ist]*norm(ph[ipoint][ist][0]);
								 den[ipoint][1] += occ[ist]*norm(ph[ipoint][ist][1]);								 
								 auto crossterm = occ[ist]*ph[ipoint][ist][0]*conj(ph[ipoint][ist][1]);
								 den[ipoint][2] += real(crossterm);
								 den[ipoint][3] += imag(crossterm);
							 }
						 });
	}
		
}

///////////////////////////////////////////////////////////////

template<class occupations_array_type, class field_set_type, class vector_field_set_type, typename VectorSpace>
void calculate_gradient_add(const occupations_array_type & occupations, field_set_type const & phi, vector_field_set_type const & gphi, basis::field<typename vector_field_set_type::basis_type, vector3<double, VectorSpace>> & gdensity){

	CALI_CXX_MARK_SCOPE("density::calculate_gradient");

	gpu::run(phi.basis().part().local_size(),
					 [nst = phi.set_part().local_size(), occs = begin(occupations),
            phip = begin(phi.matrix()), gphip = begin(gphi.matrix()), gdensityp = begin(gdensity.linear())] GPU_LAMBDA (auto ip){
						 for(int ist = 0; ist < nst; ist++) gdensityp[ip] += occs[ist]*real(conj(gphip[ip][ist])*phip[ip][ist] + conj(phip[ip][ist])*gphip[ip][ist]);
					 });

}

///////////////////////////////////////////////////////////////

template <typename ElecType>
basis::field_set<basis::real_space, double> calculate(ElecType & elec){
	
	basis::field_set<basis::real_space, double> density(elec.density_basis_, elec.states().num_density_components());

	density.fill(0.0);

	int iphi = 0;
	for(auto & phi : elec.lot()) {
		density::calculate_add(elec.occupations()[iphi], phi, density);
		iphi++;
	}

	density.all_reduce(elec.lot_states_comm_);

	return density;
}

///////////////////////////////////////////////////////////////

template <class FieldType>
void normalize(FieldType & density, const double & total_charge){

	CALI_CXX_MARK_FUNCTION;
	
	auto qq = operations::integral_sum(density);
	assert(fabs(qq) > 1e-16);

	gpu::run(density.local_set_size(), density.basis().local_size(),
					 [den = begin(density.matrix()), factor = total_charge/qq] GPU_LAMBDA (auto ist, auto ip){ 
						 den[ip][ist] *= factor;
					 });
}

///////////////////////////////////////////////////////////////

template <class BasisType, class ElementType>
basis::field<BasisType, ElementType> total(basis::field_set<BasisType, ElementType> const & spin_density){

	CALI_CXX_MARK_FUNCTION;

	assert(spin_density.set_size() == 1 or spin_density.set_size() == 2 or spin_density.set_size() == 4);
	assert(spin_density.set_size() == spin_density.local_set_size());
	
	basis::field<BasisType, ElementType> total_density(spin_density.basis());

	gpu::run(spin_density.basis().local_size(),
					 [spi = begin(spin_density.matrix()), tot = begin(total_density.linear()), nspin = spin_density.set_size()] GPU_LAMBDA (auto ip){
						 if(nspin == 1) tot[ip] = spi[ip][0];
						 else tot[ip] = spi[ip][0] + spi[ip][1];
					 });

	return total_density;
}

}
}
}
#endif

#ifdef INQ_OBSERVABLES_DENSITY_UNIT_TEST
#undef INQ_OBSERVABLES_DENSITY_UNIT_TEST

#include <basis/trivial.hpp>
#include <math/complex.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

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
		
		states::orbital_set<basis::trivial, double> aa(bas, nvec, 1, vector3<double, covariant>{0.0, 0.0, 0.0}, 0, cart_comm);

		math::array<double, 1> occ(aa.set_part().local_size());
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii).value())*(aa.set_part().local_to_global(jj).value() + 1);
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj).value() + 1);

		basis::field_set<basis::trivial, double> dd(bas, 1);
		dd.fill(0.0);
		
		observables::density::calculate_add(occ, aa, dd);

		dd.all_reduce(aa.set_comm());
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(dd.matrix()[ii][0] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));

		auto tdd = observables::density::total(dd);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(tdd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));
		
	}
	
	SECTION("complex"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, 1, vector3<double, covariant>{0.0, 0.0, 0.0}, 0, cart_comm);

		math::array<double, 1> occ(nvec);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				aa.matrix()[ii][jj] = sqrt(bas.part().local_to_global(ii).value())*(aa.set_part().local_to_global(jj).value() + 1)*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii).value()));
			}
		}

		for(int jj = 0; jj < aa.set_part().local_size(); jj++) occ[jj] = 1.0/(aa.set_part().local_to_global(jj).value() + 1);

		basis::field_set<basis::trivial, double> dd(bas, 1);
		dd.fill(0.0);
		
		observables::density::calculate_add(occ, aa, dd);

		dd.all_reduce(aa.set_comm());
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(dd.matrix()[ii][0] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));

		auto tdd = observables::density::total(dd);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) CHECK(tdd.linear()[ii] == Approx(0.5*bas.part().local_to_global(ii).value()*nvec*(nvec + 1)));
		
	}
	
	SECTION("spinor"){
		
		states::orbital_set<basis::trivial, complex> aa(bas, nvec, 2, vector3<double, covariant>{0.0, 0.0, 0.0}, 0, cart_comm);

		math::array<double, 1> occ(nvec);
		
		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			for(int jj = 0; jj < aa.spinor_local_set_size(); jj++){
				auto iig = bas.part().local_to_global(ii).value();
				auto jjg = aa.spinor_set_part().local_to_global(jj).value();
				aa.spinor_matrix()[ii][jj][0] = sqrt(iig)*(jjg + 1)*exp(complex(0.0, M_PI/65.0*iig));
				aa.spinor_matrix()[ii][jj][1] = sqrt(iig)*(jjg + 1)*exp(complex(0.0, M_PI/65.0*iig));				
			}
		}

		for(int jj = 0; jj < aa.spinor_local_set_size(); jj++) {
			auto jjg = aa.spinor_set_part().local_to_global(jj).value();			
			occ[jj] = 1.0/(jjg + 1);
		}

		basis::field_set<basis::trivial, double> dd(bas, 4);
		dd.fill(0.0);
		
		observables::density::calculate_add(occ, aa, dd);

		dd.all_reduce(aa.set_comm());

		for(int ii = 0; ii < dd.basis().part().local_size(); ii++) {
			auto iig = bas.part().local_to_global(ii).value();			
			CHECK(dd.matrix()[ii][0] == Approx(0.5*iig*nvec*(nvec + 1)));
			CHECK(dd.matrix()[ii][1] == Approx(0.5*iig*nvec*(nvec + 1)));
			CHECK(dd.matrix()[ii][2] == Approx(0.5*iig*nvec*(nvec + 1)));
			CHECK(fabs(dd.matrix()[ii][3]) < 1e-12);
		}

		
		auto tdd = observables::density::total(dd);

		for(int ii = 0; ii < dd.basis().part().local_size(); ii++) {
			auto iig = bas.part().local_to_global(ii).value();
			CHECK(tdd.linear()[ii] == Approx(2.0*0.5*iig*nvec*(nvec + 1)));
		}
		
	}
}

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 100;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	basis::trivial bas(npoint, comm);
	
	SECTION("double"){
		
		basis::field_set<basis::trivial, double> aa(bas, 1);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) aa.matrix()[ii][0] = sqrt(bas.part().local_to_global(ii).value());

		observables::density::normalize(aa, 33.3);

		CHECK(operations::integral_sum(aa) == 33.3_a);
		
	}
	
	SECTION("complex"){
		
		basis::field_set<basis::trivial, complex> aa(bas, 1);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			aa.matrix()[ii][0] = sqrt(bas.part().local_to_global(ii).value())*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii).value()));
		}

		observables::density::normalize(aa, 19.2354);

		CHECK(real(operations::integral_sum(aa)) == 19.2354_a);
		
	}
	
}
#endif
