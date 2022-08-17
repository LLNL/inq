/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__EXPONENTIAL
#define INQ__OPERATIONS__EXPONENTIAL

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

#include <operations/shift.hpp>

#include <utils/profiling.hpp>

namespace inq {
namespace operations {

	template <class operator_type, class field_set_type>
	void exponential_in_place(const operator_type & ham, typename field_set_type::element_type const & factor, field_set_type & phi, int const order = 4){

    CALI_CXX_MARK_FUNCTION;

		auto hnphi = ham(phi);
		
		typename field_set_type::element_type coeff = 1.0;
		for(int iter = 1; iter <= order; iter++){
			if(iter > 1) hnphi = ham(hnphi);
			coeff *= factor/typename complex::value_type(iter);
			shift(coeff, hnphi, phi);
		}
  }

	///////////////////////////////////////////////////////////////////////////
	
	template <class operator_type, class field_set_type>
	auto exponential(const operator_type & ham, typename field_set_type::element_type const & factor, field_set_type const & phi, int const order = 4){

		CALI_CXX_MARK_FUNCTION;
		
    auto expphi = phi;
		exponential_in_place(ham, factor, expphi, order);
		return expphi;
  }

	///////////////////////////////////////////////////////////////////////////
	
	template <class operator_type, class field_set_type>
	auto exponential_2_for_1(const operator_type & ham, typename field_set_type::element_type const & factor1, typename field_set_type::element_type const & factor2, field_set_type & phi, int const order = 4){

		CALI_CXX_MARK_FUNCTION;

    auto expphi = phi;
		auto hnphi = ham(phi);
		
		typename field_set_type::element_type coeff1 = 1.0;
		typename field_set_type::element_type coeff2 = 1.0;
		for(int iter = 1; iter <= order; iter++){
			if(iter > 1) hnphi = ham(hnphi);
			coeff1 *= factor1/typename complex::value_type(iter);
			coeff2 *= factor2/typename complex::value_type(iter);
			shift(coeff1, hnphi, expphi);
			shift(coeff2, hnphi, phi);
		}

		return expphi;		
  }

}
}

#ifdef INQ_OPERATIONS_EXPONENTIAL_UNIT_TEST
#undef INQ_OPERATIONS_EXPONENTIAL_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>
#include <ions/unit_cell.hpp>
#include <basis/trivial.hpp>
#include <operations/matrix_operator.hpp>

TEST_CASE("operations::exponential", "[operations::exponential]") {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  const int npoint = 100;
  const int nvec = 12;
  
  basis::trivial bas(npoint);

	SECTION("Diagonal double"){
		
		math::array<double, 2> diagonal_matrix({npoint, npoint});
		
		for(int ip = 0; ip < npoint; ip++){
			for(int jp = 0; jp < npoint; jp++){
				diagonal_matrix[ip][jp] = 0.0;
				if(ip == jp) diagonal_matrix[ip][jp] = ip;
			}
		}
		
		operations::matrix_operator<double> diagonal(std::move(diagonal_matrix));
		
		basis::field_set<basis::trivial, double> phi(bas, nvec);
		
		phi = 0.0;
		
		for(int ivec = 0; ivec < nvec; ivec++) phi.matrix()[ivec][ivec] = 1.0;
		
		{
			auto expphi = operations::exponential(diagonal, -0.1, phi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
						CHECK(expphi.matrix()[ivec][ivec] == Approx(exp(-0.1*ivec)));
					} else {
					CHECK(expphi.matrix()[ip][ivec] == 0.0_a);
					}
				}
			}
		}

		{
			auto expphi = phi;
			operations::exponential_in_place(diagonal, -0.1, expphi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
					CHECK(expphi.matrix()[ivec][ivec] == Approx(exp(-0.1*ivec)));
					} else {
						CHECK(expphi.matrix()[ip][ivec] == 0.0_a);
					}
				}
			}
		}

		{
			auto expphi = operations::exponential_2_for_1(diagonal, -0.1, -0.2, phi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
						CHECK(phi.matrix()[ivec][ivec] == Approx(exp(-0.2*ivec)));
						CHECK(expphi.matrix()[ivec][ivec] == Approx(exp(-0.1*ivec)));
					} else {
						CHECK(phi.matrix()[ip][ivec] == 0.0_a);
						CHECK(expphi.matrix()[ip][ivec] == 0.0_a);
					}
				}
			}
		}
		
	}
	
	SECTION("Diagonal complex"){
		
		math::array<complex, 2> diagonal_matrix({npoint, npoint});
		
		for(int ip = 0; ip < npoint; ip++){
			for(int jp = 0; jp < npoint; jp++){
				diagonal_matrix[ip][jp] = 0.0;
				if(ip == jp) diagonal_matrix[ip][jp] = ip;
			}
		}
		
		operations::matrix_operator<complex> diagonal(std::move(diagonal_matrix));
		
		basis::field_set<basis::trivial, complex> phi(bas, nvec);
		
		phi = 0.0;
		
		for(int ivec = 0; ivec < nvec; ivec++) phi.matrix()[ivec][ivec] = 1.0;

		{
			auto expphi = operations::exponential(diagonal, complex(0.0, -0.1), phi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
						CHECK(real(expphi.matrix()[ivec][ivec]) == Approx(real(exp(complex(0.0, -0.1*ivec)))));
						CHECK(imag(expphi.matrix()[ivec][ivec]) == Approx(imag(exp(complex(0.0, -0.1*ivec)))));
					} else {
						CHECK(real(expphi.matrix()[ip][ivec]) == 0.0_a);
						CHECK(imag(expphi.matrix()[ip][ivec]) == 0.0_a);
					}
				}
			}
		}

		{
			auto expphi = phi;
			operations::exponential_in_place(diagonal, complex(0.0, -0.1), expphi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
						CHECK(real(expphi.matrix()[ivec][ivec]) == Approx(real(exp(complex(0.0, -0.1*ivec)))));
						CHECK(imag(expphi.matrix()[ivec][ivec]) == Approx(imag(exp(complex(0.0, -0.1*ivec)))));
					} else {
						CHECK(real(expphi.matrix()[ip][ivec]) == 0.0_a);
						CHECK(imag(expphi.matrix()[ip][ivec]) == 0.0_a);
					}
				}
			}
			
		}
		
		{
			auto expphi = operations::exponential_2_for_1(diagonal, complex(0.0, -0.1), complex(0.0, -0.2), phi, 16);
			
			for(int ivec = 0; ivec < nvec; ivec++) {
				for(int ip = 0; ip < npoint; ip++){
					if(ip == ivec){
						CHECK(real(phi.matrix()[ivec][ivec]) == Approx(real(exp(complex(0.0, -0.2*ivec)))));
						CHECK(imag(phi.matrix()[ivec][ivec]) == Approx(imag(exp(complex(0.0, -0.2*ivec)))));
						CHECK(real(expphi.matrix()[ivec][ivec]) == Approx(real(exp(complex(0.0, -0.1*ivec)))));
						CHECK(imag(expphi.matrix()[ivec][ivec]) == Approx(imag(exp(complex(0.0, -0.1*ivec)))));
					} else {
						CHECK(real(phi.matrix()[ip][ivec]) == 0.0_a);
						CHECK(imag(phi.matrix()[ip][ivec]) == 0.0_a);
						CHECK(real(expphi.matrix()[ip][ivec]) == 0.0_a);
						CHECK(imag(expphi.matrix()[ip][ivec]) == 0.0_a);						
					}
				}
			}
		}

	}
	
}

#endif
#endif
