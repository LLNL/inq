/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__DENSITY__NORMALIZE
#define INQ__DENSITY__NORMALIZE

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

#include <operations/integral.hpp>

#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <math/complex.hpp>

namespace inq {
namespace density {

template <class FieldType>
void normalize(FieldType & density, const double & total_charge){
	
	auto qq = operations::integral(density);
	assert(fabs(qq) > 1e-16);
	for(int i = 0; i < density.basis().part().local_size(); i++) density.linear()[i] *= total_charge/qq;
	
}

}
}

#ifdef INQ_density_normalize_UNIT_TEST
#define INQ_TESTED

#include <basis/trivial.hpp>

#include <catch2/catch.hpp>

TEST_CASE("function density::normalize", "[density::normalize]") {

	using namespace inq;
	using namespace Catch::literals;

	const int npoint = 100;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	basis::trivial bas(npoint, comm);
	
	SECTION("double"){
		
		basis::field<basis::trivial, double> aa(bas);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++) aa.linear()[ii] = sqrt(bas.part().local_to_global(ii));

		density::normalize(aa, 33.3);

		CHECK(operations::integral(aa) == 33.3_a);
		
	}
	
	SECTION("complex"){
		
		basis::field<basis::trivial, complex> aa(bas);

		for(int ii = 0; ii < aa.basis().part().local_size(); ii++){
			aa.linear()[ii] = sqrt(bas.part().local_to_global(ii))*exp(complex(0.0, M_PI/65.0*bas.part().local_to_global(ii)));
		}

		density::normalize(aa, 19.2354);

		CHECK(real(operations::integral(aa)) == 19.2354_a);
		
	}
	
}


#endif

#endif
