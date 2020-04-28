/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__TRANSFER
#define OPERATIONS__TRANSFER

/*
 Copyright (C) 2020 Xavier Andrade, Alfredo A. Correa.

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

#include <gpu/run.hpp>
#include <basis/field_set.hpp>

#include <multi/adaptors/fftw.hpp>

#ifdef HAVE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

namespace operations {
	namespace transfer {

		template <class FieldType>
		auto enlarge(FieldType const & source, typename FieldType::basis_type const & new_basis) {

			FieldType destination(new_basis);

			destination = 0.0;
			
			for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
				for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
					for(int iz = 0; iz < source.basis().sizes()[2]; iz++){
						auto i2x = ix;
						auto i2y = iy;
						auto i2z = iz;
						if(ix >= source.basis().sizes()[0]/2) i2x += source.basis().sizes()[0];
						if(iy >= source.basis().sizes()[1]/2) i2y += source.basis().sizes()[1];
						if(iz >= source.basis().sizes()[2]/2) i2z += source.basis().sizes()[2];
						destination.cubic()[i2x][i2y][i2z] = source.cubic()[ix][iy][iz];
					}
				}
			}

			return destination;			
		}

		
		template <class FieldType>
		auto reduce(FieldType const & source, typename FieldType::basis_type const & new_basis) {

			FieldType destination(new_basis);
			
			destination = 0.0;
			
			for(int ix = 0; ix < destination.basis().sizes()[0]; ix++){
				for(int iy = 0; iy < destination.basis().sizes()[1]; iy++){
					for(int iz = 0; iz < destination.basis().sizes()[2]; iz++){	
						auto i2x = ix;
						auto i2y = iy;
						auto i2z = iz;
						if(ix >= destination.basis().sizes()[0]/2) i2x += destination.basis().sizes()[0];
						if(iy >= destination.basis().sizes()[1]/2) i2y += destination.basis().sizes()[1];
						if(iz >= destination.basis().sizes()[2]/2) i2z += destination.basis().sizes()[2];
						destination.cubic()[ix][iy][iz] = source.cubic()[i2x][i2y][i2z];
					}
				}
			}

			return destination;
		}
		
	}
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

namespace transfer_unit_test {
	auto outside(double xx, double ll){
		if(xx < 0.0){
			return xx - -0.5*ll < 1e-10;
		} else {
			return xx - 0.5*ll >= -1e-10;
		}
	}
}

TEST_CASE("function operations::transfer", "[operations::transfer]") {
	using namespace transfer_unit_test;
	using namespace Catch::literals;
	using math::vec3d;
	
	double ecut = 23.0;

	vec3d ll{6.66, 7.77, 9.99};

	
	ions::geometry geo;
	ions::UnitCell cell(vec3d(ll[0], 0.0, 0.0), vec3d(0.0, ll[1], 0.0), vec3d(0.0, 0.0, ll[2]));
	basis::real_space grid(cell, input::basis::cutoff_energy(ecut));

	basis::field<basis::real_space, double> small(grid);

	CHECK(small.basis().rlength()[0] == Approx(ll[0]));
	CHECK(small.basis().rlength()[1] == Approx(ll[1]));
	CHECK(small.basis().rlength()[2] == Approx(ll[2]));
	
	for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
				
				auto ixg = small.basis().cubic_dist(0).local_to_global(ix);
				auto iyg = small.basis().cubic_dist(1).local_to_global(iy);
				auto izg = small.basis().cubic_dist(2).local_to_global(iz);						
				auto rr = small.basis().rvector(ixg, iyg, izg);
				small.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
			}
		}
	}
 	
	auto large = operations::transfer::enlarge(small, grid.enlarge(2));

	CHECK(large.basis().rlength()[0] == Approx(2.0*ll[0]));
	CHECK(large.basis().rlength()[1] == Approx(2.0*ll[1]));
	CHECK(large.basis().rlength()[2] == Approx(2.0*ll[2]));

	long count_large = 0;
	long count_small = 0;
	for(int ix = 0; ix < large.basis().local_sizes()[0]; ix++){
		for(int iy = 0; iy < large.basis().local_sizes()[1]; iy++){
			for(int iz = 0; iz < large.basis().local_sizes()[2]; iz++){
				
				auto ixg = large.basis().cubic_dist(0).local_to_global(ix);
				auto iyg = large.basis().cubic_dist(1).local_to_global(iy);
				auto izg = large.basis().cubic_dist(2).local_to_global(iz);						
				auto rr = large.basis().rvector(ixg, iyg, izg);

				//				std::cout << rr[0] << '\t' << rr[1] << '\t' << rr[2] << '\t' << ll[0] << '\t' << ll[1] << '\t' << ll[2] << std::endl;
				//std::cout << rr[0] << '\t' << ll[0] << '\t' << outside(rr[0], ll[0]) << std::endl;
				
				if(outside(rr[0], ll[0]) or outside(rr[1], ll[1]) or outside(rr[2], ll[2])){
					CHECK(large.cubic()[ix][iy][iz] == 0.0_a);
					count_large++;
				} else {
					//std::cout << ix << '\t' << iy << '\t' << iz << '\t' << large.cubic()[ix][iy][iz] << '\t' << exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]) << std::endl;
					//					CHECK(large.cubic()[ix][iy][iz] == Approx(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])));
					count_small++;
				}
				
			}
		}
	}

	CHECK(count_small == small.basis().size());
	CHECK(count_large > count_small);
	CHECK(count_large == large.basis().size() - count_small);
	
	std::cout << "SMALL SIZE " << small.basis().sizes()[0] << '\t'  << small.basis().sizes()[1] << '\t' << small.basis().sizes()[2] << std::endl;
	std::cout << "LARGE SIZE " << large.basis().sizes()[0] << '\t'  << large.basis().sizes()[1] << '\t' << large.basis().sizes()[2] << std::endl;	
	
	auto small2 = operations::transfer::reduce(large, small.basis());
	
}


#endif

#endif

