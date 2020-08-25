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

#include <basis/fourier_space.hpp>
#include <basis/field_set.hpp>
#include <gpu/run.hpp>
#include <operations/space.hpp>

#include <multi/adaptors/fftw.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/cufft.hpp>
#endif

#include <cassert>

namespace inq {
namespace operations {
namespace transfer {

template <class FieldType>
auto enlarge(FieldType const & source, typename FieldType::basis_type const & new_basis, double const factor = 1.0) {

	assert(not source.basis().part().parallel());
			
	FieldType destination(new_basis);

	destination = 0.0;
			
	for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
		for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
			for(int iz = 0; iz < source.basis().sizes()[2]; iz++){

				auto ii = source.basis().to_symmetric_range(ix, iy, iz);
				auto idest = destination.basis().from_symmetric_range(ii);
						
				destination.cubic()[idest[0]][idest[1]][idest[2]] = factor*source.cubic()[ix][iy][iz];
			}
		}
	}

	return destination;			
}

//////////////////////////////////////////////////////////

template <class Type, class BasisType>
auto enlarge(basis::field_set<BasisType, Type> const & source, BasisType const & new_basis, double const factor = 1.0) {

	assert(not source.basis().part().parallel());
			
	basis::field_set<BasisType, Type> destination(new_basis, source.set_size());

	destination = 0.0;
			
	for(int ix = 0; ix < source.basis().sizes()[0]; ix++){
		for(int iy = 0; iy < source.basis().sizes()[1]; iy++){
			for(int iz = 0; iz < source.basis().sizes()[2]; iz++){
				for(int ist = 0; ist < source.set_part().local_size(); ist++){
							
					auto ii = source.basis().to_symmetric_range(ix, iy, iz);
					auto idest = destination.basis().from_symmetric_range(ii);
							
					destination.cubic()[idest[0]][idest[1]][idest[2]][ist] = factor*source.cubic()[ix][iy][iz][ist];
				}
			}
		}
	}

	return destination;			
}


//////////////////////////////////////////////////////////
		
template <class FieldType>
auto shrink(FieldType const & source, typename FieldType::basis_type const & new_basis, double const factor = 1.0) {

	assert(not source.basis().part().parallel());
			
	FieldType destination(new_basis);
			
	destination = 0.0;
			
	for(int ix = 0; ix < destination.basis().sizes()[0]; ix++){
		for(int iy = 0; iy < destination.basis().sizes()[1]; iy++){
			for(int iz = 0; iz < destination.basis().sizes()[2]; iz++){	

				auto ii = destination.basis().to_symmetric_range(ix, iy, iz);
				auto isource = source.basis().from_symmetric_range(ii);
				destination.cubic()[ix][iy][iz] = factor*source.cubic()[isource[0]][isource[1]][isource[2]];
						
			}
		}
	}

	return destination;
}

//////////////////////////////////////////////////////////
		
template <class Type, class BasisType>
auto shrink(basis::field_set<BasisType, Type> const & source, BasisType const & new_basis, double const factor = 1.0) {

	assert(not source.basis().part().parallel());
			
	basis::field_set<BasisType, Type> destination(new_basis, source.set_size());
			
	destination = 0.0;
			
	for(int ix = 0; ix < destination.basis().sizes()[0]; ix++){
		for(int iy = 0; iy < destination.basis().sizes()[1]; iy++){
			for(int iz = 0; iz < destination.basis().sizes()[2]; iz++){	
				for(int ist = 0; ist < source.set_part().local_size(); ist++){
							
					auto ii = destination.basis().to_symmetric_range(ix, iy, iz);
					auto isource = source.basis().from_symmetric_range(ii);
					destination.cubic()[ix][iy][iz][ist] = factor*source.cubic()[isource[0]][isource[1]][isource[2]][ist];
							
				}
			}
		}
	}

	return destination;
}
		
//////////////////////////////////////////////////////////

template <class FieldType>
auto refine(FieldType const & source, typename basis::real_space const & new_basis){

	static_assert(std::is_same<typename FieldType::basis_type, basis::real_space>::value, "Only implemented for real space");
	assert(new_basis.size() == 8*source.basis().size()); //only a factor of 2 has been tested
	assert(not source.basis().part().parallel());
			
	basis::fourier_space new_fourier_basis(new_basis);
	auto destination_fourier = enlarge(operations::space::to_fourier(source), new_fourier_basis, 1.0/source.basis().size());
	return operations::space::to_real(destination_fourier,  /*normalize = */ false);
}

//////////////////////////////////////////////////////////

template <template<class, class> class FieldType>
auto refine(FieldType<basis::real_space, double> const & source, typename basis::real_space const & new_basis){

	auto complex_refine = refine(source.complex(), new_basis);
	return complex_refine.real();
}

//////////////////////////////////////////////////////////
		
template <class FieldType>
auto coarsen(FieldType const & source, typename basis::real_space const & new_basis){

	assert(8*new_basis.size() == source.basis().size()); //only a factor of 2 has been tested		
	assert(not source.basis().part().parallel());
			
	basis::fourier_space new_fourier_basis(new_basis);
	auto destination_fourier = shrink(operations::space::to_fourier(source), new_fourier_basis, 1.0/source.basis().size());
	return operations::space::to_real(destination_fourier, /*normalize = */ false);
}

//////////////////////////////////////////////////////////

template <template<class, class> class FieldType>		
auto coarsen(FieldType<basis::real_space, double> const & source, typename basis::real_space const & new_basis){

	auto complex_coarsen = coarsen(source.complex(), new_basis);
	return complex_coarsen.real();
}

}

}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef INQ_OPERATIONS_TRANSFER_UNIT_TEST
#undef INQ_OPERATIONS_TRANSFER_UNIT_TEST

#include <operations/integral.hpp>

#include <catch2/catch.hpp>

namespace transfer_unit_test {
	
	auto outside(int ii, int ll){
		if(ii < 0.0){
			return ii < -(ll/2);
		} else {
			return ii >= (ll + 1)/2;
		}
	}
	
}

TEMPLATE_TEST_CASE("function operations::transfer", "[operations::transfer]", double, inq::complex) {

	using namespace transfer_unit_test;
	using namespace inq;
	using namespace Catch::literals;
	using math::vec3d;
	
	double ecut = 23.0;

	vec3d ll{6.66, 7.77, 9.99};

	ions::UnitCell cell(vec3d(ll[0], 0.0, 0.0), vec3d(0.0, ll[1], 0.0), vec3d(0.0, 0.0, ll[2]));
	basis::real_space grid(cell, input::basis::cutoff_energy(ecut));
	
	SECTION("Enlarge and shrink -- field"){
		
		basis::field<basis::real_space, TestType> small(grid);
		
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

					auto ii = large.basis().to_symmetric_range(ixg, iyg, izg);

					if(outside(ii[0], small.basis().sizes()[0]) or outside(ii[1], small.basis().sizes()[1]) or outside(ii[2], small.basis().sizes()[2])){
						CHECK(real(large.cubic()[ix][iy][iz]) == 0.0_a);
						CHECK(imag(large.cubic()[ix][iy][iz]) == 0.0_a);
						count_large++;
					} else {
						auto rr = large.basis().rvector(ix, iy, iz);
						CHECK(real(large.cubic()[ix][iy][iz]) == Approx(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])));
						count_small++;
					}
				
				}
			}
		}

		CHECK(count_small == small.basis().size());
		CHECK(count_large > count_small);
		CHECK(count_large == large.basis().size() - count_small);
	
		auto small2 = operations::transfer::shrink(large, small.basis());
	
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					CHECK(small2.cubic()[ix][iy][iz] == small.cubic()[ix][iy][iz]);
				}
			}
		}
	}

	SECTION("Enlarge and shrink -- field_set"){
		
		basis::field_set<basis::real_space, TestType> small(grid, 5);
		
		CHECK(small.basis().rlength()[0] == Approx(ll[0]));
		CHECK(small.basis().rlength()[1] == Approx(ll[1]));
		CHECK(small.basis().rlength()[2] == Approx(ll[2]));
		
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < small.set_part().local_size(); ist++){
							
						auto ixg = small.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = small.basis().cubic_dist(1).local_to_global(iy);
						auto izg = small.basis().cubic_dist(2).local_to_global(iz);						
						auto rr = small.basis().rvector(ixg, iyg, izg);
						small.cubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
					}
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
					for(int ist = 0; ist < large.set_part().local_size(); ist++){
						
						auto ixg = large.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = large.basis().cubic_dist(1).local_to_global(iy);
						auto izg = large.basis().cubic_dist(2).local_to_global(iz);						
						
						auto ii = large.basis().to_symmetric_range(ixg, iyg, izg);
						
						if(outside(ii[0], small.basis().sizes()[0]) or outside(ii[1], small.basis().sizes()[1]) or outside(ii[2], small.basis().sizes()[2])){
							CHECK(real(large.cubic()[ix][iy][iz][ist]) == 0.0_a);
							CHECK(imag(large.cubic()[ix][iy][iz][ist]) == 0.0_a);							
							count_large++;
						} else {
							auto rr = large.basis().rvector(ix, iy, iz);
							CHECK(real(large.cubic()[ix][iy][iz][ist]) == Approx(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])));
							count_small++;
						}
						
					}
				}
			}
		}

		CHECK(count_small == small.basis().size()*small.set_size());
		CHECK(count_large > count_small);
		CHECK(count_large == large.basis().size()*small.set_size() - count_small);
	
		auto small2 = operations::transfer::shrink(large, small.basis());
	
		for(int ix = 0; ix < small.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < small.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < small.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < small.set_part().local_size(); ist++){
						CHECK(small2.cubic()[ix][iy][iz][ist] == small.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
	}

	SECTION("Mesh refinement -- field"){

		basis::field<basis::real_space, TestType> coarse(grid); 

		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					
					auto ixg = coarse.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = coarse.basis().cubic_dist(1).local_to_global(iy);
					auto izg = coarse.basis().cubic_dist(2).local_to_global(iz);						
					auto rr = coarse.basis().rvector(ixg, iyg, izg);
					coarse.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
				}
			}
		}

		auto fine = operations::transfer::refine(coarse, grid.refine(2));

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		CHECK(real(operations::integral(coarse)) == 109.2100704359_a);
		CHECK(fabs(imag(operations::integral(coarse))) < 1e-14);
		CHECK(real(operations::integral(fine)) == 109.2100704359_a);
		CHECK(fabs(imag(operations::integral(fine))) < 1e-14);

			
		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					
					auto ixg = fine.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = fine.basis().cubic_dist(1).local_to_global(iy);
					auto izg = fine.basis().cubic_dist(2).local_to_global(iz);						
					auto rr = fine.basis().rvector(ixg, iyg, izg);
					CHECK(fabs(real(fine.cubic()[ix][iy][iz])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
					CHECK(fabs(imag(fine.cubic()[ix][iy][iz])) < 5e-3);
				}
			}
		}
		
	}
		
	SECTION("Mesh refinement -- field_set"){

		basis::field_set<basis::real_space, TestType> coarse(grid, 5); 

		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < coarse.set_part().local_size(); ist++){
						
						auto ixg = coarse.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = coarse.basis().cubic_dist(1).local_to_global(iy);
						auto izg = coarse.basis().cubic_dist(2).local_to_global(iz);						
						auto rr = coarse.basis().rvector(ixg, iyg, izg);
						coarse.cubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);

					}
				}
			}
		}

		auto fine = operations::transfer::refine(coarse, grid.refine(2));

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
							
						auto ixg = fine.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = fine.basis().cubic_dist(1).local_to_global(iy);
						auto izg = fine.basis().cubic_dist(2).local_to_global(iz);						
						auto rr = fine.basis().rvector(ixg, iyg, izg);
						CHECK(fabs(real(fine.cubic()[ix][iy][iz][ist])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
						CHECK(fabs(imag(fine.cubic()[ix][iy][iz][ist])) < 5e-3);
						
					}
				}
			}
		}
		
	}
	
	SECTION("Mesh coarsening -- field"){

		auto fine_grid = grid.refine(2);
		
		basis::field<basis::real_space, TestType> fine(fine_grid); 

		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					
					auto ixg = fine.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = fine.basis().cubic_dist(1).local_to_global(iy);
					auto izg = fine.basis().cubic_dist(2).local_to_global(iz);						
					auto rr = fine.basis().rvector(ixg, iyg, izg);
					fine.cubic()[ix][iy][iz] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
				}
			}
		}

		auto coarse = operations::transfer::coarsen(fine, grid);

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		CHECK(real(operations::integral(coarse)) == 109.3027787767_a);
		CHECK(fabs(imag(operations::integral(coarse))) < 1e-14);
		CHECK(real(operations::integral(fine)) == 109.3027787767_a);
		CHECK(fabs(imag(operations::integral(fine))) < 1e-14);

			
		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					
					auto ixg = coarse.basis().cubic_dist(0).local_to_global(ix);
					auto iyg = coarse.basis().cubic_dist(1).local_to_global(iy);
					auto izg = coarse.basis().cubic_dist(2).local_to_global(iz);						
					auto rr = coarse.basis().rvector(ixg, iyg, izg);
					CHECK(fabs(real(coarse.cubic()[ix][iy][iz])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
					CHECK(fabs(imag(coarse.cubic()[ix][iy][iz])) < 5e-3);
				}
			}
		}
		
	}

		SECTION("Mesh coarsening -- field_set"){

		auto fine_grid = grid.refine(2);
		
		basis::field_set<basis::real_space, TestType> fine(fine_grid, 5);

		for(int ix = 0; ix < fine.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fine.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fine.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
						
						auto ixg = fine.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = fine.basis().cubic_dist(1).local_to_global(iy);
						auto izg = fine.basis().cubic_dist(2).local_to_global(iz);						
						auto rr = fine.basis().rvector(ixg, iyg, izg);
						fine.cubic()[ix][iy][iz][ist] = exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2]);
					}
				}
			}
		}

		auto coarse = operations::transfer::coarsen(fine, grid);

		static_assert(std::is_same<decltype(fine.basis()), basis::real_space const &>::value, "The return value must be in real_space");

		CHECK(not (coarse.basis() == fine.basis()));
		CHECK(coarse.basis().refine(2) == fine.basis());
	
		CHECK(8*coarse.basis().size() == fine.basis().size());
		CHECK(coarse.basis().rspacing()[0] == 2.0*fine.basis().rspacing()[0]);
		CHECK(coarse.basis().rspacing()[1] == 2.0*fine.basis().rspacing()[1]);
		CHECK(coarse.basis().rspacing()[2] == 2.0*fine.basis().rspacing()[2]);
		
		for(int ix = 0; ix < coarse.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < coarse.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < coarse.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < fine.set_part().local_size(); ist++){
						
						auto ixg = coarse.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = coarse.basis().cubic_dist(1).local_to_global(iy);
						auto izg = coarse.basis().cubic_dist(2).local_to_global(iz);						
						auto rr = coarse.basis().rvector(ixg, iyg, izg);
						CHECK(fabs(real(coarse.cubic()[ix][iy][iz][ist])/(exp(-rr[0]*rr[0]/ll[0] - rr[1]*rr[1]/ll[1] - rr[2]*rr[2]/ll[2])) - 1.0) < 0.2);
						CHECK(fabs(imag(coarse.cubic()[ix][iy][iz][ist])) < 5e-3);
					}
				}
			}
		}
		
	}
}


#endif

#endif

