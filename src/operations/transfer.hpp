/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__TRANSFER
#define OPERATIONS__TRANSFER

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa.

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
		
	}
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::transfer", "[operations::transfer]") {

	using namespace Catch::literals;
	using math::vec3d;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance());

	auto basis_comm = cart_comm.axis(1);
	
	double ecut = 23.0;
	double ll = 6.66;
	
	ions::geometry geo;
	ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
	basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);
	
}


#endif

#endif

