/* -*- indent-tabs-mode: t -*- */
/*
 Copyright (C) 2022 Xavier Andrade

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

#include <catch2/catch_all.hpp>

#include <math/array.hpp>
#include <gpu/copy.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("speed_test::copy", "[speed_test::copy]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  auto comm = boost::mpi3::environment::get_world_instance();

  SECTION("2D copy performance benchmark"){

		if(comm.root()){
			
			long nn = 8000;
			
			math::array<complex, 2> src({nn, nn});
			math::array<complex, 2> dest({nn, nn});
			
			gpu::run(nn, nn, [sr = begin(src), de = begin(dest)] GPU_LAMBDA (auto i1, auto i2){
				sr[i2][i1] = 1.0;
				de[i2][i1] = -1.0;
			});
			
			double size = nn*nn*sizeof(complex)/1e9;
			double memcpy_rate;
			
			//MEMCPY
			{
				auto start_time = std::chrono::high_resolution_clock::now();
#ifdef ENABLE_CUDA
				cudaMemcpy(raw_pointer_cast(dest.data_elements()), raw_pointer_cast(src.data_elements()), nn*nn*sizeof(complex), cudaMemcpyDeviceToDevice);
#else
				memcpy(dest.data_elements(), src.data_elements(), nn*nn*sizeof(complex));
#endif
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				memcpy_rate = size/time.count();
			}

			std::cout << "memcpy    rate = " << memcpy_rate << " GB/s " << std::endl;

			//GPU::COPY
			{
				auto start_time = std::chrono::high_resolution_clock::now();				
				gpu::copy(nn, nn, src, dest);
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				double rate = size/time.count();
				double ratio = rate/memcpy_rate;
				
				std::cout << "gpu::copy rate = " << rate << " GB/s " << "(ratio = " << ratio << ")" << std::endl;
																																				 
				CHECK(ratio >= 0.25);
			}

		}
		
  }
  
}
