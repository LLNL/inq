/* -*- indent-tabs-mode: t -*- */
/*
 Copyright (C) 2022 Xavier Andrade, Alfredo A. Correa

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

	double const threshold = 0.1;
	
  SECTION("2D copy performance benchmark"){

		if(comm.root()){
			
			long nn = 8000;

			// WARMUP POOL MEMORY, IF ANY
			std::vector<math::array<complex, 2>> warmup(10, math::array<complex, 2>({nn, nn}, complex{1.0, 2.0}));
			CHECK( warmup[5][nn-1][nn-1] == complex{1.0, 2.0} );

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

			std::cout << "memcpy rate        = " << memcpy_rate << " GB/s " << std::endl;

			{ //GPU::COPY
				auto start_time = std::chrono::high_resolution_clock::now();				
				gpu::copy(nn, nn, src, dest);
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				double rate = size/time.count();
				double ratio = rate/memcpy_rate;
				
				std::cout << "gpu::copy rate     = " << rate << " GB/s " << "(ratio = " << ratio << ")" << std::endl;
				CHECK(ratio >= threshold);
			}

			{ //MULTI COPY CONSTRUCTOR
				{
					auto warmup = src;
				}
				auto start_time = std::chrono::high_resolution_clock::now();
				auto dest2 = src;
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				double rate = size/time.count();
				double ratio = rate/memcpy_rate;
				
				std::cout << "constructor rate   = " << rate << " GB/s " << "(ratio = " << ratio << ")" << std::endl;
				CHECK(ratio >= threshold);
			}
			
			{ //MULTI ASSIGNMENT
				auto start_time = std::chrono::high_resolution_clock::now();				
				dest = src;
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				double rate = size/time.count();
				double ratio = rate/memcpy_rate;
				
				std::cout << "assignment rate    = " << rate << " GB/s " << "(ratio = " << ratio << ")" << std::endl;
				CHECK(ratio >= threshold);
			}

			{ //MULTI SUB ARRAY ASSIGNMENT
				auto start_time = std::chrono::high_resolution_clock::now();				
				dest({0, nn - 2}, {0, nn - 2}) = src({2, nn}, {2, nn});
				std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
				double rate = size/time.count();
				double ratio = rate/memcpy_rate;
				
				std::cout << "sub array rate     = " << rate << " GB/s " << "(ratio = " << ratio << ")" << std::endl;
				CHECK(ratio >= threshold);
			}
		}
		
  }
  
}
