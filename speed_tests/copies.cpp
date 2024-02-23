/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <catch2/catch_all.hpp>

#include <gpu/array.hpp>
#include <gpu/copy.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("speed_test::copy", "[speed_test::copy]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;

  auto comm = boost::mpi3::environment::get_world_instance();

	double const threshold = 0.1;
	
  SECTION("2D copy performance benchmark"){

		// only run in series
		if(comm.size() == 1){
			
			long nn = 7500;

			// WARMUP POOL MEMORY, IF ANY
			std::vector<gpu::array<complex, 2>> warmup(10, gpu::array<complex, 2>({nn, nn}, complex{1.0, 2.0}));
			CHECK( warmup[5][nn-1][nn-1] == complex{1.0, 2.0} );

			gpu::array<complex, 2> src({nn, nn});
			gpu::array<complex, 2> dest({nn, nn});
			
			gpu::run(nn, nn, [sr = begin(src), de = begin(dest)] GPU_LAMBDA (auto i1, auto i2){
				sr[i2][i1] = 1.0;
				de[i2][i1] = -1.0;
			});
			
			double size = nn*nn*sizeof(complex)/1e9;
			double memcpy_rate;
			
			//MEMCPY
			{
				auto start_time = std::chrono::high_resolution_clock::now();
#if   defined(ENABLE_CUDA)
				cudaMemcpy(raw_pointer_cast(dest.data_elements()), raw_pointer_cast(src.data_elements()), nn*nn*sizeof(complex), cudaMemcpyDeviceToDevice);
#elif defined(ENABLE_HIP)
				hipMemcpy (raw_pointer_cast(dest.data_elements()), raw_pointer_cast(src.data_elements()), nn*nn*sizeof(complex), hipMemcpyDeviceToDevice)==HIP_SUCCESS?void():assert(0);
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
