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
			
			std::chrono::duration<double> memcpy_time, gpucopy_time;
			
			{
				auto iter_start_time = std::chrono::high_resolution_clock::now();
				gpu::copy(nn, nn, src, dest);
				gpucopy_time = std::chrono::high_resolution_clock::now() - iter_start_time;
			}
			
			{
				auto iter_start_time = std::chrono::high_resolution_clock::now();
#ifdef ENABLE_CUDA
				cudaMemcpy(raw_pointer_cast(dest.data_elements()), raw_pointer_cast(src.data_elements()), nn*nn*sizeof(complex), cudaMemcpyDeviceToDevice);
#else
				memcpy(dest.data_elements(), src.data_elements(), nn*nn*sizeof(complex));
#endif
				memcpy_time = std::chrono::high_resolution_clock::now() - iter_start_time;
			}
			
			double size = nn*nn*sizeof(complex)/1e9;
			double ratio = memcpy_time.count()/gpucopy_time.count();
			
			std::cout << "memcpy transfer rate = " << size/memcpy_time.count() << "GB/s " << std::endl;
			std::cout << "gpu::copy transfer rate = " << size/gpucopy_time.count() << "GB/s " << std::endl;
			std::cout << "ratio = " << ratio << std::endl;
			
			CHECK(ratio >= 0.25);

		}
		
  }
  
}
