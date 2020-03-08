/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SPACE
#define OPERATIONS__SPACE

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

#include <gpu/run.hpp>
#include <basis/field_set.hpp>

#include <multi/adaptors/fftw.hpp>

#include <cassert>
#include <fftw3.h>

#ifdef HAVE_CUDA
#include <cufft.h>
#endif

namespace operations {

  namespace space {

#ifdef HAVE_CUDA
		template <class field_type>
		auto cuda_fft_plan(const field_type & phi){

			// the information about the layout can be found here:
			//
			//   https://docs.nvidia.com/cuda/cufft/index.html#advanced-data-layout
			//
			// Essentially the access is:
			//
			//   input[b*idist + ((x*inembed[1] + y)*inembed[2] + z)*istride]
		
			int nn[3] = {phi.basis().sizes()[0], phi.basis().sizes()[1], phi.basis().sizes()[2]};
			
			cufftHandle plan;
			
			auto res = cufftPlanMany(/* plan = */ &plan,
															 /* rank = */ 3,
															 /* n = */ nn,
															 /* inembed = */ nn,
															 /* istride = */ phi.set_part().local_size(),
															 /* idist = */ 1,
															 /* onembed = */ nn,
															 /* ostride = */ phi.set_part().local_size(),
															 /* odist =*/ 1,
															 /* type = */ CUFFT_Z2Z,
															 /* batch = */ phi.set_part().local_size());

			assert(res == CUFFT_SUCCESS);

			return plan;
			
		}
		
#endif

		void zero_outside_sphere(const basis::field_set<basis::fourier_space, complex> & fphi){

			//DATAOPERATIONS GPU::RUN 4D
			gpu::run(fphi.set_part().local_size(), fphi.basis().sizes()[2], fphi.basis().sizes()[1], fphi.basis().sizes()[0],
							 [fphicub = begin(fphi.cubic()), bas = fphi.basis()] GPU_LAMBDA
							 (auto ist, auto iz, auto iy, auto ix){
								 if(bas.outside_sphere(bas.g2(ix, iy, iz))) fphicub[ix][iy][iz][ist] = complex(0.0);
							 });
		}
		
		///////////////////////////////////////////////////////////////
		
    basis::field_set<basis::fourier_space, complex> to_fourier(const basis::field_set<basis::real_space, complex> & phi){
      
			basis::field_set<basis::fourier_space, complex> fphi(phi.basis(), phi.set_size(), phi.full_comm());


			//DATAOPERATIONS FFT
#ifdef HAVE_CUDA

			auto plan = cuda_fft_plan(phi);
			
			auto res = cufftExecZ2Z(plan, (cufftDoubleComplex *) raw_pointer_cast(phi.data()),
															(cufftDoubleComplex *) raw_pointer_cast(fphi.data()), CUFFT_FORWARD);

			assert(res == CUFFT_SUCCESS);
			
			cudaDeviceSynchronize();
			
			cufftDestroy(plan);
						
#else
			boost::multi::fftw::dft({true, true, true, false}, phi.cubic(), fphi.cubic(), boost::multi::fftw::forward);
#endif
			
			if(fphi.basis().spherical()) zero_outside_sphere(fphi);
      
      return fphi;
    }

		///////////////////////////////////////////////////////////////
		    
		basis::field_set<basis::real_space, complex> to_real(const basis::field_set<basis::fourier_space, complex> & fphi){

			basis::field_set<basis::real_space, complex> phi(fphi.basis(), fphi.set_size(), fphi.full_comm());

			//DATAOPERATIONS FFT
#ifdef HAVE_CUDA

			auto plan = cuda_fft_plan(phi);
			
			auto res = cufftExecZ2Z(plan, (cufftDoubleComplex *) raw_pointer_cast(fphi.data()),
															(cufftDoubleComplex *) raw_pointer_cast(phi.data()), CUFFT_INVERSE);

			assert(res == CUFFT_SUCCESS);
			
			cudaDeviceSynchronize();
			
			cufftDestroy(plan);
			
#else
			boost::multi::fftw::dft({true, true, true, false}, fphi.cubic(), phi.cubic(), boost::multi::fftw::backward);
#endif
			
			double norm_factor = phi.basis().size();

			//DATAOPERATIONS LOOP + GPU::RUN 1D
#ifdef HAVE_CUDA
			auto phip = raw_pointer_cast(phi.data());
			
			gpu::run(fphi.basis().size()*phi.set_part().local_size(),
							 [=] __device__ (auto ii){
								 phip[ii] = phip[ii]/norm_factor;
							 });
#else
			for(long ii = 0; ii < fphi.basis().size()*phi.set_part().local_size(); ii++) phi.data()[ii] /= norm_factor;
#endif
			return phi;
    }

		///////////////////////////////////////////////////////////////

		basis::field<basis::fourier_space, complex> to_fourier(const basis::field<basis::real_space, complex> & phi){
			namespace fftw = boost::multi::fftw;

			auto & real_basis = phi.basis();
			basis::fourier_space fourier_basis(real_basis, phi.basis_comm());
			
			basis::field<basis::fourier_space, complex> fphi(fourier_basis, phi.basis_comm());
			
			if(not real_basis.part().parallel()) {
				
				fftw::dft(phi.cubic(), fphi.cubic(),fftw::forward);
				
			} else {

				int xblock = real_basis.cubic_dist(0).block_size();
				int zblock = fourier_basis.cubic_dist(2).block_size();
				assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);

				math::array<complex, 3> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size()});

				fftw::dft({false, true, true}, phi.cubic(), tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), fftw::forward);
				
				math::array<complex, 4> buffer = tmp.unrotated().partitioned(phi.basis_comm().size()).transposed().rotated();
				
				tmp.clear();
				
				MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, &phi.basis_comm());
				
				fftw::dft({true, false, false}, buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fphi.cubic(), fftw::forward);
				
			}

			return fphi;
			
		}

		///////////////////////////////////////////////////////////////			
				
		basis::field<basis::real_space, complex> to_real(const basis::field<basis::fourier_space, complex> & fphi){
			namespace fftw = boost::multi::fftw;

			auto & fourier_basis = fphi.basis();
			basis::real_space real_basis(fourier_basis, fphi.basis_comm());

			basis::field<basis::real_space, complex> phi(real_basis, fphi.basis_comm());
		
			if(not real_basis.part().parallel()) {
				
				fftw::dft(fphi.cubic(), phi.cubic(), fftw::backward);
				
			} else {

				int xblock = real_basis.cubic_dist(0).block_size();
				int zblock = fourier_basis.cubic_dist(2).block_size();
				
				math::array<complex, 4> buffer({fphi.basis_comm().size(), xblock, real_basis.local_sizes()[1], zblock});
				
				fftw::dft({true, true, false}, fphi.cubic(), buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fftw::backward);
									
				MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, &fphi.basis_comm());

				math::array<complex, 3> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size()});

				tmp.unrotated().partitioned(phi.basis_comm().size()).transposed().rotated() = buffer;
				
				fftw::dft({false, false, true}, tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), 	phi.cubic(), fftw::backward);

			}
			
			return phi;
		}
	}
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::space", "[operations::space]") {

	using namespace Catch::literals;
	using math::vec3d;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance());

	auto set_comm = cart_comm.axis(0);
	auto basis_comm = cart_comm.axis(1);
	
	double ecut = 23.0;
	double ll = 6.66;
	
	ions::geometry geo;
	ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
	basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);
	
	basis::field_set<basis::real_space, complex> phi(rs, 7, cart_comm);
	
	SECTION("Zero"){
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++) phi.cubic()[ix][iy][iz][ist] = 0.0;
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(fphi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		diff /= fphi.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++)	diff += fabs(phi.cubic()[ix][iy][iz][ist]);
				}
			}
		}

		diff /= phi2.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
	}
	
	SECTION("Gaussian"){
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					double r2 = rs.r2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						double sigma = 0.5*(ist + 1);
						phi.cubic()[ix][iy][iz][ist] = exp(-sigma*r2);
					}
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().sizes()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().sizes()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().sizes()[2]; iz++){
					double g2 = fphi.basis().g2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						double sigma = 0.5*(ist + 1);
						diff += fabs(fphi.cubic()[ix][iy][iz][ist] - pow(M_PI/sigma, 3.0/2.0)*exp(-0.25*g2/sigma));
					}
				}
			}
		}
		
		diff /= fphi.cubic().num_elements();

		//not sure what is wrong here
		std::cout << "DIFF1 " << diff << std::endl;

		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(phi.cubic()[ix][iy][iz][ist] - phi2.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= phi2.cubic().num_elements();
		
		REQUIRE(diff < 1e-15);
		
	}
	
}


#endif

#endif
