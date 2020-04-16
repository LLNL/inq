/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__SPACE
#define OPERATIONS__SPACE

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
namespace space {

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
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis, phi.basis_comm());
	
	basis::field_set<basis::fourier_space, complex> fphi(fourier_basis, phi.set_size(), phi.full_comm());

	if(not real_basis.part().parallel()) {
			
		//DATAOPERATIONS FFT
		fft::dft({true, true, true, false}, phi.cubic(), fphi.cubic(), boost::multi::fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	} else {

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);


		math::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size(), phi.set_part().local_size()});

		namespace multi = boost::multi;
		namespace fft = multi::fft;

		auto const real_x = real_basis.local_sizes();
		fft::dft({false, true, true, false}, phi.cubic(), tmp({0, real_x[0]}, {0, real_x[1]}, {0, real_x[2]}), fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		math::array<complex, 5> buffer = tmp.unrotated(2).partitioned(phi.basis_comm().size()).transposed().rotated().transposed().rotated();

		assert(std::get<4>(sizes(buffer)) == phi.set_part().local_size());
		
		tmp.clear();
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, phi.basis_comm().get());

		auto const fourier_x = fourier_basis.local_sizes();
		fft::dft({true, false, false, false}, buffer.flatted()({0, fourier_x[0]}, {0, fourier_x[1]}, {0, fourier_x[2]}), fphi.cubic(), fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
	}
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
	return fphi;
}

///////////////////////////////////////////////////////////////

basis::field_set<basis::real_space, complex> to_real(const basis::field_set<basis::fourier_space, complex> & fphi){
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis, fphi.basis_comm());
	
	basis::field_set<basis::real_space, complex> phi(real_basis, fphi.set_size(), fphi.full_comm());

	if(not real_basis.part().parallel()) {

		//DATAOPERATIONS FFT
		fft::dft({true, true, true, false}, fphi.cubic(), phi.cubic(), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
	} else {

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		
		math::array<complex, 5> buffer({fphi.basis_comm().size(), xblock, real_basis.local_sizes()[1], zblock, fphi.set_part().local_size()});
		
		namespace multi = boost::multi;
		namespace fft = multi::fft;
		fft::dft({true, true, false, false}, fphi.cubic(), buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, fphi.basis_comm().get());

		math::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size(), fphi.set_part().local_size()});

		tmp.unrotated(2).partitioned(phi.basis_comm().size()).transposed().rotated().transposed().rotated() = buffer;
		
		fft::dft({false, false, true, false}, tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), phi.cubic(), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
	}
	
	//DATAOPERATIONS GPU::RUN 1D
	gpu::run(fphi.basis().part().local_size()*phi.set_part().local_size(),
					 [phip = (complex *) phi.data(), norm_factor = (double) phi.basis().size()] GPU_LAMBDA (auto ii){
						 phip[ii] = phip[ii]/norm_factor;
					 });
	
	return phi;
}

///////////////////////////////////////////////////////////////

basis::field<basis::fourier_space, complex> to_fourier(const basis::field<basis::real_space, complex> & phi){
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis, phi.basis_comm());
	
	basis::field<basis::fourier_space, complex> fphi(fourier_basis, phi.basis_comm());
	
	if(not real_basis.part().parallel()) {
		
		fft::dft(phi.cubic(), fphi.cubic(),fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	} else {

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);

		math::array<complex, 3> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size()});

		namespace multi = boost::multi;
		namespace fft = multi::fft;
		fft::dft({false, true, true}, phi.cubic(), tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		math::array<complex, 4> buffer = tmp.unrotated().partitioned(phi.basis_comm().size()).transposed().rotated();
		
		tmp.clear();
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, phi.basis_comm().get());
		
		fft::dft({true, false, false}, buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fphi.cubic(), fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	}

	return fphi;
	
}

///////////////////////////////////////////////////////////////			
	
basis::field<basis::real_space, complex> to_real(const basis::field<basis::fourier_space, complex> & fphi){
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis, fphi.basis_comm());

	basis::field<basis::real_space, complex> phi(real_basis, fphi.basis_comm());

	if(not real_basis.part().parallel()) {
		
		fft::dft(fphi.cubic(), phi.cubic(), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	} else {

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		
		math::array<complex, 4> buffer({fphi.basis_comm().size(), xblock, real_basis.local_sizes()[1], zblock});

		namespace multi = boost::multi;
		namespace fft = multi::fft;				
		fft::dft({true, true, false}, fphi.cubic(), buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, fphi.basis_comm().get());

		math::array<complex, 3> tmp({xblock, real_basis.local_sizes()[1], zblock*phi.basis_comm().size()});

		tmp.unrotated().partitioned(phi.basis_comm().size()).transposed().rotated() = buffer;
		
		fft::dft({false, false, true}, tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), 	phi.cubic(), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	}
	
	return phi;
}

}}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::space", "[operations::space]") {

	using namespace Catch::literals;
	using math::vec3d;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance());

	auto basis_comm = cart_comm.axis(1);
	
	double ecut = 23.0;
	double ll = 6.66;
	
	ions::geometry geo;
	ions::UnitCell cell(vec3d(ll, 0.0, 0.0), vec3d(0.0, ll, 0.0), vec3d(0.0, 0.0, ll));
	basis::real_space rs(cell, input::basis::cutoff_energy(ecut), basis_comm);
	
	basis::field_set<basis::real_space, complex> phi(rs, 7, cart_comm);
	
	SECTION("Zero"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++) phi.cubic()[ix][iy][iz][ist] = 0.0;
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(fphi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= fphi.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++)	diff += fabs(phi.cubic()[ix][iy][iz][ist]);
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= phi2.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
	}
	
	SECTION("Gaussian"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
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
		for(int ix = 0; ix < fphi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().local_sizes()[2]; iz++){
					double g2 = fphi.basis().g2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						double sigma = 0.5*(ist + 1);
						diff += fabs(fphi.cubic()[ix][iy][iz][ist] - pow(M_PI/sigma, 3.0/2.0)*exp(-0.25*g2/sigma));
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= fphi.cubic().num_elements();

		//not sure what is wrong here
		std::cout << "DIFF1 " << diff << std::endl;

		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(phi.cubic()[ix][iy][iz][ist] - phi2.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});		

		diff /= phi2.cubic().num_elements();
		
		REQUIRE(diff < 1e-15);
		
	}
	
}


#endif

#endif

