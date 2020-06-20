/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SPACE
#define INQ__OPERATIONS__SPACE

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

namespace inq {
namespace operations {
namespace space {

void zero_outside_sphere(const basis::field<basis::fourier_space, complex> & fphi){
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), bas = fphi.basis()] GPU_LAMBDA
					 (auto iz, auto iy, auto ix){
						 if(bas.outside_sphere(bas.g2(ix, iy, iz))) fphicub[ix][iy][iz] = complex(0.0);
					 });
}

///////////////////////////////////////////////////////////////
		
void zero_outside_sphere(const basis::field_set<basis::fourier_space, complex> & fphi){
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(fphi.set_part().local_size(), fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), bas = fphi.basis()] GPU_LAMBDA
					 (auto ist, auto iz, auto iy, auto ix){
						 if(bas.outside_sphere(bas.g2(ix, iy, iz))) fphicub[ix][iy][iz][ist] = complex(0.0);
					 });
}

///////////////////////////////////////////////////////////////

template <class Comm, class InArray4D, class OutArray4D>
void to_fourier(basis::real_space const & real_basis, basis::fourier_space const & fourier_basis, Comm & comm, InArray4D const & array_rs, OutArray4D && array_fs) {
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	if(not real_basis.part().parallel()) {
		//DATAOPERATIONS FFT
		fft::dft({true, true, true, false}, array_rs, array_fs, boost::multi::fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
	} else {
		
		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);
 		auto last_dim = std::get<3>(sizes(array_rs));

		math::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*comm.size(), last_dim});

		namespace multi = boost::multi;
		namespace fft = multi::fft;

		auto const real_x = real_basis.local_sizes();
		fft::dft({false, true, true, false}, array_rs, tmp({0, real_x[0]}, {0, real_x[1]}, {0, real_x[2]}), fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		math::array<complex, 5> buffer = tmp.unrotated(2).partitioned(comm.size()).transposed().rotated().transposed().rotated();

		assert(std::get<4>(sizes(buffer)) == last_dim);
		
		tmp.clear();
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, comm.get());

		auto const fourier_x = fourier_basis.local_sizes();
		fft::dft({true, false, false, false}, buffer.flatted()({0, fourier_x[0]}, {0, fourier_x[1]}, {0, fourier_x[2]}), array_fs, fft::forward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	}
	
}
		
///////////////////////////////////////////////////////////////

template <class Comm, class InArray4D, class OutArray4D>
void to_real(basis::fourier_space const & fourier_basis, basis::real_space const & real_basis, Comm & comm, InArray4D const & array_fs, OutArray4D && array_rs) {
	namespace multi = boost::multi;
	namespace fft = multi::fft;

	if(not real_basis.part().parallel()) {

		//DATAOPERATIONS FFT
		fft::dft({true, true, true, false}, array_fs, array_rs, fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
	} else {

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		auto last_dim = std::get<3>(sizes(array_fs));
		
		math::array<complex, 5> buffer({comm.size(), xblock, real_basis.local_sizes()[1], zblock, last_dim});
		
		namespace multi = boost::multi;
		namespace fft = multi::fft;
		fft::dft({true, true, false, false}, array_fs, buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
		
		MPI_Alltoall(MPI_IN_PLACE, buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, static_cast<complex *>(buffer.data()), buffer[0].num_elements(), MPI_CXX_DOUBLE_COMPLEX, comm.get());

		math::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*comm.size(), last_dim});

		tmp.unrotated(2).partitioned(comm.size()).transposed().rotated().transposed().rotated() = buffer;
		
		fft::dft({false, false, true, false}, tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]}), array_rs, fft::backward);
#ifdef HAVE_CUDA
		cudaDeviceSynchronize();
#endif
	}

}

///////////////////////////////////////////////////////////////
		
basis::field_set<basis::fourier_space, complex> to_fourier(const basis::field_set<basis::real_space, complex> & phi){

	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis, phi.basis_comm());
	
	basis::field_set<basis::fourier_space, complex> fphi(fourier_basis, phi.set_size(), phi.full_comm());

	to_fourier(real_basis, fourier_basis, phi.basis_comm(), phi.cubic(), fphi.cubic());
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
	
	return fphi;
}

///////////////////////////////////////////////////////////////

basis::field_set<basis::real_space, complex> to_real(const basis::field_set<basis::fourier_space, complex> & fphi, bool const normalize = true){

	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis, fphi.basis_comm());
	
	basis::field_set<basis::real_space, complex> phi(real_basis, fphi.set_size(), fphi.full_comm());

	to_real(fourier_basis, real_basis, phi.basis_comm(), fphi.cubic(), phi.cubic());
 
	if(normalize){
		//DATAOPERATIONS GPU::RUN 1D
		gpu::run(fphi.basis().part().local_size()*phi.set_part().local_size(),
						 [phip = (complex *) phi.data(), norm_factor = (double) phi.basis().size()] GPU_LAMBDA (auto ii){
							 phip[ii] = phip[ii]/norm_factor;
						 });
	}
	return phi;
}

///////////////////////////////////////////////////////////////

basis::field<basis::fourier_space, complex> to_fourier(const basis::field<basis::real_space, complex> & phi){

	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis, phi.basis_comm());
	
	basis::field<basis::fourier_space, complex> fphi(fourier_basis, phi.basis_comm());

	to_fourier(real_basis, fourier_basis, phi.basis_comm(), phi.hypercubic(), fphi.hypercubic());
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
			
	return fphi;
	
}

///////////////////////////////////////////////////////////////			
	
basis::field<basis::real_space, complex> to_real(const basis::field<basis::fourier_space, complex> & fphi, bool normalize = true){

	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis, fphi.basis_comm());

	basis::field<basis::real_space, complex> phi(real_basis, fphi.basis_comm());

	to_real(fourier_basis, real_basis, phi.basis_comm(), fphi.hypercubic(), phi.hypercubic());
 
	if(normalize){
		gpu::run(phi.linear().size(),
						 [phil = begin(phi.linear()), factor = 1.0/phi.basis().size()] GPU_LAMBDA (auto ip){
							 phil[ip] = factor*phil[ip];
						 });
	}
			
	return phi;
}

}
}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::space", "[operations::space]") {

	using namespace inq;
	using namespace Catch::literals;
	using math::vec3d;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

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

		CHECK(diff < 1e-15);
		
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

		CHECK(diff < 1e-15);
		
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
		
		CHECK(diff < 1e-15);
		
	}
	
}


#endif

#endif

