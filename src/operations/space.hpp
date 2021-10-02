/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SPACE
#define INQ__OPERATIONS__SPACE

/*
 Copyright (C) 2019-2021 Xavier Andrade, Alfredo A. Correa.

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

#include <inq_config.h> //for ENABLE_HEFFTE

#include <gpu/run.hpp>
#include <gpu/alltoall.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/fourier_space.hpp>
#include <states/orbital_set.hpp>

#ifdef ENABLE_CUDA
#include <multi/adaptors/fft.hpp>
#else
#include <multi/adaptors/fftw.hpp>
#endif

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

#ifdef ENABLE_HEFFTE
#include <heffte.h>
#endif

#include <cassert>

namespace inq {
namespace operations {
namespace space {

void zero_outside_sphere(basis::field<basis::fourier_space, complex> & fphi){
		CALI_CXX_MARK_FUNCTION;
		
	//DATAOPERATIONS GPU::RUN 3D
	gpu::run(fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), point_op = fphi.basis().point_op()] GPU_LAMBDA
					 (auto iz, auto iy, auto ix){
						 if(point_op.outside_sphere(point_op.g2(ix, iy, iz))) fphicub[ix][iy][iz] = complex(0.0);
					 });
}

///////////////////////////////////////////////////////////////
		
void zero_outside_sphere(basis::field<basis::fourier_space, math::vector3<complex>> & fphi){
		CALI_CXX_MARK_FUNCTION;
		
	//DATAOPERATIONS GPU::RUN 3D
	gpu::run(fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), point_op = fphi.basis().point_op()] GPU_LAMBDA
					 (auto iz, auto iy, auto ix){
						 if(point_op.outside_sphere(point_op.g2(ix, iy, iz))) fphicub[ix][iy][iz] = {0.0, 0.0, 0.0};
					 });
}

///////////////////////////////////////////////////////////////

void zero_outside_sphere(basis::field_set<basis::fourier_space, complex>& fphi){
	CALI_CXX_MARK_FUNCTION;
	
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(fphi.set_part().local_size(), fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), point_op = fphi.basis().point_op()] GPU_LAMBDA
					 (auto ist, auto iz, auto iy, auto ix){
						 if(point_op.outside_sphere(point_op.g2(ix, iy, iz))) fphicub[ix][iy][iz][ist] = complex(0.0);
					 });
}

///////////////////////////////////////////////////////////////

void zero_outside_sphere(states::orbital_set<basis::fourier_space, complex>& fphi){
	CALI_CXX_MARK_FUNCTION;
	
	//DATAOPERATIONS GPU::RUN 4D
	gpu::run(fphi.set_part().local_size(), fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.cubic()), point_op = fphi.basis().point_op()] GPU_LAMBDA
					 (auto ist, auto iz, auto iy, auto ix){
						 if(point_op.outside_sphere(point_op.g2(ix, iy, iz))) fphicub[ix][iy][iz][ist] = complex(0.0);
					 });
}

///////////////////////////////////////////////////////////////

template <class InArray4D, class OutArray4D>
void to_fourier(basis::real_space const & real_basis, basis::fourier_space const & fourier_basis, InArray4D const & array_rs, OutArray4D && array_fs) {

	CALI_CXX_MARK_FUNCTION;

#ifdef ENABLE_HEFFTE

	CALI_MARK_BEGIN("heffte_initialization");
 
	heffte::box3d<> const rs_box = {{int(real_basis.cubic_dist(2).start()), int(real_basis.cubic_dist(1).start()), int(real_basis.cubic_dist(0).start())},
																	{int(real_basis.cubic_dist(2).end()) - 1, int(real_basis.cubic_dist(1).end()) - 1, int(real_basis.cubic_dist(0).end()) - 1}};
	
	heffte::box3d<> const fs_box = {{int(fourier_basis.cubic_dist(2).start()), int(fourier_basis.cubic_dist(1).start()), int(fourier_basis.cubic_dist(0).start())},
																	{int(fourier_basis.cubic_dist(2).end()) - 1, int(fourier_basis.cubic_dist(1).end()) - 1, int(fourier_basis.cubic_dist(0).end()) - 1}};

#ifdef ENABLE_CUDA
	heffte::fft3d<heffte::backend::cufft>
#else
	heffte::fft3d<heffte::backend::fftw>
#endif
		fft(rs_box, fs_box, real_basis.comm().get());

	CALI_MARK_END("heffte_initialization");
	
	math::array<complex, 1> input(fft.size_inbox());
	math::prefetch(input);
	math::array<complex, 1> output(fft.size_outbox());	
	math::prefetch(output);
	
	for(int ist = 0; ist < size(array_rs[0][0][0]); ist++){

		input({0, real_basis.local_size()}) = array_rs.flatted().flatted().transposed()[ist];
		
		{
			CALI_CXX_MARK_SCOPE("heffte_forward");
			fft.forward(raw_pointer_cast(input.data_elements()), static_cast<complex *>(output.data_elements()));
		}

		array_fs.flatted().flatted().transposed()[ist] = output({0, fourier_basis.local_size()});
	}
		
#else
	
	namespace multi = boost::multi;
#ifdef ENABLE_CUDA
	namespace fft = multi::fft;
#else
	namespace fft = multi::fftw;
#endif
	
	if(not real_basis.part().parallel()) {
		CALI_CXX_MARK_SCOPE("fft_forward_3d");

		fft::dft({true, true, true, false}, array_rs, array_fs, fft::forward);
		gpu::sync();		

	} else {

		auto & comm = real_basis.comm();
		
		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);
 		auto last_dim = std::get<3>(sizes(array_rs));

		math::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*comm.size(), last_dim});
		
		math::prefetch(tmp);
		
		{
			CALI_CXX_MARK_SCOPE("fft_forward_2d");

			auto const real_x = real_basis.local_sizes();
			fft::dft({false, true, true, false}, array_rs, tmp({0, real_x[0]}, {0, real_x[1]}, {0, real_x[2]}), fft::forward);
			gpu::sync();
		}
		
		CALI_MARK_BEGIN("fft_forward_transpose");		
		math::array<complex, 5> buffer({comm.size(), xblock, real_basis.local_sizes()[1], zblock, last_dim});

		for(int i4 = 0; i4 < comm.size(); i4++){
			gpu::run(last_dim, zblock, real_basis.local_sizes()[1], xblock, 
							 [i4,
								buf = begin(buffer),
								rot = begin(tmp.unrotated(2).partitioned(comm.size()).transposed().rotated().transposed().rotated())]
							 GPU_LAMBDA (auto i0, auto i1, auto i2, auto i3){
								 buf[i4][i3][i2][i1][i0] = rot[i4][i3][i2][i1][i0];
							 });
		}
		CALI_MARK_END("fft_forward_transpose");

		assert(std::get<4>(sizes(buffer)) == last_dim);
		
		tmp.clear();

		{
			CALI_CXX_MARK_SCOPE("fft_forward_alltoall");
			gpu::alltoall(buffer, comm);
		}

		{
			CALI_CXX_MARK_SCOPE("fft_forward_1d");
			
			auto const fourier_x = fourier_basis.local_sizes();
			fft::dft({true, false, false, false}, buffer.flatted()({0, fourier_x[0]}, {0, fourier_x[1]}, {0, fourier_x[2]}), array_fs, fft::forward);
			gpu::sync();
		}

	}
	
#endif
	
}
		
///////////////////////////////////////////////////////////////

template <class InArray4D, class OutArray4D>
void to_real(basis::fourier_space const & fourier_basis, basis::real_space const & real_basis, InArray4D const & array_fs, OutArray4D && array_rs, bool normalize) {

	CALI_CXX_MARK_FUNCTION;

#ifdef ENABLE_HEFFTE
	
	CALI_MARK_BEGIN("heffte_initialization");
	
	heffte::box3d<> const rs_box = {{int(real_basis.cubic_dist(2).start()), int(real_basis.cubic_dist(1).start()), int(real_basis.cubic_dist(0).start())},
																	{int(real_basis.cubic_dist(2).end()) - 1, int(real_basis.cubic_dist(1).end()) - 1, int(real_basis.cubic_dist(0).end()) - 1}};
	
	heffte::box3d<> const fs_box = {{int(fourier_basis.cubic_dist(2).start()), int(fourier_basis.cubic_dist(1).start()), int(fourier_basis.cubic_dist(0).start())},
																	{int(fourier_basis.cubic_dist(2).end()) - 1, int(fourier_basis.cubic_dist(1).end()) - 1, int(fourier_basis.cubic_dist(0).end()) - 1}};

#ifdef ENABLE_CUDA
	heffte::fft3d<heffte::backend::cufft>
#else
	heffte::fft3d<heffte::backend::fftw>
#endif
		fft(rs_box, fs_box, real_basis.comm().get());

	CALI_MARK_END("heffte_initialization");

	auto scaling = heffte::scale::none;
	if(normalize) scaling = heffte::scale::full;
	
	math::array<complex, 1> input(fft.size_inbox());
	math::prefetch(input);	
	math::array<complex, 1> output(fft.size_outbox());	
	math::prefetch(output);
	
	for(int ist = 0; ist < size(array_rs[0][0][0]); ist++){

		input({0, fourier_basis.local_size()}) = array_fs.flatted().flatted().transposed()[ist];
		
		{
			CALI_CXX_MARK_SCOPE("heffte_backward");
			fft.backward(raw_pointer_cast(input.data_elements()), static_cast<complex *>(output.data_elements()), scaling);
		}

		array_rs.flatted().flatted().transposed()[ist] = output({0, real_basis.local_size()});
		
	}
		
#else //Heftte_FOUND
	
	namespace multi = boost::multi;
#ifdef ENABLE_CUDA
	namespace fft = multi::fft;
#else
	namespace fft = multi::fftw;
#endif
	
	if(not real_basis.part().parallel()) {
		CALI_CXX_MARK_SCOPE("fft_backward_3d");
		
		fft::dft({true, true, true, false}, array_fs, array_rs, fft::backward);
		gpu::sync();
		
	} else {

		auto & comm = fourier_basis.comm();

		int xblock = real_basis.cubic_dist(0).block_size();
		int zblock = fourier_basis.cubic_dist(2).block_size();
		auto last_dim = std::get<3>(sizes(array_fs));
		
		math::array<complex, 5> buffer({comm.size(), xblock, real_basis.local_sizes()[1], zblock, last_dim});
		math::prefetch(buffer);
		
		{
			CALI_CXX_MARK_SCOPE("fft_backward_2d");
			
			fft::dft({true, true, false, false}, array_fs, buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fft::backward);
			gpu::sync();
		}

		{
			CALI_CXX_MARK_SCOPE("fft_backward_alltoall");
			gpu::alltoall(buffer, comm);
		}
		
		math::array<complex, 4> tmp({real_basis.local_sizes()[0], real_basis.local_sizes()[1], zblock*comm.size(), last_dim});
		math::prefetch(tmp);
		
		{
			CALI_CXX_MARK_SCOPE("fft_backward_transpose");
			for(int i4 = 0; i4 < comm.size(); i4++){
				gpu::run(last_dim, zblock, real_basis.local_sizes()[1], real_basis.local_sizes()[0],
								 [i4, 
									buf = begin(buffer),
									rot = begin(tmp.unrotated(2).partitioned(comm.size()).transposed().rotated().transposed().rotated())]
								 GPU_LAMBDA (auto i0, auto i1, auto i2, auto i3){
									 rot[i4][i3][i2][i1][i0] = buf[i4][i3][i2][i1][i0];
								 });
			}
		}

		{
			CALI_CXX_MARK_SCOPE("fft_backward_1d");

			auto tmpsub = tmp({0, real_basis.local_sizes()[0]}, {0, real_basis.local_sizes()[1]}, {0, real_basis.local_sizes()[2]});

#ifdef ENABLE_CUDA
			//when using cuda, do a loop explicitly over the last dimension (n1), otherwise multi makes a loop over the 2 first ones (n0 and n1). And we know that n3 << n0*n1.
			for(auto ii : extension(tmpsub.unrotated())) {
				fft::dft({false, false, true}, tmpsub.unrotated()[ii], array_rs.unrotated()[ii], fft::backward);
			}
#else
			fft::dft({false, false, true, false}, tmpsub, array_rs, fft::backward);
#endif
			
			gpu::sync();

		}
	}

	if(normalize){
		CALI_CXX_MARK_SCOPE("fft_normalize");
		gpu::run(size(array_rs[0][0][0]), real_basis.local_size(), 
						 [ar = begin(array_rs.flatted().flatted()), factor = 1.0/real_basis.size()] GPU_LAMBDA (auto ist, auto ip){
							 ar[ip][ist] = factor*ar[ip][ist];
						 });
	}

#endif //Heftte_FOUND

}

///////////////////////////////////////////////////////////////
		
basis::field_set<basis::fourier_space, complex> to_fourier(const basis::field_set<basis::real_space, complex> & phi){

	CALI_CXX_MARK_SCOPE("to_fourier(field_set)");
		
	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis);
	
	basis::field_set<basis::fourier_space, complex> fphi(fourier_basis, phi.set_size(), phi.full_comm());

	to_fourier(real_basis, fourier_basis, phi.cubic(), fphi.cubic());
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
	
	return fphi;
}

///////////////////////////////////////////////////////////////
		
states::orbital_set<basis::fourier_space, complex> to_fourier(const states::orbital_set<basis::real_space, complex> & phi){

	CALI_CXX_MARK_SCOPE("to_fourier(field_set)");
		
	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis);
	
	states::orbital_set<basis::fourier_space, complex> fphi(fourier_basis, phi.set_size(), phi.full_comm());

	to_fourier(real_basis, fourier_basis, phi.cubic(), fphi.cubic());
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
	
	return fphi;
}

///////////////////////////////////////////////////////////////

basis::field_set<basis::real_space, complex> to_real(const basis::field_set<basis::fourier_space, complex> & fphi, bool const normalize = true){

	CALI_CXX_MARK_SCOPE("to_real(field_set)");
	
	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis);
	
	basis::field_set<basis::real_space, complex> phi(real_basis, fphi.set_size(), fphi.full_comm());

	to_real(fourier_basis, real_basis, fphi.cubic(), phi.cubic(), normalize);

	return phi;
}


///////////////////////////////////////////////////////////////

states::orbital_set<basis::real_space, complex> to_real(const states::orbital_set<basis::fourier_space, complex> & fphi, bool const normalize = true){

	CALI_CXX_MARK_SCOPE("to_real(orbital_set)");
	
	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis);
	
	states::orbital_set<basis::real_space, complex> phi(real_basis, fphi.set_size(), fphi.full_comm());

	to_real(fourier_basis, real_basis, fphi.cubic(), phi.cubic(), normalize);

	return phi;
}


///////////////////////////////////////////////////////////////

basis::field<basis::fourier_space, complex> to_fourier(const basis::field<basis::real_space, complex> & phi){

	CALI_CXX_MARK_SCOPE("to_fourier(field)");
	
	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis);
	
	basis::field<basis::fourier_space, complex> fphi(fourier_basis);

	to_fourier(real_basis, fourier_basis, phi.hypercubic(), fphi.hypercubic());
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
			
	return fphi;
	
}

///////////////////////////////////////////////////////////////			
	
basis::field<basis::real_space, complex> to_real(const basis::field<basis::fourier_space, complex> & fphi, bool normalize = true){

	CALI_CXX_MARK_SCOPE("to_real(field)");
	
	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis);

	basis::field<basis::real_space, complex> phi(real_basis);

	to_real(fourier_basis, real_basis, fphi.hypercubic(), phi.hypercubic(), normalize);
			
	return phi;
}

///////////////////////////////////////////////////////////////

auto to_fourier(const basis::field<basis::real_space, math::vector3<complex>> & phi){

	CALI_CXX_MARK_SCOPE("to_fourier(vector_field)");
	
	auto & real_basis = phi.basis();
	basis::fourier_space fourier_basis(real_basis);
	
	basis::field<basis::fourier_space, math::vector3<complex>> fphi(fourier_basis);
	
	to_fourier(
		real_basis, 
		fourier_basis, 
		phi .cubic().reinterpret_array_cast<complex const>(3), 
		fphi.cubic().reinterpret_array_cast<complex      >(3)
	);

	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
			
	return fphi;
	
}

///////////////////////////////////////////////////////////////

basis::field<basis::real_space, math::vector3<complex>> 
to_real(const basis::field<basis::fourier_space, math::vector3<complex>> & fphi, bool normalize = true){

	CALI_CXX_MARK_SCOPE("to_real(vector_field)");
	
	auto & fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis);

	basis::field<basis::real_space, math::vector3<complex>> phi(real_basis);

	to_real(
		fourier_basis, real_basis, 
		fphi.cubic().reinterpret_array_cast<complex const>(3), 
		phi .cubic().reinterpret_array_cast<complex      >(3), 
		normalize
	);

	return phi;
}

///////////////////////////////////////////////////////////////

basis::field_set<basis::real_space, math::vector3<complex>> 
to_real(basis::field_set<basis::fourier_space, math::vector3<complex>> const& fphi, bool normalize = true){

	CALI_CXX_MARK_SCOPE("to_real(vector_field_set)");

	auto const& fourier_basis = fphi.basis();
	basis::real_space real_basis(fourier_basis);

	basis::field_set<basis::real_space, math::vector3<complex>> phi(real_basis, fphi.set_size(), fphi.full_comm());

	auto const& fphi_as_scalar = fphi.cubic().reinterpret_array_cast<complex const>(3).rotated(3).flatted().rotated();
	auto &&     phi_as_scalar  = phi .cubic().reinterpret_array_cast<complex      >(3).rotated(3).flatted().rotated();

	to_real(fourier_basis, real_basis, fphi_as_scalar, phi_as_scalar, normalize);

	return phi;
}

///////////////////////////////////////////////////////////////

}
}
}

#ifdef INQ_OPERATIONS_SPACE_UNIT_TEST
#undef INQ_OPERATIONS_SPACE_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("function operations::space", "[operations::space]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using math::vector3;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto basis_comm = cart_comm.axis(1);
	
	auto ecut = 23.0_Ha;
	double ll = 6.66;
	
	ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
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
					double r2 = rs.point_op().r2(ix, iy, iz);
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
					double g2 = fphi.basis().point_op().g2(ix, iy, iz);
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

