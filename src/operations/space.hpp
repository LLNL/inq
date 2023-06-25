/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OPERATIONS__SPACE
#define INQ__OPERATIONS__SPACE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h> //for ENABLE_HEFFTE

#include <gpu/run.hpp>
#include <parallel/alltoall.hpp>
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

///////////////////////////////////////////////////////////////

template <typename FieldSetType>
void zero_outside_sphere(FieldSetType & fphi){
	CALI_CXX_MARK_FUNCTION;
	
	gpu::run(fphi.local_set_size(), fphi.basis().local_sizes()[2], fphi.basis().local_sizes()[1], fphi.basis().local_sizes()[0],
					 [fphicub = begin(fphi.hypercubic()), point_op = fphi.basis().point_op()] GPU_LAMBDA
					 (auto ist, auto iz, auto iy, auto ix){
						 if(point_op.outside_sphere(ix, iy, iz)) fphicub[ix][iy][iz][ist] = zero<typename FieldSetType::element_type>();
					 });
}

///////////////////////////////////////////////////////////////

#ifdef ENABLE_HEFFTE
template <class InArray4D, class OutArray4D>
void to_fourier_array(basis::real_space const & real_basis, basis::fourier_space const & fourier_basis, InArray4D const & array_rs, OutArray4D && array_fs) {

	CALI_CXX_MARK_FUNCTION;

	assert(std::get<3>(sizes(array_rs)) == std::get<3>(sizes(array_fs)));
	
	CALI_MARK_BEGIN("heffte_initialization");
 
	heffte::box3d<> const rs_box = {{int(real_basis.cubic_part(2).start()), int(real_basis.cubic_part(1).start()), int(real_basis.cubic_part(0).start())},
																	{int(real_basis.cubic_part(2).end()) - 1, int(real_basis.cubic_part(1).end()) - 1, int(real_basis.cubic_part(0).end()) - 1}};
	
	heffte::box3d<> const fs_box = {{int(fourier_basis.cubic_part(2).start()), int(fourier_basis.cubic_part(1).start()), int(fourier_basis.cubic_part(0).start())},
																	{int(fourier_basis.cubic_part(2).end()) - 1, int(fourier_basis.cubic_part(1).end()) - 1, int(fourier_basis.cubic_part(0).end()) - 1}};

#ifdef ENABLE_CUDA
	heffte::fft3d<heffte::backend::cufft>
#else
	heffte::fft3d<heffte::backend::fftw>
#endif
		fft(rs_box, fs_box, real_basis.comm().get());

	CALI_MARK_END("heffte_initialization");

	// we don't need a copy when there is just one field
	if(size(array_rs[0][0][0]) == 1) {
		CALI_CXX_MARK_SCOPE("heffte_forward_1");
		fft.forward((const std::complex<double> *) raw_pointer_cast(array_rs.base()), (std::complex<double> *) raw_pointer_cast(array_fs.base()));
		return;
	}
	
	gpu::array<complex, 1> input(fft.size_inbox());
	gpu::prefetch(input);
	gpu::array<complex, 1> output(fft.size_outbox());	
	gpu::prefetch(output);

	for(int ist = 0; ist < size(array_rs[0][0][0]); ist++){

		{
			CALI_CXX_MARK_SCOPE("heffte_forward_copy_1");
			input({0, real_basis.local_size()}) = array_rs.flatted().flatted().transposed()[ist];
		}
		
		{
			CALI_CXX_MARK_SCOPE("heffte_forward");
			fft.forward((const std::complex<double> *) raw_pointer_cast(input.data_elements()), (std::complex<double> *) raw_pointer_cast(output.data_elements()));
		}
		{
			CALI_CXX_MARK_SCOPE("heffte_forward_copy_2");
			array_fs.flatted().flatted().transposed()[ist] = output({0, fourier_basis.local_size()});
		}
	}

}

#else // no HEFFTE

///////////////////////////////////////////////////////////////

template <class InArray4D, class OutArray4D>
void to_fourier_array(basis::real_space const & real_basis, basis::fourier_space const & fourier_basis, InArray4D const & array_rs, OutArray4D && array_fs) {

	CALI_CXX_MARK_FUNCTION;

	assert(std::get<3>(sizes(array_rs)) == std::get<3>(sizes(array_fs)));
	
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
		
		int xblock = real_basis.cubic_part(0).max_local_size();
		int zblock = fourier_basis.cubic_part(2).max_local_size();
		assert(real_basis.local_sizes()[1] == fourier_basis.local_sizes()[1]);
 		auto last_dim = std::get<3>(sizes(array_rs));

		gpu::array<complex, 4> tmp({xblock, real_basis.local_sizes()[1], zblock*comm.size(), last_dim});
		
		gpu::prefetch(tmp);
		
		{
			CALI_CXX_MARK_SCOPE("fft_forward_2d");

			auto const real_x = real_basis.local_sizes();
			fft::dft({false, true, true, false}, array_rs, tmp({0, real_x[0]}, {0, real_x[1]}, {0, real_x[2]}), fft::forward);
			gpu::sync();
		}
		
		CALI_MARK_BEGIN("fft_forward_transpose");		
		gpu::array<complex, 5> buffer({comm.size(), xblock, real_basis.local_sizes()[1], zblock, last_dim});

		for(int i4 = 0; i4 < comm.size(); i4++){
			gpu::run(last_dim, zblock, real_basis.local_sizes()[1], xblock, 
							 [i4,
								buf = begin(buffer),
								rot = begin(tmp.unrotated().unrotated().partitioned(comm.size()).transposed().rotated().transposed().rotated())]
							 GPU_LAMBDA (auto i0, auto i1, auto i2, auto i3){
								 buf[i4][i3][i2][i1][i0] = rot[i4][i3][i2][i1][i0];
							 });
		}
		CALI_MARK_END("fft_forward_transpose");

		assert(std::get<4>(sizes(buffer)) == last_dim);
		
		tmp.clear();

		{
			CALI_CXX_MARK_SCOPE("fft_forward_alltoall");
			parallel::alltoall(buffer, comm);
		}

		{
			CALI_CXX_MARK_SCOPE("fft_forward_1d");
			
			auto const fourier_x = fourier_basis.local_sizes();
			fft::dft({true, false, false, false}, buffer.flatted()({0, fourier_x[0]}, {0, fourier_x[1]}, {0, fourier_x[2]}), array_fs, fft::forward);
			gpu::sync();
		}

	}
}

#endif

///////////////////////////////////////////////////////////////

#ifdef ENABLE_HEFFTE
template <class InArray4D, class OutArray4D>
void to_real_array(basis::fourier_space const & fourier_basis, basis::real_space const & real_basis, InArray4D const & array_fs, OutArray4D && array_rs, bool normalize) {

	CALI_CXX_MARK_FUNCTION;

	CALI_MARK_BEGIN("heffte_initialization");
	
	heffte::box3d<> const rs_box = {{int(real_basis.cubic_part(2).start()), int(real_basis.cubic_part(1).start()), int(real_basis.cubic_part(0).start())},
																	{int(real_basis.cubic_part(2).end()) - 1, int(real_basis.cubic_part(1).end()) - 1, int(real_basis.cubic_part(0).end()) - 1}};
	
	heffte::box3d<> const fs_box = {{int(fourier_basis.cubic_part(2).start()), int(fourier_basis.cubic_part(1).start()), int(fourier_basis.cubic_part(0).start())},
																	{int(fourier_basis.cubic_part(2).end()) - 1, int(fourier_basis.cubic_part(1).end()) - 1, int(fourier_basis.cubic_part(0).end()) - 1}};

#ifdef ENABLE_CUDA
	heffte::fft3d<heffte::backend::cufft>
#else
	heffte::fft3d<heffte::backend::fftw>
#endif
		fft(rs_box, fs_box, real_basis.comm().get());

	CALI_MARK_END("heffte_initialization");

	auto scaling = heffte::scale::none;
	if(normalize) scaling = heffte::scale::full;

	// we don't need a copy when there is just one field
	if(size(array_rs[0][0][0]) == 1) {
		CALI_CXX_MARK_SCOPE("heffte_backward_1");
		fft.backward((const std::complex<double> *) raw_pointer_cast(array_fs.base()), (std::complex<double> *) raw_pointer_cast(array_rs.base()), scaling);
		return;
	}
	
	gpu::array<complex, 1> input(fft.size_inbox());
	gpu::prefetch(input);	
	gpu::array<complex, 1> output(fft.size_outbox());	
	gpu::prefetch(output);
	
	for(int ist = 0; ist < size(array_rs[0][0][0]); ist++){

		input({0, fourier_basis.local_size()}) = array_fs.flatted().flatted().transposed()[ist];
		
		{
			CALI_CXX_MARK_SCOPE("heffte_backward");
			fft.backward((const std::complex<double> *) raw_pointer_cast(input.data_elements()), (std::complex<double> *) raw_pointer_cast(output.data_elements()), scaling);
		}

		array_rs.flatted().flatted().transposed()[ist] = output({0, real_basis.local_size()});
		
	}
}
	
#else //NO HEFTTE

///////////////////////////////////////////////////////////////

template <class InArray4D, class OutArray4D>
void to_real_array(basis::fourier_space const & fourier_basis, basis::real_space const & real_basis, InArray4D const & array_fs, OutArray4D && array_rs, bool normalize) {

	CALI_CXX_MARK_FUNCTION;
	
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

		int xblock = real_basis.cubic_part(0).max_local_size();
		int zblock = fourier_basis.cubic_part(2).max_local_size();
		auto last_dim = std::get<3>(sizes(array_fs));
		
		gpu::array<complex, 5> buffer({comm.size(), xblock, real_basis.local_sizes()[1], zblock, last_dim});
		gpu::prefetch(buffer);
		
		{
			CALI_CXX_MARK_SCOPE("fft_backward_2d");
			
			fft::dft({true, true, false, false}, array_fs, buffer.flatted()({0, fourier_basis.local_sizes()[0]}, {0, fourier_basis.local_sizes()[1]}, {0, fourier_basis.local_sizes()[2]}), fft::backward);
			gpu::sync();
		}

		{
			CALI_CXX_MARK_SCOPE("fft_backward_alltoall");
			parallel::alltoall(buffer, comm);
		}
		
		gpu::array<complex, 4> tmp({real_basis.local_sizes()[0], real_basis.local_sizes()[1], zblock*comm.size(), last_dim});
		gpu::prefetch(tmp);
		
		{
			CALI_CXX_MARK_SCOPE("fft_backward_transpose");
			for(int i4 = 0; i4 < comm.size(); i4++){
				gpu::run(last_dim, zblock, real_basis.local_sizes()[1], real_basis.local_sizes()[0],
								 [i4, 
									buf = begin(buffer),
									rot = begin(tmp.unrotated().unrotated().partitioned(comm.size()).transposed().rotated().transposed().rotated())]
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
		gpu::run(size(array_rs[0][0][0])*real_basis.local_size(), 
						 [ar = begin(array_rs.flatted().flatted().flatted()), factor = 1.0/real_basis.size()] GPU_LAMBDA (auto ip){
							 ar[ip] = factor*ar[ip];
						 });
	}
}
#endif

///////////////////////////////////////////////////////////////

template <class FieldSetType>
auto to_fourier(const FieldSetType & phi){

	CALI_CXX_MARK_SCOPE("to_fourier");

	auto fphi = FieldSetType::reciprocal(phi.skeleton());

	using type = typename FieldSetType::element_type;
	
	if constexpr (not is_vector3<type>::value){
		static_assert(std::is_same<type, complex>::value, "Only implemented for complex");
		
		to_fourier_array(phi.basis(), fphi.basis(), phi.hypercubic(), fphi.hypercubic());

	} else {

		static_assert(std::is_same<typename type::element_type, complex>::value, "Only implemented for complex vector3");
		
		auto &&    fphi_as_scalar = fphi.hypercubic().template reinterpret_array_cast<complex      >(3).rotated().rotated().rotated().flatted().rotated();
		auto const& phi_as_scalar = phi .hypercubic().template reinterpret_array_cast<complex const>(3).rotated().rotated().rotated().flatted().rotated();

		to_fourier_array(phi.basis(), fphi.basis(), phi_as_scalar, fphi_as_scalar);
	}
	
	if(fphi.basis().spherical()) zero_outside_sphere(fphi);
	
	return fphi;
}

///////////////////////////////////////////////////////////////

template <class FieldSetType>
auto to_real(const FieldSetType & fphi, bool const normalize = true){

	CALI_CXX_MARK_SCOPE("to_real");

	auto phi = FieldSetType::reciprocal(fphi.skeleton());
	
	using type = typename FieldSetType::element_type;
	
	if constexpr (not is_vector3<type>::value){

		static_assert(std::is_same<type, complex>::value, "Only valid for complex");
		
		to_real_array(fphi.basis(), phi.basis(), fphi.hypercubic(), phi.hypercubic(), normalize);

	} else {

		static_assert(std::is_same<typename type::element_type, complex>::value, "Only valid for complex vector3");
		
		auto const& fphi_as_scalar = fphi.hypercubic().template reinterpret_array_cast<complex const>(3).rotated().rotated().rotated().flatted().rotated();
		auto &&     phi_as_scalar  = phi .hypercubic().template reinterpret_array_cast<complex      >(3).rotated().rotated().rotated().flatted().rotated();
		
		to_real_array(fphi.basis(), phi.basis(), fphi_as_scalar, phi_as_scalar, normalize);

	}
	
	return phi;
}

///////////////////////////////////////////////////////////////

}
}
}
#endif

#ifdef INQ_OPERATIONS_SPACE_UNIT_TEST
#undef INQ_OPERATIONS_SPACE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	auto basis_comm = basis::basis_subcomm(cart_comm);
	basis::real_space rs(ions::unit_cell::cubic(6.66_b), /*spacing =*/ 0.46320257, basis_comm);
	
	basis::field_set<basis::real_space, complex> phi(rs, 7, cart_comm);
	
	SECTION("Zero"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++) phi.hypercubic()[ix][iy][iz][ist] = 0.0;
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().local_sizes()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().local_sizes()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(fphi.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= fphi.hypercubic().num_elements();

		CHECK(diff < 1e-15);
		
		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++)	diff += fabs(phi.hypercubic()[ix][iy][iz][ist]);
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= phi2.hypercubic().num_elements();

		CHECK(diff < 1e-15);
		
	}
	
	SECTION("Gaussian"){
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					double r2 = rs.point_op().r2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						double sigma = 0.5*(ist + 1);
						phi.hypercubic()[ix][iy][iz][ist] = exp(-sigma*r2);
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
						diff += fabs(fphi.hypercubic()[ix][iy][iz][ist] - pow(M_PI/sigma, 3.0/2.0)*exp(-0.25*g2/sigma));
					}
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});
		
		diff /= fphi.hypercubic().num_elements();

		//not sure what is wrong here
		std::cout << "DIFF1 " << diff << std::endl;

		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					for(int ist = 0; ist < phi.set_part().local_size(); ist++){
						diff += fabs(phi.hypercubic()[ix][iy][iz][ist] - phi2.hypercubic()[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		cart_comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});		

		diff /= phi2.hypercubic().num_elements();
		
		CHECK(diff < 1e-15);
		
	}
	
}
#endif

