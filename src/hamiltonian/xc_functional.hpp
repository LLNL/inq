/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_FUNCTIONAL
#define INQ__HAMILTONIAN__XC_FUNCTIONAL

/*
 Copyright (C) 2019 Xavier Andrade, Alexey Kartsev

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

#include <xc.h>
#include <operations/gradient.hpp>
#include <operations/integral.hpp>
#include <operations/divergence.hpp>
#include <operations/laplacian.hpp>
#include <basis/field.hpp>

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {
	class xc_functional {

		public:
			
		xc_functional(const int functional_id){
			true_functional_ = functional_id > 0;
			if(true_functional_ and xc_func_init(&func_, functional_id, XC_UNPOLARIZED) != 0){
				fprintf(stderr, "Functional '%d' not found\n", functional_id);
				exit(1);
			}
		}

		~xc_functional(){
			if(true_functional_) xc_func_end(&func_);
		}

		auto true_functional() const {
			return true_functional_;
		}
		
		template <class field_type>
		void operator()(field_type const & density, double & xc_energy, field_type & vxc) const {

			assert(true_functional_);
			
			CALI_CXX_MARK_SCOPE("xc_functional");
			
			field_type exc(vxc.skeleton());
			unpolarized(density.linear().size(), density, exc, vxc);

			xc_energy = operations::integral_product(density, exc);

		}

		auto & libxc_func() const {
			return func_;
		}
		
		auto libxc_func_ptr() const {
			return &func_;
		}

		auto exx_coefficient() const {
			if(xc_hyb_type(&func_) == XC_HYB_HYBRID) return xc_hyb_exx_coef(&func_);
			return 0.0;
		}
		
		//this function has to be public because of cuda limitations
		template <class density_type, class exc_type, class vxc_type>
		void ggafunctional(long size, density_type const & density, exc_type & exc, vxc_type & vxc) const { 

			CALI_CXX_MARK_FUNCTION;

			auto grad_real = operations::gradient(density);
			basis::field<basis::real_space, double> sigma(vxc.basis());

			gpu::run(vxc.basis().local_size(),
							 [sig = begin(sigma.linear()), grad = begin(grad_real.linear()), metric = density.basis().cell().metric()] GPU_LAMBDA (auto ip){
								 sig[ip] = metric.norm(grad[ip]);
							 });

			basis::field<basis::real_space, double> vsigma(vxc.basis());
			
			xc_gga_exc_vxc(&func_, size, density.data(), sigma.data(), exc.data(), vxc.data(), vsigma.data());
			gpu::sync();
					
			basis::field<basis::real_space, math::vector3<double, math::covariant>> vxc_extra(vxc.basis());

			gpu::run(vxc.basis().local_size(),
							 [vex = begin(vxc_extra.linear()), vsig = begin(vsigma.linear()), grad = begin(grad_real.linear())] GPU_LAMBDA (auto ip){
								 vex[ip] = vsig[ip]*grad[ip];
							 });

			auto divvxcextra = operations::divergence(vxc_extra);

			gpu::run(vxc.basis().local_size(),
							 [vv = begin(vxc.linear()), div = begin(divvxcextra.linear())] GPU_LAMBDA (auto ip){
								 vv[ip] -= 2.0*div[ip];
							 });
		}

	private:
		
		template <class density_type, class exc_type, class vxc_type>
		void unpolarized(long size, density_type const & density, exc_type & exc, vxc_type & vxc) const{

			CALI_CXX_MARK_FUNCTION;
			
			switch(func_.info->family) {
				case XC_FAMILY_LDA:{
					xc_lda_exc_vxc(&func_, size, density.data(), exc.data(), vxc.data());
					gpu::sync();
					break;
				}
				case XC_FAMILY_GGA:{
					ggafunctional(size, density, exc, vxc);
					break;
				}	
				default:{
					assert(false);
					break;
				}
			}
		}

	private:

		bool true_functional_;
		xc_func_type func_;
			
	};

}
}


///////////////////////////////////////////////////////////////////

#ifdef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <operations/randomize.hpp>
#include <operations/gradient.hpp>
#include <basis/field.hpp>
#include <math/array.hpp>
#include <utils/finite_difference.hpp>
#include <ions/geometry.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>

GPU_FUNCTION auto sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	return  sin(n*dot(kvec, rr))*sin(n*dot(kvec, rr));
}

auto laplacian_sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	auto ff = 0.0;
	for(int idir = 0; idir < 3 ; idir++) ff += 2*n*n*kvec[idir]*kvec[idir]*cos(2*n*dot(kvec, rr));
	return ff;
}

auto gradient_sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	inq::math::vector3<double> ff;
	for(int idir = 0; idir < 3 ; idir++) ff[idir] = 2*n*kvec[idir]*cos(n*dot(kvec, rr))*sin(n*dot(kvec, rr));
	return ff;
}

auto gaussiansigma(inq::math::vector3<double> rr, double sigma){
	return  pow(2*M_PI, -1.5)*pow(sigma, -3.0)*exp(-0.5*norm(rr)*pow(sigma, -2.0));
}

auto gaussian(inq::math::vector3<double> rr){  // sigma = 1/sqrt(2)
	return  pow(M_PI, -1.5)*exp(-norm(rr));
}

auto dgaussian(inq::math::vector3<double> rr){
	inq::math::vector3<double> ff;
	for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian(rr);
	return ff;
}

auto laplacian_gaussian(inq::math::vector3<double> rr){
	return 4.0*dot(rr, rr)*gaussian(rr) - 6.0*gaussian(rr);
}

TEST_CASE("function hamiltonian::xc_functional", "[hamiltonian::xc_functional]") {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	using namespace operations;
	using namespace inq::utils;
	using math::vector3;

	double lx = 9;
	double ly = 12;
	double lz = 10;

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	SECTION("LDA"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		basis::real_space rs(box, cart_comm);

		basis::field<basis::real_space, double> gaussian_field(rs);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					gaussian_field.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}

		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussianVxc(rs);
		double gaussianExc;

		CHECK(ldafunctional.exx_coefficient() == 0.0);
		
		ldafunctional(gaussian_field, gaussianExc, gaussianVxc);
		CHECK(gaussianExc == -0.270646_a);

		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					math::array<double, 1> local_density{gaussian(vec)};
					math::array<double, 1> local_exc{1};
					math::array<double, 1> local_vxc{1};
					xc_lda_exc_vxc(ldafunctional.libxc_func_ptr(), 1, raw_pointer_cast(local_density.data_elements()), raw_pointer_cast(local_exc.data_elements()), raw_pointer_cast(local_vxc.data_elements()));
					gpu::sync();
					
					CHECK(Approx(local_vxc[0]) == gaussianVxc.cubic()[ix][iy][iz]);

					int_xc_energy += local_exc[0]*local_density[0]*rs.volume_element();
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&int_xc_energy, 1, std::plus<>{});
		
		CHECK(Approx(gaussianExc) == int_xc_energy);

		if(rs.part().contains(1)) CHECK(gaussianVxc.linear()[rs.part().global_to_local(utils::global_index(1))] == -0.5111609291_a);
		if(rs.part().contains(8233)) CHECK(gaussianVxc.linear()[rs.part().global_to_local(utils::global_index(8233))] == -0.00406881_a);
		if(rs.part().contains(233)) CHECK(fabs(gaussianVxc.linear()[rs.part().global_to_local(utils::global_index(233))]) < 1e-10);
		if(rs.part().contains(rs.size() - 1)) CHECK(gaussianVxc.linear()[rs.part().global_to_local(utils::global_index(rs.size() - 1))] == -0.4326883849_a);

	}

	SECTION("GGA"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(90.0_Ha);
		basis::real_space rs(box, cart_comm);
		
		basis::field<basis::real_space, double> field(rs);
		gpu::run(rs.local_sizes()[2], rs.local_sizes()[1], rs.local_sizes()[0],
						 [pop = rs.point_op(), fie = begin(field.cubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 auto vec = pop.rvector_cartesian(ix, iy, iz);
							 fie[ix][iy][iz] = sqwave(vec, 3);
						 });
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		basis::field<basis::real_space, double> Vxc(rs);
		basis::field<basis::real_space, double> local_vsigma_output(rs);

		CHECK(ggafunctional.exx_coefficient() == 0.0);
		
		double Exc = 0.0;
		ggafunctional(field, Exc, Vxc);

		CHECK(Exc == -393.4604748792_a);
		
		if(rs.part().contains(1)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(1))] == -0.5607887985_a);
		if(rs.part().contains(33)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(33))] == -1.1329131862_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(rs.size() - 1))] == -1.1461742979_a);
	}

	SECTION("Uniform"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		basis::real_space rs(box, cart_comm);
		
		basis::field<basis::real_space, double> gaussian_field(rs);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					gaussian_field.cubic()[ix][iy][iz] = 0.393;
				}
			}
		}
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussianVxcLDA(rs) , gaussianVxcGGA(rs);
		double gaussianExcLDA, gaussianExcGGA;

		ggafunctional(gaussian_field, gaussianExcLDA, gaussianVxcLDA);
		ldafunctional(gaussian_field, gaussianExcGGA, gaussianVxcGGA);
		CHECK(gaussianExcLDA == gaussianExcGGA);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					CHECK(Approx(gaussianVxcLDA.cubic()[ix][iy][iz]) == gaussianVxcGGA.cubic()[ix][iy][iz]);
				}
			}
		}

	}

	SECTION("HYBRIDS"){
		inq::hamiltonian::xc_functional b3lyp(XC_HYB_GGA_XC_B3LYP);
		CHECK(b3lyp.exx_coefficient() == 0.2_a);

		inq::hamiltonian::xc_functional pbeh(XC_HYB_GGA_XC_PBEH);
		CHECK(pbeh.exx_coefficient() == 0.25_a);
	}
	
}

#endif
#endif
