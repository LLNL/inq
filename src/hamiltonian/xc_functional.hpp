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

namespace inq {
namespace hamiltonian {
	class xc_functional {

		public:
			
		xc_functional(const int functional_id){
			if(xc_func_init(&func_, functional_id, XC_UNPOLARIZED) != 0){
				fprintf(stderr, "Functional '%d' not found\n", functional_id);
				exit(1);
			}
		}

		~xc_functional(){
			xc_func_end(&func_);
		}

		template <class field_type>
		void operator()(field_type const & density, double & xc_energy, field_type & vxc) const {

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
		
		//this function has to be public because of cuda limitations
		template <class density_type, class exc_type, class vxc_type>
		void ggafunctional(long size, density_type const & density, exc_type & exc, vxc_type & vxc) const { 

			CALI_CXX_MARK_FUNCTION;

			auto grad_real = operations::gradient(density);
			basis::field<basis::real_space, double> sigma(vxc.basis());

			gpu::run(vxc.basis().local_size(),
							 [sig = begin(sigma.linear()), grad = begin(grad_real.linear())] GPU_LAMBDA (auto ip){
								 sig[ip] = norm(grad[ip]);
							 });

			basis::field<basis::real_space, double> vsigma(vxc.basis());
			
			xc_gga_exc_vxc(&func_, size, density.data(), sigma.data(), exc.data(), vxc.data(), vsigma.data());
			gpu::sync();
					
			basis::field<basis::real_space, math::vector3<double>> vxc_extra(vxc.basis());

			gpu::run(vxc.basis().local_size(),
							 [vex = begin(vxc_extra.linear()), vsig = begin(vsigma.linear()), grad = begin(grad_real.linear())] GPU_LAMBDA (auto ip){
								 vex[ip][0] = vsig[ip]*grad[ip][0];
								 vex[ip][1] = vsig[ip]*grad[ip][1];
								 vex[ip][2] = vsig[ip]*grad[ip][2];
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

			xc_func_type func_;
			
	};

}
}


///////////////////////////////////////////////////////////////////

#ifdef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST

#include <catch2/catch.hpp>
#include <operations/randomize.hpp>
#include <operations/gradient.hpp>
#include <basis/field.hpp>
#include <math/array.hpp>
#include <utils/finite_difference.hpp>
#include <ions/geometry.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>

	auto sqwave(inq::math::vector3<double> rr, int n){
		inq::math::vector3<double> kvec = 2.0 * M_PI * inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
		return  sin(n*dot(kvec, rr))*sin(n*dot(kvec, rr));
	}

	auto laplacian_sqwave(inq::math::vector3<double> rr, int n){
		inq::math::vector3<double> kvec = 2.0 * M_PI * inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
		auto ff = 0.0;
		for(int idir = 0; idir < 3 ; idir++) ff += 2*n*n*kvec[idir]*kvec[idir]*cos(2*n*dot(kvec, rr));
		return ff;
	}

	auto gradient_sqwave(inq::math::vector3<double> rr, int n){
		inq::math::vector3<double> kvec = 2.0 * M_PI * inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
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
	using namespace Catch::literals;
	using namespace operations;
	using namespace inq::utils;
	using math::vector3;

	double lx = 9;
	double ly = 12;
	double lz = 10;

	ions::UnitCell cell(vector3<double>(lx, 0.0, 0.0), vector3<double>(0.0, ly, 0.0), vector3<double>(0.0, 0.0, lz));

	boost::mpi3::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});
	
	SECTION("LDA"){
		
		basis::real_space rs(cell, input::basis::cutoff_energy(20.0), cart_comm);

		basis::field<basis::real_space, double> gaussian_field(rs);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					gaussian_field.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}

		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussianVxc(rs);
		double gaussianExc;

		ldafunctional(gaussian_field, gaussianExc, gaussianVxc);
		CHECK(gaussianExc == -0.270646_a);

		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					math::array<double, 1> local_density{{gaussian(vec)}};
					math::array<double, 1> local_exc{1};
					math::array<double, 1> local_vxc{1};
					xc_lda_exc_vxc(ldafunctional.libxc_func_ptr(), 1, static_cast<double *>(local_density.data_elements()), static_cast<double *>(local_exc.data_elements()), static_cast<double *>(local_vxc.data_elements()));
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
		
		basis::real_space rs(cell, input::basis::cutoff_energy(90.0), cart_comm);
		basis::field<basis::real_space, double> field(rs);
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					field.cubic()[ix][iy][iz] = sqwave(vec, 3);
				}
			}
		}
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		basis::field<basis::real_space, double> Vxc(rs);
		basis::field<basis::real_space, double> local_vsigma_output(rs);

		double Exc = 0.0;
		ggafunctional(field, Exc, Vxc);

		CHECK(Exc == -393.4604748792_a);
		
		std::ofstream fout("fout.dat");

		double diff_ways = 0.0;
		double diff = 0.0;
		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){

					auto vec = rs.rvector(ix, iy, iz);
					math::array<double, 1> local_exc{1};
					math::array<double, 1> local_vxc{1};
					math::array<double, 1> local_vsigma{1};
					math::array<double, 1> local_density{{sqwave(vec, 3)}};
					math::array<double, 1> local_sigma{{dot(gradient_sqwave(vec, 3), gradient_sqwave(vec, 3))}};
					
					xc_gga_exc_vxc(ggafunctional.libxc_func_ptr(), 1, static_cast<double *>(local_density.data_elements()), static_cast<double *>(local_sigma.data_elements()),
												 static_cast<double *>(local_exc.data_elements()), static_cast<double *>(local_vxc.data_elements()), static_cast<double *>(local_vsigma.data_elements()));
					gpu::sync();
					
					auto calc_vsigma = [func = ggafunctional.libxc_func_ptr()] (auto point){
						math::array<double, 1> local_density{{sqwave(point, 3)}};
						math::array<double, 1> local_sigma{{dot(gradient_sqwave(point, 3), gradient_sqwave(point, 3))}};
						math::array<double, 1> local_exc{1};
						math::array<double, 1> local_vxc{1};
						math::array<double, 1> local_vsigma{1};
						xc_gga_exc_vxc(func, 1, static_cast<double *>(local_density.data_elements()), static_cast<double *>(local_sigma.data_elements()),
													 static_cast<double *>(local_exc.data_elements()), static_cast<double *>(local_vxc.data_elements()), static_cast<double *>(local_vsigma.data_elements()));
						gpu::sync();
						return local_vsigma[0];
					};
					
					auto grad_vsigma = finite_difference_gradient5p(calc_vsigma, vec);

					auto calc_vsigmadn = [calc_vsigma] (auto point){
						return gradient_sqwave(point, 3)*calc_vsigma(point);
					};
					
					auto vxc_extra = finite_difference_divergence5p(calc_vsigmadn, vec);

					diff_ways = std::max(diff_ways, fabs(dot(grad_vsigma, gradient_sqwave(vec, 3)) + local_vsigma[0]*laplacian_sqwave(vec, 3) - vxc_extra));

					local_vxc[0] -= 2.0*vxc_extra;
					int_xc_energy += local_exc[0]*local_density[0]*rs.volume_element();
					
					diff = std::max(diff, fabs(local_vxc[0] - Vxc.cubic()[ix][iy][iz])*local_density[0]);
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&int_xc_energy, 1, std::plus<>{});
		
		CHECK(diff == Approx(0).margin(0.033132713));
		CHECK(diff_ways == Approx(0).margin(0.000001359));
		CHECK(Approx(Exc) == int_xc_energy);
		
		if(rs.part().contains(1)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(1))] == -0.5607887985_a);
		if(rs.part().contains(33)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(33))] == -1.1329131862_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(Vxc.linear()[rs.part().global_to_local(utils::global_index(rs.size() - 1))] == -1.1461742979_a);
	}

	SECTION("Uniform"){
		
		basis::real_space rs(cell, input::basis::cutoff_energy(20.0), cart_comm);
		
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
}

#endif
#endif
