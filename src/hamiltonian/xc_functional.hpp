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
#include <operations/divergence.hpp>
#include <operations/laplacian.hpp>
#include <basis/field.hpp>



namespace inq {
namespace hamiltonian {
	class xc_functional {

		public:
			
	//	static constexpr double sigma_threshold = 1.0e-15;
	//	static constexpr double dens_threshold = 1.0e-8;	

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

			field_type exc(vxc.skeleton());
			unpolarized(density.linear().size(), density, exc, vxc);

			xc_energy = operations::integral_product(density, exc);
				
		}

		auto & libxc_func() const {
			return func_;
		}

	private:
		
		template <class density_type, class exc_type, class vxc_type>
		void unpolarized(long size, density_type const & density, exc_type & exc, vxc_type & vxc) const{
			switch(func_.info->family) {
				case XC_FAMILY_LDA:{
					auto param_lda = func_;
		//			xc_func_set_dens_threshold(&param_lda, inq::hamiltonian::xc_functional::dens_threshold);
					xc_lda_exc_vxc(&param_lda, size, density.data(), exc.data(), vxc.data());
					break;
				}
				case XC_FAMILY_GGA:{
					// How to compute Vxc terms for GGA http://mbpt-domiprod.wikidot.com/calculation-of-gga-kernel
					auto grad_real = operations::gradient(density);
					// Compute sigma as a square of the gradient in the Real space
					basis::field<basis::real_space, double> sigma(vxc.basis());
					//std::ofstream generic("genericvsigma.dat");
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								sigma.cubic()[ix][iy][iz] = 0.0;
								for(int idir = 0; idir < 3 ; idir++) sigma.cubic()[ix][iy][iz] += grad_real.cubic()[ix][iy][iz][idir]*grad_real.cubic()[ix][iy][iz][idir];
					//			if (fabs(sigma.cubic()[ix][iy][iz]) < inq::hamiltonian::xc_functional::sigma_threshold) sigma.cubic()[ix][iy][iz] = 0.0;
							}
						}
					}
					//Call libxc to computer vxc and vsigma
					basis::field<basis::real_space, double> vsigma(vxc.basis());
					auto param = func_;	
	//				xc_func_set_dens_threshold(&param, inq::hamiltonian::xc_functional::dens_threshold);
					xc_gga_exc_vxc(&param, size, density.data(), sigma.data(), exc.data(), vxc.data(), vsigma.data());
					//Compute extra term for Vxc using diverdence: Vxc_extra=2 nabla *[vsigma*grad(n)]
					basis::field_set<basis::real_space, double> vxc_extra(vxc.basis(), 3);
					//Compute field-set as a product between vsigma(0) and gradient field-set
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){	// Iterating over x-,y- and z- components of the each gradient field-set
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								// Iterating over each vectorial components of the grad field-set at each (ix,iy,iz) point in the Real space
								for(int idir = 0; idir < 3 ; idir++) vxc_extra.cubic()[ix][iy][iz][idir] = vsigma.cubic()[ix][iy][iz]*grad_real.cubic()[ix][iy][iz][idir];
							}
						}
					}
					//Taking diverdence of [vsigma * Grad(n)]
					auto divvxcextra = operations::divergence(vxc_extra);
					// Add extra component to Vxc
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){	// Iterating over x-,y- and z- components of the each gradient field-set
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								// Iterating over each (ix,iy,iz) point in the Real space
								vxc.cubic()[ix][iy][iz] -= 2.0*divvxcextra.cubic()[ix][iy][iz];
							}
						}
					}
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

#ifdef INQ_UNIT_TEST

#include <catch2/catch.hpp>
#include <operations/randomize.hpp>
#include <operations/gradient.hpp>
#include <basis/field.hpp>
#include <math/array.hpp>
#include <utils/finite_difference.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>

	auto sqwave(inq::math::vec3d rr, int n){
		inq::math::vec3d kvec = 2.0 * M_PI * inq::math::vec3d(1.0/9.0, 1.0/12.0, 1.0/10.0);
		return  sin(n*(kvec | rr))*sin(n*(kvec | rr));
	}

	auto laplacian_sqwave(inq::math::vec3d rr, int n){
		inq::math::vec3d kvec = 2.0 * M_PI * inq::math::vec3d(1.0/9.0, 1.0/12.0, 1.0/10.0);
		auto ff = 0.0;
		for(int idir = 0; idir < 3 ; idir++) ff += 2*n*n*kvec[idir]*kvec[idir]*cos(2*n*(kvec|rr));
		return ff;
	}

	auto gradient_sqwave(inq::math::vec3d rr, int n){
		inq::math::vec3d kvec = 2.0 * M_PI * inq::math::vec3d(1.0/9.0, 1.0/12.0, 1.0/10.0);
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = 2*n*kvec[idir]*cos(n*(kvec | rr))*sin(n*(kvec | rr));
		return ff;
	}

	auto gaussiansigma(inq::math::vec3d rr, double sigma){
		return  pow(2*M_PI, -1.5)*pow(sigma, -3.0)*exp(-0.5*norm(rr)*pow(sigma, -2.0));
	}
	

	auto gaussian(inq::math::vec3d rr){  // sigma = 1/sqrt(2)
		return  pow(M_PI, -1.5)*exp(-norm(rr));
	}

	auto dgaussian(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian(rr);
		return ff;
	}

	auto laplacian_gaussian(inq::math::vec3d rr){
		return 4.0*(rr|rr)*gaussian(rr) - 6.0*gaussian(rr);
	}

TEST_CASE("function hamiltonian::xc_functional", "[hamiltonian::xc_functional]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace operations;
	using namespace inq::utils;
	using math::vec3d;

	//UnitCell size
	double lx = 9;
	double ly = 12;
	double lz = 10;

	ions::geometry geo;
	ions::UnitCell cell(vec3d(lx, 0.0, 0.0), vec3d(0.0, ly, 0.0), vec3d(0.0, 0.0, lz));

	SECTION("LDA"){
		basis::real_space rslda(cell, input::basis::cutoff_energy(20.0));
		basis::field<basis::real_space, double> gaussian_field(rslda);
		for(int ix = 0; ix < rslda.sizes()[0]; ix++){
			for(int iy = 0; iy < rslda.sizes()[1]; iy++){
				for(int iz = 0; iz < rslda.sizes()[2]; iz++){
					auto vec = rslda.rvector(ix, iy, iz);
					gaussian_field.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}

		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussianVxc(rslda);
		double gaussianExc;

		ldafunctional(gaussian_field, gaussianExc, gaussianVxc);
		CHECK(gaussianExc == -0.270646_a);
		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rslda.sizes()[0]; ix++){
			for(int iy = 0; iy < rslda.sizes()[1]; iy++){
				for(int iz = 0; iz < rslda.sizes()[2]; iz++){
					auto vec = rslda.rvector(ix, iy, iz);
					auto local_density = gaussian(vec);
					double local_exc, local_vxc;
					auto param_lda = ldafunctional.libxc_func();
		//			xc_func_set_dens_threshold(&param_lda, inq::hamiltonian::xc_functional::dens_threshold);
					xc_lda_exc_vxc(&param_lda, 1, &local_density, &local_exc, &local_vxc);
					CHECK(Approx(local_vxc) == gaussianVxc.cubic()[ix][iy][iz]);
					int_xc_energy += local_exc*local_density*rslda.volume_element();
				}
			}
		}
	CHECK(Approx(gaussianExc) == int_xc_energy);
	CHECK(gaussianVxc.linear()[1] == -0.5111609291_a);
	CHECK(gaussianVxc.linear()[8233] == -0.00406881_a);
	CHECK(Approx(gaussianVxc.linear()[233]).margin(1.0e-10) == -0.00000039110);
	CHECK(gaussianVxc.linear()[rslda.size()-1] == -0.4326883849_a);
	}

	SECTION("GGA"){
		basis::real_space rsgga(cell, input::basis::cutoff_energy(90.0));
		basis::field<basis::real_space, double> field(rsgga);
		for(int ix = 0; ix < rsgga.sizes()[0]; ix++){
			for(int iy = 0; iy < rsgga.sizes()[1]; iy++){
				for(int iz = 0; iz < rsgga.sizes()[2]; iz++){
					auto vec = rsgga.rvector(ix, iy, iz);
					field.cubic()[ix][iy][iz] = sqwave(vec, 3);
				}
			}
		}
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		basis::field<basis::real_space, double> Vxc(rsgga);
		basis::field<basis::real_space, double> local_vsigma_output(rsgga);

		double Exc = 0.0;
		ggafunctional(field, Exc, Vxc);

		CHECK(Exc == -393.4604748792_a);
		std::ofstream fout("fout.dat");
		double diff_ways = 0.0;
		double diff = 0.0;
		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rsgga.sizes()[0]; ix++){
			for(int iy = 0; iy < rsgga.sizes()[1]; iy++){
				for(int iz = 0; iz < rsgga.sizes()[2]; iz++){
					auto vec = rsgga.rvector(ix, iy, iz);
					double local_exc = 0.0;
					double local_vxc = 0.0;
					double local_vsigma = 0.0;
					auto local_density = sqwave(vec, 3);
					auto local_sigma = gradient_sqwave(vec, 3) | gradient_sqwave(vec, 3);
		//			if (fabs(local_sigma) < inq::hamiltonian::xc_functional::sigma_threshold) local_sigma = inq::hamiltonian::xc_functional::sigma_threshold*1.1e0;
					auto param = ggafunctional.libxc_func();
		//			xc_func_set_dens_threshold(&param, inq::hamiltonian::xc_functional::dens_threshold);
					xc_gga_exc_vxc(&param, 1, &local_density, &local_sigma, &local_exc, &local_vxc, &local_vsigma);

					auto calc_vsigma = [func = ggafunctional.libxc_func()] (auto point){
													auto local_density = sqwave(point, 3);
													auto local_sigma = gradient_sqwave(point, 3) | gradient_sqwave(point, 3);
													double local_exc, local_vxc, local_vsigma;
													auto param = func;
							//						if (fabs(local_sigma) < inq::hamiltonian::xc_functional::sigma_threshold) local_sigma = inq::hamiltonian::xc_functional::sigma_threshold*1.1e0;
							//						xc_func_set_dens_threshold(&param, inq::hamiltonian::xc_functional::dens_threshold);
													xc_gga_exc_vxc(&param, 1, &local_density, &local_sigma, &local_exc, &local_vxc, &local_vsigma);
													return local_vsigma;
													};
					auto grad_vsigma = finite_difference_gradient5p(calc_vsigma, vec);

					auto calc_vsigmadn = [func = ggafunctional.libxc_func()] (auto point){
													auto local_density = sqwave(point, 3);
													auto local_sigma = gradient_sqwave(point, 3) | gradient_sqwave(point, 3);
													double local_exc, local_vxc, local_vsigma;
													auto param = func;
						//							if (fabs(local_sigma) < inq::hamiltonian::xc_functional::sigma_threshold) local_sigma = inq::hamiltonian::xc_functional::sigma_threshold*1.1e0;
						//							xc_func_set_dens_threshold(&param, inq::hamiltonian::xc_functional::dens_threshold);
													xc_gga_exc_vxc(&param, 1, &local_density, &local_sigma, &local_exc, &local_vxc, &local_vsigma);
													return local_vsigma*gradient_sqwave(point, 3);
													};
					auto vxc_extra = finite_difference_divergence5p(calc_vsigmadn, vec);

					diff_ways = std::max(diff_ways, fabs((grad_vsigma | gradient_sqwave(vec, 3)) + local_vsigma*laplacian_sqwave(vec, 3) - vxc_extra));

					local_vxc -= 2.0*vxc_extra;
					int_xc_energy += local_exc*local_density*rsgga.volume_element();

					diff = std::max(diff, fabs(local_vxc - Vxc.cubic()[ix][iy][iz])*local_density);
			}
		}
	}
	CHECK(Approx(diff_ways).margin(1.0e-10) == 0.0000013589);
	CHECK(Approx(diff) == 0.033132712);
	CHECK(Approx(Exc) == int_xc_energy);
	CHECK(Vxc.linear()[1] == -0.5607887985_a);
	CHECK(Vxc.linear()[33] == -1.1329131862_a);
	CHECK(Vxc.linear()[rsgga.size()-1] == -1.1461742979_a);
	}

	SECTION("UNIFORM"){ //Check LDA==GGA for unifrom electronic density
		basis::real_space rs(cell, input::basis::cutoff_energy(20.0));
		basis::field<basis::real_space, double> gaussian_field(rs);
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
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
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					CHECK(Approx(gaussianVxcLDA.cubic()[ix][iy][iz]) == gaussianVxcGGA.cubic()[ix][iy][iz]);
				}
			}
		}

	}
}

#endif
#endif
