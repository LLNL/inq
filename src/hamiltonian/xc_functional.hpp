/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_FUNCTIONAL
#define INQ__HAMILTONIAN__XC_FUNCTIONAL

/*
 Copyright (C) 2019 Xavier Andrade

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
#include <basis/field.hpp>

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
				case XC_FAMILY_LDA:
					xc_lda_exc_vxc(&func_, size, density.data(), exc.data(), vxc.data());
					break;
				case XC_FAMILY_GGA:{
					// How to compute Vxc terms for GGA http://mbpt-domiprod.wikidot.com/calculation-of-gga-kernel
					auto grad_real = operations::gradient(density);
					// Compute sigma as a square of the gradient in the Real space
					basis::field<basis::real_space, double> sigma(vxc.basis());
					
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								sigma.cubic()[ix][iy][iz] = 0.0;
								for(int idir = 0; idir < 3 ; idir++) sigma.cubic()[ix][iy][iz] += grad_real.cubic()[ix][iy][iz][idir]*grad_real.cubic()[ix][iy][iz][idir];
							}
						}
					}
					// Initialize derivative of xc energy vsigma = d Exc/d sigma as a scalar field
					basis::field<basis::real_space, double> vsigma(vxc.basis());
					//Call libxc to computer vxc and vsigma
					xc_gga_exc_vxc(&func_, size, density.data(), sigma.data(), exc.data(), vxc.data(), vsigma.data());
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

					auto grad_vsigma = operations::gradient(vsigma);
					
					//Taking diverdence of [vsigma * Grad(n)]
					auto divvxcextra = operations::divergence(vxc_extra);
					// Add extra component to Vxc
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){	// Iterating over x-,y- and z- components of the each gradient field-set
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								// Iterating over each (ix,iy,iz) point in the Real space
								vxc.cubic()[ix][iy][iz] -= 2.0*divvxcextra.cubic()[ix][iy][iz];
								//								vxc.cubic()[ix][iy][iz] = vsigma.cubic()[ix][iy][iz];
								//vxc.cubic()[ix][iy][iz] = grad_vsigma.cubic()[ix][iy][iz][0];
							}
						}
					}
					break;
				}
				case XC_FAMILY_HYB_GGA:
					assert(false);
					break;
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
#include <math/array.hpp>

	auto gaussian(inq::math::vec3d rr){
		return pow(M_PI, -1.5)*exp(-norm(rr)); // sigma = 1/sqrt(2)
	}

	auto dgaussian(inq::math::vec3d rr){
		inq::math::vec3d ff;
		for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian(rr);
		return ff;
	}

auto laplacian_gaussian(inq::math::vec3d rr){
	return 4.0*(rr|rr)*gaussian(rr) - 2.0*gaussian(rr);
}

template <class FunctionType, class VecType>
auto finite_difference_gradient(FunctionType const & function, VecType const & point){
	double const delta = 0.1;

	VecType gradient = {0.0, 0.0, 0.0};
						
	for(int idir = 0; idir < 3; idir++){
		for(auto isign : {-1, 1}){
			auto point_plus_delta = point;
			point_plus_delta[idir] += isign*delta;
			gradient[idir] += isign*function(point_plus_delta)/(2.0*delta);
		}
	}

	return gradient;
}

TEST_CASE("function hamiltonian::xc_functional", "[hamiltonian::xc_functional]") {

	using namespace inq;
	using namespace Catch::literals;
	using namespace operations;
	using math::vec3d;

	//UnitCell size
	double lx = 9;
	double ly = 12;
	double lz = 10;

	ions::geometry geo;
	ions::UnitCell cell(vec3d(lx, 0.0, 0.0), vec3d(0.0, ly, 0.0), vec3d(0.0, 0.0, lz));

	basis::real_space rs(cell, input::basis::cutoff_energy(20.0));

	SECTION("LDA"){

		basis::field<basis::real_space, double> gaussian_field(rs);

		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					gaussian_field.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}
	
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussian_vxc(rs);
		double gaussian_exc;

		ldafunctional(gaussian_field, gaussian_exc, gaussian_vxc);
		
		CHECK(gaussian_exc == -0.270646_a);
		double int_xc_energy = 0.0;

		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					auto local_density = gaussian(vec);
					double local_exc, local_vxc;
					xc_lda_exc_vxc(&ldafunctional.libxc_func(), 1, &local_density, &local_exc, &local_vxc);
					CHECK(Approx(local_vxc) == gaussian_vxc.cubic()[ix][iy][iz]);
					int_xc_energy += local_exc*local_density*rs.volume_element();
				}
			}
		}
		
		CHECK(Approx(gaussian_exc) == int_xc_energy);
		CHECK(gaussian_vxc.linear()[1] == -0.5111609291_a);
		CHECK(gaussian_vxc.linear()[333] == -0.0000452004_a);
		CHECK(gaussian_vxc.linear()[rs.size()-1] == -0.4326883849_a);

	}

	SECTION("GGA"){
		
		basis::field<basis::real_space, double> gaussian_density(rs);
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					gaussian_density.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}
		
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		basis::field<basis::real_space, double> gaussian_vxc(rs);
		double gaussian_exc;

		ggafunctional(gaussian_density, gaussian_exc, gaussian_vxc);

		CHECK(gaussian_exc == -0.3139862364_a);

		double int_xc_energy = 0.0;
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					auto local_density = gaussian(vec);
					auto local_sigma = dgaussian(vec) | dgaussian(vec);
					double local_exc,local_vxc, local_vsigma;
					xc_gga_exc_vxc(&ggafunctional.libxc_func(), 1, &local_density, &local_sigma, &local_exc, &local_vxc, &local_vsigma);

					auto calc_vsigma = [func = &ggafunctional.libxc_func()] (auto point){
															 auto local_density = gaussian(point);
															 auto local_sigma = dgaussian(point) | dgaussian(point);
															 double local_exc, local_vxc,  local_vsigma;
															 
															 xc_gga_exc_vxc(func, 1, &local_density, &local_sigma, &local_exc, &local_vxc, &local_vsigma);

															 return local_vsigma;
														 };
					
					auto grad_vsigma = finite_difference_gradient(calc_vsigma, vec);

					local_vxc -= 2.0*((grad_vsigma|dgaussian(vec)) + local_vsigma*laplacian_gaussian(vec));


					//std::cout << local_density << '\t' << local_vxc << '\t' << gaussian_vxc.cubic()[ix][iy][iz] << std::endl;
					//std::cout << local_density << '\t' << local_vsigma << '\t' << gaussian_vxc.cubic()[ix][iy][iz] << std::endl;
					//std::cout << local_density << '\t' << grad_vsigma[0] << '\t' << gaussian_vxc.cubic()[ix][iy][iz] << std::endl;
					// Local divergence ???
					//CHECK(Approx(local_vxc) == gaussian_vxc.cubic()[ix][iy][iz]);

					int_xc_energy += local_exc*local_density*rs.volume_element();
				}
			}
		}
		
		CHECK(Approx(gaussian_exc) == int_xc_energy);
		CHECK(gaussian_vxc.linear()[1] == -0.5931073473_a);
		CHECK(gaussian_vxc.linear()[33] == -0.0166346166_a);
		CHECK(gaussian_vxc.linear()[rs.size()-1] == -0.5021400086_a);

	}
	
	SECTION("UNIFORM"){ //Check LDA==GGA for unifrom electronic density
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
		basis::field<basis::real_space, double> gaussian_vxcLDA(rs) , gaussian_vxcGGA(rs);
		double gaussian_excLDA, gaussian_excGGA;

		ggafunctional(gaussian_field, gaussian_excLDA, gaussian_vxcLDA);
		ldafunctional(gaussian_field, gaussian_excGGA, gaussian_vxcGGA);
		CHECK(gaussian_excLDA == gaussian_excGGA);
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					CHECK(Approx(gaussian_vxcLDA.cubic()[ix][iy][iz]) == gaussian_vxcGGA.cubic()[ix][iy][iz]);
				}
			}
		}
	}
/*
	SECTION("NONUNIFORM"){ //Check LDA==GGA for nonunifrom electronic density in the case grad[n]=0
		basis::field<basis::real_space, double> gaussian_field(rs);
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					auto vec = rs.rvector(ix, iy, iz);
					gaussian_field.cubic()[ix][iy][iz] = gaussian(vec);
				}
			}
		}
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE);
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X);
		basis::field<basis::real_space, double> gaussian_vxcLDA(rs) , gaussian_vxcGGA(rs);
		double gaussian_excLDA, gaussian_excGGA;

		ggafunctional(gaussian_field, gaussian_excLDA, gaussian_vxcLDA);
		ldafunctional(gaussian_field, gaussian_excGGA, gaussian_vxcGGA);
		std::cout << gaussian_excLDA << "  " << gaussian_excGGA << "\n";
		
		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					std::cout << gaussian_vxcLDA.cubic()[ix][iy][iz]  << "  " << double(gaussian_vxcGGA.cubic()[ix][iy][iz]) << '\n';
				}
			}
		}
	}
*/
/*
	SECTION("QE_DENSITY"){
		std::ifstream cubefile("starting_density.cube");
		xsize << cubefile;
		ysize << cubefile;
		zsize << cubefile;

		//Back reverse engineering for the basis-set and system

		for(int ix = 0; ix < rs.sizes()[0]; ix++){
			for(int iy = 0; iy < rs.sizes()[1]; iy++){
				for(int iz = 0; iz < rs.sizes()[2]; iz++){
					cubefile >> density.cubic()[ix][iy][iz]
				}
			}
		}
	}	
*/
}


#endif
#endif
