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
					sigma = 0.0;
					for(int ix = 0; ix < vxc.basis().sizes()[0]; ix++){
						for(int iy = 0; iy < vxc.basis().sizes()[1]; iy++){
							for(int iz = 0; iz < vxc.basis().sizes()[2]; iz++){
								// Iterating over each vectorial components of the grad field-set at (ix,iy,iz) point in the space
								for(int idir = 0; idir < 3 ; idir++) sigma.cubic()[ix][iy][iz] += norm(grad_real.cubic()[ix][iy][iz][idir]);
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

#endif
