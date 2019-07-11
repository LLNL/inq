/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

#include <multi/adaptors/fftw.hpp>
#include <states/coefficients.hpp>
#include <cassert>

namespace operations {

  namespace space {

    auto to_fourier(const basis::coefficients_set<basis::real_space, complex> & phi){
      
      namespace multi = boost::multi;
      namespace fftw = boost::multi::fftw;
	
			basis::coefficients_set<basis::fourier_space, complex> fphi(phi.basis(), phi.set_size());
			
      //DATAOPERATIONS
      
      multi::array<complex, 3> fftgrid(phi.basis().rsize());
      
      for(int ist = 0; ist < phi.set_size(); ist++){
				
				// for the moment we have to copy to single grid, since the
				// fft interfaces assumes the transform is over the last indices					
				for(int ix = 0; ix < phi.basis().rsize()[0]; ix++){
					for(int iy = 0; iy < phi.basis().rsize()[1]; iy++){
						for(int iz = 0; iz < phi.basis().rsize()[2]; iz++){
							fftgrid[ix][iy][iz] = phi.cubic[ix][iy][iz][ist];
						}
					}
				}
				
				fftw::dft_inplace(fftgrid, fftw::forward);
				
				for(int ix = 0; ix < fphi.basis().gsize()[0]; ix++){
					for(int iy = 0; iy < fphi.basis().gsize()[1]; iy++){
						for(int iz = 0; iz < fphi.basis().gsize()[2]; iz++){
							fphi.cubic[ix][iy][iz][ist] = fftgrid[ix][iy][iz];
						}
					}
				}
				
      }
      
      return fphi;    
    }
    
    auto to_real(const basis::coefficients_set<basis::fourier_space, complex> & fphi){
      
      namespace multi = boost::multi;
      namespace fftw = boost::multi::fftw;

			basis::coefficients_set<basis::real_space, complex> phi(fphi.basis(), fphi.set_size());

      //DATAOPERATIONS
			
      multi::array<complex, 3> fftgrid(fphi.basis().rsize());
      
      for(int ist = 0; ist < fphi.set_size(); ist++){
				
				// for the moment we have to copy to single grid, since the
				// fft interfaces assumes the transform is over the last indices					
				for(int ix = 0; ix < fphi.basis().gsize()[0]; ix++){
					for(int iy = 0; iy < fphi.basis().gsize()[1]; iy++){
						for(int iz = 0; iz < fphi.basis().gsize()[2]; iz++){
							fftgrid[ix][iy][iz] = phi.cubic[ix][iy][iz][ist];
						}
					}
				}
	
				fftw::dft_inplace(fftgrid, fftw::backward);

				double norm_factor = phi.basis().num_points();
	
				for(int ix = 0; ix < phi.basis().rsize()[0]; ix++){
					for(int iy = 0; iy < phi.basis().rsize()[1]; iy++){
						for(int iz = 0; iz < phi.basis().rsize()[2]; iz++){
							phi.cubic[ix][iy][iz][ist] = fftgrid[ix][iy][iz]/norm_factor;
						}
					}
				}
				
      }

			return phi;
    }

	}
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::space", "[space]") {

	using namespace Catch::literals;
}


#endif

#endif
