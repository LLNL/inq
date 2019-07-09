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
  
    auto to_fourier(const states::ks_states st, const basis::real_space & basis, const states::coefficients & phi) {
      
      namespace multi = boost::multi;
      namespace fftw = boost::multi::fftw;
      
      states::coefficients fphi(st, basis);
			
      //DATAOPERATIONS
      
      multi::array<complex, 3> fftgrid(basis.rsize());
      
      for(int ist = 0; ist < st.num_states(); ist++){
				
				// for the moment we have to copy to single grid, since the
				// fft interfaces assumes the transform is over the last indices					
				for(int ix = 0; ix < basis.rsize()[0]; ix++){
					for(int iy = 0; iy < basis.rsize()[1]; iy++){
						for(int iz = 0; iz < basis.rsize()[2]; iz++){
							fftgrid[ix][iy][iz] = phi.cubic[ix][iy][iz][ist];
						}
					}
				}
				
				fftw::dft_inplace(fftgrid, fftw::forward);
				
				for(int ix = 0; ix < basis.gsize()[0]; ix++){
					for(int iy = 0; iy < basis.gsize()[1]; iy++){
						for(int iz = 0; iz < basis.gsize()[2]; iz++){
							fphi.cubic[ix][iy][iz][ist] = fftgrid[ix][iy][iz];
						}
					}
				}
				
      }
      
      return fphi;    
    }
    
    void to_real_inplace(const states::ks_states st, const basis::real_space & basis, states::coefficients & phi){
      
      namespace multi = boost::multi;
      namespace fftw = boost::multi::fftw;
      
      //DATAOPERATIONS
			
      multi::array<complex, 3> fftgrid(basis.rsize());
      
      for(int ist = 0; ist < st.num_states(); ist++){
				
				// for the moment we have to copy to single grid, since the
				// fft interfaces assumes the transform is over the last indices					
				for(int ix = 0; ix < basis.gsize()[0]; ix++){
					for(int iy = 0; iy < basis.gsize()[1]; iy++){
						for(int iz = 0; iz < basis.gsize()[2]; iz++){
							fftgrid[ix][iy][iz] = phi.cubic[ix][iy][iz][ist];
						}
					}
				}
	
				fftw::dft_inplace(fftgrid, fftw::backward);

				double norm_factor = basis.num_points();
	
				for(int ix = 0; ix < basis.rsize()[0]; ix++){
					for(int iy = 0; iy < basis.rsize()[1]; iy++){
						for(int iz = 0; iz < basis.rsize()[2]; iz++){
							phi.cubic[ix][iy][iz][ist] = fftgrid[ix][iy][iz]/norm_factor;
						}
					}
				}
				
      }
      
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
