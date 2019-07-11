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
	using math::d3vector;

	double ecut = 23.0;
	double ll = 6.66;
	
	ions::geometry geo;
	ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
	basis::real_space rs(cell, ecut);
	
	basis::coefficients_set<basis::real_space, complex> phi(rs, 7);
	
	SECTION("Zero"){
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++) phi.cubic[ix][iy][iz][ist] = 0.0;
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().gsize()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().gsize()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().gsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(fphi.cubic[ix][iy][iz][ist]);
					}
				}
			}
		}
		
		diff /= fphi.cubic.num_elements();

		REQUIRE(diff < 1e-15);
		
		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++)	diff += fabs(phi.cubic[ix][iy][iz][ist]);
				}
			}
		}

		diff /= phi2.cubic.num_elements();

		REQUIRE(diff < 1e-15);
		
	}
	
	SECTION("Gaussian"){
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double r2 = rs.r2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_size(); ist++){
						double sigma = 0.5*(ist + 1);
						sigma = 1.0;
						phi.cubic[ix][iy][iz][ist] = exp(-sigma*r2);
					}
				}
			}
		}
		
		auto fphi = operations::space::to_fourier(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < fphi.basis().gsize()[0]; ix++){
			for(int iy = 0; iy < fphi.basis().gsize()[1]; iy++){
				for(int iz = 0; iz < fphi.basis().gsize()[2]; iz++){
					double g2 = fphi.basis().g2(ix, iy, iz);
					for(int ist = 0; ist < phi.set_size(); ist++){
						double sigma = 0.5*(ist + 1);
						sigma = 1.0;
						diff += fabs(fphi.cubic[ix][iy][iz][ist] - sqrt(M_PI/sigma)*exp(-0.25*g2/sigma));
					}
				}
			}
		}
		
		diff /= fphi.cubic.num_elements();

		std::cout << diff << std::endl;

		auto phi2 = operations::space::to_real(fphi);

		diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(phi.cubic[ix][iy][iz][ist] - phi2.cubic[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= phi2.cubic.num_elements();

		REQUIRE(diff < 1e-15);
		
		
	}
	
}


#endif

#endif
