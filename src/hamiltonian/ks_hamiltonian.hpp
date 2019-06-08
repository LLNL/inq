#ifndef HAMILTONIAN_KS_HAMILTONIAN
#define HAMILTONIAN_KS_HAMILTONIAN

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

#include <states/ks_states.hpp>
#include <multi/adaptors/fftw.hpp>

namespace hamiltonian {
  template <class basis_type>
  class ks_hamiltonian {
		
  public:

    ks_hamiltonian(const basis_type & basis):
			scalar_potential(basis.rsize()){

			for(int ix = 0; ix < basis.rsize()[0]; ix++){
				for(int iy = 0; iy < basis.rsize()[1]; iy++){
					for(int iz = 0; iz < basis.rsize()[2]; iz++){
						scalar_potential[ix][iy][iz] = 0.0;
					}
				}
			}
    }

		boost::multi::array<double, 3> scalar_potential;

		template <class array_dim5>
		void apply(const basis_type & basis, const states::ks_states st, const array_dim5 && phi, array_dim5 && hphi){

			namespace multi = boost::multi;
			namespace fftw = boost::multi::fftw;
			
			multi::array<complex, 3> fftgrid(basis.rsize());

			//the Laplacian in Fourier space
			for(int ist = 0; ist < st.num_states(); ist++){
				for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++){

					// for the moment we have to copy to single grid, since the
					// fft interfaces assumes the transform is over the last indices					
					for(int ix = 0; ix < basis.rsize()[0]; ix++){
						for(int iy = 0; iy < basis.rsize()[1]; iy++){
							for(int iz = 0; iz < basis.rsize()[2]; iz++){
								fftgrid[ix][iy][iz] = phi[ix][iy][iz][ist][ispinor];
							}
						}
					}

					fftw::dft_inplace(fftgrid, fftw::forward);

					double scal = -0.5/basis.rtotalsize();
					for(int ix = 0; ix < basis.gsize()[0]; ix++){
						for(int iy = 0; iy < basis.gsize()[1]; iy++){
							for(int iz = 0; iz < basis.gsize()[2]; iz++){
								fftgrid[ix][iy][iz] *= -scal*basis.g2(ix, iy, iz);
							}
						}
					}
					
					fftw::dft_inplace(fftgrid, fftw::backward);
					
					for(int ix = 0; ix < basis.rsize()[0]; ix++){
						for(int iy = 0; iy < basis.rsize()[1]; iy++){
							for(int iz = 0; iz < basis.rsize()[2]; iz++){
								hphi[ix][iy][iz][ist][ispinor] = fftgrid[ix][iy][iz];
							}
						}
					}
					
				}
			}

			//the scalar local potential in real space
			for(int ix = 0; ix < basis.rsize()[0]; ix++){
				for(int iy = 0; iy < basis.rsize()[1]; iy++){
					for(int iz = 0; iz < basis.rsize()[2]; iz++){

						double vv  = scalar_potential[ix][iy][iz];
						
						for(int ist = 0; ist < st.num_states(); ist++){
							for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++) {
								hphi[ix][iy][iz][ist][ispinor] += vv*phi[ix][iy][iz][ist][ispinor];
							}
						}
						
					}
				}
			}
			
		}
		
  private:

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/plane_wave.hpp>

TEST_CASE("Class hamiltonian::ks_hamiltonian", "[ks_hamiltonian]"){

  using namespace Catch::literals;
  using math::d3vector;
  
  double ecut = 20.0;
  double ll = 10.0;
  
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::plane_wave pw(cell, ecut);
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);  

	states::ks_states::coeff phi(st.coeff_dimensions(pw.rsize()));
	states::ks_states::coeff hphi(st.coeff_dimensions(pw.rsize()));

	hamiltonian::ks_hamiltonian<basis::plane_wave> ham(pw);


	SECTION("Constant function"){
		
		
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					for(int ist = 0; ist < st.num_states(); ist++){
						for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++) {
							phi[0][ix][iy][iz][ist][ispinor] = 1.0;
						}
					}
				}
			}
		}
		
		ham.apply(pw, st, phi[0], hphi[0]);
		
		double diff = 0.0;
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					for(int ist = 0; ist < st.num_states(); ist++){
						for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++) {
							diff += fabs(hphi[0][ix][iy][iz][ist][ispinor] - 0.0);
						}
					}
				}
			}
		}

		diff /= hphi.num_elements();
		
		REQUIRE(diff == 0.0_a);
		
	}
	
	SECTION("Plane wave"){
		
		double kk = 2.0*M_PI/pw.rsize()[0];
		
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					for(int ist = 0; ist < st.num_states(); ist++){
						for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++) {
							phi[0][ix][iy][iz][ist][ispinor] = complex(cos(ist*kk*ix), sin(ist*kk*ix));
						}
					}
				}
			}
		}
		
		ham.apply(pw, st, phi[0], hphi[0]);
		
		double diff = 0.0;
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					for(int ist = 0; ist < st.num_states(); ist++){
						for(int ispinor = 0; ispinor < st.num_spinors(); ispinor++) {
							diff += fabs(hphi[0][ix][iy][iz][ist][ispinor] - 0.5*ist*kk*ist*kk*phi[0][ix][iy][iz][ist][ispinor]);
						}
					}
				}
			}
		}

		diff /= hphi.num_elements();
		
		REQUIRE(diff == 5.18263e-16_a);
		
	}
}

#endif

#endif
