/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SOLVERS_POISSON
#define SOLVERS_POISSON

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

#include <math/complex.hpp>
#include <math/d3vector.hpp>
#include <multi/array.hpp>
#include <multi/adaptors/fftw.hpp>

namespace solvers {

  class poisson {

	public:

		template <class basis_type>
		auto solve(const basis_type & basis, const boost::multi::array<double, 3> & density){

			//For the moment we copy to a complex array.
			
			boost::multi::array<complex, 3> complex_density(extensions(density));

			for(int ix = 0; ix < basis.rsize()[0]; ix++){
				for(int iy = 0; iy < basis.rsize()[1]; iy++){
					for(int iz = 0; iz < basis.rsize()[2]; iz++){
						complex_density[ix][iy][iz] = density[ix][iy][iz];
					}
				}
			}

			auto complex_potential = solve(basis, complex_density);
			boost::multi::array<double, 3> potential(extensions(density));
			
			for(int ix = 0; ix < basis.rsize()[0]; ix++){
				for(int iy = 0; iy < basis.rsize()[1]; iy++){
					for(int iz = 0; iz < basis.rsize()[2]; iz++){
						potential[ix][iy][iz] = real(complex_potential[ix][iy][iz]);
					}
				}
			}

			return potential;			
		}
		
		template <class basis_type>
		auto solve(const basis_type & basis, const boost::multi::array<complex, 3> & density){
			namespace fftw = boost::multi::fftw;

			auto potential = fftw::dft(density, fftw::forward);

			const double scal = (-4.0*M_PI)/basis.rtotalsize();
			
			for(int ix = 0; ix < basis.gsize()[0]; ix++){
				for(int iy = 0; iy < basis.gsize()[1]; iy++){
					for(int iz = 0; iz < basis.gsize()[2]; iz++){

						if(basis.g_is_zero(ix, iy, iz)){
							potential[0][0][0] = 0;
							continue;
						}
						potential[ix][iy][iz] *= -scal/basis.g2(ix, iy, iz);
					}
				}
			}
			
			fftw::dft_inplace(potential, fftw::backward);

			return potential;
		}
		
	private:
		
  };    
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/plane_wave.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class solvers::poisson", "[poisson]") {

	using namespace Catch::literals;
	namespace multi = boost::multi;
	
  double ll = 10.0;

  ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll});
  basis::plane_wave pw(cell, {100, 100, 100} );

	multi::array<complex, 3> density(pw.rsize());
	solvers::poisson psolver;

	SECTION("Point charge"){
		
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					density[ix][iy][iz] = 0.0;
				}
			}
		}

		density[0][0][0] = -1.0;
		
		auto potential = psolver.solve(pw, density);
		
		double sumreal = 0.0;
		double sumimag = 0.0;
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					sumreal += fabs(real(potential[ix][iy][iz]));
					sumimag += fabs(imag(potential[ix][iy][iz]));
				}
			}
		}

		// These values haven't been validated against anything, they are
		// just for consistency. Of course the imaginary part has to be
		// zero, since the density is real.
		
		REQUIRE(sumreal == 59.7758543176_a);
		REQUIRE(sumimag == 3.87333e-13_a);
		
		REQUIRE(real(potential[0][0][0]) == -0.0241426581_a);
	}

	SECTION("Plane wave"){

		double kk = 2.0*M_PI/pw.rlength()[0];
		
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					double xx = pw.rvector(ix, iy, iz)[0];
					density[ix][iy][iz] = complex(cos(kk*xx), sin(kk*xx));
					
				}
			}
		}

		auto potential = psolver.solve(pw, density);

		double diff = 0.0;
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					diff += fabs(potential[ix][iy][iz] - 4*M_PI/kk/kk*density[ix][iy][iz]);
				}
			}
		}

		diff /= pw.size();
		
		REQUIRE(diff == 7.33009e-15_a);
	
	}


	SECTION("Real plane wave"){

		multi::array<complex, 3> rdensity(pw.rsize());

		double kk = 8.0*M_PI/pw.rlength()[1];
		
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					double yy = pw.rvector(ix, iy, iz)[1];
					rdensity[ix][iy][iz] = cos(kk*yy);
				}
			}
		}

		auto rpotential = psolver.solve(pw, rdensity);

		double diff = 0.0;
		for(int ix = 0; ix < pw.rsize()[0]; ix++){
			for(int iy = 0; iy < pw.rsize()[1]; iy++){
				for(int iz = 0; iz < pw.rsize()[2]; iz++){
					diff += fabs(rpotential[ix][iy][iz] - 4*M_PI/kk/kk*rdensity[ix][iy][iz]);
				}
			}
		}

		diff /= pw.size();
		
		REQUIRE(diff < 1e-8);

	}
	

}


#endif


#endif
