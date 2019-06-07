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
		void solve(const basis_type & basis, const boost::multi::array<complex, 3> & density, boost::multi::array<complex, 3> & potential){

			namespace fftw = boost::multi::fftw;
			
			potential = fftw::dft(density, fftw::forward);

			const double scal = (-4.0*M_PI)/basis.rtotalsize();
			
			for(int ix = 0; ix < basis.rsize()[0]; ix++){
				for(int iy = 0; iy < basis.rsize()[1]; iy++){
					for(int iz = 0; iz < basis.rsize()[2]; iz++){

						math::d3vector g{ix*basis.gspacing()[0], iy*basis.gspacing()[1], iz*basis.gspacing()[2]};

						for(int idx = 0; idx != 3; ++idx) if(g[idx] > basis.glength()[idx]/2) g[idx] -= basis.glength()[idx];
					
						if(ix == 0 and iy == 0 and iz == 0){
							potential[0][0][0] = 0;
							continue;
						}
						
						potential[ix][iy][iz] *= -scal/norm(g);
					}
				}
			}
			
			fftw::dft_inplace(potential, fftw::backward);
		}
		
  public:

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

	for(int ix = 0; ix < pw.rsize()[0]; ix++){
		for(int iy = 0; iy < pw.rsize()[1]; iy++){
			for(int iz = 0; iz < pw.rsize()[2]; iz++){
				density[ix][iy][iz] = 0.0;
			}
		}
	}
	
  density[0][0][0] = -1.0;
  
  namespace fftw = multi::fftw;
  
  multi::array<complex, 3> potential(extensions(density));

	solvers::poisson psolver;

	psolver.solve(pw, density, potential);

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
	
	REQUIRE(sumreal == 239.1034172704_a);
	REQUIRE(sumimag == 1.54933e-12_a);

	REQUIRE(real(potential[0][0][0]) == -0.0965706326_a);
}


#endif


#endif

// Local Variables:
// mode: c++
// End:
