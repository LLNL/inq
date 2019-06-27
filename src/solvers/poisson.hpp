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
#include <basis/coefficients.hpp>

namespace solvers {

	template<class basis_type>
  class poisson {

	public:

		auto operator()(const basis::coefficients<basis_type, complex> & density){
			namespace fftw = boost::multi::fftw;

			basis::coefficients<basis_type, complex> potential(density.basis());
			
			potential.cubic = fftw::dft(density.cubic, fftw::forward);

			const double scal = (-4.0*M_PI)/potential.basis().size();
			
			for(int ix = 0; ix < potential.basis().gsize()[0]; ix++){
				for(int iy = 0; iy < potential.basis().gsize()[1]; iy++){
					for(int iz = 0; iz < potential.basis().gsize()[2]; iz++){

						if(potential.basis().g_is_zero(ix, iy, iz)){
							potential.cubic[0][0][0] = 0;
							continue;
						}
						potential.cubic[ix][iy][iz] *= -scal/potential.basis().g2(ix, iy, iz);
					}
				}
			}
			
			fftw::dft_inplace(potential.cubic, fftw::backward);

			return potential;
		}
		
		auto operator()(const basis::coefficients<basis_type, double> & density){

			//For the moment we copy to a complex array.
			
			basis::coefficients<basis_type, complex> complex_density(density.basis());

			//DATAOPERATIONS
			for(long ic = 0; ic < density.basis().size(); ic++) complex_density.linear[ic] = density.linear[ic];

			auto complex_potential = operator()(complex_density);

			basis::coefficients<basis_type, double> potential(density.basis());

			//DATAOPERATIONS
			for(long ic = 0; ic < potential.basis().size(); ic++) potential.linear[ic] = std::real(complex_potential.linear[ic]);

			return potential;
		}
		
	private:
		
  };    
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class solvers::poisson", "[poisson]") {

	using namespace Catch::literals;
	namespace multi = boost::multi;
	using namespace basis;
	
  double ll = 10.0;

  ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll});
  basis::real_space rs(cell, 493.48);

	REQUIRE(rs.rsize()[0] == 100);
	REQUIRE(rs.rsize()[1] == 100);
	REQUIRE(rs.rsize()[2] == 100);
	
	
	coefficients<real_space, complex> density(rs);
	solvers::poisson<basis::real_space> psolver;

	SECTION("Point charge"){
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					density.cubic[ix][iy][iz] = 0.0;
				}
			}
		}

		density.cubic[0][0][0] = -1.0;
		
		auto potential = psolver(density);
		
		double sumreal = 0.0;
		double sumimag = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					sumreal += fabs(real(potential.cubic[ix][iy][iz]));
					sumimag += fabs(imag(potential.cubic[ix][iy][iz]));
				}
			}
		}

		// These values haven't been validated against anything, they are
		// just for consistency. Of course the imaginary part has to be
		// zero, since the density is real.
		
		REQUIRE(sumreal == 59.7758543176_a);
		REQUIRE(sumimag == 3.87333e-13_a);
		
		REQUIRE(real(potential.cubic[0][0][0]) == -0.0241426581_a);
	}

	SECTION("Plane wave"){

		double kk = 2.0*M_PI/rs.rlength()[0];
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double xx = rs.rvector(ix, iy, iz)[0];
					density.cubic[ix][iy][iz] = complex(cos(kk*xx), sin(kk*xx));
					
				}
			}
		}

		auto potential = psolver(density);

		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					diff += fabs(potential.cubic[ix][iy][iz] - 4*M_PI/kk/kk*density.cubic[ix][iy][iz]);
				}
			}
		}

		diff /= rs.size();
		
		REQUIRE(diff == 7.33009e-15_a);
	
	}


	SECTION("Real plane wave"){

		coefficients<real_space, double> rdensity(rs);

		double kk = 8.0*M_PI/rs.rlength()[1];
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double yy = rs.rvector(ix, iy, iz)[1];
					rdensity.cubic[ix][iy][iz] = cos(kk*yy);
				}
			}
		}

		auto rpotential = psolver(rdensity);

		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					diff += fabs(rpotential.cubic[ix][iy][iz] - 4*M_PI/kk/kk*rdensity.cubic[ix][iy][iz]);
				}
			}
		}

		diff /= rs.size();
		
		REQUIRE(diff < 1e-8);

	}
	

}


#endif


#endif
