/* -*- indent-tabs-mode: t -*- */

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
#include <basis/field.hpp>
#include <basis/fourier_space.hpp>

namespace solvers {

	template<class basis_type>
  class poisson {

	public:

		auto operator()(const basis::field<basis_type, complex> & density){
			if(density.basis().periodic_dimensions() == 3){
				return solve_periodic(density);
			} else {
				return solve_finite(density);
			}
		}

		auto solve_periodic(const basis::field<basis_type, complex> & density){
			namespace fftw = boost::multi::fftw;
			
			basis::fourier_space fourier_basis(density.basis());

			basis::field<basis::fourier_space, complex> potential_fs(fourier_basis);
			
			potential_fs.cubic() = fftw::dft(density.cubic(), fftw::forward);

			const double scal = (-4.0*M_PI)/potential_fs.basis().size();
			
			for(int ix = 0; ix < potential_fs.basis().gsize()[0]; ix++){
				for(int iy = 0; iy < potential_fs.basis().gsize()[1]; iy++){
					for(int iz = 0; iz < potential_fs.basis().gsize()[2]; iz++){
						
						if(potential_fs.basis().g_is_zero(ix, iy, iz)){
							potential_fs.cubic()[0][0][0] = 0;
							continue;
						}
						potential_fs.cubic()[ix][iy][iz] *= -scal/potential_fs.basis().g2(ix, iy, iz);
					}
				}
			}
			
			basis::field<basis_type, complex> potential_rs(density.basis());

			potential_rs.cubic() = fftw::dft(potential_fs.cubic(), fftw::backward);

			return potential_rs;
		}

		auto solve_finite(const basis::field<basis_type, complex> & density){
			namespace fftw = boost::multi::fftw;

			basis::field<basis_type, complex> potential2x(density.basis().enlarge(2));

			potential2x = 0.0;
			
			for(int ix = 0; ix < density.basis().rsize()[0]; ix++){
				for(int iy = 0; iy < density.basis().rsize()[1]; iy++){
					for(int iz = 0; iz < density.basis().rsize()[2]; iz++){
						potential2x.cubic()[ix][iy][iz] = density.cubic()[ix][iy][iz];
					}
				}
			}
			
			fftw::dft_inplace(potential2x.cubic(), fftw::forward);

			basis::fourier_space fourier_basis(potential2x.basis());
			
			const auto scal = (-4.0*M_PI)/fourier_basis.size();
			const auto cutoff_radius = 0.5*fourier_basis.min_rlength();
			
			for(int ix = 0; ix < fourier_basis.gsize()[0]; ix++){
				for(int iy = 0; iy < fourier_basis.gsize()[1]; iy++){
					for(int iz = 0; iz < fourier_basis.gsize()[2]; iz++){
						
						// this is the kernel of C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006).
						if(fourier_basis.g_is_zero(ix, iy, iz)){
							potential2x.cubic()[0][0][0] *= scal*cutoff_radius*cutoff_radius/2.0;
							continue;
						}
						auto g2 = fourier_basis.g2(ix, iy, iz);
						//potential2x.cubic()[ix][iy][iz] *= -scal/fourier_basis.g2(ix, iy, iz);
						potential2x.cubic()[ix][iy][iz] *= scal*(1.0 - cos(cutoff_radius*sqrt(g2)))/g2;
						//						potential2x.cubic()[ix][iy][iz] *= scal;
						
					}
				}
			}

			fftw::dft_inplace(potential2x.cubic(), fftw::backward);
			
			basis::field<basis_type, complex> potential(density.basis());

			for(int ix = 0; ix < potential.basis().rsize()[0]; ix++){
				for(int iy = 0; iy < potential.basis().rsize()[1]; iy++){
					for(int iz = 0; iz < potential.basis().rsize()[2]; iz++){
						potential.cubic()[ix][iy][iz] = potential2x.cubic()[ix][iy][iz];
					}
				}
			}
			
			return potential;
		}
		
		auto operator()(const basis::field<basis_type, double> & density){

			//For the moment we copy to a complex array.
			
			basis::field<basis_type, complex> complex_density(density.basis());

			//DATAOPERATIONS
			for(long ic = 0; ic < density.basis().size(); ic++) complex_density[ic] = density[ic];

			auto complex_potential = operator()(complex_density);

			basis::field<basis_type, double> potential(density.basis());

			//DATAOPERATIONS
			for(long ic = 0; ic < potential.basis().size(); ic++) potential[ic] = std::real(complex_potential[ic]);

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

	{
		
		double ll = 10.0;
		
		ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll});
		basis::real_space rs(cell, input::basis::spacing(0.1));

		SECTION("Grid periodic"){
		
			REQUIRE(rs.periodic_dimensions() == 3);
			
			REQUIRE(rs.rsize()[0] == 100);
			REQUIRE(rs.rsize()[1] == 100);
			REQUIRE(rs.rsize()[2] == 100);

		}
		
		field<real_space, complex> density(rs);
		solvers::poisson<basis::real_space> psolver;
		
		SECTION("Point charge"){
		
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						density.cubic()[ix][iy][iz] = 0.0;
					}
				}
			}

			density.cubic()[0][0][0] = -1.0;
		
			auto potential = psolver(density);
		
			double sumreal = 0.0;
			double sumimag = 0.0;
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						sumreal += fabs(real(potential.cubic()[ix][iy][iz]));
						sumimag += fabs(imag(potential.cubic()[ix][iy][iz]));
					}
				}
			}

			// These values haven't been validated against anything, they are
			// just for consistency. Of course the imaginary part has to be
			// zero, since the density is real.
		
			REQUIRE(sumreal == 59.7758543176_a);
			REQUIRE(sumimag == 3.87333e-13_a);
		
			REQUIRE(real(potential.cubic()[0][0][0]) == -0.0241426581_a);
		}

		SECTION("Plane wave"){

			double kk = 2.0*M_PI/rs.rlength()[0];
		
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						double xx = rs.rvector(ix, iy, iz)[0];
						density.cubic()[ix][iy][iz] = complex(cos(kk*xx), sin(kk*xx));
					
					}
				}
			}

			auto potential = psolver(density);

			double diff = 0.0;
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						diff += fabs(potential.cubic()[ix][iy][iz] - 4*M_PI/kk/kk*density.cubic()[ix][iy][iz]);
					}
				}
			}

			diff /= rs.size();
		
			REQUIRE(diff < 1.0e-14);
	
		}


		SECTION("Real plane wave"){

			field<real_space, double> rdensity(rs);

			double kk = 8.0*M_PI/rs.rlength()[1];
		
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						double yy = rs.rvector(ix, iy, iz)[1];
						rdensity.cubic()[ix][iy][iz] = cos(kk*yy);
					}
				}
			}

			auto rpotential = psolver(rdensity);

			double diff = 0.0;
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						diff += fabs(rpotential.cubic()[ix][iy][iz] - 4*M_PI/kk/kk*rdensity.cubic()[ix][iy][iz]);
					}
				}
			}

			diff /= rs.size();
		
			REQUIRE(diff < 1e-8);

		}
	}

	{

		const double ll = 9.0;
		
		ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll}, 0);
		basis::real_space rs(cell, input::basis::spacing(1.0));

		solvers::poisson<basis::real_space> psolver;

		SECTION("Grid finite"){		

			REQUIRE(rs.periodic_dimensions() == 0);
			
			REQUIRE(rs.rsize()[0] == 90);
			REQUIRE(rs.rsize()[1] == 90);
			REQUIRE(rs.rsize()[2] == 90);

		}
		
		field<real_space, complex> density(rs);

		SECTION("Point charge finite"){

			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						density.cubic()[ix][iy][iz] = 0.0;
						if(norm(rs.rvector(ix, iy, iz)) < 1e-10) density.cubic()[ix][iy][iz] = -1.0;
					}
				}
			}
			
			auto potential = psolver(density);
		
			double sumreal = 0.0;
			double sumimag = 0.0;
			for(int ix = 0; ix < rs.rsize()[0]; ix++){
				for(int iy = 0; iy < rs.rsize()[1]; iy++){
					for(int iz = 0; iz < rs.rsize()[2]; iz++){
						double r2 = norm(rs.rvector(ix, iy, iz));
						
						std::cout << ix << '\t' << iy << '\t' << iz << "\t->\t" << std::real(potential.cubic()[ix][iy][iz]) << '\t' << 1/sqrt(r2) << std::endl;
							//1/sqrt(r2) << '\t' << std::real(potential.cubic()[ix][iy][iz]) << '\t' << 1/sqrt(r2)/std::real(potential.cubic()[ix][iy][iz]) << std::endl;
					}
				}
			}

			// These values haven't been validated against anything, they are
			// just for consistency. Of course the imaginary part has to be
			// zero, since the density is real.
		
			//			REQUIRE(sumreal == 59.7758543176_a);
			//			REQUIRE(sumimag == 3.87333e-13_a);
		
			REQUIRE(real(potential.cubic()[1][0][0]) == -0.0241426581_a);
		}

	}

}


#endif


#endif
