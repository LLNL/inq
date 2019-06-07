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


#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <basis/plane_wave.hpp>
#include <ions/geometry.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <multi/array.hpp>
#include <multi/adaptors/fftw.hpp>

#include <complex>
#include <iostream>

using std::cout;

int main(){

	using math::d3vector;

	double ecut = 30.0;
  double ll = 10.0;
  
 // ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");
 // hamiltonian::atomic_potential pot(geo.number_of_atoms(), geo.atoms());

  ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll});
  basis::plane_wave pw(cell, {100, 100, 100} );

  namespace multi = boost::multi;
  using complex = std::complex<double>;
  
  multi::array<complex, 3> density(pw.rsize());
  assert( density[1][2][3] == 0. );
  
  density[0][0][ 0] = -1.0;
  density[0][0][30] = +1.0;
  
  namespace fftw = multi::fftw;
  
  multi::array<complex, 3> potential = fftw::dft(density, fftw::forward);
  //	cout << "000 " << potential[0][0][0] << " " << << '\n';
  
  //	for(int ix = 0; ix != pw.rsize()[0]; ++ix)
  //		cout<< potential[ix][2][2] <<'\n';
  
#if 1
  for(int ix = 0; ix < pw.rsize()[0]; ix++){
		for(int iy = 0; iy < pw.rsize()[1]; iy++){
			for(int iz = 0; iz < pw.rsize()[2]; iz++){
				d3vector g{ix*pw.gspacing()[0], iy*pw.gspacing()[1], iz*pw.gspacing()[2]};
				for(int idx = 0; idx != 3; ++idx){if( g[idx] > pw.glength()[idx]/2) g[idx] -= pw.glength()[idx];};
				if(ix == 0 and iy == 0 and iz == 0){potential[0][0][0] = 0; continue;}
				potential[ix][iy][iz] *= -1./norm(g);
			}
		}
	}
#endif
	fftw::dft_inplace(potential, fftw::backward);
	for(int ix = 0; ix != pw.rsize()[0]; ++ix)
		cout<< real(potential[ix][0][0]) <<'\n';

	multi::array<complex, 3> potential2(extensions(potential), 0.);
//		cout<< density[ix][0][0] <<'\n';


#if 0

	fftw_inplace(potential);

	multi::array<double, 3> density_check(extensions(potential));
	for(int ix = 0; ix < pw.rsize()[0]; ix++){
		for(int iy = 0; iy < pw.rsize()[1]; iy++){
			for(int iz = 0; iz < pw.rsize()[2]; iz++){
				d3vector g = ix*pw.gspacing()[0] + iy*pw.gspacing()[1] + iz*pw.gspacing()[2];
				potential[ix][iy][iz] *= norm(g);
			}
		}
	}
	fftw_inplace(density_check);
	assert( density == density_check );

#endif

}

// Local variables:
// eval:(setq indent-tabs-mode t)
// End:
