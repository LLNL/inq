/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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
#include <basis/real_space.hpp>
#include <ions/geometry.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <multi/array.hpp>
#include <multi/adaptors/fftw.hpp>
#include <parser/input_file.hpp>
#include <states/ks_states.hpp>
#include <states/coefficients.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/scal.hpp>
#include <operations/orthogonalization.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <solvers/steepest_descent.hpp>

#include <complex>
#include <iostream>

using std::cout;

int main(int argc, char ** argv){

  parser::input_file input(argv[1]);
  
  auto coordinates_file = input.parse<std::string>("Coordinates");
  
  ions::geometry geo(coordinates_file);
  
  geo.info(std::cout);
  
  auto lx = input.parse<double>("Lx");	
  auto ly = input.parse<double>("Ly");
  auto lz = input.parse<double>("Lz");
  
  ions::UnitCell cell({lx, 0.0, 0.0}, {0.0, ly, 0.0}, {0.0, 0.0, lz});
  
  cell.info(std::cout);
  
  auto ecut = input.parse<double>("CutoffEnergy");
  
  basis::grid pw(cell, ecut);
  
  pw.info(std::cout);
  
  hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());
  
  pot.info(std::cout);
  
  states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, pot.num_electrons());
  
  st.info(std::cout);
  
  hamiltonian::ks_hamiltonian<basis::grid> ham(pw, cell, pot, geo);
  
  ham.info(std::cout);

  states::coefficients phi(st, pw);
	states::coefficients hphi(st, pw);

	operations::randomize(st, pw, phi);

	for(int ii = 0; ii < 2000; ii++){
		
		operations::scal_invsqrt(st, pw, operations::overlap_diagonal(st, pw, phi), phi);
		
		hphi = ham.apply(st, pw, phi);
		
		auto overlap = operations::overlap_diagonal(st, pw, hphi, phi);

		std::cout << ii << '\t' << std::scientific << real(overlap[0]) << std::endl;
		
		solvers::steepest_descent(st, pw, ham, phi);

	}
	
}
