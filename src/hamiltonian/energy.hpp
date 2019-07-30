/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef HAMILTONIAN__ENERGY
#define HAMILTONIAN__ENERGY

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


namespace hamiltonian {

	struct energy {
		
		double ion;
		double eigenvalues;
		double external;
		double hartree;
		double xc;
		double nvxc;

		energy(){
			ion = 0.0;
			eigenvalues = 0.0;
			external = 0.0;
			hartree = 0.0;
			xc = 0.0;
			nvxc = 0.0;
		}

		auto kinetic() const {
			return eigenvalues - 2.0*hartree - nvxc - external;
		}
		
		auto total() const {
			return ion + kinetic() + external + hartree + xc;
		}

		template <class out_type>
		void print(out_type & out) const {
		
			out << std::endl;
			out << "  total       = " << total()     << std::endl;
			out << "  kinetic     = " << kinetic()   << std::endl;
			out << "  eigenvalues = " << eigenvalues << std::endl;
			out << "  external    = " << external    << std::endl;
			out << "  hartree     = " << hartree     << std::endl;
			out << "  xc          = " << xc          << std::endl;
			out << "  intnvxc     = " << nvxc        << std::endl;
			out << std::endl;

		}
		
	};

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::energy", "[energy]"){

  using namespace Catch::literals;
	
}

#endif

#endif
