/* -*- indent-tabs-mode: t -*- */

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
		double ion_sr_lr;
		double eigenvalues;
		double external;
		double coulomb;
		double xc;
		double nvxc;

		energy(){
			ion = 0.0;
			ion_sr_lr = 0.0;
			eigenvalues = 0.0;
			external = 0.0;
			coulomb = 0.0;
			xc = 0.0;
			nvxc = 0.0;
		}

		auto kinetic() const {
			return eigenvalues - 2.0*coulomb - nvxc - external;
		}
		
		auto total() const {
			return ion + ion_sr_lr + eigenvalues - 2.0*coulomb - nvxc + coulomb + xc;
		}

		template <class out_type>
		void print(out_type & out) const {
		
			out << std::endl;
			out << "  total       = " << total()     << std::endl;
			out << "  eigenvalues = " << eigenvalues << std::endl;
			out << "  coulomb     = " << coulomb     << std::endl;
			out << "  xc          = " << xc          << std::endl;
			out << "  intnvxc     = " << nvxc        << std::endl;
			out << "  ion         = " << ion         << std::endl;
			out << "  ion sr lr   = " << ion_sr_lr   << std::endl;
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
