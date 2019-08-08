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

#include <tinyformat/tinyformat.h>

namespace hamiltonian {

	struct energy {
		
		double ion;
		double ion_sr_lr;
		double self;
		double eigenvalues;
		double external;
		double hartree;
		double nvhartree;
		double xc;
		double nvxc;

		energy(){
			ion = 0.0;
			self = 0.0;
			eigenvalues = 0.0;
			external = 0.0;
			hartree = 0.0;
			nvhartree = 0.0;
			xc = 0.0;
			nvxc = 0.0;
		}

		auto coulomb() const {
			return hartree + ion + self;
		}
		
		auto kinetic() const {
			return eigenvalues - nvhartree - nvxc - external;
		}
		
		auto total() const {
			return kinetic() + hartree + self + ion + external + xc;
		}

		template <class out_type>
		void print(out_type & out) const {

			tfm::format(out, "\n");
			tfm::format(out, "  total       = %d16\n", total());			
			tfm::format(out, "  kinetic     = %d\n", kinetic());
			tfm::format(out, "  eigenvalues = %d\n", eigenvalues);
			tfm::format(out, "  coulomb     = %d\n", coulomb());
			tfm::format(out, "  hartree     = %d\n", hartree);
			tfm::format(out, "  nvhartree   = %d\n", nvhartree);
			tfm::format(out, "  external    = %d\n", external);
			tfm::format(out, "  xc          = %d\n", xc);
			tfm::format(out, "  intnvxc     = %d\n", nvxc);
			tfm::format(out, "  ion         = %d\n", ion);
			tfm::format(out, "  self        = %d\n", self);
			tfm::format(out, "\n");

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
