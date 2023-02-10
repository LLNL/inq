/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__ENERGY
#define INQ__HAMILTONIAN__ENERGY

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

namespace inq {
namespace hamiltonian {

	struct energy {
		
		double ion;
		double ion_sr_lr;
		double eigenvalues;
		double external;
		double nonlocal;
		double hartree;
		double xc;
		double nvxc;
		double hf_exchange;

		energy(){
			ion = 0.0;
			eigenvalues = 0.0;
			external = 0.0;
			nonlocal = 0.0;
			hartree = 0.0;
			xc = 0.0;
			nvxc = 0.0;
			hf_exchange = 0.0;
		}

		auto kinetic() const {
			return eigenvalues - 2.0*hartree - nvxc - 2.0*hf_exchange - external - nonlocal;
		}
		
		auto total() const {
			return kinetic() + hartree + external + nonlocal + xc + hf_exchange + ion;
		}

		template <class out_type>
		void print(out_type & out) const {

			tfm::format(out, "\n");
			tfm::format(out, "  total          = %20.12f\n", total());			
			tfm::format(out, "  kinetic        = %20.12f\n", kinetic());
			tfm::format(out, "  eigenvalues    = %20.12f\n", eigenvalues);
			tfm::format(out, "  hartree        = %20.12f\n", hartree);
			tfm::format(out, "  external       = %20.12f\n", external);
			tfm::format(out, "  nonlocal       = %20.12f\n", nonlocal);
			tfm::format(out, "  xc             = %20.12f\n", xc);
			tfm::format(out, "  intnvxc        = %20.12f\n", nvxc);
			tfm::format(out, "  HF exchange    = %20.12f\n", hf_exchange);
			tfm::format(out, "  ion            = %20.12f\n", ion);
			tfm::format(out, "\n");

		}
		
		template<class OStream>
		friend OStream& operator<<(OStream& os, energy const& self){
			self.print(os);
			return os;
		}
		
	};

}
}
#endif

#ifdef INQ_HAMILTONIAN_ENERGY_UNIT_TEST
#undef INQ_HAMILTONIAN_ENERGY_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::energy", "[energy]"){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif

