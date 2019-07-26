/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef HAMILTONIAN__KS_POTENTIAL
#define HAMILTONIAN__KS_POTENTIAL

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

#include <basis/field.hpp>
#include <functionals/lda.hpp>
#include <solvers/poisson.hpp>
#include <operations/sum.hpp>

namespace hamiltonian {
	template <class vexternal_type, class density_type>
	auto ks_potential(const vexternal_type & vexternal, const density_type & density) {

		solvers::poisson<basis::real_space> poisson_solver;
		
		assert(vexternal.basis() == density.basis()); //for the moment they must be equal
		
		auto vhartree = poisson_solver(density);

		vexternal_type exc(vexternal.basis());
		vexternal_type vxc(vexternal.basis());
		
		functionals::lda::xc_unpolarized(density.basis().size(), density, exc, vxc);

		auto vks = operations::sum(vexternal, vhartree, vxc);

		return vks;
	}
	
}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::ks_potential", "[ks_potential]"){

  using namespace Catch::literals;
	
}

#endif

#endif
