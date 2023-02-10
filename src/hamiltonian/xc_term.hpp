/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_TERM
#define INQ__HAMILTONIAN__XC_TERM

/*
 Copyright (C) 2023 Xavier Andrade

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
#include <solvers/poisson.hpp>
#include <observables/density.hpp>
#include <operations/add.hpp>
#include <operations/integral.hpp>
#include <input/interaction.hpp>
#include <hamiltonian/xc_functional.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <perturbations/none.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

class xc_term {
	
public:
	
	xc_term(input::interaction interaction, int const spin_components):
		exchange_(int(interaction.exchange()), spin_components),
		correlation_(int(interaction.correlation()), spin_components)
	{
	}

  ////////////////////////////////////////////////////////////////////////////////////////////
	
  template <typename SpinDensityType, typename CoreDensityType, typename VKSType>
  void operator()(SpinDensityType const & spin_density, CoreDensityType const & core_density_, VKSType & vks, double & exc, double & nvxc) const {
    
    exc = 0.0;
		nvxc = 0.0;
		
		if(exchange_.true_functional() or correlation_.true_functional()){

			auto full_density = operations::add(spin_density, core_density_);
			
			double efunc = 0.0;
			basis::field_set<basis::real_space, double> vfunc(spin_density.skeleton());

			if(exchange_.true_functional()){
				exchange_(full_density, efunc, vfunc);
				exc += efunc;
				operations::increment(vks, vfunc);
				nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
			}
				
			if(correlation_.true_functional()){
				correlation_(full_density, efunc, vfunc);
				exc += efunc;
				operations::increment(vks, vfunc);
				nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
			}
		}
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
	
	auto & exchange() const {
		return exchange_;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
private:
	hamiltonian::xc_functional exchange_;
	hamiltonian::xc_functional correlation_;
	
};
}
}
#endif

#ifdef INQ_HAMILTONIAN_XC_TERM_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_TERM_UNIT_TEST

#include <ions/unit_cell.hpp>
#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::xc_term", "[xc_term]"){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif
