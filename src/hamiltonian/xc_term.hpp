/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_TERM
#define INQ__HAMILTONIAN__XC_TERM

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

	template <typename SpinDensityType, typename CoreDensityType>
	SpinDensityType process_density(SpinDensityType const & spin_density, CoreDensityType const & core_density) const{

		SpinDensityType full_density(spin_density.basis(), std::min(2, spin_density.set_size()));

		if(spin_density.set_size() == 4) {
			//for spinors convert the density to 2 components
			gpu::run(spin_density.basis().local_size(),
							 [spi = begin(spin_density.matrix()), ful = begin(full_density.matrix()), cor = begin(core_density.linear())] GPU_LAMBDA (auto ip){
								 auto dtot = spi[ip][0] + spi[ip][1];
								 auto dd = spi[ip][0] - spi[ip][1];
								 auto dpol = sqrt(dd*dd + 4.0*(spi[ip][2]*spi[ip][2] + spi[ip][3]*spi[ip][3]));
								 ful[ip][0] = 0.5*(dtot + dpol);
								 ful[ip][1] = 0.5*(dtot - dpol);
								 for(int ispin = 0; ispin < 2; ispin++){
									 if(ful[ip][ispin] < 0.0) ful[ip][ispin] = 0.0;
									 ful[ip][ispin] += 0.5*cor[ip];
								 }
							 });
		} else {
			assert(spin_density.set_size() == 1 or spin_density.set_size() == 2);
			
			gpu::run(spin_density.basis().local_size(),
							 [spi = begin(spin_density.matrix()), ful = begin(full_density.matrix()), cor = begin(core_density.linear()), nspin = spin_density.set_size()] GPU_LAMBDA (auto ip){
								 for(int ispin = 0; ispin < nspin; ispin++){
									 ful[ip][ispin] = spi[ip][ispin];
									 if(ful[ip][ispin] < 0.0) ful[ip][ispin] = 0.0;
									 ful[ip][ispin] += cor[ip]/nspin;
								 }
							 });
		}
		
		return full_density;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
  template <typename SpinDensityType, typename CoreDensityType, typename VKSType>
  void operator()(SpinDensityType const & spin_density, CoreDensityType const & core_density, VKSType & vks, double & exc, double & nvxc) const {
    
    exc = 0.0;
		nvxc = 0.0;
		if(not exchange_.true_functional() and not correlation_.true_functional()) return;
		
		auto full_density = process_density(spin_density, core_density);
		
		double efunc = 0.0;
		basis::field_set<basis::real_space, double> vfunc(spin_density.skeleton());

		auto density_gradient = std::optional<decltype(operations::gradient(full_density))>{};
		if(exchange_.requires_gradient() or correlation_.requires_gradient()){
			density_gradient.emplace(operations::gradient(full_density));
		}
			
		if(exchange_.true_functional()){
			exchange_(full_density, density_gradient, efunc, vfunc);
			exc += efunc;
			operations::increment(vks, vfunc);
			nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
		}
		
		if(correlation_.true_functional()){
			correlation_(full_density, density_gradient, efunc, vfunc);
			exc += efunc;
			operations::increment(vks, vfunc);
			nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
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

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace Catch::literals;
	
}
#endif
