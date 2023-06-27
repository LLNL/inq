/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__CALCULATE
#define INQ__GROUND_STATE__CALCULATE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/energy.hpp>
#include <hamiltonian/forces.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/preconditioner.hpp>
#include <operations/integral.hpp>
#include <observables/density.hpp>
#include <parallel/gather.hpp>
#include <mixers/linear.hpp>
#include <mixers/broyden.hpp>
#include <eigensolvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <ions/interaction.hpp>
#include <observables/dipole.hpp>
#include <options/ground_state.hpp>
#include <systems/electrons.hpp>
#include <ground_state/eigenvalue_output.hpp>
#include <ground_state/result.hpp>
#include <ground_state/subspace_diagonalization.hpp>

#include<tinyformat/tinyformat.h>

#include<spdlog/spdlog.h>
#include<spdlog/sinks/stdout_color_sinks.h>

#include<memory>

#include <utils/profiling.hpp>

namespace inq {
namespace ground_state {

template <typename NormResType>
auto state_convergence(systems::electrons & el, NormResType const & normres) {
	auto state_conv = 0.0;
	
	for(int iphi = 0; iphi < el.kpin_size(); iphi++){
		state_conv += operations::sum(el.occupations()[iphi], normres[iphi], [](auto occ, auto nres){ return fabs(occ*nres); });
	}
	
	el.kpin_states_comm().all_reduce_n(&state_conv, 1);
	state_conv /= el.states().num_electrons();
	
	return state_conv;
}

ground_state::result calculate(const systems::ions & ions, systems::electrons & electrons, const options::theory & inter = {}, options::ground_state const & solver = {}){

	CALI_CXX_MARK_FUNCTION;

	assert(electrons.kpin()[0].full_comm() == electrons.states_basis_comm());
	
	auto console = electrons.logger();
	if(console) console->trace("calculate started");
	hamiltonian::self_consistency sc(inter, electrons.states_basis(), electrons.density_basis(), electrons.states().num_density_components());
	
	hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis(), electrons.brillouin_zone(), electrons.states(), electrons.atomic_pot(), inter.fourier_pseudo_value(), ions, electrons.states().num_states(), sc.exx_coefficient(), /* use_ace = */ true);
	
	if(electrons.full_comm().root()) ham.info(std::cout);
		
	ground_state::result res;
		
	operations::preconditioner prec;

	using mix_arr_type = std::remove_reference_t<decltype(electrons.spin_density().matrix().flatted())>;
	
	auto mixer = [&]()->std::unique_ptr<mixers::base<mix_arr_type>>{
		switch(solver.mixing_algorithm()){
		case options::ground_state::mixing_algo::LINEAR : return std::make_unique<mixers::linear <mix_arr_type>>(solver.mixing());
		case options::ground_state::mixing_algo::BROYDEN: return std::make_unique<mixers::broyden<mix_arr_type>>(4, solver.mixing(), electrons.spin_density().matrix().flatted().size(), electrons.density_basis().comm());
		} __builtin_unreachable();
	}();
	
	auto old_energy = std::numeric_limits<double>::max();
		
	sc.update_ionic_fields(electrons.states_comm(), ions, electrons.atomic_pot());
	sc.update_hamiltonian(ham, res.energy, electrons.spin_density());
		
	res.energy.ion(inq::ions::interaction_energy(ions.cell(), ions, electrons.atomic_pot()));

	double old_exe = ham.exchange.update(electrons);
	double exe_diff = fabs(old_exe);
	auto update_hf = false;

	electrons.full_comm().barrier();
	auto iter_start_time = std::chrono::high_resolution_clock::now();

	int conv_count = 0;
	for(int iiter = 0; iiter < solver.scf_steps(); iiter++){

		CALI_CXX_MARK_SCOPE("scf_iteration");

		if(solver.subspace_diag()) {
			int ilot = 0;
			for(auto & phi : electrons.kpin()) {
				electrons.eigenvalues()[ilot] = subspace_diagonalization(ham, phi);
				ilot++;
			}
			electrons.update_occupations(electrons.eigenvalues());
		}

		if(update_hf){
			auto exe = ham.exchange.update(electrons);
			exe_diff = fabs(exe - old_exe);
			old_exe = exe;
		}
		
		for(auto & phi : electrons.kpin()) {
			auto fphi = operations::space::to_fourier(std::move(phi));
				
			switch(solver.eigensolver()){
					
			case options::ground_state::scf_eigensolver::STEEPEST_DESCENT:
				eigensolvers::steepest_descent(ham, prec, fphi);
				break;
				
			default:
				assert(false);
			}

			phi = operations::space::to_real(std::move(fphi));
		}

		CALI_MARK_BEGIN("mixing");

		double density_diff = 0.0;
		{
			auto new_density = observables::density::calculate(electrons);
			density_diff = operations::integral_sum_absdiff(electrons.spin_density(), new_density);
			density_diff /= electrons.states().num_electrons();
				
			if(inter.self_consistent()) {
				auto tmp = +electrons.spin_density().matrix().flatted();
				mixer->operator()(tmp, new_density.matrix().flatted());
				electrons.spin_density().matrix().flatted() = tmp;
				observables::density::normalize(electrons.spin_density(), electrons.states().num_electrons());
			} else {
				electrons.spin_density() = std::move(new_density);
			}
		}
		
		sc.update_hamiltonian(ham, res.energy, electrons.spin_density());
		
		CALI_MARK_END("mixing");

		{
			auto normres = res.energy.calculate(ham, electrons);
			auto energy_diff = (res.energy.eigenvalues() - old_energy)/electrons.states().num_electrons();

			electrons.full_comm().barrier();
			std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - iter_start_time;

			electrons.full_comm().barrier();
			iter_start_time = std::chrono::high_resolution_clock::now();

			auto state_conv = state_convergence(electrons, normres);
			auto ev_out = eigenvalues_output(electrons, normres);
			
			if(solver.verbose_output() and console){
				console->info("\nSCF iter {} : wtime = {:5.2f}s e = {:.10f} de = {:5.0e} dexe = {:5.0e} dn = {:5.0e} dst = {:5.0e}\n{}", 
											iiter, elapsed_seconds.count(), res.energy.total(), energy_diff, exe_diff, density_diff, state_conv, ev_out);
			}
			
			if(fabs(energy_diff) < solver.energy_tolerance()){
				conv_count++;
				if(conv_count > 2 and exe_diff < solver.energy_tolerance()) break;
				if(conv_count > 2) update_hf = true;
			} else {
				conv_count = 0; 
			}

			old_energy = res.energy.eigenvalues();
		}
	}

	//make sure we have a density consistent with phi
	electrons.spin_density() = observables::density::calculate(electrons);
	sc.update_hamiltonian(ham, res.energy, electrons.spin_density());
	auto normres = res.energy.calculate(ham, electrons);
			
	if(solver.calc_forces()) res.forces = hamiltonian::calculate_forces(ions, electrons, ham);

	auto ev_out = eigenvalues_output(electrons, normres);		
	
	if(solver.verbose_output() and console) {
		console->info("\nSCF iters ended with resulting eigenvalues and energies:\n\n{}{}", ev_out.full(), res.energy);
	}

	if(ions.cell().periodicity() == 0){
		res.dipole = observables::dipole(ions, electrons);
	} else {
		res.dipole = vector3<double>(0.);
	}
	
	if(console) console->trace("calculate ended normally");
	return res;
}
}
}
#endif

#ifdef INQ_GROUND_STATE_CALCULATE_UNIT_TEST
#undef INQ_GROUND_STATE_CALCULATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif

