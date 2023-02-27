/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__CALCULATE
#define INQ__GROUND_STATE__CALCULATE

/*
 Copyright (C) 2019-2023 Xavier Andrade, Alfredo A. Correa

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
#include <mixers/pulay.hpp>
#include <mixers/broyden.hpp>
#include <eigensolvers/steepest_descent.hpp>
#include <math/complex.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <ions/interaction.hpp>
#include <input/scf.hpp>
#include <observables/dipole.hpp>
#include <systems/electrons.hpp>
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
	
	for(int iphi = 0; iphi < el.lot_size(); iphi++){
		state_conv += operations::sum(el.occupations()[iphi], normres[iphi], [](auto occ, auto nres){ return fabs(occ*nres); });
	}
	
	el.lot_states_comm_.all_reduce_n(&state_conv, 1);
	state_conv /= el.states().num_electrons();
	
	return state_conv;
}

template <typename NormResType, typename ConsoleType>
void output_eigenvalues(systems::electrons const & el, NormResType const & normres, ConsoleType & console){

	math::array<int, 2> kpoint_index({el.lot_part().local_size(), el.max_local_set_size()});
	math::array<int, 2> spin_index({el.lot_part().local_size(), el.max_local_set_size()});
	math::array<int, 2> state_index({el.lot_part().local_size(), el.max_local_set_size()});
	math::array<double, 2> occs({el.lot_part().local_size(), el.max_local_set_size()});

	auto iphi = 0;
	for(auto & phi : el.lot()){
		auto ik = el.kpoint_index(phi);
		for(int ist = 0; ist < el.max_local_set_size(); ist++){
			kpoint_index[iphi][ist] = ik;
			spin_index[iphi][ist] = phi.spin_index();
			state_index[iphi][ist] = ist;
			occs[iphi][ist] = 0.0;
			if(fabs(el.lot_weights()[iphi]) > 1e-14) occs[iphi][ist] = el.occupations()[iphi][ist]/el.lot_weights()[iphi];
		}
		iphi++;
	}
	
	auto all_kpoint_index = parallel::gather(+kpoint_index.flatted(), el.lot_states_part(), el.lot_states_comm_, 0);
	auto all_spin_index = parallel::gather(+spin_index.flatted(), el.lot_states_part(), el.lot_states_comm_, 0);
	auto all_states_index = parallel::gather(+state_index.flatted(), el.lot_states_part(), el.lot_states_comm_, 0);
	auto all_eigenvalues = parallel::gather(+el.eigenvalues().flatted(), el.lot_states_part(), el.lot_states_comm_, 0);
	auto all_occupations = parallel::gather(+occs.flatted(), el.lot_states_part(), el.lot_states_comm_, 0);
	auto all_normres = parallel::gather(+normres.flatted(), el.lot_states_part(), el.lot_states_comm_, 0);

	if(not console) return;

	auto order = math::array<int, 1>(all_eigenvalues().size());
	
	for(int iorder = 0; iorder < order.size(); iorder++) order[iorder] = iorder;

	std::sort(order.begin(), order.end(), [evs = all_eigenvalues](auto io, auto jo){ return evs[io] < evs[jo]; });
	
	for(int iorder = 0; iorder < order.size(); iorder++) {
		auto ieig = order[iorder];
		console->info(" kp = {:4d}  sp = {:2d}  st = {:4d}  occ = {:4.3f}  evalue = {:18.12f}  res = {:5.0e}",
									all_kpoint_index[ieig], all_spin_index[ieig], all_states_index[ieig], all_occupations[ieig], real(all_eigenvalues[ieig]), real(all_normres[ieig]));
	}

}

ground_state::result calculate(const systems::ions & ions, systems::electrons & electrons, const input::interaction & inter = {}, const input::scf & solver = {}){

	CALI_CXX_MARK_FUNCTION;

	assert(electrons.lot()[0].full_comm() == electrons.states_basis_comm_);
	
	auto console = electrons.logger();
	if(console) console->trace("calculate started");
	hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_, electrons.states().num_density_components());
	
	hamiltonian::ks_hamiltonian<double> ham(electrons.states_basis_, electrons.states(), electrons.atomic_pot_, inter.fourier_pseudo_value(), ions.geo(), electrons.states().num_states(), sc.exx_coefficient(), /* use_ace = */ true);
	
	if(electrons.full_comm_.root()) ham.info(std::cout);
		
	ground_state::result res;
		
	operations::preconditioner prec;

	using mix_arr_type = std::remove_reference_t<decltype(electrons.spin_density().matrix().flatted())>;
	
	auto mixer = [&]()->std::unique_ptr<mixers::base<mix_arr_type>>{
		switch(solver.mixing_algorithm()){
		case input::scf::mixing_algo::LINEAR : return std::make_unique<mixers::linear <mix_arr_type>>(solver.mixing());
		case input::scf::mixing_algo::PULAY  : return std::make_unique<mixers::pulay  <mix_arr_type>>(4, solver.mixing(), electrons.spin_density().matrix().flatted().size(), electrons.density_basis_.comm());
		case input::scf::mixing_algo::BROYDEN: return std::make_unique<mixers::broyden<mix_arr_type>>(4, solver.mixing(), electrons.spin_density().matrix().flatted().size(), electrons.density_basis_.comm());
		} __builtin_unreachable();
	}();
	
	auto old_energy = std::numeric_limits<double>::max();
		
	sc.update_ionic_fields(electrons.states_comm_, ions, electrons.atomic_pot_);
	sc.update_hamiltonian(ham, res.energy, electrons.spin_density());
		
	res.energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);

	double old_exe = ham.exchange.update(electrons);
	double exe_diff = fabs(old_exe);
	auto update_hf = false;

	electrons.full_comm_.barrier();
	auto iter_start_time = std::chrono::high_resolution_clock::now();

	int conv_count = 0;
	for(int iiter = 0; iiter < solver.scf_steps(); iiter++){

		CALI_CXX_MARK_SCOPE("scf_iteration");

		if(solver.subspace_diag()) {
			int ilot = 0;
			for(auto & phi : electrons.lot()) {
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
		
		for(auto & phi : electrons.lot()) {
			auto fphi = operations::space::to_fourier(std::move(phi));
				
			switch(solver.eigensolver()){
					
			case input::scf::scf_eigensolver::STEEPEST_DESCENT:
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
			auto energy_diff = (res.energy.eigenvalues - old_energy)/electrons.states().num_electrons();

			electrons.full_comm_.barrier();
			std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - iter_start_time;

			electrons.full_comm_.barrier();
			iter_start_time = std::chrono::high_resolution_clock::now();

			auto state_conv = state_convergence(electrons, normres);

			if(solver.verbose_output() and console){
				console->info("SCF iter {} : wtime = {:5.2f}s e = {:.10f} de = {:5.0e} dexe = {:5.0e} dn = {:5.0e} dst = {:5.0e}", 
											iiter, elapsed_seconds.count(), res.energy.total(), energy_diff, exe_diff, density_diff, state_conv);
			}

			if(solver.verbose_output()) output_eigenvalues(electrons, normres, console);
				
			if(fabs(energy_diff) < solver.energy_tolerance()){
				conv_count++;
				if(conv_count > 2 and exe_diff < solver.energy_tolerance()) break;
				if(conv_count > 2) update_hf = true;
			} else {
				conv_count = 0; 
			}

			old_energy = res.energy.eigenvalues;
		}
	}

	//make sure we have a density consistent with phi
	electrons.spin_density() = observables::density::calculate(electrons);
	sc.update_hamiltonian(ham, res.energy, electrons.spin_density());
	res.energy.calculate(ham, electrons);
	
	if(solver.calc_forces()) res.forces = hamiltonian::calculate_forces(ions, electrons, ham);

	if(solver.verbose_output() and console) console->info("SCF iters ended with result energies {}", res.energy);

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

