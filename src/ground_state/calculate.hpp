/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__GROUND_STATE__CALCULATE
#define INQ__GROUND_STATE__CALCULATE

#include <cfloat>

#include <systems/ions.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/calculate_energy.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/energy.hpp>
#include <hamiltonian/forces.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/overlap.hpp>
#include <operations/scal.hpp>
#include <operations/orthogonalize.hpp>
#include <operations/preconditioner.hpp>
#include <operations/integral.hpp>
#include <density/calculate.hpp>
#include <density/normalize.hpp>
#include <mixers/linear.hpp>
#include <mixers/pulay.hpp>
#include <mixers/broyden.hpp>
#include <eigensolvers/conjugate_gradient.hpp>
#include <eigensolvers/steepest_descent.hpp>
#include <eigensolvers/davidson.hpp>
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

ground_state::result calculate(const systems::ions & ions, systems::electrons & electrons, const input::interaction & inter = {}, const input::scf & solver = {}){

	CALI_CXX_MARK_FUNCTION;

	auto console = electrons.logger();
	if(console) console->trace("calculate started");
		
	hamiltonian::ks_hamiltonian<basis::real_space> ham(electrons.states_basis_, ions.cell(), electrons.atomic_pot_, inter.fourier_pseudo_value(), ions.geo(),
																										 electrons.states_.num_states(), inter.exchange_coefficient(), electrons.states_basis_comm_, /* use_ace = */ true);
		
	if(electrons.full_comm_.root()) ham.info(std::cout);
		
	hamiltonian::self_consistency sc(inter, electrons.states_basis_, electrons.density_basis_);

	ground_state::result res;
		
	operations::preconditioner prec;
		
	auto mixer = [&]()->std::unique_ptr<mixers::base<double>>{
		switch(solver.mixing_algorithm()){
		case input::scf::mixing_algo::LINEAR : return std::make_unique<mixers::linear <double>>(solver.mixing());
		case input::scf::mixing_algo::PULAY  : return std::make_unique<mixers::pulay  <double>>(4, solver.mixing(), electrons.states_basis_.part().local_size(), electrons.density_basis_.comm());
		case input::scf::mixing_algo::BROYDEN: return std::make_unique<mixers::broyden<double>>(4, solver.mixing(), electrons.states_basis_.part().local_size(), electrons.density_basis_.comm());
		} __builtin_unreachable();
	}();
		
	auto old_energy = std::numeric_limits<double>::max();
		
	sc.update_ionic_fields(ions, electrons.atomic_pot_);
		
	ham.scalar_potential = sc.ks_potential(electrons.density_, res.energy);
		
	res.energy.ion = inq::ions::interaction_energy(ions.cell(), ions.geo(), electrons.atomic_pot_);

	
	double old_exe = ham.exchange.update(electrons);
	double exe_diff = fabs(old_exe);
	auto update_hf = false;

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
			auto fphi = operations::space::to_fourier(phi);
				
			switch(solver.eigensolver()){
					
			case input::scf::scf_eigensolver::STEEPEST_DESCENT:
				eigensolvers::steepest_descent(ham, prec, fphi);
				break;
					
			case input::scf::scf_eigensolver::CONJUGATE_GRADIENT:
				eigensolvers::conjugate_gradient(ham, prec, fphi);
				break;
					
			case input::scf::scf_eigensolver::DAVIDSON:
				eigensolvers::davidson(ham, prec, fphi);
				break;
					
			default:
				assert(false);
			}

			//This fails, I don't know why. XA
			// phi = operations::space::to_real(fphi);			

			phi.fields() = operations::space::to_real(fphi.fields());
		}

		CALI_MARK_BEGIN("mixing");

		double density_diff = 0.0;
		{
			auto new_density = density::calculate(electrons);
			density_diff = operations::integral_absdiff(electrons.density_, new_density);				
			density_diff /= electrons.states_.num_electrons();
				
			if(inter.self_consistent() and solver.mix_density()) {
				mixer->operator()(electrons.density_.linear(), new_density.linear());
				density::normalize(electrons.density_, electrons.states_.num_electrons());
			} else {
				electrons.density_ = std::move(new_density);
			}
		}
			
		auto vks = sc.ks_potential(electrons.density_, res.energy);
			
		if(inter.self_consistent() and solver.mix_potential()) {
			mixer->operator()(ham.scalar_potential.linear(), vks.linear());
		} else {
			ham.scalar_potential = std::move(vks);
		}

		CALI_MARK_END("mixing");

		{
			CALI_CXX_MARK_SCOPE("energy_calculation");

			auto ecalc = hamiltonian::calculate_energy(ham, electrons);
			
			res.energy.eigenvalues = ecalc.sum_eigenvalues_;
			res.energy.nonlocal = ecalc.nonlocal_;
			res.energy.hf_exchange = ecalc.hf_exchange_;

			auto energy_diff = (res.energy.eigenvalues - old_energy)/ions.geo().num_atoms();
				
			if(solver.verbose_output() and console){
				console->info("SCF iter {} : e = {:.10f} de = {:5.0e} dexe = {:5.0e} dn = {:5.0e} dst = {:5.0e}", 
											iiter, res.energy.total(), energy_diff, exe_diff, density_diff, ecalc.state_conv_);
					
				for(int istate = 0; istate < electrons.states_.num_states(); istate++){
					for(int ilot = 0; ilot < electrons.lot_size(); ilot++){
						console->info("	k-point {:4d} state {:4d}  occ = {:4.3f}  evalue = {:18.12f}  res = {:5.0e}",
													ilot + 1, istate + 1, electrons.occupations()[ilot][istate]/electrons.lot_weights()[ilot], real(ecalc.eigenvalues_[ilot][istate]), real(ecalc.normres_[ilot][istate])
													);
					}
				}
			}
				
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
	electrons.density_ = density::calculate(electrons);

	if(solver.calc_forces()) res.forces = hamiltonian::calculate_forces(ions, electrons, ham);

	if(solver.verbose_output() and console) console->info("SCF iters ended with result energies {}", res.energy);

	if(ions.cell().periodic_dimensions() == 0){
		res.dipole = observables::dipole(ions, electrons);
	} else {
		res.dipole = math::vector3<double>(0.);
	}


	if(console) console->trace("calculate ended normally");
	return res;
}
}
}

#ifdef INQ_GROUND_STATE_CALCULATE_UNIT_TEST
#undef INQ_GROUND_STATE_CALCULATE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("ground_state::calculate", "[ground_state::calculate]") {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}

#endif
#endif

