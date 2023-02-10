/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2020-2021 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__REAL_TIME__VIEWABLES
#define INQ__REAL_TIME__VIEWABLES

#include <systems/ions.hpp>
#include <hamiltonian/calculate_energy.hpp>
#include <hamiltonian/self_consistency.hpp>
#include <hamiltonian/forces.hpp>
#include <operations/overlap_diagonal.hpp>
#include <observables/dipole.hpp>
#include <observables/current.hpp>
#include <perturbations/none.hpp>
#include <ions/propagator.hpp>
#include <systems/electrons.hpp>
#include <real_time/crank_nicolson.hpp>
#include <real_time/etrs.hpp>
#include <utils/profiling.hpp>

#include <chrono>

namespace inq {
namespace real_time {

template <class ForcesType, class HamiltonianType, class Perturbation>
class viewables {
	bool last_iter_;
	int iter_;
	double time_;
	systems::ions const & ions_;
	systems::electrons const & electrons_;
	hamiltonian::energy const & energy_;
	ForcesType const & forces_;
  HamiltonianType const & ham_;
	Perturbation const & pert_;
	
public:

	viewables(bool last_iter, int iter, double time, systems::ions const & ions, systems::electrons const & electrons, hamiltonian::energy const & energy, ForcesType const & forces, HamiltonianType const & ham, Perturbation const & pert)
		:last_iter_(last_iter), iter_(iter), time_(time), ions_(ions), electrons_(electrons), energy_(energy), forces_(forces), ham_(ham), pert_(pert){
	}

	auto iter() const {
		return iter_;
	}

	auto last_iter() const {
		return last_iter_;
	}

	auto every(int every_iter) const {
		if(iter() == 0) return false;
		return (iter()%every_iter == 0) or last_iter(); 
	}

	auto root() const {
		return electrons_.full_comm_.root();
	}
	
	auto time() const {
		return time_;
	}
	
	auto coordinates(int iatom) const {
		return ions_.geo().coordinates()[iatom];
	}
	
	auto velocities(int iatom) const {
		return ions_.geo().velocities()[iatom];
	}

	auto forces(int iatom) const {
		return forces_[iatom];
	}

	auto energy() const {
		return energy_.total();
	}

	auto dipole() const {
		return observables::dipole(ions_, electrons_);
	}

	auto laser_field() const {
		return pert_.uniform_electric_field(time_);
	}

	auto vector_field() const{
		return ham_.uniform_vector_potential();
	}

	auto num_electrons() const {
		return operations::integral(electrons_.density());
	}

  auto current() const {
    return ions_.cell().metric().to_cartesian(observables::current(ions_, electrons_, ham_));
  }
	
};

}
}
#endif

#ifdef INQ_REAL_TIME_VIEWABLES_UNIT_TEST
#undef INQ_REAL_TIME_VIEWABLES_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
