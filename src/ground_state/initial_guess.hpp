/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__INITIAL_GUESS
#define INQ__GROUND_STATE__INITIAL_GUESS

#include <cfloat>

#include <systems/ions.hpp>
#include <operations/randomize.hpp>
#include <operations/orthogonalize.hpp>
#include <observables/density.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace ground_state {
	
void initial_guess(const systems::ions & ions, systems::electrons & electrons){

	int iphi = 0;
	for(auto & phi : electrons.lot()) {
		operations::randomize(phi, iphi + electrons.lot_part().start());
		operations::orthogonalize(phi);
		for(long ist = 0; ist < phi.local_set_size(); ist++) electrons.eigenvalues()[iphi][ist] = ist + phi.set_part().start() + (iphi + electrons.lot_part().start())/double(electrons.lot_size());

		iphi++;
	}
	
	electrons.update_occupations(electrons.eigenvalues());
	
	if(ions.geo().num_atoms() > 0){
		electrons.spin_density() = electrons.atomic_pot_.atomic_electronic_density(electrons.states_comm_, electrons.density_basis_, ions.geo(), electrons.states());
	} else {
		electrons.spin_density() = observables::density::calculate(electrons);
	}

	assert(fabs(operations::integral_sum(electrons.spin_density())) > 1e-16);
	
  observables::density::normalize(electrons.spin_density(), electrons.states().num_electrons());

}
}
}
#endif

#ifdef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST
#undef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif

