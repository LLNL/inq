/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__INITIAL_GUESS
#define INQ__GROUND_STATE__INITIAL_GUESS

#include <cfloat>

#include <systems/ions.hpp>
#include <operations/randomize.hpp>
#include <operations/orthogonalize.hpp>
#include <density/normalize.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace ground_state {
	
void initial_guess(const systems::ions & ions, systems::electrons & electrons){
  electrons.density_ = electrons.atomic_pot_.atomic_electronic_density(electrons.density_basis_, ions.cell(), ions.geo());

  density::normalize(electrons.density_, electrons.states_.num_electrons());
  
  operations::randomize(electrons.phi_);
  operations::orthogonalize(electrons.phi_);

	math::array<double, 1> eigenvalues(electrons.phi_.local_set_size());
	
	for(long ist = 0; ist < electrons.phi_.local_set_size(); ist++) eigenvalues[ist] = ist + electrons.phi_.set_part().start();
	
	electrons.update_occupations(eigenvalues);
	
}
}
}

#ifdef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST
#undef INQ_GROUND_STATE_INITIAL_GUESS_UNIT_TEST
#endif

#endif

