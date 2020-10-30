/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__INITIALIZE
#define INQ__GROUND_STATE__INITIALIZE

#include <cfloat>

#include <systems/ions.hpp>
#include <operations/randomize.hpp>
#include <operations/orthogonalize.hpp>
#include <density/normalize.hpp>
#include <systems/electrons.hpp>

namespace inq {
namespace ground_state {
	
void initialize(const systems::ions & ions, systems::electrons & electrons){
  electrons.density_ = electrons.atomic_pot_.atomic_electronic_density(electrons.density_basis_, ions.cell(), ions.geo());

  density::normalize(electrons.density_, electrons.states_.total_charge());
  
  operations::randomize(electrons.phi_);
  operations::orthogonalize(electrons.phi_);

}
}
}

#ifdef INQ_GROUND_STATE_INITIALIZE_UNIT_TEST
#undef INQ_GROUND_STATE_INITIALIZE_UNIT_TEST
#endif

#endif

