/* -*- indent-tabs-mode: t -*- */

#include <inq/inq.hpp>

using namespace inq;
using namespace inq::input;
using namespace inq::math;

#include <iostream>
#include <vector>

using namespace std;

int main(){

	boost::mpi3::environment env;
	
	double const distance = 2.0;
	
	vector<atom> geo;
	geo.push_back( "N" | vec3d(0.0, 0.0, -distance/2.0));
	geo.push_back( "N" | vec3d(0.0, 0.0,  distance/2.0));

	cell super = cell::cubic(2.5, 2.5, 5.0) | cell::periodic();

	systems::ions ions(super, geo);

	systems::electrons electrons(ions, input::basis::cutoff_energy(30.0));

  auto result = ground_state::calculate(ions,
																				electrons,
																				interaction::dft(),
																				scf::conjugate_gradient() | scf::mixing(0.1));
	
  cout << "\nCALCULATED ENERGY " << result.energy.total() << "\n\n";

}
