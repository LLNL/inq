/* -*- indent-tabs-mode: t -*- */

#include <iostream>
#include <ground_state/calculate.hpp>

int main(){

	boost::mpi3::environment env;
	
	using namespace std;
	using namespace inq;
	using namespace inq::input;
	using namespace inq::math;

	double const distance = 2.0;
	
	std::vector<atom> geo;
	geo.push_back( "N" | vec3d(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | vec3d(0.0, 0.0,  0.5*distance));

	double const hbox = 2.5;

	systems::ions ions(cell::cubic(hbox, hbox, 2*hbox) | cell::periodic(),
										 geo);

	systems::electrons electrons(ions, input::basis::cutoff_energy(30.0));

  auto const result = ground_state::calculate(ions,
																							electrons,
																							input::interaction::dft(),
																							input::scf::conjugate_gradient() | input::scf::mixing(0.1));


	
  cout << "\nCALCULATED ENERGY " << result.energy.total() << "\n\n";

}
