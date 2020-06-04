/* -*- indent-tabs-mode: t -*- */

#include <iostream>

#include <ground_state/calculate.hpp>

using namespace std;
using namespace inq;
using namespace inq::input;
using namespace inq::math;

auto nitrogen_energy(double distance){

	std::vector<atom> geo;
	geo.push_back( "N" | vec3d(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | vec3d(0.0, 0.0,  0.5*distance));

	double const hbox = 3.0;

	systems::ions ions(cell::cubic(hbox, hbox, 2*hbox) | cell::periodic(),
										 geo);

    systems::electrons electrons(ions, input::basis::cutoff_energy(30.0));

    auto const result = ground_state::calculate(
        ions,
        electrons,
        input::interaction::dft(),
        input::scf::conjugate_gradient() | input::scf::mixing(0.05)
    );

    cout << "\nDISTANCE " << distance << " CALCULATED ENERGY " << result.energy.total() << "\n\n";

	return result.energy.total();
}

int main(){

	boost::mpi3::environment env;

	double const d_0   = 1.50;
	double const del   = 0.05;
	int    const n_max = 20;

	std::vector<std::pair<double, double>> e_vs_d;

	for(int n = 0; n != n_max; ++n){
		double const distance = d_0 + n*del;
		e_vs_d.emplace_back(distance, nitrogen_energy(distance));
	}
	
	{
		std::ofstream ofs{"nitrogen_e_vs_d.dat"};
		
		ofs << "#nitrogen energy" << std::endl;
		for(int n = 0; n != n_max; ++n){
			ofs<< e_vs_d[n].first <<" "<< e_vs_d[n].second <<'\n';
		}
	}

}


