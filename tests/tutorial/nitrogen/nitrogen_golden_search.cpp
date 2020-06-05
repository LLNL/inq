/* -*- indent-tabs-mode: t -*- */

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa, Alexey Kartsev

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
...
*/

#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <config/path.hpp>
#include <input/atom.hpp>
#include <utils/match.hpp>
#include <ground_state/calculate.hpp>

#include <iostream>
#include<fstream>
#include<iomanip>

namespace mpi3 = boost::mpi3;
using std::cout;

double nitrogen_energy(double distance){

	cout<<"calculating ground state for distance "<< distance <<'\n';

	std::vector<input::atom> geo;
	geo.push_back( "N" | math::vec3d(0.0, 0.0, -0.5*distance));
	geo.push_back( "N" | math::vec3d(0.0, 0.0,  0.5*distance));

	double const hbox = 3.;
	systems::ions ions{
		input::cell::cubic(hbox, hbox, 2*hbox) | input::cell::periodic(),
		geo
	};

	input::config conf;
	conf.extra_states = 4;

	systems::electrons electrons{
		ions, 
		input::basis::cutoff_energy(/*ecut*/30.) | input::basis::density_factor(2), 
		conf
	};

	auto const energies = ground_state::calculate(
		electrons, 
		input::interaction::dft(), 
		input::scf::conjugate_gradient() | input::scf::mixing(/*beta*/0.1)
	);
	cout<<"calculated energy "<< energies.total() <<'\n';
	return energies.total();

}

std::ofstream ofs{"nitrogen_golden_search.dat"};

template<class F>
std::pair<double, double> golden_search(F f, double a, double m, double b, double tol = 1e-3){
	std::pair<double, double> fxmin; // return variable
	ofs << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	static constexpr double R = 0.61803399, C = 1. - R; // golden ratios
	double x1 = m, x2 = m;
	double f1 = f(x1);
	double f2 = f1;
	while(b - a > tol*(std::abs(x1)+std::abs(x2))){
		fxmin = f1<f2?std::make_pair(f1, x1):std::make_pair(f2, x2);
		ofs<< fxmin.second <<" "<< fxmin.first <<std::endl;
		if(f2 < f1){
			a = x1; 
			x1 = x2; x2 = R*x2 + C*b;
			f1 = f2; f2 = f(x2);
		}else{
			b = x2;
			x2 = x1; x1 = R*x1 + C*a;
			f2 = f1; f1 = f(x1);
		}
	}
	return fxmin;
}

int main(){

	mpi3::environment env;

	auto min = golden_search(nitrogen_energy, 1.50, 2.25, 3.50, 1e-3);

	auto x = min.second;
	auto f = min.first;

	std::cout<< "min found at x = "<< x << " with value f = "<< f <<'\n';

}

