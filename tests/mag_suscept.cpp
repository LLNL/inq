/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>
#include <perturbations/magnetic.hpp>

using namespace inq;
using namespace inq::magnitude;

int main(int argc, char ** argv){
    auto env = input::environment{};

	utils::match energy_match(3.0e-5);

    // metallic Na
    // bcc : a = 0.428 nm

    auto a = 4.28_A;
    auto cell = systems::cell::lattice({-a/2.0, a/2.0, a/2.0}, {a/2.0, -a/2.0, a/2.0}, {a/2.0, a/2.0, -a/2.0});
    auto ions = systems::ions(cell);
    ions.insert_fractional("Na", {0.0, 0.0, 0.0});
    ions.insert_fractional("Na", {0.5, 0.5, 0.5});
	
	systems::electrons electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).spin_non_collinear().extra_states(5), input::kpoints::grid({2, 2, 2}));
	ground_state::initial_guess(ions, electrons);

    auto scf_options = options::ground_state{}.energy_tolerance(1.0e-9_Ha).max_steps(500).broyden_mixing();
    std::vector<double> mz;
    std::vector<double> bz;

    for (auto i=0; i<10; i++) {
        vector3 bvec = {0.0_amu, 0.0_amu, 0.1_amu};
        bvec = bvec * double(i);
        std::cout << bvec[0] << " " << bvec[1] << " " << bvec[2] << std::endl;
        perturbations::magnetic B{bvec};
        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), scf_options, B);
        auto mag = observables::total_magnetization(electrons.spin_density());
        mz.push_back(mag[2]/cell.volume());                   // bohr^-3
        bz.push_back(bvec[2].in_atomic_units()*2.3505e+05);   // tesla
    }

    if (electrons.kpin_states_comm().rank() == 0) {
        std::ofstream outfile("chi_al.txt");
        if (outfile.is_open()) {
            for (auto i=0; i<mz.size(); i++){
                outfile << bz[i] << "          " << mz[i] << std::endl;
            }
            outfile.close();
            std::cout << "File written successfully" << std::endl;
        }
        else {
            std::cerr << "Error opening file" << std::endl;
        }
    }
	
	return energy_match.fail();
}