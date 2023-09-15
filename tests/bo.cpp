/* -*- indent-tabs-mode: t -*- */

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq/inq.hpp>

int main(int argc, char ** argv){

	using namespace inq;
	using namespace inq::magnitude;
	
	input::environment env{};
		
	utils::match match(2.0e-5);

  //SINGLE ATOM SANITY CHECK
  {
    systems::ions ions(systems::cell::cubic(8.0_b).finite());
    ions.insert(input::species("He").nofilter(), {0.0_b, 0.0_b, 0.0_b});
    
    systems::electrons electrons(env.par(), ions, options::electrons{}.extra_states(3).cutoff(30.0_Ha));
    ground_state::initial_guess(ions, electrons);

    std::vector<double> energy;
		auto output = [&energy](auto data){
			energy.push_back(data.energy().total());
		};
    
    bomd::propagate(ions, electrons, output, options::theory{}.pbe(), options::real_time{}.num_steps(10).dt(1.0_fs));
    
    match.check("BOMD: energy step   0", energy[0],   -2.861530896904);
    match.check("BOMD: energy step   1", energy[1],   -2.861530967190);
    match.check("BOMD: energy step   2", energy[2],   -2.861530964221);
    match.check("BOMD: energy step   3", energy[3],   -2.861530914085);
    match.check("BOMD: energy step   4", energy[4],   -2.861530961752);
    match.check("BOMD: energy step   5", energy[5],   -2.861530904608);
    match.check("BOMD: energy step   6", energy[6],   -2.861530844482);
    match.check("BOMD: energy step   7", energy[7],   -2.861530756606);
    match.check("BOMD: energy step   8", energy[8],   -2.861530594045);
    match.check("BOMD: energy step   9", energy[9],   -2.861530393755);
    match.check("BOMD: energy step  10", energy[10],  -2.861530113884);    

  }
  
}
