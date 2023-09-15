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
    std::vector<systems::ions::positions_type> positions;
		auto output = [&](auto data){
			energy.emplace_back(data.energy().total());
      positions.emplace_back(data.positions());
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

    match.check("BOMD: position step   0", positions[0][0],   {0.000000000000, 0.000000000000, 0.000000000000});
    match.check("BOMD: position step   1", positions[1][0],   {0.000004716990, 0.000004726350, 0.000004741760});
    match.check("BOMD: position step   2", positions[2][0],   {0.000018854944, 0.000018867598, 0.000018898803});
    match.check("BOMD: position step   3", positions[3][0],   {0.000042330832, 0.000042329281, 0.000042377634});
    match.check("BOMD: position step   4", positions[4][0],   {0.000075012237, 0.000074967862, 0.000075035464});
    match.check("BOMD: position step   5", positions[5][0],   {0.000116707921, 0.000116581132, 0.000116670889});
    match.check("BOMD: position step   6", positions[6][0],   {0.000167161049, 0.000166900859, 0.000167016421});
    match.check("BOMD: position step   7", positions[7][0],   {0.000226005295, 0.000225547773, 0.000225693618});
    match.check("BOMD: position step   8", positions[8][0],   {0.000292687099, 0.000291960315, 0.000292141254});
    match.check("BOMD: position step   9", positions[9][0],   {0.000366358531, 0.000365294393, 0.000365514492});
    match.check("BOMD: position step  10", positions[10][0],  {0.000445997470, 0.000444524872, 0.000444788213});

  }
  
}
