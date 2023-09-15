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
		
	utils::match energy_match(2.0e-5);

  //SINGLE ATOM SANITY CHECK
  {
    systems::ions ions(systems::cell::cubic(8.0_b).finite());
    ions.insert(input::species("He").nofilter(), {0.0_b, 0.0_b, 0.0_b});
    
    systems::electrons electrons(env.par(), ions, options::electrons{}.extra_states(3).cutoff(30.0_Ha));
    ground_state::initial_guess(ions, electrons);
    
    bomd::propagate(ions, electrons, options::theory{}.pbe(), options::real_time{}.num_steps(10).dt(1.0_fs));
  }
  
}
