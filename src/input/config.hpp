/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__CONFIG
#define INPUT__CONFIG

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <input/electronic_theory.hpp>
#include <cassert>

namespace input {

  struct config {

		config(){
			extra_states = 0;
			excess_charge = 0.0;
			theory = electronic_theory::DENSITY_FUNCTIONAL;
		}

		int extra_states;
		double excess_charge;
		electronic_theory theory;
		
  };
}

////////////////////////////////////////////////////////

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#endif
   
#endif
