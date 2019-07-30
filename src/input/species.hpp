/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef INPUT__SPECIES
#define INPUT__SPECIES

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

#include <pseudo/element.hpp>
#include <vector>
#include <cmath>

namespace input {
  class species : public pseudo::element {

  public:
		
		species(const pseudo::element & arg_el, std::string arg_pseudo_file = ""):
			pseudo::element(arg_el),
			pseudo_file_(arg_pseudo_file){
		}
		
	private:

		std::string pseudo_file_;		
		
  };
  
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class ions::species", "[species]") {
  
  using namespace Catch::literals;

	input::species s(pseudo::element("Xe"));

	REQUIRE(s.atomic_number() == 54);
	
}


#endif


#endif

// Local Variables:
// eval:(setq indent-tabs-mode: t tab-width: 2)
// mode: c++
// coding: utf-8
// End:
