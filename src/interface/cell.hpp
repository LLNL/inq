/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__CELL
#define INQ__INTERFACE__CELL

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace interface {

struct {

	std::string name() const {
		return "cell";
	}

	std::string one_line() const {
		return "Defines the simulation cell.";
	}
	
	void help() const {
		
		std::cout << R""""(

The 'cell' command defines the inq simulation cell. There are four
variations of the command depending on the type of cell you want to
define.

In all cases you must include the units for the cell. For the moment
note you can use 'bohr' (or 'b') and "angstrom" (or 'A').  See `inq
help units` for other supported units.

An optional last argument defines the periodicity of the cell. The
periodicity can be specified as a number from '0d' to '3d', or as
'finite', 'wire', 'slab' and 'periodic'. If the periodicity is not
given, the system is considered periodic.

Note that defining the cell is not necessary if you are going to read
it from a file (using the `inq ion file` command).

The following are the accepted forms of the cell command:

- inq cell cubic <a> <units> [periodicity]

  Defines a cubic cell of side <a>.

  For example 'inq cell 5.0 A finite'.


- inq cell orthorhombic <a> <b> <c> <units> [periodicity]

  Defines a orthorhombic (parallelepipedic) cell of sides a, b and c.

  For example 'inq cell orthorhombic 10.0 10.0 12.0 bohr'.


- inq cell  <a1> <a2> <a3>  <b1> <b2> <b3>  <c1> <c2> <c3>  <units> [periodicity]

  Creates a general cell defined by the lattice vectors a, b, and c
  (given by components). The units are applied to all the vector
  components.

  For example 'inq cell  4.6478 0 0  -2.3239 4.02512 0  0 0 10.0 b 2d'.


- inq cell  <a1> <a2> <a3>  <b1> <b2> <b3>  <c1> <c2> <c3>  scale <s> <units> [periodicity]

  Like the previous case, creates a general cell defined by the lattice vectors a, b, and c
  (given by components). However in this case a general scale factor is applied to all the vectors.

  For example 'inq cell  0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0  scale 3.57 angstrom'.

)"""";

		exit(0);
	}

	
	void operator()() const {
		auto cell = systems::ions::load(".inq/default_ions").cell();
		if(input::environment::global().comm().root()) std::cout << cell;
	}
	
	void cubic(quantity<magnitude::length> const aa, int periodicity = 3) const {
		systems::ions ions(systems::cell::cubic(aa).periodicity(periodicity));
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	void orthorhombic(quantity<magnitude::length> const aa, quantity<magnitude::length> const bb, quantity<magnitude::length> const cc, int periodicity = 3) const {
		systems::ions ions(systems::cell::orthorhombic(aa, bb, cc).periodicity(periodicity));
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

	void operator()(quantity<magnitude::length> const aa0, quantity<magnitude::length> const aa1, quantity<magnitude::length> const aa2,
									quantity<magnitude::length> const bb0, quantity<magnitude::length> const bb1, quantity<magnitude::length> const bb2,
									quantity<magnitude::length> const cc0, quantity<magnitude::length> const cc1, quantity<magnitude::length> const cc2,
									int periodicity = 3) const {
		systems::ions ions(systems::cell::lattice({aa0, aa1, aa2}, {bb0, bb1, bb2}, {cc0, cc1, cc2}).periodicity(periodicity));
		ions.save(input::environment::global().comm(), ".inq/default_ions");
	}

private:

	static auto parse_periodicity(std::string per_string){
		
		if(per_string == "finite" or per_string == "0" or per_string == "0d") {
			return 0;
		} else if (per_string == "wire" or per_string == "1" or per_string == "1d"){
			return 1;
		} else if (per_string == "slab" or per_string == "2" or per_string == "2d"){
			return 2;
		} else if (per_string == "periodic" or per_string == "3" or per_string == "3d"){
			return 3;
		} else {
			throw std::runtime_error("Error: unknown periodicity '" + per_string + "'.");
		}
	}

public:
	
	template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {
		
		if(args.size() == 0) {
			operator()();
			exit(0);
		}
		
		if(args[0] == "cubic"){
			if(args.size() != 3 and args.size() != 4) {
				std::cerr << "Error: Wrong arguments for a cubic cell definition.\nUse: inq cell cubic <lattice_parameter> <units> [periodicity]" << std::endl;
				exit(1);
			}
			
			auto aa = atof(args[1].c_str())*magnitude::length::parse(args[2]);

			int per = 3;
			if(args.size() == 4) per = parse_periodicity(args[3]);
			
			cubic(aa, per);
			if(not quiet) operator()();
			exit(0);
		}
		
		if(args[0] == "orthorhombic"){
			if(args.size() != 5 and args.size() != 6) {
				std::cerr << "Error: Wrong arguments for an orthorhombic cell definition.\nUse: inq cell orthorhombic <a> <b> <c> <units> [periodicity]" << std::endl;
				exit(1);
			}

			auto unit = magnitude::length::parse(args[4]);
			auto aa = atof(args[1].c_str())*unit;
			auto bb = atof(args[2].c_str())*unit;
			auto cc = atof(args[3].c_str())*unit;

			int per = 3;
			if(args.size() == 6) per = parse_periodicity(args[5]);
			
			orthorhombic(aa, bb, cc, per);
			if(not quiet) operator()();
			exit(0);
		}
		
		
		if(args.size() >= 10 and args.size() <= 13) {

			auto scale = 1.0;
			auto units_name = args[9];
			auto per_location = 10;
			
			if(units_name == "scale"){
				scale = atof(args[10].c_str());
				units_name = args[11];
				per_location = 12;
			}

			auto unit = scale*magnitude::length::parse(units_name);
				
			auto aa0 = atof(args[0].c_str())*unit;
			auto aa1 = atof(args[1].c_str())*unit;
			auto aa2 = atof(args[2].c_str())*unit;
			
			auto bb0 = atof(args[3].c_str())*unit;
			auto bb1 = atof(args[4].c_str())*unit;
			auto bb2 = atof(args[5].c_str())*unit;

			auto cc0 = atof(args[6].c_str())*unit;
			auto cc1 = atof(args[7].c_str())*unit;
			auto cc2 = atof(args[8].c_str())*unit;

			int per = 3;
			if(long(args.size()) > per_location) per = parse_periodicity(args[per_location]);
			
			operator()(aa0, aa1, aa2, bb0, bb1, bb2, cc0, cc1, cc2, per);
			
			if(not quiet) operator()();
			exit(0);
		}

		std::cerr << "Error: Invalid syntax in the 'cell' command" << std::endl;
		exit(1);
	}
		
} const cell ;

}
}
#endif

#ifdef INQ_INTERFACE_CELL_UNIT_TEST
#undef INQ_INTERFACE_CELL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
