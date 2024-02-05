/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__UNITS
#define INQ__INTERFACE__UNITS

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
		return "units";
	}
	
	void help() const {
		
		std::cout << R""""(

Units in inq
============

In inq all physical quantities given to the code must also include the
units they are in. There are two main objectives for this. First, it
allows users to use the units directly from the source without having
to waste time converting. Second, it avoids errors when users assume
the code is expecting different units.

In the command line interface units are given as an argument after a
value, or a set of values. There are no default input units in
inq. Not putting the units is a syntax error and the command most
likely fail or produce unexpected behavior.

This page lists the supported units in inq for different
magnitudes. All units are case unsensitive, have abbreviations and can
be given in plural form.

- Energy
  * 'Hartree' or 'Ha'
  * 'Electronvolt' or 'eV': 0.0367493221756544 Ha
  * 'Rydberg' or 'Ry': 0.5 Ha
  * 'Kelvin' or 'K' (multiplied by the Boltzmann constant): 3.16681156345556e-06 Ha

- Length
  * 'Bohr' or 'b'
  * 'Angstrom' or 'A': 1.88972612462938 b
  * 'nanometer' or 'nm': 10 A
  * 'picometer' or 'pm': 0.01 A

- Time
  * 'Atomictime' or 'atu'
  * 'attosecond' or 'as': 0.0413413733352975 atu
  * 'femtosecond' or 'fs': 41.3413733352975 atu
  * 'nanosecond' or 'ns': 1000 fs
  * 'picosecond' or 'ps': 1000 ns


)"""";

		exit(0);
	}
		
} const units ;

}
}
#endif

#ifdef INQ_INTERFACE_UNITS_UNIT_TEST
#undef INQ_INTERFACE_UNITS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
