/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__THEORY
#define INQ__INTERFACE__THEORY

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <options/theory.hpp>

namespace inq {
namespace interface {

struct {		

	std::string name() const {
		return "theory";
	}

	std::string one_line() const {
		return "Defines the theory used to represent the electrons-electron interaction.";
	}

	void operator()() const {
		auto theo = options::theory::load(".inq/default_theory");
		std::cout << theo;
	}
	
	void non_interacting() const{
		auto theo = options::theory::load(".inq/default_theory").non_interacting();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}

	void hartree() const{
		auto theo = options::theory::load(".inq/default_theory").hartree();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}

	void hartree_fock() const{
		auto theo = options::theory::load(".inq/default_theory").hartree_fock();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}

	void lda() const{
		auto theo = options::theory::load(".inq/default_theory").lda();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}
	
	void pbe() const{
		auto theo = options::theory::load(".inq/default_theory").pbe();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}

	void rpbe() const{
		auto theo = options::theory::load(".inq/default_theory").rpbe();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}
	
	void pbe0() const{
		auto theo = options::theory::load(".inq/default_theory").pbe0();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}
	
	void b3lyp() const{
		auto theo = options::theory::load(".inq/default_theory").b3lyp();
		theo.save(input::environment::global().comm(), ".inq/default_theory");
	}
	
	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0){
			operator()();
			exit(0);
			
		} else if(args.size() == 1 and args[0] == "non-interacting") {
			non_interacting();
		} else if(args.size() == 1 and args[0] == "hartree") {
			hartree();
		} else if(args.size() == 1 and args[0] == "hartree-fock") {
			hartree_fock();
		} else if(args.size() == 1 and args[0] == "lda") {
			lda();
		} else if(args.size() == 1 and args[0] == "pbe") {
			pbe();
		} else if(args.size() == 1 and args[0] == "rpbe") {
			rpbe();
		} else if(args.size() == 1 and args[0] == "pbe0") {
			pbe0();
		} else if(args.size() == 1 and args[0] == "b3lyp") {
			b3lyp();
		} else {				
			std::cerr << "Error: Invalid syntax in 'theory' command" << std::endl;
			exit(1);
		}

		if(not quiet) operator()();
		exit(0);
	}
	
} const theory;

}
}
#endif

#ifdef INQ_INTERFACE_THEORY_UNIT_TEST
#undef INQ_INTERFACE_THEORY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
