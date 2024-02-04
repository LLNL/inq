/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__ENERGY
#define INQ__INTERFACE__ENERGY

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <hamiltonian/energy.hpp>

namespace inq {
namespace interface {

struct {		

	std::string name() const {
		return "energy";
	}

	std::string one_line() const {
		return "Get information about the energy obtained by a ground-state calculation.";
	}

	void operator()() const {
		auto ener = hamiltonian::energy::load(".inq/default_energy");
		std::cout << ener;
	}

  double total() const{
    return hamiltonian::energy::load(".inq/default_energy").total();
  }
	
  double kinetic() const{
    return hamiltonian::energy::load(".inq/default_energy").kinetic();
  }

  double eigenvalues() const{
    return hamiltonian::energy::load(".inq/default_energy").eigenvalues();
  }

  double external() const{
    return hamiltonian::energy::load(".inq/default_energy").external();
  }
  
  double non_local() const{
    return hamiltonian::energy::load(".inq/default_energy").nonlocal();
  }
  
  double hartree() const{
    return hamiltonian::energy::load(".inq/default_energy").hartree();
  }
  
  double xc() const{
    return hamiltonian::energy::load(".inq/default_energy").xc();
  }

  double nvxc() const{
    return hamiltonian::energy::load(".inq/default_energy").nvxc();
  }

  double exact_exchange() const{
    return hamiltonian::energy::load(".inq/default_energy").exact_exchange();
  }
  
  double ion() const{
    return hamiltonian::energy::load(".inq/default_energy").ion();
  }

	template <typename ArgsType>
	void command(ArgsType args, bool quiet) const {

		if(args.size() == 0) {
			operator()();
			exit(0);
    }			

    if(args.size() == 1 and args[0] == "total"){
      std::cout << total() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "kinetic"){
      std::cout << kinetic() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "eigenvalues"){
      std::cout << eigenvalues() << std::endl;
      exit(0);
    }
    
    if(args.size() == 1 and args[0] == "external"){
      std::cout << external() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "non-local"){
      std::cout << non_local() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "hartree"){
      std::cout << hartree() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "xc"){
      std::cout << xc() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "nvxc"){
      std::cout << nvxc() << std::endl;
      exit(0);
    }

    if(args.size() == 1 and args[0] == "exact-exchange"){
      std::cout << exact_exchange() << std::endl;
      exit(0);
    }
        
    if(args.size() == 1 and args[0] == "ion"){
      std::cout << ion() << std::endl;
      exit(0);
    }
      
		std::cerr << "Error: Invalid syntax in the 'energy' command" << std::endl;
		exit(1);
    
	}
	
} const energy;

}
}
#endif

#ifdef INQ_INTERFACE_ENERGY_UNIT_TEST
#undef INQ_INTERFACE_ENERGY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
