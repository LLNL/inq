/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__SPECTRUM
#define INQ__INTERFACE__SPECTRUM

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/environment.hpp>
#include <ground_state/initial_guess.hpp>
#include <ground_state/calculate.hpp>
#include <real_time/propagate.hpp>
#include <real_time/results.hpp>

namespace inq {
namespace interface {

struct {

	constexpr auto name() const {
		return "spectrum";
	}

	constexpr auto one_line() const {
		return "Calculates different types of spectra based on calculation results.";
	}
	
	constexpr auto help() const {
		return R""""(

The 'spectrum' command
==================

This command is used to calculate spectra of different types using the
results of calculations.

These are the options available:

-  Shell:  `spectrum file <filename>`
   Python: `spectrum.file(filename)`

   Calculates an absorption spectra based on a dipole series. The
   spectrum data is printed to the terminal as two data columns, the
   first one is the frequency/energy and the second one is the
   spectrum for all input columns..

   Shell examples:  `inq spectrum file data.txt`
   Python examples: `pinq.spectrum.file("data.txt")`


)"""";
	}

  static void file(std::string const & filename) {
    std::ifstream inp(filename);

    std::vector<double> values;
    
    auto ncolumns = 0;
    while(true) {
      auto line = std::string();
      getline(inp, line);
      if(inp.eof()) break;
      if(line[0] == '#') continue;

      //      std::cout << line << std::endl;

      if(ncolumns == 0) { //count the columns
        std::istringstream iss(line + " ");
        
        std::string data_item;
        while(iss) {
          iss >> data_item;
          if(iss.eof()) break;
          ncolumns++;
        }
      }
      
      std::istringstream iss(line + " ");
      
      double val;
      for(int icol = 0; icol < ncolumns; icol++) {
        iss >> val;
        values.push_back(val);
      }

    }

    assert(values.size()%ncolumns == 0);

    auto ntime = (long) (values.size())/ncolumns;
    
    gpu::array<double, 1> time(ntime);
    gpu::array<double, 2> time_series({ncolumns - 1, ntime});
    

    for(auto itime = 0l; itime < ntime; itime++){
      time[itime] = values[itime*ncolumns];
      for(auto icol = 1; icol < ncolumns; icol++) time_series[icol - 1][itime] = values[itime*ncolumns + icol];
    }

    using namespace magnitude;
    auto maxw = 20.0_eV;
    auto dw = 0.01_eV;

    long nfreq = maxw/dw + 1;

    gpu::array<complex, 2> freq_series({ncolumns - 1, nfreq});
    
    for(auto icol = 0; icol < ncolumns - 1; icol++){
      freq_series[icol] = observables::spectrum(maxw, dw, time, time_series[icol]);
    }

    for(auto ifreq = 0; ifreq < nfreq; ifreq++){
      std::cout << std::scientific << ifreq*dw;
      for(int icol = 0; icol < ncolumns - 1; icol++) std::cout << '\t' << real(freq_series[icol][ifreq]) << '\t' << imag(freq_series[icol][ifreq]);
      std::cout << '\n';
    }
    
  }

  template <typename ArgsType>
	void command(ArgsType const & args, bool quiet) const {

		if(args.size() == 2 and args[0] == "file") {
      file(args[1]);
			actions::normal_exit();
		}
  }
	
#ifdef INQ_PYTHON_INTERFACE
	template <class PythonModule>
	void python_interface(PythonModule & module) const {
		namespace py = pybind11;
		using namespace pybind11::literals;
 
		auto sub = module.def_submodule(name(), help());
	}
#endif

} const spectrum;

}
}
#endif

#ifdef INQ_INTERFACE_SPECTRUM_UNIT_TEST
#undef INQ_INTERFACE_SPECTRUM_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
