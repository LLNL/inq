/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__SYSTEMS__IONS
#define INQ__SYSTEMS__IONS

#include <cfloat>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <input/cell.hpp>
#include <mpi3/environment.hpp>

#include<spdlog/spdlog.h>
#include<spdlog/sinks/stdout_color_sinks.h>
#include<spdlog/fmt/ostr.h> // print user defined types

namespace inq {
namespace systems {

class ions {

public:

	ions(const input::cell & arg_cell_input, const inq::ions::geometry & geo_arg = inq::ions::geometry()):
		cell_(arg_cell_input, arg_cell_input.periodic_dimensions()),
		geo_(geo_arg){

		if(boost::mpi3::environment::get_world_instance().root()){
			auto console = spdlog::stdout_color_mt("ions");
			console->info("geometry {}", geo_);
			console->info("cell {}", cell_);
		}
	}

	auto & geo() const {
		return geo_;
	}

	auto & cell() const {
		return cell_;
	}
    
	inq::ions::UnitCell cell_;
	inq::ions::geometry geo_;

};
  
}
}

#ifdef INQ_SYSTEMS_IONS_UNIT_TEST
#undef INQ_SYSTEMS_IONS_UNIT_TEST

#endif
#endif

