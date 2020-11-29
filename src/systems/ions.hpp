/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__SYSTEMS__IONS
#define INQ__SYSTEMS__IONS

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <input/cell.hpp>

#include <math/array.hpp>

namespace inq {
namespace systems {

class ions {

public:

	ions(const input::cell & arg_cell_input, const inq::ions::geometry & geo_arg = inq::ions::geometry()):
		cell_(arg_cell_input, arg_cell_input.periodic_dimensions()),
		geo_(geo_arg), velocities_(geo_arg.num_atoms(), math::vector3<double>{0., 0., 0.})
	{}

	auto & geo() const {
		return geo_;
	}
	auto & geo() {
		return geo_;
	}
	auto& velocities() const{return velocities_;}
	auto& velocities()      {return velocities_;}
	
	auto & cell() const {
		return cell_;
	}

	inq::ions::UnitCell cell_;
	inq::ions::geometry geo_;
	math::array<math::vector3<double>, 1> velocities_;
};

}
}

#ifdef INQ_SYSTEMS_IONS_UNIT_TEST
#undef INQ_SYSTEMS_IONS_UNIT_TEST

#endif
#endif

