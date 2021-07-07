/* -*- indent-tabs-mode: t -*- */

//  Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

#ifndef INQ__SYSTEMS__IONS
#define INQ__SYSTEMS__IONS

#include <spglib.h>

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
		geo_(geo_arg)
	{
		char symbol[11];

		std::vector<int> types(geo_.num_atoms());

		for(int iatom = 0; iatom < geo_.num_atoms(); iatom++){
			types[iatom] = geo_.atoms()[iatom].atomic_number();
		}
		
		auto symnum = spg_get_international(symbol, reinterpret_cast<double (*)[3] >(const_cast<double *>(cell_.amat())), reinterpret_cast<double (*)[3] >(const_cast<math::vector3<double> *>(geo_.coordinates().data())),
																				types.data(), geo_.num_atoms(), 1e-4);
		
		std::cout << "symmetry: " << symbol << '\t' << symnum << std::endl;
	}

	auto & geo() const {
		return geo_;
	}
	auto & geo() {
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

