/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SYSTEMS__IONS
#define INQ__SYSTEMS__IONS

#include <cfloat>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <input/cell.hpp>

namespace inq {
namespace systems {

class ions {

public:

	ions(const input::cell & arg_cell_input, const inq::ions::geometry & geo_arg = inq::ions::geometry()):
		cell_(arg_cell_input, arg_cell_input.periodic_dimensions()),
		geo_(geo_arg){
      
		geo_.info(std::cout);
		cell_.info(std::cout);

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

