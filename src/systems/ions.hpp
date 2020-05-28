/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__IONS
#define SYSTEMS__IONS

#include <cfloat>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <input/cell.hpp>

namespace inq {
namespace systems {

class ions {

public:

	ions(const input::cell & arg_cell_input, const ::ions::geometry & geo_arg = ::ions::geometry()):
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
    
private:
    
	::ions::UnitCell cell_;
	::ions::geometry geo_;

};
  
}
}
#endif

