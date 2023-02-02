/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__BERKELEYGW
#define INQ__OBSERVABLES__BERKELEYGW

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

#include <hdf5/serial/hdf5.h>

namespace inq {
namespace observables {
namespace berkeleygw {

void save(systems::ions const & ions, systems::electrons const & electrons){
	auto h5file = H5Pcreate(H5P_FILE_ACCESS);
	H5Pclose(h5file);
}

}
}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_BERKELEYGW_UNIT_TEST
#undef INQ_OBSERVABLES_BERKELEYGW_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("observables::berkeleygw", "[observables::berkeleygw]") {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	using math::vector3;


}

#endif
#endif
