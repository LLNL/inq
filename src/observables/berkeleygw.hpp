/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__BERKELEYGW
#define INQ__OBSERVABLES__BERKELEYGW

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/field.hpp>
#include <systems/ions.hpp>
#include <systems/electrons.hpp>
#include <physics/constants.hpp>

#include <hdf5.h>

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

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	

}

#endif
#endif
