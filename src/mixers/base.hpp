/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MIXERS__BASE
#define INQ__MIXERS__BASE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/array.hpp>

namespace inq {
namespace mixers {

template <class FieldType>
class base {
	
public:
	virtual ~base(){};
	virtual void operator()(FieldType & input_value, FieldType const & output_value) = 0;
	
};

}
}
#endif

#ifdef INQ_MIXERS_BASE_UNIT_TEST
#undef INQ_MIXERS_BASE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif
