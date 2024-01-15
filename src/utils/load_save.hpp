/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__LOAD_SAVE
#define INQ__UTILS__LOAD_SAVE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>
#include <array>
#include <cmath>

namespace inq {
namespace utils {

template <typename Type>
void save_value(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	if(not value.has_value()) return;
	
	auto file = std::ofstream(filename);
	file.precision(25);
	
	if(not file) {
		auto exception_happened = true;
		comm.broadcast_value(exception_happened);
		throw std::runtime_error(error_message);
	}
	file << *value << std::endl;
}

template <typename Type>
static void load_value(std::string const & filename, std::optional<Type> & value) {
	auto file = std::ifstream(filename);
	if(file){
		Type readval;
		file >> readval;
		value = readval;
	}
}

}
}
#endif

#ifdef INQ_UTILS_LOAD_SAVE_UNIT_TEST
#undef INQ_UTILS_LOAD_SAVE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {


}
#endif
