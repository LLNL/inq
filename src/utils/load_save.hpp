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
#include <fstream>
#include <filesystem>

namespace inq {
namespace utils {

template <typename Type>
void save_value(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	auto file = std::ofstream(filename);

	if(not file) {
		auto exception_happened = true;
		comm.broadcast_value(exception_happened);
		throw std::runtime_error(error_message);
	}

	file.precision(25);

	file << value << std::endl;
}

template <typename Type>
void save_optional(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	if(not value.has_value()) {
		std::filesystem::remove(filename);
		return;
	}

	save_value(comm, filename, *value, error_message);
}

template <typename Type>
void save_optional_enum(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	if(not value.has_value()) {
		std::filesystem::remove(filename);
		return;
	}

	save_value(comm, filename, static_cast<int>(*value), error_message);

}
template <typename Type>
void save_array(parallel::communicator & comm, std::string const & filename, Type const & array, std::string const & error_message) {
	auto file = std::ofstream(filename);

	if(not file) {
		auto exception_happened = true;
		comm.broadcast_value(exception_happened);
		throw std::runtime_error(error_message);
	}

	file.precision(25);
	for(int ii = 0; ii < long(array.size()); ii++) file << array[ii] << '\n';
	
}

template <typename Type>
static void load_value(std::string const & filename, Type & value, std::string const & error_message){
	auto file = std::ifstream(filename);
	if(not file) throw std::runtime_error(error_message);
	file >> value;
}
	
template <typename Type>
static void load_optional(std::string const & filename, std::optional<Type> & value) {
	auto file = std::ifstream(filename);
	if(file){
		Type readval;
		file >> readval;
		value = readval;
	}
}

template <typename Type>
static void load_optional_enum(std::string const & filename, std::optional<Type> & value) {
	auto file = std::ifstream(filename);
	if(file){
		int readval;
		file >> readval;
		value = static_cast<Type>(readval);
	}
}

template <typename Type>
static void load_array(std::string const & filename, Type & array, std::string const & error_message){
	auto file = std::ifstream(filename);

	if(not file) throw std::runtime_error(error_message);

	for(int ii = 0; ii < long(array.size()); ii++) file >> array[ii];
	
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
