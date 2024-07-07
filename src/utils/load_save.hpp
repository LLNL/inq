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

template <typename CommunicatorType>
void create_directory(CommunicatorType & comm, std::string const & dirname) {
	comm.barrier();
	
	auto exception_happened = true;
	if(comm.root()) {
		
		try { std::filesystem::create_directories(dirname); }
		catch(...) {
			comm.broadcast_value(exception_happened);
			throw std::runtime_error("INQ Error: cannot create directory '" + dirname + "'.");
		}

		exception_happened = false;
		comm.broadcast_value(exception_happened);
	} else {
		comm.broadcast_value(exception_happened);
		if(exception_happened) throw std::runtime_error("INQ Error: cannot create directory '" + dirname + "'.");
	}
		
	comm.barrier();
}

template <typename Type>
void save_value(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {

	comm.barrier();

	auto exception_happened = true;
	if(comm.root()) {
		auto file = std::ofstream(filename);
		
		if(not file) {
			exception_happened = true;
			comm.broadcast_value(exception_happened);
			throw std::runtime_error(error_message);
		}
		
		file.precision(25);
		
		file << value << std::endl;

		exception_happened = false;
		comm.broadcast_value(exception_happened);
	} else {
		comm.broadcast_value(exception_happened);
		if(exception_happened) throw std::runtime_error(error_message);
	}

	comm.barrier();
}

template <typename Type>
void save_optional(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	if(not value.has_value()) {
		comm.barrier();
		if(comm.root()) std::filesystem::remove(filename);
		comm.barrier();
		return;
	}

	save_value(comm, filename, *value, error_message);
}

template <typename Type>
void save_optional_enum(parallel::communicator & comm, std::string const & filename, Type const & value, std::string const & error_message) {
	if(not value.has_value()) {
		comm.barrier();
		if(comm.root()) std::filesystem::remove(filename);
		comm.barrier();
		return;
	}

	save_value(comm, filename, static_cast<int>(*value), error_message);
}

template <typename Type, typename Function>
void save_container(parallel::communicator & comm, std::string const & filename, Type const & container, std::string const & error_message, Function const & func) {
	
	comm.barrier();

	if(container.empty()) {
		if(comm.root()) std::filesystem::remove(filename);
		comm.barrier();
		return;
	}
	
	auto exception_happened = true;
	if(comm.root()) {
		
		auto file = std::ofstream(filename);
		
		if(not file) {
			auto exception_happened = true;
			comm.broadcast_value(exception_happened);
			throw std::runtime_error(error_message);
		}
		
		file.precision(25);
		for(auto const & el : container) file << func(el) << '\n';

 		exception_happened = false;
		comm.broadcast_value(exception_happened);
	} else {
		comm.broadcast_value(exception_happened);
		if(exception_happened) throw std::runtime_error(error_message);
	}

	comm.barrier();
}

template <typename Type>
void save_container(parallel::communicator & comm, std::string const & filename, Type const & container, std::string const & error_message) {
	save_container(comm, filename, container, error_message, [](auto val){return val;});
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
	if(array.size() == 0) return;
	
	auto file = std::ifstream(filename);

	if(not file) throw std::runtime_error(error_message);

	for(int ii = 0; ii < long(array.size()); ii++) file >> array[ii];
	
}

template <typename Type>
static void load_container(std::string const & filename, Type & container){
	auto file = std::ifstream(filename);
	if(not file) return;

	while(true) {
		std::string line;
		std::getline(file, line);
		if(file.eof()) break;

		auto el = typename Type::value_type{};
		std::stringstream ss{line};
		ss >> el;
		container.emplace(std::move(el));
	}
}

template <typename Type>
static void load_vector(std::string const & filename, Type & vec){
	vec.clear();
	
	auto file = std::ifstream(filename);
	if(not file) return;

	while(true) {
		std::string line;
		std::getline(file, line);
		if(file.eof()) break;

		auto el = typename Type::value_type{};
		std::stringstream ss{line};
		ss >> el;
		vec.emplace_back(std::move(el));
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
