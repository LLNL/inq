/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__IO
#define OPERATIONS__IO

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <inq_config.h>

#include <gpu/array.hpp>

#include <parallel/communicator.hpp>
#include <mpi3/detail/datatype.hpp>

#include <utils/num_str.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>

#include <string>
#include <cstdio>
#include <iostream>

#include <filesystem>

namespace inq {
namespace operations {
namespace io {

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class ArrayType, class PartType, class CommType>
void save_array(std::string const & filename, CommType & comm, PartType const & part, ArrayType const & array){
	CALI_CXX_MARK_FUNCTION;

	MPI_File fh;

	assert(array.num_elements() == part.local_size());

	using Type = typename ArrayType::element_type;
	auto mpi_type = boost::mpi3::detail::basic_datatype<Type>();

	auto mpi_err = MPI_File_open(comm.get(), filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	
	if(mpi_err != MPI_SUCCESS){
		std::cerr << "Error: cannot create restart file '" << filename << "'." << std::endl;
		exit(1);
	}

	MPI_Status status;
	mpi_err = MPI_File_write_at(fh, sizeof(Type)*part.start(), raw_pointer_cast(array.data_elements()), part.local_size(), mpi_type, &status);
	
	if(mpi_err != MPI_SUCCESS){
		std::cerr << "Error: cannot write restart file '" << filename << "'." << std::endl;
		exit(1);
	}
	
	int data_written;
	MPI_Get_count(&status, mpi_type, &data_written);
	assert(data_written == long(part.local_size()));
	
	MPI_File_close(&fh);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class ArrayType, class PartType, class CommType>
auto load_array(std::string const & filename, CommType & comm, PartType const & part, ArrayType & array){
	CALI_CXX_MARK_FUNCTION;

	using Type = typename ArrayType::element_type;
	auto mpi_type = boost::mpi3::detail::basic_datatype<Type>();
	
	assert(array.num_elements() == part.local_size());
	
	MPI_File fh;
	auto mpi_err = MPI_File_open(comm.get(), filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	
	if(mpi_err != MPI_SUCCESS){
		return false;
	}
	
	MPI_Status status;
	mpi_err = MPI_File_read_at(fh, sizeof(Type)*part.start(), raw_pointer_cast(array.data_elements()), part.local_size(), mpi_type, &status);
	
	if(mpi_err != MPI_SUCCESS){
		MPI_File_close(&fh);
		return false;
	}
	
	int data_read;
	MPI_Get_count(&status, mpi_type, &data_read);
	assert(data_read == long(part.local_size()));
	
	MPI_File_close(&fh);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class ArrayType, class PartType, class CommType>
void save(std::string const & dirname, CommType & comm, PartType const & part, ArrayType const & array){
	CALI_CXX_MARK_SCOPE("save(array)");
	utils::create_directory(comm, dirname);
	save_array(dirname + "/array", comm, part, array);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class ArrayType, class PartType, class CommType>
auto load(std::string const & dirname, CommType & comm, PartType const & part, ArrayType & array){
	CALI_CXX_MARK_SCOPE("load(array)");
	return load_array(dirname + "/array", comm, part, array);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
void save(std::string const & dirname, basis::field_set<Basis, Type> const & phi){

	CALI_CXX_MARK_SCOPE("save(field_set)");
	
	gpu::array<Type, 1> buffer(phi.basis().part().local_size());

	utils::create_directory(phi.full_comm(), dirname);
	
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){
		auto filename = dirname + "/" + utils::num_to_str(ist + phi.set_part().start());
		buffer = +phi.matrix().rotated()[ist];
		save_array(filename, phi.basis().comm(), phi.basis().part(), buffer);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
auto load(std::string const & dirname, basis::field_set<Basis, Type> & phi){

	CALI_CXX_MARK_FUNCTION;

	gpu::array<Type, 1> buffer(phi.basis().part().local_size());

	DIR* dir = opendir(dirname.c_str());
	if (!dir) {
		return false;
	}
	closedir(dir);
			
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){
		auto filename = dirname + "/" + utils::num_to_str(ist + phi.set_part().start());
		auto success = load_array(filename, phi.basis().comm(), phi.basis().part(), buffer);
		if(not success) return false;
		phi.matrix().rotated()[ist] = buffer;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
void save(std::string const & dirname, states::orbital_set<Basis, Type> const & phi){

	CALI_CXX_MARK_SCOPE("save(field_set)");
	
	gpu::array<Type, 1> buffer(phi.basis().part().local_size());

	utils::create_directory(phi.full_comm(), dirname);
	
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){
		auto filename = dirname + "/" + utils::num_to_str(ist + phi.set_part().start());
		buffer = +phi.matrix().rotated()[ist];
		save_array(filename, phi.basis().comm(), phi.basis().part(), buffer);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class Basis, class Type>
auto load(std::string const & dirname, states::orbital_set<Basis, Type> & phi){

	CALI_CXX_MARK_FUNCTION;

	gpu::array<Type, 1> buffer(phi.basis().part().local_size());

	DIR* dir = opendir(dirname.c_str());
	if (!dir) {
		return false;
	}
	closedir(dir);
			
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){
		auto filename = dirname + "/" + utils::num_to_str(ist + phi.set_part().start());
		auto success = load_array(filename, phi.basis().comm(), phi.basis().part(), buffer);
		if(not success) return false;
		phi.matrix().rotated()[ist] = buffer;
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
		
}
}
}
#endif

#ifdef INQ_OPERATIONS_IO_UNIT_TEST
#undef INQ_OPERATIONS_IO_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	
	using namespace inq;
	using namespace Catch::literals;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	SECTION("array"){

		int const size = 12345;

		parallel::partition part(size, comm);

		gpu::array<int, 1> arr(part.local_size());

		for(int ii = 0; ii < part.local_size(); ii++){
			arr[ii] = part.local_to_global(ii).value();
		}
		
		operations::io::save("array_restart", comm, part, arr);
		
		gpu::array<int, 1> arr2(part.local_size());

		CHECK(operations::io::load("array_restart", comm, part, arr2));

		for(int ii = 0; ii < part.local_size(); ii++){
			CHECK(arr[ii] == arr2[ii]);
		}

		CHECK(not operations::io::load("directory_that_doesnt_exist", comm, part, arr2));
	}

	SECTION("field_set"){
		
		const int npoint = 100;
		const int nvec = 12;
		
		parallel::cartesian_communicator<2> cart_comm(comm, {});
		
		auto basis_comm = basis::basis_subcomm(cart_comm);
		
		basis::trivial bas(npoint, basis_comm);
		
		basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
		basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				auto jjg = aa.set_part().local_to_global(jj);
				auto iig = bas.part().local_to_global(ii);
				aa.matrix()[ii][jj] = 20.0*(iig.value() + 1)*sqrt(jjg.value());
			}
		}
		
		operations::io::save("restart/", aa);
		
		CHECK(operations::io::load("restart/", bb));
		
		for(int ii = 0; ii < bas.part().local_size(); ii++){
			for(int jj = 0; jj < aa.set_part().local_size(); jj++){
				CHECK(aa.matrix()[ii][jj] == bb.matrix()[ii][jj]);
			}
		}
		
		CHECK(not operations::io::load("directory_that_doesnt_exist", bb));
		
	}
}
#endif
