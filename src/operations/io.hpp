/* -*- indent-tabs-mode: t -*- */

#ifndef OPERATIONS__IO
#define OPERATIONS__IO

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <math/array.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>

#include <string>
#include <cstdio>
#include <iostream>

#include <mpi3/detail/datatype.hpp>

#include <caliper/cali.h>

namespace inq {
namespace operations {
namespace io {

auto numstr(long num){
	char numcstr[12]; 
	snprintf(numcstr, 11, "%010ld", num);
	return std::string(numcstr);				
}

auto createdir(std::string const & dirname){

	DIR* dir = opendir(dirname.c_str());
	if (dir) {
		// The directory exists, we just use it
		closedir(dir);

	} else {
		// The directory does not exist
		auto status = mkdir(dirname.c_str(), 0777);
		if(status != 0){
			std::cerr << "Error: cannot create restart directory '" << dirname << "/'." << std::endl;
			exit(1);
		}
	}

}

template <class FieldSet>
void save(std::string const & dirname, FieldSet const & phi){

	CALI_CXX_MARK_FUNCTION;

	using Type = typename FieldSet::element_type;
	auto mpi_type = boost::mpi3::detail::basic_datatype<Type>();
	
	math::array<Type, 1> buffer(phi.basis().part().local_size());

	if(phi.full_comm().rank() == 0) createdir(dirname);
	phi.full_comm().barrier();
				 
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){

		auto filename = dirname + "/" + numstr(ist + phi.set_part().start());

		buffer = phi.matrix().rotated()[ist];

		MPI_File fh;

		auto mpi_err = MPI_File_open(phi.basis().comm().get(), filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		if(mpi_err != MPI_SUCCESS){
			std::cerr << "Error: cannot create restart file '" << filename << "'." << std::endl;
			exit(1);
		}

		MPI_Status status;
		mpi_err = MPI_File_write_at_all(fh, sizeof(Type)*phi.basis().part().start(), buffer.data(), buffer.size(), mpi_type, &status);
		
		if(mpi_err != MPI_SUCCESS){
			std::cerr << "Error: cannot write restart file '" << filename << "'." << std::endl;
			exit(1);
		}

		int data_written;
		MPI_Get_count(&status, mpi_type, &data_written);
		assert(data_written == long(buffer.size()));

		MPI_File_close(&fh);
		
	}

}
		
template <class FieldSet>
void load(std::string const & dirname, FieldSet & phi){

	CALI_CXX_MARK_FUNCTION;

	using Type = typename FieldSet::element_type;
	auto mpi_type = boost::mpi3::detail::basic_datatype<Type>();
	
	math::array<Type, 1> buffer(phi.basis().part().local_size());

	DIR* dir = opendir(dirname.c_str());
	if (!dir) {
		std::cerr << "Error: cannot open restart directory '" << dirname << "/'." << std::endl;
		exit(1);
	}
	closedir(dir);
			
	for(int ist = 0; ist < phi.set_part().local_size(); ist++){

		auto filename = dirname + "/" + numstr(ist + phi.set_part().start());

		MPI_File fh;

		auto mpi_err = MPI_File_open(phi.basis().comm().get(), filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

		if(mpi_err != MPI_SUCCESS){
			std::cerr << "Error: cannot open restart file '" << filename << "'." << std::endl;
			exit(1);
		}
		
		MPI_Status status;
		mpi_err = MPI_File_read_at_all(fh, sizeof(Type)*phi.basis().part().start(), buffer.data(), buffer.size(), mpi_type, &status);
		
		if(mpi_err != MPI_SUCCESS){
			std::cerr << "Error: cannot read restart file '" << filename << "'." << std::endl;
			exit(1);
		}

		int data_read;
		MPI_Get_count(&status, mpi_type, &data_read);
		assert(data_read == long(buffer.size()));

		MPI_File_close(&fh);

		phi.matrix().rotated()[ist] = buffer;
				
	}
}
		
}
}
}

#ifdef INQ_OPERATIONS_IO_UNIT_TEST
#undef INQ_OPERATIONS_IO_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::io", "[operations::io]") {
	
	using namespace inq;
	using namespace Catch::literals;
	
	const int npoint = 100;
	const int nvec = 12;
	
	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {});
	
	auto basis_comm = cart_comm.axis(1);
	
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

	operations::io::load("restart/", bb);
	
	for(int ii = 0; ii < bas.part().local_size(); ii++){
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			CHECK(aa.matrix()[ii][jj] == bb.matrix()[ii][jj]);
		}
	}
	
}


#endif
#endif
