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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <cstdio>


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

			using Type = typename FieldSet::element_type;
			
			assert(not phi.basis().part().parallel());

			boost::multi::array<Type, 1> buffer(phi.basis().part().local_size());

			if(phi.full_comm().rank() == 0) createdir(dirname);
			phi.full_comm().barrier();
				 
			for(int ist = 0; ist < phi.set_part().local_size(); ist++){

				auto filename = dirname + "/" + numstr(ist + phi.set_part().start());

				buffer = phi.matrix().rotated()[ist];

				auto fd = open(filename.data(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
				if(fd == -1){
					std::cerr << "Error: cannot create restart file '" << filename << "'." << std::endl;
					exit(1);
				}
				
				[[maybe_unused]] auto data_written = write(fd, buffer.data(), buffer.size()*sizeof(Type));
				assert(data_written == long(buffer.size()*sizeof(Type)));

				close(fd);
				
			}
		}
		
		template <class FieldSet>
		void load(std::string const & dirname, FieldSet & phi){

			using Type = typename FieldSet::element_type;
			
			assert(not phi.basis().part().parallel());

			boost::multi::array<Type, 1> buffer(phi.basis().part().local_size());

			DIR* dir = opendir(dirname.c_str());
			if (!dir) {
				std::cerr << "Error: cannot open restart directory '" << dirname << "/'." << std::endl;
				exit(1);
			}
			
			for(int ist = 0; ist < phi.set_part().local_size(); ist++){

				auto filename = dirname + "/" + numstr(ist + phi.set_part().start());
				
				auto fd = open(filename.c_str(), O_RDONLY);
				if(fd == -1){
					std::cerr << "Error: cannot open restart file '" << filename << "'." << std::endl;
					exit(1);
				}
				
				[[maybe_unused]] auto data_read = read(fd, buffer.data(), buffer.size()*sizeof(Type));
				assert(data_read == long(buffer.size()*sizeof(Type)));

				close(fd);

				phi.matrix().rotated()[ist] = buffer;
				
			}
		}
		
	}
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/field.hpp>
#include <basis/field_set.hpp>
#include <basis/trivial.hpp>

TEST_CASE("function operations::io", "[operations::io]") {
	
	using namespace Catch::literals;
	
	const int npoint = 100;
	const int nvec = 12;
	
	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {comm.size(), 1});
	
	auto basis_comm = cart_comm.axis(1);
	
	CHECK(basis_comm.size() == 1);
		
	basis::trivial bas(npoint, basis_comm);

	basis::field_set<basis::trivial, double> aa(bas, nvec, cart_comm);
	basis::field_set<basis::trivial, double> bb(bas, nvec, cart_comm);
	
	for(int ii = 0; ii < bas.part().local_size(); ii++){
		for(int jj = 0; jj < aa.set_part().local_size(); jj++){
			auto jjg = aa.set_part().local_to_global(jj);
			auto iig = bas.part().local_to_global(ii);
			aa.matrix()[ii][jj] = 20.0*(iig + 1)*sqrt(jjg);
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
