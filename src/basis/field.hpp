/* -*- indent-tabs-mode: t -*- */

#ifndef BASIS_FIELD
#define BASIS_FIELD

/*
 Copyright (C) 2019 Xavier Andrade

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

#include <math/array.hpp>
#include <tinyformat/tinyformat.h>
#include <algorithm>

namespace basis {
	
	template<class b_type, class type>
  class field {

  public:

		typedef math::array<type, 1> internal_array_type;
		typedef b_type basis_type;
		typedef type element_type;
		
    field(const basis_type & basis, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			basis_comm_(comm),
			linear_(basis.dist().local_size()),
			basis_(basis){

			assert(basis_.dist().comm_size() == basis_comm_.size());
    }

		field(const field & coeff) = delete;
		field(field && coeff) = default;
		field & operator=(const field & coeff) = default;
		field & operator=(field && coeff) = default;

		//set to a scalar value
		field & operator=(const type value) {

			//DATAOPERATIONS GPU::RUN FILL
			gpu::run(linear_.size(),
							 [lin = begin(linear_), value] GPU_LAMBDA (auto ii){
								 lin[ii] = value;
							 });

			return *this;
		}

		auto size() const {
			return basis_.size();
		}
		
		const auto & basis() const {
			return basis_;
		}

		auto cubic() const {
			return linear_.partitioned(basis_.sizes()[1]*basis_.sizes()[0]).partitioned(basis_.sizes()[0]);
		}

		auto cubic() {
			return linear_.partitioned(basis_.sizes()[1]*basis_.sizes()[0]).partitioned(basis_.sizes()[0]);
		}

		auto data() {
			return linear_.data();
		}

		auto data() const {
			return linear_.data();
		}

		auto & linear() const {
			return linear_;
		}

		auto & linear() {
			return linear_;
		}

		template <int dir = 2>
		friend void print_debug(const field & fld, const std::string & filename){

			std::ofstream file(filename);

			std::array<int, 3> point = {0, 0, 0};

			auto size = std::get<dir>(sizes(fld.cubic()));
			
			for(int ii = 0; ii < size; ii++){
				auto ip = ii + size/2;
				if(ip >= size) ip -= size;
				point[dir] = ip;
				auto rr = fld.basis().rvector(point);
				tfm::format(file, "%f %e %e\n", rr[dir], real(fld.cubic()[point[0]][point[1]][point[2]]), imag(fld.cubic()[point[0]][point[1]][point[2]]));
			}
		}

		auto & basis_comm() const {
			return basis_comm_;
		}
		
				
	private:
		mutable boost::mpi3::communicator basis_comm_;
		internal_array_type linear_;
		basis_type basis_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <multi/adaptors/fftw.hpp>

TEST_CASE("Class basis::field", "[basis::field]"){

  using namespace Catch::literals;
  using math::vec3d;
  
  double ecut = 40.0;

  ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 4.0, 0.0), vec3d(0.0, 0.0, 7.0));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	basis::field<basis::real_space, double> ff(rs);

	namespace fftw = boost::multi::fftw;

	REQUIRE(sizes(rs)[0] == 28);
	REQUIRE(sizes(rs)[1] == 11);
	REQUIRE(sizes(rs)[2] == 20);	

	REQUIRE(std::get<0>(sizes(ff.linear())) == 6160);

	REQUIRE(std::get<0>(sizes(ff.cubic())) == 28);
	REQUIRE(std::get<1>(sizes(ff.cubic())) == 11);
	REQUIRE(std::get<2>(sizes(ff.cubic())) == 20);

	ff = 12.2244;

	for(int ii = 0; ii < rs.size(); ii++) REQUIRE(ff.linear()[ii] == 12.2244_a);	

}

#endif

#endif
