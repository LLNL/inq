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

#include <multi/array.hpp>
#include <tinyformat/tinyformat.h>

namespace basis {
	
	template<class b_type, class type>
  class field : public boost::multi::array<type, 1> {

  public:

		typedef type value_type;
		typedef b_type basis_type;
		
    field(const basis_type & basis):
			boost::multi::array<type, 1>(basis.size()),
			basis_(basis){
    }

		field(const field & coeff) = delete;
		field(field && coeff) = default;
		field & operator=(const field & coeff) = default;
		field & operator=(field && coeff) = default;

		//set to a scalar value
		field & operator=(const value_type value) {
			//DATAOPERATIONS
			for(int ii = 0; ii < basis_.size(); ii++) (*this)[ii] = value;

			return *this;
		}
		
		const auto & basis() const {
			return basis_;
		}

		auto cubic() const {
			return this->partitioned(basis_.sizes()[1]*basis_.sizes()[0]).partitioned(basis_.sizes()[0]);
		}

		auto cubic() {
			return this->partitioned(basis_.sizes()[1]*basis_.sizes()[0]).partitioned(basis_.sizes()[0]);
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
				
	private:
		
		basis_type basis_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <multi/adaptors/fftw.hpp>

TEST_CASE("Class basis::field", "[basis::field]"){

  using namespace Catch::literals;
  using math::d3vector;
  
  double ecut = 40.0;

  ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 4.0, 0.0), d3vector(0.0, 0.0, 7.0));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	basis::field<basis::real_space, double> ff(rs);


	namespace fftw = boost::multi::fftw;
	using boost::multi::array_ref;

	REQUIRE(sizes(rs)[0] == 28);
	REQUIRE(sizes(rs)[1] == 11);
	REQUIRE(sizes(rs)[2] == 20);	

	REQUIRE(std::get<0>(sizes(ff)) == 6160);

	REQUIRE(std::get<0>(sizes(ff.cubic())) == 28);
	REQUIRE(std::get<1>(sizes(ff.cubic())) == 11);
	REQUIRE(std::get<2>(sizes(ff.cubic())) == 20);
	
}

#endif

#endif
