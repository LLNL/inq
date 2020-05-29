/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__TRIVIAL
#define INQ__BASIS__TRIVIAL

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

#include <cassert>
#include <array>

#include <basis/base.hpp>

namespace inq {
namespace basis {

	/*
		This is a class that implements a very simple basis object. Useful for testing.
	*/
	
  class trivial : public base {

  public:

		const static int dimension = 1;
		
		trivial(const long size, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			base(size, comm),
      size_(size){
		}
		
    long size() const {
      return size_;
    }

		friend auto sizes(const trivial & ss){
      return std::array<long, 1>{ss.size_};
		}

		auto sizes() const {
			return std::array<long, 1>{size_};
		}

    double volume_element() const {
      return 1.0/size_;
    }    

    friend bool operator==(trivial const & b1, trivial const & b2){
      return b1.size_ == b2.size_;
    }
    
	protected:

    long size_;
		
  };

}
}

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::trivial", "[basis::trivial]") {
  
	using namespace inq;
	using namespace Catch::literals;
  using math::vec3d;
  
}
#endif

    
#endif
