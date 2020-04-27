/* -*- indent-tabs-mode: t -*- */

#ifndef UTILS__MATCH
#define UTILS__MATCH

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

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

namespace utils {

  class match {

  public:
		
		match(double arg_tol)
      :tol_(arg_tol),
       ok_(true)
		{
		}

    template <class Type>
    auto check(const std::string & match_name, const Type & reference, const Type & value){

      auto diff = fabs(reference - value);
      
      if(diff > tol_){
        std::cout << std::endl;        
        std::cout << "Match '" + match_name + "': FAILED" << std::endl;
        std::cout << "  reference value  = " << reference << std::endl;
        std::cout << "  calculated value = " << value << std::endl;
        std::cout << "  difference       = " << diff << std::endl;
        std::cout << "  tolerance        = " << tol_ << std::endl;
        std::cout << std::endl;
        ok_ = false;
        return false;
      } else {
        std::cout << "Match '" + match_name + "': SUCCESS" << std::endl;
        return true;
      }
    }

    auto ok() const {
      return ok_;
    }

    auto fail() const {
      return not ok();
    }
    
	protected:

    double tol_;
    bool ok_;   
    
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

#include <mpi3/environment.hpp>

TEST_CASE("class utils::match", "[utils::match]") {

  utils::match mtc(1e-7);

  CHECK(mtc.ok());
  CHECK(not mtc.fail());  
  
  CHECK(mtc.check("test true", 10.0, 10.0 + 1e-8));

  CHECK(mtc.ok());
  CHECK(not mtc.fail());  

  CHECK(not mtc.check("test false", 3.0, 4.0));

  CHECK(not mtc.ok());
  CHECK(mtc.fail());
  
}

#endif

    
#endif
