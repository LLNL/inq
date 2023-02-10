/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__MATCH
#define INQ__UTILS__MATCH

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
#include <cmath>

#include <tinyformat/tinyformat.h>

namespace inq {
namespace utils {

  class match {

  public:
		
		match(double arg_tol)
      :tol_(arg_tol),
       ok_(true)
		{
		}

    template <class Type>
    auto check(const std::string & match_name, const Type & value, const Type & reference){

      auto diff = fabs(reference - value);
      
      if(diff > tol_){//tfm::format(std::cout, "step %9d :  t =  %9.3f e = %.12f\n", 0, 0.0, energy.total());
				
        tfm::format(std::cout, "\nMatch '%s': FAILED\n", match_name);
        tfm::format(std::cout, "  reference value  = %.12f\n", reference);
        tfm::format(std::cout, "  calculated value = %.12f\n", value);
        tfm::format(std::cout, "  difference       = %.1e\n", diff);
				tfm::format(std::cout, "  tolerance        = %.1e\n\n", tol_);
        ok_ = false;
        return false;
      } else {
        tfm::format(std::cout, "Match '%s': SUCCESS\n", match_name);
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
}
#endif

#ifdef INQ_UTILS_MATCH_UNIT_TEST
#undef INQ_UTILS_MATCH_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class utils::match", "[utils::match]") {

	inq::utils::match mtc(1e-7);

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
