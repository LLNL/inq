#ifndef BASIS_HPP
#define BASIS_HPP

#include <math/d3vector.hpp>

namespace basis {
  class basis {

  public:

    basis(ions::UnitCell & cell, const double ecut):nr_(3), ng_(3) {
      ecut_ = ecut;
      rspacing_ = math::d3vector(sqrt(0.5*M_PI/ecut));

      std::cout << "SPACING1 " << rspacing_[0] << std::endl;
      
      //make the spacing conmensurate with the grid
      for(int idir = 0; idir < 3; idir++){
	nr_[idir] = ceil(length(cell[idir])/rspacing_[idir]);
	rspacing_[idir] = length(cell[idir])/nr_[idir];
	gspacing_[idir] = M_PI/rspacing_[idir];
	ng_[idir] = nr_[idir];
      }

      std::cout << "SPACING2 " << rspacing_[0] << std::endl;

    }

    const std::vector<int> & real_space_size() const{
      return nr_;
    }

  private:
    double ecut_;
    std::vector<int> nr_;
    std::vector<int> ng_;
    math::d3vector rspacing_;
    math::d3vector gspacing_;
    
  };
}

#ifdef UNIT_TEST
#include <catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::basis", "[basis]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));
    
    SECTION("Cubic cell 0"){

      basis::basis bas(cell, 20.0);

      REQUIRE(bas.real_space_size()[0] == 36);
      REQUIRE(bas.real_space_size()[1] == 36);
      REQUIRE(bas.real_space_size()[2] == 36);
    }

  }
}
#endif

    
#endif

