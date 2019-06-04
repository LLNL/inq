#ifndef BASIS_HPP
#define BASIS_HPP

#include <math/d3vector.hpp>
#include <cassert>

namespace basis {
  class basis {

  public:

    basis(ions::UnitCell & cell, const double & ecut):nr_(3), ng_(3) {
      ecut_ = ecut;
      rspacing_ = math::d3vector(M_PI*sqrt(0.5/ecut));

      // make the spacing conmensurate with the grid
      // OPTIMIZATION: we can select a good size here for the FFT
      for(int idir = 0; idir < 3; idir++){
	rlength_[idir] = length(cell[idir]);
	
	nr_[idir] = round(rlength_[idir]/rspacing_[idir]);
	ng_[idir] = nr_[idir];

	rspacing_[idir] = rlength_[idir]/nr_[idir];
	glength_[idir] = M_PI/rspacing_[idir];
	gspacing_[idir] = glength_[idir]/ng_[idir];
      }

    }

    const double & ecut() const {
      return ecut_;
    }
    
    const std::vector<int> & rsize() const{
      return nr_;
    }

    const std::vector<int> & gsize() const{
      return ng_;
    }

    const math::d3vector & rspacing() const{
      return rspacing_;
    }

    const math::d3vector & gspacing() const{
      return gspacing_;
    }
    
    const math::d3vector & rlength() const{
      return rlength_;
    }

    const math::d3vector & glength() const{
      return glength_;
    }

    int rtotalsize() const {
      return nr_[0]*nr_[1]*nr_[2];
    }

    int gtotalsize() const {
      return ng_[0]*ng_[1]*ng_[2];
    }
    
  private:
    double ecut_;
    std::vector<int> nr_;
    std::vector<int> ng_;

    math::d3vector rspacing_;
    math::d3vector gspacing_;
    
    math::d3vector rlength_;
    math::d3vector glength_;
    
    
  };
}

#ifdef UNIT_TEST
#include <catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::basis", "[basis]") {
  
  using namespace Catch::literals;
  using math::d3vector;
  
  {
    
    SECTION("Cubic cell"){

      ions::UnitCell cell(d3vector(10.0, 0.0, 0.0), d3vector(0.0, 10.0, 0.0), d3vector(0.0, 0.0, 10.0));

      double ecut = 20.0;
      
      basis::basis bas(cell, ecut);

      REQUIRE(bas.ecut() == Approx(ecut));

      REQUIRE(bas.rtotalsize() == 8000);
      REQUIRE(bas.gtotalsize() == 8000);
      
      REQUIRE(bas.rspacing()[0] == 0.5_a);
      REQUIRE(bas.rspacing()[1] == 0.5_a);
      REQUIRE(bas.rspacing()[2] == 0.5_a);
      
      REQUIRE(bas.gspacing()[0] == 0.3141592654_a);
      REQUIRE(bas.gspacing()[1] == 0.3141592654_a);
      REQUIRE(bas.gspacing()[2] == 0.3141592654_a);
      
      REQUIRE(bas.rsize()[0] == 20);
      REQUIRE(bas.rsize()[1] == 20);
      REQUIRE(bas.rsize()[2] == 20);

      REQUIRE(bas.gsize()[0] == 20);
      REQUIRE(bas.gsize()[1] == 20);
      REQUIRE(bas.gsize()[2] == 20);

    }

    SECTION("Parallelepipedic cell"){

      ions::UnitCell cell(d3vector(77.7, 0.0, 0.0), d3vector(0.0, 14.14, 0.0), d3vector(0.0, 0.0, 23.25));

      double ecut = 37.9423091;
      
      basis::basis bas(cell, ecut);

      REQUIRE(bas.ecut() == Approx(ecut));

      REQUIRE(bas.rtotalsize() == 536640);
      REQUIRE(bas.gtotalsize() == 536640);
	    
      REQUIRE(bas.rspacing()[0] == 0.3613953488_a);
      REQUIRE(bas.rspacing()[1] == 0.3625641026_a);
      REQUIRE(bas.rspacing()[2] == 0.36328125_a);
      
      REQUIRE(bas.gspacing()[0] == 0.0404323379_a);
      REQUIRE(bas.gspacing()[1] == 0.2221776983_a);
      REQUIRE(bas.gspacing()[2] == 0.1351222647_a);
      
      REQUIRE(bas.rsize()[0] == 215);
      REQUIRE(bas.rsize()[1] == 39);
      REQUIRE(bas.rsize()[2] == 64);

      REQUIRE(bas.gsize()[0] == 215);
      REQUIRE(bas.gsize()[1] == 39);
      REQUIRE(bas.gsize()[2] == 64);

    }

  }
}
#endif

    
#endif

