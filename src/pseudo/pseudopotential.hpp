#ifndef PSEUDOPOTENTIAL_HPP
#define PSEUDOPOTENTIAL_HPP

/*
 Copyright (C) 2018-2019 Xavier Andrade

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

#include "../pseudo/psml.hpp"
#include "../pseudo/qso.hpp"
#include "../pseudo/upf1.hpp"
#include "../pseudo/upf2.hpp"
#include "../pseudo/psp8.hpp"
#include "../pseudo/detect_format.hpp"

#include "../math/spline.hpp"

namespace pseudo {

  class pseudopotential {

  public:

    enum class error {
      FILE_NOT_FOUND,
      UNKNOWN_FORMAT,
      UNSUPPORTED_FORMAT,
      UNSUPPORTED_TYPE
    };
    
    pseudopotential(const std::string & filename){

      //PARSE THE FILE
      base * pseudo;    
      pseudo::format format = pseudo::detect_format(filename);
      
      if(format == pseudo::format::FILE_NOT_FOUND) throw error::FILE_NOT_FOUND;
      if(format == pseudo::format::UNKNOWN) throw error::UNKNOWN_FORMAT;
      
      switch(format){
      case pseudo::format::QSO:
	pseudo = new pseudo::qso(filename);
	break;
      case pseudo::format::UPF1:
	pseudo = new pseudo::upf1(filename);
	break;
      case pseudo::format::UPF2:
	pseudo = new pseudo::upf2(filename);
	break;
      case pseudo::format::PSML:
	pseudo = new pseudo::psml(filename);
	break;
      case pseudo::format::PSP8:
	pseudo = new pseudo::psp8(filename);
	break;
      default:
	delete pseudo;
	throw error::UNSUPPORTED_FORMAT;
      }

      if(pseudo->type() != pseudo::type::KLEINMAN_BYLANDER) {
	delete pseudo;
	throw error::UNSUPPORTED_TYPE;
      }
      
      std::vector<double> grid, local_potential;

      pseudo->grid(grid);

      valence_charge_ = pseudo->valence_charge();

      //SEPARATE THE LOCAL PART
      
      sigma_erf_ = 0.625; // the constant to separate the pseudo
      
      pseudo->local_potential(local_potential);

      for(unsigned ii = 0; ii < local_potential.size(); ii++){
	local_potential[ii] -= long_range_potential(grid[ii]);
      }

      short_range_.fit(grid.data(), local_potential.data(), local_potential.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC);

      delete pseudo;
    }

    const double & valence_charge() const {
      return valence_charge_;
    }

    double long_range_potential(double rr) const {
      if(rr < 1e-8) return -valence_charge_*2.0/(sqrt(2.0*M_PI)*sigma_erf_);
      return -valence_charge_*erf(rr/(sigma_erf_*sqrt(2.0)))/rr;
    }

    //outside this, the long_range_density is zero
    double long_range_density_radius() const {
      const double tol = 1e-15;
      return sqrt(-2.0*log(tol))*sigma_erf_;
    }    
    
    double long_range_density(double rr) const {
      const double exp_arg = -0.5*pow(rr/sigma_erf_, 2);
      return -valence_charge_*exp(exp_arg)/pow(sigma_erf_*sqrt(2.0*M_PI), 3);
    }
    
    const math::spline & short_range_potential() const {
      return short_range_;
    }
    
  private:

    double sigma_erf_;
    math::spline short_range_;
    double valence_charge_;
    
  };
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class pseudo::pseudopotential", "[pseudopotential]") {
  
  using namespace Catch::literals;

  SECTION("Non-existing file"){
    REQUIRE_THROWS(pseudo::pseudopotential("/this_file_doesnt_exists"));
  }

  SECTION("Non-pseudopotential file"){
    REQUIRE_THROWS(pseudo::pseudopotential(config::path::unit_tests_data() + "benzene.xyz"));
  }

  SECTION("Non-supported format QSO/Qbox pseudopotential file"){
    REQUIRE_THROWS(pseudo::pseudopotential(config::path::unit_tests_data() + "I_HSCV_LDA-1.0.xml"));
  }

  SECTION("UPF2 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "W_ONCV_PBE-1.0.upf");

    REQUIRE(ps.valence_charge() == 28.0_a);

    REQUIRE(ps.long_range_density_radius() == 5.1945566758_a);
    
    //values validated with Octopus
    REQUIRE(ps.long_range_potential(0.00000000E+00) == -3.57452283E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-04) == -3.57452282E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-02) == -3.57437033E+01_a);
    REQUIRE(ps.long_range_potential(5.00000000E-02) == -3.57071367E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-01) == -3.55932992E+01_a);
    REQUIRE(ps.long_range_potential(5.00000000E-01) == -3.22721954E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-00) == -2.49312397E+01_a);
    REQUIRE(ps.long_range_potential(5.00000000E-00) == -5.60000000E+00_a);
    REQUIRE(ps.long_range_potential(6.01000000E-00) == -4.65890183E+00_a);

    REQUIRE(ps.long_range_density(0.00000000E+00) == -7.28195812E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-04) == -7.28195802E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-02) == -7.28102609E+00_a);
    REQUIRE(ps.long_range_density(5.00000000E-02) == -7.25869310E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-01) == -7.18934306E+00_a);
    REQUIRE(ps.long_range_density(5.00000000E-01) == -5.28778688E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-00) == -2.02465598E+00_a);
    REQUIRE(ps.long_range_density(5.00000000E-00) == -9.22199231E-14_a);
    REQUIRE(ps.long_range_density(6.01000000E-00) == -6.07009137E-20_a);
    
    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == 4.99765777E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == 4.99014665E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == 4.84335957E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == 4.68743986E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == 3.30185087E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == 7.75571234E-01_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-00) == -2.17300001E-06_a);
    REQUIRE(ps.short_range_potential().value(6.01000000E-00) == -1.05361714E-06_a);

  }

  SECTION("UPF1 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "F.UPF");

    REQUIRE(ps.valence_charge() == 7.0_a);

    REQUIRE(ps.long_range_density_radius() == 5.1945566758_a);
  }

  SECTION("PSP8 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "78_Pt_r.oncvpsp.psp8");

    REQUIRE(ps.valence_charge() == 18.0_a);

    REQUIRE(ps.long_range_density_radius() == 5.1945566758_a);
    
    //values validated with Octopus
    REQUIRE(ps.long_range_potential(0.00000000E+00) == -2.29790754E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-04) == -2.29790753E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-02) == -2.29780949E+01_a);
    REQUIRE(ps.long_range_potential(5.00000000E-02) == -2.29545879E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-01) == -2.28814066E+01_a);
    REQUIRE(ps.long_range_potential(5.00000000E-01) == -2.07464113E+01_a);
    REQUIRE(ps.long_range_potential(1.00000000E-00) == -1.60272255E+01_a);
    REQUIRE(ps.long_range_potential(4.99000000E+00) == -3.60721443E+00_a);

    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == 5.77774525E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == 5.77700949E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == 5.75937414E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == 5.70458232E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == 4.18518028E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == 1.52278621E+00_a);
    REQUIRE(ps.short_range_potential().value(4.99000000E+00) == 8.23104510E-07_a);
    
  }

  SECTION("QSO/Qbox pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml");

    REQUIRE(ps.valence_charge() == 4.0_a);

    REQUIRE(ps.long_range_density_radius() == 5.1945566758_a);
    
    //values validated with Octopus
    REQUIRE(ps.long_range_potential(0.00000000E+00) == -5.10646119E+00_a);
    REQUIRE(ps.long_range_potential(1.00000000E-04) == -5.10646116E+00_a);
    REQUIRE(ps.long_range_potential(1.00000000E-02) == -5.10624332E+00_a);
    REQUIRE(ps.long_range_potential(5.00000000E-02) == -5.10101952E+00_a);
    REQUIRE(ps.long_range_potential(1.00000000E-01) == -5.08475703E+00_a);
    REQUIRE(ps.long_range_potential(5.00000000E-01) == -4.61031362E+00_a);
    REQUIRE(ps.long_range_potential(1.00000000E-00) == -3.56160567E+00_a);
    REQUIRE(ps.long_range_potential(5.00000000E-00) == -8.00000000E-01_a);
    REQUIRE(ps.long_range_potential(6.01000000E-00) == -6.65557404E-01_a);

    REQUIRE(ps.long_range_density(0.00000000E+00) == -1.04027973E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-02) == -1.04014658E+00_a);
    REQUIRE(ps.long_range_density(5.00000000E-02) == -1.03695616E+00_a);
    REQUIRE(ps.long_range_density(1.00000000E-01) == -1.02704901E+00_a);
    REQUIRE(ps.long_range_density(5.00000000E-01) == -7.55398125E-01_a);
    REQUIRE(ps.long_range_density(1.00000000E-00) == -2.89236568E-01_a);
    REQUIRE(ps.long_range_density(5.00000000E-00) == -1.31742747E-14_a);
    REQUIRE(ps.long_range_density(6.01000000E-00) == -8.67155910E-21_a);
    
    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == -5.42158341E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == -5.41538787E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == -5.27845780E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == -4.99611728E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == -2.56127957E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == -4.36986709E-01_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-00) == -9.28740001E-07_a);
    REQUIRE(ps.short_range_potential().value(6.01000000E-00) == -7.59993877E-07_a);

  }

}
#endif
 
#endif
