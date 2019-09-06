/* -*- indent-tabs-mode: t -*- */

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

#include <pseudo/psml.hpp>
#include <pseudo/qso.hpp>
#include <pseudo/upf1.hpp>
#include <pseudo/upf2.hpp>
#include <pseudo/psp8.hpp>
#include <pseudo/detect_format.hpp>
#include <math/spline.hpp>
#include <math/erf_range_separation.hpp>

namespace pseudo {

  class pseudopotential {

  public:

    enum class error {
      FILE_NOT_FOUND,
      UNKNOWN_FORMAT,
      UNSUPPORTED_FORMAT,
      UNSUPPORTED_TYPE,
      PROJECTOR_INDEX_OUT_OF_RANGE
    };
    
    pseudopotential(const std::string & filename, const math::erf_range_separation & sep){

      //PARSE THE FILE
      base * pseudo_file;    
      pseudo::format format = pseudo::detect_format(filename);
      
      if(format == pseudo::format::FILE_NOT_FOUND) throw error::FILE_NOT_FOUND;
      if(format == pseudo::format::UNKNOWN) throw error::UNKNOWN_FORMAT;
      
      switch(format){
      case pseudo::format::QSO:
				pseudo_file = new pseudo::qso(filename);
				break;
      case pseudo::format::UPF1:
				pseudo_file = new pseudo::upf1(filename);
				break;
      case pseudo::format::UPF2:
				pseudo_file = new pseudo::upf2(filename);
				break;
      case pseudo::format::PSML:
				pseudo_file = new pseudo::psml(filename);
				break;
      case pseudo::format::PSP8:
				pseudo_file = new pseudo::psp8(filename);
				break;
      default:
				delete pseudo_file;
				throw error::UNSUPPORTED_FORMAT;
      }

      if(pseudo_file->type() != pseudo::type::KLEINMAN_BYLANDER) {
				delete pseudo_file;
				throw error::UNSUPPORTED_TYPE;
      }
      
      std::vector<double> local_potential;

      pseudo_file->grid(grid_);

      valence_charge_ = pseudo_file->valence_charge();

      //SEPARATE THE LOCAL PART
      
      pseudo_file->local_potential(local_potential);

      for(unsigned ii = 0; ii < local_potential.size(); ii++){
				local_potential[ii] -= valence_charge_*sep.long_range_potential(grid_[ii]);
      }

      short_range_.fit(grid_.data(), local_potential.data(), local_potential.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC);

      //THE PROJECTORS

      std::vector<double> proj;
      
      nproj_lm_ = 0;
      for(int ll = 0; ll <= pseudo_file->lmax(); ll++){
				for(int ichan = 0; ichan < pseudo_file->nchannels(); ichan++){
					if(ll == pseudo_file->llocal()) continue;
					
					pseudo_file->projector(ll, ichan, proj);
					
					if(proj.size() == 0) continue;
					
					assert(proj.size() <= grid_.size());
					
					projectors_.push_back(math::spline(grid_.data(), proj.data(), proj.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC));
					projectors_l_.push_back(ll);
					kb_coeff_.push_back(pseudo_file->d_ij(ll, ichan, ichan));
					nproj_lm_ += 2*ll + 1;
				}
      }
      
      delete pseudo_file;
    }
		
    const double & valence_charge() const {
      return valence_charge_;
    }
      
    const math::spline & short_range_potential() const {
      return short_range_;
    }

		double short_range_potential_radius() const {
			/*
				Precision in energy seems to be 1 to 2 orders of magnitude
				lower than the threshold. This value gives around 1e-7 error
				in the energy.

				Since this value does not affect the cost of the application
				of the hamiltonian we can use a small threshold.
				
				It makes the calculation of the local potential more expensive, so we
				should tweak if that becomes a bottleneck.

			*/
			const double threshold = 1e-5;
      return short_range_.cutoff_radius(threshold);
    }
		
    double projector_radius() const {
      const double threshold = 0.001;
      double radius = 0.0;
      for(unsigned iproj = 0; iproj < projectors_.size(); iproj++){
				radius = std::max(radius, projectors_[iproj].cutoff_radius(threshold));
      }
      return radius;
    }

    // the number of projectors with different l
    int num_projectors_l() const {
      return projectors_.size();
    }

    //the number of projectors with different l and m
    int num_projectors_lm() const {
      return nproj_lm_;
    }

    const math::spline & projector(int iproj) const {
      if(iproj < 0 or unsigned(iproj) >= projectors_.size()) throw error::PROJECTOR_INDEX_OUT_OF_RANGE;
      return projectors_[iproj];
    }

    int projector_l(int iproj) const {
      if(iproj < 0 or unsigned(iproj) >= projectors_l_.size()) throw error::PROJECTOR_INDEX_OUT_OF_RANGE;
      return projectors_l_[iproj];
    }

    double kb_coeff(int iproj){
      if(iproj < 0 or unsigned(iproj) >= kb_coeff_.size()) throw error::PROJECTOR_INDEX_OUT_OF_RANGE;
      return kb_coeff_[iproj];
    }
    
  private:
    
    std::vector<double> grid_;
    double sigma_erf_;
    math::spline short_range_;
    double valence_charge_;
    std::vector<math::spline> projectors_;
    std::vector<int> projectors_l_;
    std::vector<double> kb_coeff_;
    int nproj_lm_;
    
  };
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class pseudo::pseudopotential", "[pseudopotential]") {
  
  using namespace Catch::literals;

	const math::erf_range_separation sep(0.625);
	
  SECTION("Non-existing file"){
    REQUIRE_THROWS(pseudo::pseudopotential("/this_file_doesnt_exists", sep));
  }

  SECTION("Non-pseudopotential file"){
    REQUIRE_THROWS(pseudo::pseudopotential(config::path::unit_tests_data() + "benzene.xyz", sep));
  }

  SECTION("Non-supported format QSO/Qbox pseudopotential file"){
    REQUIRE_THROWS(pseudo::pseudopotential(config::path::unit_tests_data() + "I_HSCV_LDA-1.0.xml", sep));
  }

  SECTION("UPF2 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "W_ONCV_PBE-1.0.upf", sep);

    REQUIRE(ps.valence_charge() == 28.0_a);

    REQUIRE(sep.long_range_density_radius() == 5.0922713869_a);
    
    //values validated with Octopus
    REQUIRE(ps.valence_charge()*sep.long_range_potential(0.00000000E+00) == -3.57452283E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-04) == -3.57452282E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-02) == -3.57437033E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-02) == -3.57071367E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-01) == -3.55932992E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-01) == -3.22721954E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-00) == -2.49312397E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-00) == -5.60000000E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(6.01000000E-00) == -4.65890183E+00_a);

    REQUIRE(ps.valence_charge()*sep.long_range_density(0.00000000E+00) == -7.28195812E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-04) == -7.28195802E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-02) == -7.28102609E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-02) == -7.25869310E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-01) == -7.18934306E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-01) == -5.28778688E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-00) == -2.02465598E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-00) == -9.22199231E-14_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(6.01000000E-00) == -6.07009137E-20_a);

    REQUIRE(ps.short_range_potential_radius() == 3.37_a);
		
    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == 4.99765777E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == 4.99014665E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == 4.84335957E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == 4.68743986E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == 3.30185087E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == 7.75571234E-01_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-00) == -2.17300001E-06_a);
    REQUIRE(ps.short_range_potential().value(6.01000000E-00) == -1.05361714E-06_a);

    REQUIRE(ps.num_projectors_l() == 8);
    REQUIRE(ps.num_projectors_lm() == 32);
    
    REQUIRE(ps.projector(0).value(0.00000000E+00) == -6.16416761E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-02) == -6.16245244E+00_a);
    REQUIRE(ps.projector(0).value(5.00000000E-02) == -6.12149052E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-01) == -5.99451036E+00_a);
    REQUIRE(ps.projector(0).value(5.00000000E-01) == -2.72974026E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-00) ==  4.73465778E-01_a);
    REQUIRE(ps.projector(0).value(5.00000000E-00) ==  0.00000000E+00_a);
    REQUIRE(ps.projector(0).value(6.01000000E-00) ==  0.00000000E+00_a);

    REQUIRE_THROWS(ps.projector(-1));
    REQUIRE_THROWS(ps.projector(8));

    REQUIRE(ps.projector_l(0) == 0);
    REQUIRE(ps.projector_l(1) == 0);
    REQUIRE(ps.projector_l(2) == 1);
    REQUIRE(ps.projector_l(3) == 1);
    REQUIRE(ps.projector_l(4) == 2);
    REQUIRE(ps.projector_l(5) == 2);
    REQUIRE(ps.projector_l(6) == 3);
    REQUIRE(ps.projector_l(7) == 3);

		REQUIRE(ps.kb_coeff(0) == 8.243498986_a);
    REQUIRE(ps.kb_coeff(1) == 0.9729808132_a);
    REQUIRE(ps.kb_coeff(2) == 3.8116529178_a);
    REQUIRE(ps.kb_coeff(3) == 0.4990369016_a);
    REQUIRE(ps.kb_coeff(4) == 1.4609009471_a);
    REQUIRE(ps.kb_coeff(5) == 0.2713485048_a);
    REQUIRE(ps.kb_coeff(6) == -5.0048898015_a);
		REQUIRE(ps.kb_coeff(7) == -2.9667606006_a);

    REQUIRE_THROWS(ps.projector_l(-100));    
    REQUIRE_THROWS(ps.projector_l(8));    
    
    REQUIRE(ps.projector_radius() == 2.88_a);
    
  }

  SECTION("UPF1 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "F.UPF", sep);

    REQUIRE(ps.valence_charge() == 7.0_a);

    REQUIRE(sep.long_range_density_radius() == 5.0922713869_a);

    REQUIRE(ps.short_range_potential_radius() == 2.8262095048_a);
		
    REQUIRE(ps.num_projectors_l() == 1);
    REQUIRE(ps.num_projectors_lm() == 1);
    
    REQUIRE(ps.projector_radius() == 0.8138477088_a);
  }

  SECTION("PSP8 pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "78_Pt_r.oncvpsp.psp8", sep);

    REQUIRE(ps.valence_charge() == 18.0_a);

    REQUIRE(sep.long_range_density_radius() == 5.0922713869_a);

    REQUIRE(ps.short_range_potential_radius() == 3.09_a);
		
    //values validated with Octopus
    REQUIRE(ps.valence_charge()*sep.long_range_potential(0.00000000E+00) == -2.29790754E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-04) == -2.29790753E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-02) == -2.29780949E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-02) == -2.29545879E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-01) == -2.28814066E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-01) == -2.07464113E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-00) == -1.60272255E+01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(4.99000000E+00) == -3.60721443E+00_a);

    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == 5.77774525E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == 5.77700949E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == 5.75937414E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == 5.70458232E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == 4.18518028E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == 1.52278621E+00_a);
    REQUIRE(ps.short_range_potential().value(4.99000000E+00) == 8.23104510E-07_a);

    REQUIRE(ps.projector_radius() == 3.03_a);
  }

  SECTION("QSO/Qbox pseudopotential file"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "C_ONCV_PBE-1.2.xml", sep);

    REQUIRE(ps.valence_charge() == 4.0_a);

    REQUIRE(sep.long_range_density_radius() == 5.0922713869_a);

    REQUIRE(ps.short_range_potential_radius() == 2.86_a);
		
    //values validated with Octopus
    REQUIRE(ps.valence_charge()*sep.long_range_potential(0.00000000E+00) == -5.10646119E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-04) == -5.10646116E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-02) == -5.10624332E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-02) == -5.10101952E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-01) == -5.08475703E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-01) == -4.61031362E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(1.00000000E-00) == -3.56160567E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(5.00000000E-00) == -8.00000000E-01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_potential(6.01000000E-00) == -6.65557404E-01_a);

    REQUIRE(ps.valence_charge()*sep.long_range_density(0.00000000E+00) == -1.04027973E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-02) == -1.04014658E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-02) == -1.03695616E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-01) == -1.02704901E+00_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-01) == -7.55398125E-01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(1.00000000E-00) == -2.89236568E-01_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(5.00000000E-00) == -1.31742747E-14_a);
    REQUIRE(ps.valence_charge()*sep.long_range_density(6.01000000E-00) == -8.67155910E-21_a);
    
    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == -5.42158341E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == -5.41538787E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == -5.27845780E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == -4.99611728E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == -2.56127957E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == -4.36986709E-01_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-00) == -9.28740001E-07_a);
    REQUIRE(ps.short_range_potential().value(6.01000000E-00) == -7.59993877E-07_a);

    REQUIRE(ps.projector_radius() == 1.33_a);

  }
	
  SECTION("UPF nitrogen"){
    pseudo::pseudopotential ps(config::path::unit_tests_data() + "N.upf", sep);

    REQUIRE(ps.valence_charge() == 5.0_a);

    REQUIRE(sep.long_range_density_radius() == 5.0922713869_a);

    REQUIRE(ps.short_range_potential_radius() == 2.9_a);

		REQUIRE(ps.num_projectors_l() == 4);
    REQUIRE(ps.num_projectors_lm() == 8);

		REQUIRE(ps.projector_l(0) == 0);
		REQUIRE(ps.projector_l(1) == 0);
		REQUIRE(ps.projector_l(2) == 1);
		REQUIRE(ps.projector_l(3) == 1);

		REQUIRE(ps.kb_coeff(0) == 7.494508815_a);
		REQUIRE(ps.kb_coeff(1) == 0.6363049519_a);
		REQUIRE(ps.kb_coeff(2) == -4.2939052122_a);
		REQUIRE(ps.kb_coeff(3) == -1.0069878791_a);

    //values validated with Octopus
    REQUIRE(ps.short_range_potential().value(0.00000000E+00) == -2.18602017E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-02) == -2.18565897E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-02) == -2.17726246E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-01) == -2.15397947E+00_a);
    REQUIRE(ps.short_range_potential().value(5.00000000E-01) == -1.79645565E+00_a);
    REQUIRE(ps.short_range_potential().value(1.00000000E-00) == -5.32185535E-01_a);
		REQUIRE(ps.short_range_potential().value(2.00000000E-00) == -3.45198884E-03_a);
		REQUIRE(ps.short_range_potential().value(3.00000000E-00) == -6.18686051E-06_a);

    REQUIRE(ps.projector_radius() == 1.35_a);
    
    REQUIRE(ps.projector(0).value(0.00000000E+00) ==  9.45961545E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-02) ==  9.45376646E+00_a);
    REQUIRE(ps.projector(0).value(5.00000000E-02) ==  9.31684484E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-01) ==  8.90580060E+00_a);
    REQUIRE(ps.projector(0).value(5.00000000E-01) ==  2.21693494E+00_a);
    REQUIRE(ps.projector(0).value(1.00000000E-00) == -1.24268183E-01_a);
    REQUIRE(ps.projector(0).value(1.20000000E+00) == -1.28767214E-04_a);

		REQUIRE(ps.projector(1).value(0.00000000E+00) == -6.07818791E+00_a);
		REQUIRE(ps.projector(1).value(1.00000000E-02) == -6.06581452E+00_a);
    REQUIRE(ps.projector(1).value(5.00000000E-02) == -5.78271364E+00_a);
    REQUIRE(ps.projector(1).value(1.00000000E-01) == -4.96486497E+00_a);
    REQUIRE(ps.projector(1).value(5.00000000E-01) ==  1.62029039E+00_a);
    REQUIRE(ps.projector(1).value(1.00000000E-00) == -5.84285308E-03_a);
    REQUIRE(ps.projector(1).value(1.20000000E+00) == -1.62497839E-04_a);
		
		REQUIRE(ps.projector(2).value(0.00000000E+00) == -1.12401448E-03_a);
		REQUIRE(ps.projector(2).value(1.00000000E-02) ==  2.65800505E-01_a);
    REQUIRE(ps.projector(2).value(5.00000000E-02) ==  1.30653138E+00_a);
    REQUIRE(ps.projector(2).value(1.00000000E-01) ==  2.47873383E+00_a);
    REQUIRE(ps.projector(2).value(5.00000000E-01) ==  3.17520408E+00_a);
    REQUIRE(ps.projector(2).value(1.00000000E-00) == -4.57974975E-01_a);
    REQUIRE(ps.projector(2).value(1.20000000E+00) ==  4.69107012E-02_a);
		
		REQUIRE(ps.projector(3).value(0.00000000E+00) ==  1.78824837E-03_a);
		REQUIRE(ps.projector(3).value(1.00000000E-02) == -1.93734531E-01_a);
    REQUIRE(ps.projector(3).value(5.00000000E-02) == -9.32927457E-01_a);
    REQUIRE(ps.projector(3).value(1.00000000E-01) == -1.65577954E+00_a);
    REQUIRE(ps.projector(3).value(5.00000000E-01) ==  8.85876846E-01_a);
    REQUIRE(ps.projector(3).value(1.00000000E-00) ==  6.61170461E-01_a);
    REQUIRE(ps.projector(3).value(1.20000000E+00) == -3.55726746E-01_a);
		
  }
}
#endif
 
#endif
