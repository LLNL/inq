/* -*- indent-tabs-mode: t -*- */

#ifndef KS_STATES
#define KS_STATES

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

#include <math/complex.hpp>
#include <multi/array.hpp>
#include <basis/real_space.hpp>

namespace states {
  class ks_states {

  public:

    typedef complex coeff_type;
    
    enum class spin_config {
      UNPOLARIZED,
      POLARIZED,
      NON_COLLINEAR
    };
        
    ks_states(const spin_config spin, const double nelectrons, const int extra_states = 0){

      if(spin == spin_config::NON_COLLINEAR){
				nstates_ = ceil(nelectrons);
      } else {
				nstates_ = ceil(0.5*nelectrons);
      }

			nstates_ += extra_states;
			
      nquantumnumbers_ = 1;
      if(spin == spin_config::POLARIZED) nquantumnumbers_ = 2;

			occs_.reextent({nstates_});

			auto rem_electrons = nelectrons;
			for(int ist = 0; ist < nstates_; ist++){
				occs_[ist] = std::min(2.0, rem_electrons);
				rem_electrons -= occs_[ist];
			}
			
    }

    int num_states() const {
      return nstates_;
    }

    int num_quantum_numbers() const {
      return nquantumnumbers_;
    }

    template <class array_type>
    std::array<long int, 4> cubic_dims(const array_type & basis_dims) const {
      return {basis_dims[0], basis_dims[1], basis_dims[2], nstates_};
    }

    template <class array_type>
    std::array<long int, 2> linear_dims(const array_type & basis_dims) const {
      return {basis_dims[0]*basis_dims[1]*basis_dims[2], nstates_};
    }
    
    template <class output_stream>
    void info(output_stream & out) const {
      out << "KOHN-SHAM STATES:" << std::endl;
      out << "  Number of states = " << num_states() << std::endl;
      out << std::endl;
    }

		auto & occupations() const {
			return occs_;
		}
		
  private:

    int nstates_;
    int nquantumnumbers_;
		boost::multi::array<double, 1> occs_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>

TEST_CASE("Class states::ks_states", "[ks_states]"){

  using math::d3vector;
  
  double ecut = 30.0;
  double ll = 10.0;
  
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::real_space pw(cell, input::basis::cutoff_energy(ecut));
  
  SECTION("Spin unpolarized"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);
    
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 1);
		REQUIRE(st.cubic_dims(pw.rsize())[0] == pw.rsize()[0]);
		REQUIRE(st.cubic_dims(pw.rsize())[1] == pw.rsize()[1]);
		REQUIRE(st.cubic_dims(pw.rsize())[2] == pw.rsize()[2]);
		REQUIRE(st.cubic_dims(pw.rsize())[3] == st.num_states());

		REQUIRE(st.occupations()[0] == 2.0);
		REQUIRE(st.occupations()[1] == 2.0);
		REQUIRE(st.occupations()[2] == 2.0);
		REQUIRE(st.occupations()[3] == 2.0);
		REQUIRE(st.occupations()[4] == 2.0);
		REQUIRE(st.occupations()[5] == 1.0);
		
  }

  SECTION("Spin polarized"){
    
    states::ks_states st(states::ks_states::spin_config::POLARIZED, 11.0);
    
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 2);
    REQUIRE(st.cubic_dims(pw.rsize())[0] == pw.rsize()[0]);
    REQUIRE(st.cubic_dims(pw.rsize())[1] == pw.rsize()[1]);
    REQUIRE(st.cubic_dims(pw.rsize())[2] == pw.rsize()[2]);
    REQUIRE(st.cubic_dims(pw.rsize())[3] == st.num_states());
  }

  SECTION("Non-collinear spin"){
    
    states::ks_states st(states::ks_states::spin_config::NON_COLLINEAR, 11.0);
    
    REQUIRE(st.num_states() == 11);
    REQUIRE(st.num_quantum_numbers() == 1);
    REQUIRE(st.cubic_dims(pw.rsize())[0] == pw.rsize()[0]);
    REQUIRE(st.cubic_dims(pw.rsize())[1] == pw.rsize()[1]);
    REQUIRE(st.cubic_dims(pw.rsize())[2] == pw.rsize()[2]);
    REQUIRE(st.cubic_dims(pw.rsize())[3] == st.num_states());
  }

  
}

#endif

#endif
