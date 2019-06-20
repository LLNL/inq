/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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
#include <basis/plane_wave.hpp>

namespace states {
  class ks_states {

  public:

    typedef complex coeff_type;
    typedef boost::multi::array<complex, 4> coeff;
    
    enum class spin_config {
      UNPOLARIZED,
      POLARIZED,
      NON_COLLINEAR
    };
        
    ks_states(const spin_config spin, const double nelectrons){

      if(spin == spin_config::NON_COLLINEAR){
				nstates_ = ceil(nelectrons);
      } else {
				nstates_ = ceil(0.5*nelectrons);
      }

      nquantumnumbers_ = 1;
      if(spin == spin_config::POLARIZED) nquantumnumbers_ = 2;
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
    
  private:

    int nspinor_;
    int nstates_;
    int nquantumnumbers_;

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
  basis::plane_wave pw(cell, ecut);
  
  SECTION("Spin unpolarized"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);
    
    REQUIRE(st.num_states() == 6);
    REQUIRE(st.num_quantum_numbers() == 1);
		REQUIRE(st.cubic_dims(pw.rsize())[0] == pw.rsize()[0]);
		REQUIRE(st.cubic_dims(pw.rsize())[1] == pw.rsize()[1]);
		REQUIRE(st.cubic_dims(pw.rsize())[2] == pw.rsize()[2]);
		REQUIRE(st.cubic_dims(pw.rsize())[3] == st.num_states());

		states::ks_states::coeff(st.cubic_dims(pw.rsize()));
		
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

// Local variables:
// eval: (setq indent-tabs-mode: t tab-width: 2)
// End:
