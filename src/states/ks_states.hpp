/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__STATES__KS_STATES
#define INQ__STATES__KS_STATES

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa

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
#include <math/array.hpp>
#include <operations/sum.hpp>

namespace inq {
namespace states {

class ks_states {

public:

	typedef complex coeff_type;
    
	enum class spin_config {
													UNPOLARIZED,
													POLARIZED,
													NON_COLLINEAR
	};
        
	ks_states(const spin_config spin, const double nelectrons, const int extra_states = 0, double temperature = 0.0):
		temperature_(temperature)
	{
		
		max_occ_ = 1.0;
		if(spin == spin_config::UNPOLARIZED) max_occ_ = 2.0;
			
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

		num_electrons_ = nelectrons;

	}

	int num_states() const {
		return nstates_;
	}

	int num_quantum_numbers() const {
		return nquantumnumbers_;
	}

	template <class output_stream>
	void info(output_stream & out) const {
		out << "KOHN-SHAM STATES:" << std::endl;
		out << "  Number of states = " << num_states() << std::endl;
		out << std::endl;
	}
	
	template<class OStream>
	friend OStream& operator<<(OStream& os, ks_states const& self){
		self.info(os);
		return os;
	}
	
	auto & occupations() const {
		return occs_;
	}

	auto num_electrons() const {
		return num_electrons_;
	}

	auto smear_function(double xx) const {

		double const maxarg = 200.0;
		
		if(xx > maxarg){
			return 1.0;
		} else if(xx > -maxarg) {
			return 1.0/(1.0 + exp(-xx));
		} else {
			return 0.0;
		}
		
	}
	
	template <typename ArrayType>
	void update_occupations(ArrayType const eigenval) {

		assert(nstates_ == eigenval.size());
		
		double const tol = 1e-10;
		double efermi;

		if(temperature_ == 0.0){

			auto rem_electrons = num_electrons_;
			for(int ist = 0; ist < nstates_; ist++){
				occs_[ist] = std::min(2.0, rem_electrons);
				rem_electrons -= occs_[ist];
			}
			
		} else {

			int const nitmax = 200;

			auto dsmear = std::max(1e-14, temperature_);
			auto drange = dsmear*sqrt(-log(tol*0.01));

			auto emin = real(eigenval[0]) - drange;
			auto emax = real(eigenval[nstates_ - 1]) + drange;

			//check that the eigenvalues are sorted
			for(int ist = 1; ist < nstates_; ist++){
				assert(real(eigenval[ist]) >= real(eigenval[ist - 1]));
			}
		
			double sumq;
			for(int iter = 0; iter < nitmax; iter++){
				efermi = 0.5*(emin + emax);
				sumq = 0.0;

				for(int ist = 0; ist < nstates_; ist++){
					auto xx = (efermi - real(eigenval[ist]))/dsmear;
					sumq = sumq + max_occ_*smear_function(xx);
				}

				if(fabs(sumq - num_electrons_) <= tol) break;
				if(sumq <= num_electrons_) emin = efermi;
				if(sumq >= num_electrons_) emax = efermi;
				
			}

			for(int ist = 0; ist < nstates_; ist++){
				auto xx = (efermi - real(eigenval[ist]))/dsmear;
				occs_[ist] = max_occ_*smear_function(xx);
			}

			assert(fabs(operations::sum(occs_) - num_electrons_) <= tol);
			
		}
		
	}
	
private:

	double num_electrons_;
	int nstates_;
	int nquantumnumbers_;
	math::array<double, 1> occs_;
	double temperature_;
	double max_occ_;

};

}
}

#ifdef INQ_STATES_KS_STATES_UNIT_TEST
#undef INQ_STATES_KS_STATES_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("Class states::ks_states", "[ks_states]"){

	using namespace Catch::literals;
	using namespace inq;
  
  SECTION("Spin unpolarized"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

		CHECK(st.num_electrons() == 11.0);
    CHECK(st.num_states() == 6);
    CHECK(st.num_quantum_numbers() == 1);

		st.update_occupations(math::array<double, 1>{0.1, 0.2, 0.3, 0.3, 0.4, 1.0});
		
		CHECK(st.occupations()[0] == 2.0);
		CHECK(st.occupations()[1] == 2.0);
		CHECK(st.occupations()[2] == 2.0);
		CHECK(st.occupations()[3] == 2.0);
		CHECK(st.occupations()[4] == 2.0);
		CHECK(st.occupations()[5] == 1.0);

  }
	
  SECTION("Spin unpolarized with temperature"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 4.0, 4, 0.01);

		CHECK(st.num_electrons() == 4.0);    
    CHECK(st.num_states() == 6);
    CHECK(st.num_quantum_numbers() == 1);

		st.update_occupations(math::array<double, 1>{0.0, 0.2, 0.3, 0.3, 0.4, 1.0});

		CHECK(st.occupations()[0] == 2.0_a);
		CHECK(st.occupations()[1] == 1.9810772793_a);
		CHECK(st.occupations()[2] == 0.0094611446_a);
		CHECK(st.occupations()[3] == 0.0094611446_a);
		CHECK(st.occupations()[4] == 4.315768121820668e-07_a);
		CHECK(st.occupations()[5] == 3.779107816290222e-33_a);

  }

  SECTION("Spin polarized"){
    
    states::ks_states st(states::ks_states::spin_config::POLARIZED, 11.0);

		CHECK(st.num_electrons() == 11.0);    
    CHECK(st.num_states() == 6);
    CHECK(st.num_quantum_numbers() == 2);

  }

  SECTION("Non-collinear spin"){
    
    states::ks_states st(states::ks_states::spin_config::NON_COLLINEAR, 11.0);

		CHECK(st.num_electrons() == 11.0);
    CHECK(st.num_states() == 11);
    CHECK(st.num_quantum_numbers() == 1);

  }

  
}

#endif

#endif
