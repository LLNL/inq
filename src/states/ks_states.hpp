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
#include <utils/partition.hpp>

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
        
	ks_states(const spin_config spin, const double nelectrons, const int extra_states = 0, double temperature = 0.0, int nkpoints = 1):
		temperature_(temperature),
		nkpoints_(nkpoints)
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
	
	auto num_electrons() const {
		return num_electrons_;
	}

	static auto smear_function(double xx)  {

		double const maxarg = 200.0;
		
		if(xx > maxarg){
			return 1.0;
		} else if(xx > -maxarg) {
			return 1.0/(1.0 + exp(-xx));
		} else {
			return 0.0;
		}
		
	}

	template <typename EigType, typename FunctionType>
	auto get_efermi(boost::mpi3::communicator & comm, double nelec, EigType const & eig, FunctionType function){

			
		int const nitmax = 200;
		double const tol = 1e-10;
		auto drange = sqrt(-log(tol*0.01));

		double emin = std::numeric_limits<double>::max();
		double emax = std::numeric_limits<double>::min();
		
		for(long ie = 0; ie < eig.size(); ie++){
			emin = std::min(emin, eig[ie]);
			emax = std::max(emax, eig[ie]);
		}
		
		emin = comm.all_reduce_value(emin, boost::mpi3::min<>{}) - drange;
		emax = comm.all_reduce_value(emax, boost::mpi3::max<>{}) + drange;
		
		double efermi;
		
		for(int iter = 0; iter < nitmax; iter++){
			efermi = 0.5*(emin + emax);

			double sumq = 0.0;
			for(long ie = 0; ie < eig.size(); ie++){
				sumq = sumq + max_occ_*function(efermi, real(eig[ie]));
			}
			comm.all_reduce_in_place_n(&sumq, 1, std::plus<>{});
			
			if(fabs(sumq - nelec) <= tol) break;
			if(sumq <= nelec) emin = efermi;
			if(sumq >= nelec) emax = efermi;
			
		}

		return efermi;		
	}
	
	template <typename EigenvalType, typename OccsType>
	void update_occupations(boost::mpi3::communicator & comm, utils::partition const & part, EigenvalType const & eigenval, OccsType & occs) {

		assert(part.local_size() == eigenval[0].size());
		assert(part.local_size() == occs[0].size());
		assert(sizes(eigenval) == sizes(occs));
		
		double efermi;
		
		auto feig = eigenval.flatted();
		auto focc = occs.flatted();
		
		if(temperature_ == 0.0){

			auto rem_electrons = num_electrons_;
			for(int ist = 0; ist < nstates_; ist++){
				auto occ = std::min(max_occ_, rem_electrons);
				if(part.contains(ist)) occs[0][ist - part.start()] = occ;
				rem_electrons -= occ;
			}
			
		} else {

			auto dsmear = std::max(1e-14, temperature_);

			auto func = [dsmear] (auto efermi, auto eig){
				return smear_function((efermi - eig)/dsmear);
			};
			
			efermi = get_efermi(comm, num_electrons_, feig, func);
			
			for(long ie = 0; ie < feig.size(); ie++){
				focc[ie] = max_occ_*func(efermi, feig[ie]);
			}
			
		}

		assert(fabs(comm.all_reduce_value(operations::sum(focc)) - num_electrons_) <= 1e-10);

		for(long ie = 0; ie < feig.size(); ie++) focc[ie] /= nkpoints_;
		
	}
	
private:

	double num_electrons_;
	int nstates_;
	int nquantumnumbers_;
	double temperature_;
	double max_occ_;
	int nkpoints_;

};

}
}

#ifdef INQ_STATES_KS_STATES_UNIT_TEST
#undef INQ_STATES_KS_STATES_UNIT_TEST

#include <catch2/catch.hpp>
#include <mpi3/environment.hpp>

TEST_CASE("Class states::ks_states", "[ks_states]"){

	using namespace Catch::literals;
	using namespace inq;

	auto comm = boost::mpi3::environment::get_world_instance();

  SECTION("Spin unpolarized"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

		CHECK(st.num_electrons() == 11.0);
    CHECK(st.num_states() == 6);
    CHECK(st.num_quantum_numbers() == 1);

		utils::partition part(st.num_states(), comm);

		math::array<double, 2> eigenvalues({1, part.local_size()});
		math::array<double, 2> occupations({1, part.local_size()});

		if(part.contains(0)) eigenvalues[0][part.global_to_local(utils::global_index(0))] = 0.1;
		if(part.contains(1)) eigenvalues[0][part.global_to_local(utils::global_index(1))] = 0.2;
		if(part.contains(2)) eigenvalues[0][part.global_to_local(utils::global_index(2))] = 0.3;
		if(part.contains(3)) eigenvalues[0][part.global_to_local(utils::global_index(3))] = 0.3;
		if(part.contains(4)) eigenvalues[0][part.global_to_local(utils::global_index(4))] = 0.4;
		if(part.contains(5)) eigenvalues[0][part.global_to_local(utils::global_index(5))] = 1.0;
		
		st.update_occupations(comm, part, eigenvalues, occupations);
		
		if(part.contains(0)) CHECK(occupations[0][part.global_to_local(utils::global_index(0))] == 2.0);
		if(part.contains(1)) CHECK(occupations[0][part.global_to_local(utils::global_index(1))] == 2.0);
		if(part.contains(2)) CHECK(occupations[0][part.global_to_local(utils::global_index(2))] == 2.0);
		if(part.contains(3)) CHECK(occupations[0][part.global_to_local(utils::global_index(3))] == 2.0);
		if(part.contains(4)) CHECK(occupations[0][part.global_to_local(utils::global_index(4))] == 2.0);
		if(part.contains(5)) CHECK(occupations[0][part.global_to_local(utils::global_index(5))] == 1.0);

  }
	
  SECTION("Spin unpolarized with temperature"){
    
    states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 4.0, 4, 0.01);

		CHECK(st.num_electrons() == 4.0);    
    CHECK(st.num_states() == 6);
    CHECK(st.num_quantum_numbers() == 1);

		comm.barrier();
		
		utils::partition part(st.num_states(), comm);
		
		math::array<double, 2> eigenvalues({1, part.local_size()});
		math::array<double, 2> occupations({1, part.local_size()});

		if(part.contains(0)) eigenvalues[0][part.global_to_local(utils::global_index(0))] = 0.0;
		if(part.contains(1)) eigenvalues[0][part.global_to_local(utils::global_index(1))] = 0.2;
		if(part.contains(2)) eigenvalues[0][part.global_to_local(utils::global_index(2))] = 0.3;
		if(part.contains(3)) eigenvalues[0][part.global_to_local(utils::global_index(3))] = 0.3;
		if(part.contains(4)) eigenvalues[0][part.global_to_local(utils::global_index(4))] = 0.4;
		if(part.contains(5)) eigenvalues[0][part.global_to_local(utils::global_index(5))] = 1.0;
		
		st.update_occupations(comm, part, eigenvalues, occupations);
		
		if(part.contains(0)) CHECK(occupations[0][part.global_to_local(utils::global_index(0))] == 2.0_a);
		if(part.contains(1)) CHECK(occupations[0][part.global_to_local(utils::global_index(1))] == 1.9810772793_a);
		if(part.contains(2)) CHECK(occupations[0][part.global_to_local(utils::global_index(2))] == 0.0094611446_a);
		if(part.contains(3)) CHECK(occupations[0][part.global_to_local(utils::global_index(3))] == 0.0094611446_a);
		if(part.contains(4)) CHECK(occupations[0][part.global_to_local(utils::global_index(4))] == 4.315768121820668e-07_a);
		if(part.contains(5)) CHECK(occupations[0][part.global_to_local(utils::global_index(5))] == 3.779107816290222e-33_a);
		
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
