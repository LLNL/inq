/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__GROUND_STATE__EIGENVALUE_OUTPUT
#define INQ__GROUND_STATE__EIGENVALUE_OUTPUT

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cfloat>
#include <systems/electrons.hpp>

#include<tinyformat/tinyformat.h>

#include <utils/profiling.hpp>

namespace inq {
namespace ground_state {

class eigenvalues_output {

	int nspin_;
	int nkpoints_;
	gpu::array<int, 1> all_kpoint_index;
	gpu::array<int, 1> all_spin_index;
	gpu::array<int, 1> all_states_index;
	gpu::array<double, 1> all_eigenvalues;
	gpu::array<double, 1> all_occupations;
	gpu::array<complex, 1> all_normres;

public:
	
	template <typename NormResType>
	eigenvalues_output(systems::electrons const & el, NormResType const & normres):
		nspin_(el.states().num_spin_indices()),
		nkpoints_(el.brillouin_zone().size())
	{
		
		gpu::array<int, 2> kpoint_index({el.kpin_part().local_size(), el.max_local_set_size()});
		gpu::array<int, 2> spin_index({el.kpin_part().local_size(), el.max_local_set_size()});
		gpu::array<int, 2> state_index({el.kpin_part().local_size(), el.max_local_set_size()});
		gpu::array<double, 2> occs({el.kpin_part().local_size(), el.max_local_set_size()});
		
		auto iphi = 0;
		for(auto & phi : el.kpin()){
			auto ik = el.kpoint_index(phi);
			for(int ist = 0; ist < el.max_local_set_size(); ist++){
				kpoint_index[iphi][ist] = ik;
				spin_index[iphi][ist] = phi.spin_index();
				state_index[iphi][ist] = phi.set_part().local_to_global(ist).value();
			occs[iphi][ist] = 0.0;
			if(fabs(el.kpin_weights()[iphi]) > 1e-14) occs[iphi][ist] = el.occupations()[iphi][ist]/el.kpin_weights()[iphi];
			}
			iphi++;
		}
		
		all_kpoint_index = parallel::gather(+kpoint_index.flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		all_spin_index = parallel::gather(+spin_index.flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		all_states_index = parallel::gather(+state_index.flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		all_eigenvalues = parallel::gather(+el.eigenvalues().flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		all_occupations = parallel::gather(+occs.flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		all_normres = parallel::gather(+normres.flatted(), el.kpin_states_part(), el.kpin_states_comm(), 0);
		
	}

	static std::string spin_string(int index){
		if(index == 0) return "\u21D1";
		return "\u21D3";
	}

	template <class OStream>
	friend OStream& operator<<(OStream & out, eigenvalues_output const & self){

		auto order = gpu::array<int, 1>(self.all_eigenvalues().size());
		
		for(int iorder = 0; iorder < order.size(); iorder++) order[iorder] = iorder;
		
		std::sort(order.begin(), order.end(), [evs = self.all_eigenvalues](auto io, auto jo){ return evs[io] < evs[jo]; });
		
		auto const print_range = 8;
		
		//define the LUMO as the first state with ocupation below 0.1 (this is only for output purposes)
		auto lumo_index = order.size() - 1;
		for(int iorder = 0; iorder < order.size(); iorder++) {
			if(self.all_occupations[order[iorder]] >= 0.1) continue;
			lumo_index = iorder;
			break;
		}
		
		int skipped = 0;
		double minres = 1000.0;
		double maxres = 0.0;
		for(int iorder = 0; iorder < order.size(); iorder++) {
			auto ieig = order[iorder];
			auto print = (iorder < print_range) or (abs(iorder - lumo_index) < print_range) or (order.size() - 1 - iorder < print_range);
			
			if(not print) {
				skipped++;
				minres = std::min(minres, fabs(self.all_normres[ieig]));
				maxres = std::max(minres, fabs(self.all_normres[ieig]));
				continue;
			}
			
			if(skipped > 0) {
				tfm::format(out, "    [output of %d eigenvalues suppressed,  minres = %5.0e  maxres = %5.0e]\n", skipped, minres, maxres);
				skipped = 0;
				minres = 1000.0;
				maxres = 0.0;			
			}
			
			if(self.nkpoints_ > 1) tfm::format(out, "  kpt = %4d", self.all_kpoint_index[ieig] + 1);
			if(self.nspin_    > 1) tfm::format(out, "  spin = %s", spin_string(self.all_spin_index[ieig]));
			tfm::format(out, "  st = %4d  occ = %4.3f  evalue = %18.12f  res = %5.0e\n",
									self.all_states_index[ieig] + 1, self.all_occupations[ieig], real(self.all_eigenvalues[ieig]), fabs(self.all_normres[ieig]));
		}

		return out;
	}
};

}
}
#endif

#ifdef INQ_GROUND_STATE_EIGENVALUE_OUTPUT_UNIT_TEST
#undef INQ_GROUND_STATE_EIGENVALUE_OUTPUT_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
}
#endif

