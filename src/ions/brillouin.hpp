/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__BRILLOUIN
#define INQ__IONS__BRILLOUIN

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <input/kpoints.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace ions {

class brillouin {

	std::vector<vector3<double, covariant>> kpoints_;
	std::vector<double> weights_;	
	
public:
  
  brillouin(inq::systems::ions const & ions, input::kpoints::grid const & kpts):
		kpoints_(kpts.size()),
		weights_(kpts.size())
  {
		auto num_atoms = std::max(1, ions.geo().num_atoms());
		
		std::vector<int> types(num_atoms);
		std::vector<double> positions(3*num_atoms);
		
		//add a dummy atom, since spg doesn't work without atoms
		types[0] = 0;
		positions[0] = 0.0;
		positions[1] = 0.0;
		positions[2] = 0.0;		
		
		for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++){
			types[iatom] = ions.geo().atoms()[iatom].atomic_number();
			auto pos = ions.cell().metric().to_contravariant(ions.cell().position_in_cell(ions.geo().coordinates()[iatom]));
			positions[3*iatom + 0] = pos[0];
			positions[3*iatom + 1] = pos[1];
			positions[3*iatom + 2] = pos[2];
		}

		double amat[9];
		amat[0] = ions.cell().lattice(0)[0];
		amat[1] = ions.cell().lattice(0)[1];
		amat[2] = ions.cell().lattice(0)[2];
		amat[3] = ions.cell().lattice(1)[0];
		amat[4] = ions.cell().lattice(1)[1];
		amat[5] = ions.cell().lattice(1)[2];
		amat[6] = ions.cell().lattice(2)[0];
		amat[7] = ions.cell().lattice(2)[1];
		amat[8] = ions.cell().lattice(2)[2];

		auto is_shifted = kpts.is_shifted();
		auto grid_address = std::vector<int>(3*kpts.size());
    auto map = std::vector<int>(kpts.size());
		
		spg_get_ir_reciprocal_mesh(reinterpret_cast<int (*)[3]>(grid_address.data()), map.data(), (int const *) &kpts.dims(), (int const *) &is_shifted, 0,
															 reinterpret_cast<double (*)[3]>(amat), reinterpret_cast<double (*)[3]>(positions.data()), types.data(), num_atoms, 1e-4);
		
		for(int ik = 0; ik < kpts.size(); ik++){
			auto kpr = vector3<double, covariant>{
				(grid_address[3*ik + 0] + 0.5*is_shifted[0])/kpts.dims()[0],
				(grid_address[3*ik + 1] + 0.5*is_shifted[1])/kpts.dims()[1],
				(grid_address[3*ik + 2] + 0.5*is_shifted[2])/kpts.dims()[2]};
			kpr.transform([](auto xx){ return (xx >= 0.5) ? xx - 1.0 : xx; });
			kpoints_[ik] = 2.0*M_PI*kpr;
			weights_[ik] = 1.0/kpts.size();
		}
		
  }
	
  brillouin(inq::systems::ions const & ions, input::kpoints::list const & kpts):
		kpoints_(kpts.size()),
		weights_(kpts.size())
  {

		double totalw = 0.0;
		for(int ik = 0; ik < kpts.size(); ik++){
			auto kpr = kpts.kpoint(ik);
			kpr.transform([](auto xx){ return (xx >= 0.5) ? xx - 1.0 : xx; });
			kpoints_[ik] = 2.0*M_PI*kpr;
			weights_[ik] = kpts.weight(ik);
			totalw += weights_[ik];
		}

		if(totalw < 1.0e-15) throw std::runtime_error("inq error: the k-points have zero total weight");

		for(int ik = 0; ik < kpts.size(); ik++) weights_[ik] /= totalw;

  }
	
  auto size() const {
    return (long) kpoints_.size();
  }

  auto kpoint(int ik) const {
		return kpoints_[ik];
  }
  
  auto kpoint_weight(int ik) const {
    return weights_[ik];
  }
	
	template<class OStream>
	friend OStream& operator<<(OStream& os, brillouin const & self){
		os << std::endl;		
		os << "  Number of kpoints = " << self.size() << std::endl;
		for(int ikpt = 0; ikpt < self.size(); ikpt++){
			auto kk = self.kpoint(ikpt)/(2.0*M_PI);
			tfm::format(os, "  k-point %7d = %7.3f %7.3f %7.3f   weight = %5.3f\n", ikpt, kk[0], kk[1], kk[2], self.kpoint_weight(ikpt));
		}
		os << std::endl;
		return os;
	}
};

}
}
#endif

#ifdef INQ_IONS_BRILLOUIN_UNIT_TEST
#undef INQ_IONS_BRILLOUIN_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using Catch::Approx;
	using namespace Catch::literals;

	SECTION("Diamond"){

		auto a =  3.567095_A;

		auto box = systems::box::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b});
		auto ions = systems::ions(box);
		
		ions.insert_fractional("C", {0.0 , 0.0 , 0.0 });
		ions.insert_fractional("C", {0.25, 0.25, 0.25});

		auto bz1 = ions::brillouin(ions, input::kpoints::gamma());

		CHECK(bz1.size() == 1);
		CHECK(bz1.kpoint(0)[0] == 0.0_a);
		CHECK(bz1.kpoint(0)[1] == 0.0_a);
		CHECK(bz1.kpoint(0)[2] == 0.0_a);

		auto bz2 = ions::brillouin(ions, input::kpoints::grid({1, 2, 3}));

		CHECK(bz2.size() == 6);

		CHECK(bz2.kpoint(0)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(0)[1]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(0)[2]/(2*M_PI) == 0.0_a);
		
		CHECK(bz2.kpoint(1)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(1)[1]/(2*M_PI) == -0.5_a);
		CHECK(bz2.kpoint(1)[2]/(2*M_PI) == 0.0_a);
		
		CHECK(bz2.kpoint(2)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(2)[1]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(2)[2]/(2*M_PI) == 0.3333333333_a);
		
		CHECK(bz2.kpoint(3)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(3)[1]/(2*M_PI) == -0.5_a);
		CHECK(bz2.kpoint(3)[2]/(2*M_PI) == 0.3333333333_a);
		
		CHECK(bz2.kpoint(4)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(4)[1]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(4)[2]/(2*M_PI) == -0.3333333333_a);
		
		CHECK(bz2.kpoint(5)[0]/(2*M_PI) == 0.0_a);
		CHECK(bz2.kpoint(5)[1]/(2*M_PI) == -0.5_a);
		CHECK(bz2.kpoint(5)[2]/(2*M_PI) == -0.3333333333_a);

		auto bz3 = ions::brillouin(ions, input::kpoints::grid({2, 2, 2}, true));

		CHECK(bz3.size() == 8);

		CHECK(bz3.kpoint(0)[0]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(0)[1]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(0)[2]/(2*M_PI) ==  0.25_a);
		
		CHECK(bz3.kpoint(1)[0]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(1)[1]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(1)[2]/(2*M_PI) ==  0.25_a);
		
		CHECK(bz3.kpoint(2)[0]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(2)[1]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(2)[2]/(2*M_PI) ==  0.25_a);
		
		CHECK(bz3.kpoint(3)[0]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(3)[1]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(3)[2]/(2*M_PI) ==  0.25_a);
		
		CHECK(bz3.kpoint(4)[0]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(4)[1]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(4)[2]/(2*M_PI) == -0.25_a);
		
		CHECK(bz3.kpoint(5)[0]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(5)[1]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(5)[2]/(2*M_PI) == -0.25_a);
		
		CHECK(bz3.kpoint(6)[0]/(2*M_PI) ==  0.25_a);
		CHECK(bz3.kpoint(6)[1]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(6)[2]/(2*M_PI) == -0.25_a);
		
		CHECK(bz3.kpoint(7)[0]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(7)[1]/(2*M_PI) == -0.25_a);
		CHECK(bz3.kpoint(7)[2]/(2*M_PI) == -0.25_a);

		auto kpts = input::kpoints::list();
		
		kpts.insert({0.5, 0.5, -0.5}, 1.0);
		kpts.insert({0.1, 0.2,  0.3}, 2.0);
		kpts.insert({1.0, -0.01, 0.0}, 0.0);
		
		auto bz4 = ions::brillouin(ions, kpts);

		CHECK(bz4.size() == 3);

		CHECK(bz4.kpoint(0)[0]/(2*M_PI) == -0.5_a);
		CHECK(bz4.kpoint(0)[1]/(2*M_PI) == -0.5_a);
		CHECK(bz4.kpoint(0)[2]/(2*M_PI) == -0.5_a);
		
		CHECK(bz4.kpoint(1)[0]/(2*M_PI) ==  0.1_a);
		CHECK(bz4.kpoint(1)[1]/(2*M_PI) ==  0.2_a);
		CHECK(bz4.kpoint(1)[2]/(2*M_PI) ==  0.3_a);
		
		CHECK(bz4.kpoint(2)[0]/(2*M_PI) ==  0.0_a);
		CHECK(bz4.kpoint(2)[1]/(2*M_PI) == -0.01_a);
		CHECK(bz4.kpoint(2)[2]/(2*M_PI) ==  0.0_a);

		CHECK(bz4.kpoint_weight(0) ==  0.33333333_a);
		CHECK(bz4.kpoint_weight(1) ==  0.66666666_a);
		CHECK(bz4.kpoint_weight(2) ==  0.0_a);
		
	}

}
#endif
