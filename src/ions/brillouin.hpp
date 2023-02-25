/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__BRILLOUIN
#define INQ__IONS__BRILLOUIN

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa.

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

#include <input/kpoints.hpp>
#include <systems/ions.hpp>

namespace inq {
namespace ions {

class brillouin {

  input::kpoints grid_;
  std::vector<int> grid_address_;
  std::vector<int> map_;
  vector3<int> is_shifted_;

	auto kpoint_raw(int ik) const {
    return vector3<double, covariant>{
			(grid_address_[3*ik + 0] + 0.5*is_shifted_[0])/grid_.dims()[0],
      (grid_address_[3*ik + 1] + 0.5*is_shifted_[1])/grid_.dims()[1],
      (grid_address_[3*ik + 2] + 0.5*is_shifted_[2])/grid_.dims()[2]};
	}
	
public:
  
  brillouin(inq::systems::ions const & ions, input::kpoints const & kpts):
    grid_(kpts),
    grid_address_(3*grid_.num()),
    map_(grid_.num()),
    is_shifted_(grid_.is_shifted())
  {

		std::vector<int> types(ions.geo().num_atoms());
		std::vector<double> positions(3*ions.geo().num_atoms());
		
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
		
		spg_get_ir_reciprocal_mesh(reinterpret_cast<int (*)[3]>(grid_address_.data()), map_.data(), (int const *) &grid_.dims(), (int const *) &is_shifted_, 0,
															 reinterpret_cast<double (*)[3]>(amat), reinterpret_cast<double (*)[3]>(positions.data()), types.data(), ions.geo().num_atoms(), 1e-4);

  }

  auto size() const {
    return grid_.num();
  }

  auto kpoint(int ik) const {
		auto kpr = kpoint_raw(ik);
		kpr.transform([](auto xx){ return (xx >= 0.5) ? xx - 1.0 : xx; });
    return 2.0*M_PI*kpr;
  }
  
  auto kpoint_weight(int ik) const {
    return 1.0/grid_.num();
  }
	
	template<class OStream>
	friend OStream& operator<<(OStream& os, brillouin const & self){
		os << std::endl;		
		os << "  Number of kpoints = " << self.size() << std::endl;
		for(int ikpt = 0; ikpt < self.size(); ikpt++){
			os << "  k-point\t" << ikpt << " =\t" << self.kpoint(ikpt)/(2.0*M_PI) << "\tweight =\t" << self.kpoint_weight(ikpt) << std::endl;
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

		auto box = systems::box::lattice({0.0_b, a/2.0, a/2.0}, {a/2, 0.0_b, a/2.0}, {a/2.0, a/2.0, 0.0_b}).cutoff_energy(35.0_Ha);
		auto ions = systems::ions(box);
		
		ions.insert("C", {0.0_crys,  0.0_crys,  0.0_crys });
		ions.insert("C", {0.25_crys, 0.25_crys, 0.25_crys});

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
		
	}

}
#endif
