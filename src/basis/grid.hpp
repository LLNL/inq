/* -*- indent-tabs-mode: t -*- */

#ifndef PLANE_WAVE_HPP
#define PLANE_WAVE_HPP

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

#include "../math/vec3d.hpp"
#include <cassert>
#include <array>

#include <basis/base.hpp>

namespace basis {

  class grid : public base {

  public:

		const static int dimension = 3;
		
		grid(const ions::UnitCell & cell, std::array<int, 3> nr, bool spherical_grid, int periodic_dimensions, boost::mpi3::communicator & comm) :
			base(nr[0], comm),
			cubic_dist_({base::dist_, utils::distribution(nr[1]), utils::distribution(nr[2])}),
			cell_(cell),
			nr_(nr),
			spherical_g_grid_(spherical_grid),
			periodic_dimensions_(periodic_dimensions){
			
			for(int idir = 0; idir < 3; idir++){
				rlength_[idir] = length(cell[idir]);
				ng_[idir] = nr_[idir];
				rspacing_[idir] = rlength_[idir]/nr_[idir];
				glength_[idir] = 2.0*M_PI/rspacing_[idir];
				gspacing_[idir] = glength_[idir]/ng_[idir];
			}

			base::dist_ *= nr_[1]*long(nr_[2]);
			
			npoints_ = nr_[0]*long(nr_[1])*nr_[2];

			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();
			
		}

    GPU_FUNCTION const math::vec3d & rspacing() const{
      return rspacing_;
    }

		double diagonal_length() const {
			return length(rlength_);
		}
		
    GPU_FUNCTION const math::vec3d & rlength() const{
      return rlength_;
    }

    double min_rlength() const{
      return std::min(rlength_[0], std::min(rlength_[1], rlength_[2]));
    }

		
    long size() const {
      return npoints_;
    }

		long num_points() const {
			return npoints_;
		}

		GPU_FUNCTION
		friend auto sizes(const grid & gr){
			return gr.nr_;
		}

		GPU_FUNCTION
		auto & sizes() const {
			return nr_;
		}

		auto & local_sizes() const {
			return nr_local_;
		}
		
		auto periodic_dimensions() const {
			return periodic_dimensions_;
		}

		template <class output_stream>
    void info(output_stream & out) const {
      out << "PLANE WAVE BASIS SET:" << std::endl;
      out << "  Grid size   = " << sizes()[0] << " x " << sizes()[1] << " x " << sizes()[2] << std::endl;
			out << "  Spacing [b] = " << rspacing() << std::endl;
			out << std::endl;
    }

		auto & cubic_dist(int dim) const {
			return cubic_dist_[dim];
		}
		
	protected:

		std::array<utils::distribution, 3> cubic_dist_;
		ions::UnitCell cell_;

    std::array<int, 3> nr_;

		std::array<int, 3> nr_local_;
		
    std::array<int, 3> ng_;

    math::vec3d rspacing_;
    math::vec3d gspacing_;
    
    math::vec3d rlength_;
    math::vec3d glength_;

		long npoints_;

		bool spherical_g_grid_;

		int periodic_dimensions_;
		
  };
}

#ifdef UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>

TEST_CASE("class basis::grid", "[basis::grid]") {
  
  using namespace Catch::literals;
  using math::vec3d;

  ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 4.0, 0.0), vec3d(0.0, 0.0, 7.0));

	auto comm = boost::mpi3::environment::get_world_instance();
	
	basis::grid gr(cell, {120, 45, 77}, true, 3, comm);

	REQUIRE(gr.sizes()[0] == 120);
	REQUIRE(gr.sizes()[1] == 45);
	REQUIRE(gr.sizes()[2] == 77);

	if(comm.size() == 1) REQUIRE(gr.cubic_dist(0).local_size() == 120);
	if(comm.size() == 2) REQUIRE(gr.cubic_dist(0).local_size() == 60);
	if(comm.size() == 3) REQUIRE(gr.cubic_dist(0).local_size() == 40);
	if(comm.size() == 4) REQUIRE(gr.cubic_dist(0).local_size() == 30);
	if(comm.size() == 5) REQUIRE(gr.cubic_dist(0).local_size() == 24);
	if(comm.size() == 6) REQUIRE(gr.cubic_dist(0).local_size() == 20);
	if(comm.size() == 8) REQUIRE(gr.cubic_dist(0).local_size() == 15);
	if(comm.size() == 10) REQUIRE(gr.cubic_dist(0).local_size() == 12);
	if(comm.size() == 12) REQUIRE(gr.cubic_dist(0).local_size() == 10);
	REQUIRE(gr.cubic_dist(1).local_size() == 45);
	REQUIRE(gr.cubic_dist(2).local_size() == 77);
	
}
#endif

    
#endif
