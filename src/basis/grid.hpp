/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__GRID
#define INQ__BASIS__GRID

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa

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


#include <cassert>
#include <array>

#include <basis/base.hpp>
#include <basis/double_grid.hpp>
#include <math/vector3.hpp>
#include <ions/unit_cell.hpp>

namespace inq {
namespace basis {

  class grid : public base {

  public:

		const static int dimension = 3;
		
		grid(const ions::unit_cell & cell, std::array<int, 3> nr, bool spherical_grid, bool double_grid, int periodic_dimensions, boost::mpi3::communicator & comm) :
			base(nr[0], comm),
			cubic_dist_({base::part_, inq::parallel::partition(nr[1]), inq::parallel::partition(nr[2])}),
			cell_(cell),
			nr_(nr),
			spherical_g_grid_(spherical_grid),
			periodic_dimensions_(periodic_dimensions),
			double_grid_(double_grid){

			if(base::part_.local_size() == 0){
				std::cerr << "\n  Partition " << comm.rank() << " has 0 points. Please change the number of processors.\n" << std::endl;
				comm.abort(1);
			}
			
			for(int idir = 0; idir < 3; idir++){
				rlength_[idir] = length(cell[idir]);
				ng_[idir] = nr_[idir];
				rspacing_[idir] = rlength_[idir]/nr_[idir];
				conspacing_[idir] = 1.0/nr_[idir];
				covspacing_[idir] = 2.0*M_PI;				
				glength_[idir] = 2.0*M_PI/rspacing_[idir];
				gspacing_[idir] = glength_[idir]/ng_[idir];
			}

			base::part_ *= nr_[1]*long(nr_[2]);
			
			npoints_ = nr_[0]*long(nr_[1])*nr_[2];

			for(int idir = 0; idir < 3; idir++) nr_local_[idir] = cubic_dist_[idir].local_size();
			
		}
		
		grid(grid && old, boost::mpi3::communicator & new_comm):
			grid(old.cell_, old.nr_, old.spherical_g_grid_, old.double_grid_.enabled(), old.periodic_dimensions_, new_comm)
		{
		}
		
    GPU_FUNCTION const math::vector3<double> & rspacing() const{
      return rspacing_;
    }

		double diagonal_length() const {
			return length(rlength_);
		}
		
    GPU_FUNCTION const math::vector3<double> & rlength() const{
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
    
		template<class OStream>
		friend OStream& operator<<(OStream& os, grid const& self){
			self.info(os);
			return os;
		}

		constexpr auto & cubic_dist(int dim) const {
			return cubic_dist_[dim];
		}

		GPU_FUNCTION static auto to_symmetric_range(std::array<int, 3> const & nr, const int ix, const int iy, const int iz) {
			math::vector3<int> ii{ix, iy, iz};
			for(int idir = 0; idir < 3; idir++) {
				if(ii[idir] >= (nr[idir] + 1)/2) ii[idir] -= nr[idir];
			}
			return ii;
		}

		GPU_FUNCTION static auto to_symmetric_range(std::array<int, 3> const & nr, parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) {
			return to_symmetric_range(nr, ix.value(), iy.value(), iz.value());
		}

		GPU_FUNCTION static auto from_symmetric_range(std::array<int, 3> const & nr, math::vector3<int> ii) {
			for(int idir = 0; idir < 3; idir++) {
				if(ii[idir] < 0) ii[idir] += nr[idir];
			}
			return ii;
		}

		GPU_FUNCTION auto to_symmetric_range(const int ix, const int iy, const int iz) const {
			return to_symmetric_range(nr_, ix, iy, iz);
		}

		GPU_FUNCTION auto to_symmetric_range(parallel::global_index ix, parallel::global_index iy, parallel::global_index iz) const {
			return to_symmetric_range(nr_, ix.value(), iy.value(), iz.value());
		}
			
		GPU_FUNCTION auto from_symmetric_range(math::vector3<int> ii) const {
			return from_symmetric_range(nr_, ii);
		}

		auto symmetric_range_begin(int idir) const {
			return -nr_[idir]/2;
		}

		auto symmetric_range_end(int idir) const {
			return nr_[idir]/2 + nr_[idir]%2;
		}

		auto linear_index(int const ix, int const iy, const int iz) const {
			return (ix*nr_[1] + iy)*nr_[2] + iz;
		}

		auto local_contains(math::vector3<int> const & ii) const {
			bool contains = true;
			for(int idir = 0; idir < 3; idir++){
				contains = contains and cubic_dist_[idir].contains(ii[idir]);
			}
			return contains;
		}

		auto & double_grid() const {
			return double_grid_;
		}

		auto & cell() const {
			return cell_;
		}		

	protected:

		std::array<inq::parallel::partition, 3> cubic_dist_;
		ions::unit_cell cell_;

    std::array<int, 3> nr_;

		std::array<int, 3> nr_local_;
		
    std::array<int, 3> ng_;

    math::vector3<double> rspacing_;
    math::vector3<double, math::contravariant> conspacing_;
    math::vector3<double> gspacing_;
    math::vector3<double, math::covariant> covspacing_;
		
    math::vector3<double> rlength_;
    math::vector3<double> glength_;

		long npoints_;

		bool spherical_g_grid_;

		int periodic_dimensions_;

		basis::double_grid double_grid_;

  };

}
}

#ifdef INQ_BASIS_GRID_UNIT_TEST
#undef INQ_BASIS_GRID_UNIT_TEST

#include <ions/geometry.hpp>
#include <catch2/catch_all.hpp>

TEST_CASE("class basis::grid", "[basis::grid]") {

	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;

  ions::unit_cell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 4.0, 0.0), vector3<double>(0.0, 0.0, 7.0));

	auto comm = boost::mpi3::environment::get_world_instance();

	basis::grid gr(cell, {120, 45, 77}, true, false, 3, comm);

	CHECK(gr.sizes()[0] == 120);
	CHECK(gr.sizes()[1] == 45);
	CHECK(gr.sizes()[2] == 77);

	if(comm.size() == 1) CHECK(gr.cubic_dist(0).local_size() == 120);
	if(comm.size() == 2) CHECK(gr.cubic_dist(0).local_size() == 60);
	if(comm.size() == 3) CHECK(gr.cubic_dist(0).local_size() == 40);
	if(comm.size() == 4) CHECK(gr.cubic_dist(0).local_size() == 30);
	if(comm.size() == 5) CHECK(gr.cubic_dist(0).local_size() == 24);
	if(comm.size() == 6) CHECK(gr.cubic_dist(0).local_size() == 20);
	if(comm.size() == 8) CHECK(gr.cubic_dist(0).local_size() == 15);
	if(comm.size() == 10) CHECK(gr.cubic_dist(0).local_size() == 12);
	if(comm.size() == 12) CHECK(gr.cubic_dist(0).local_size() == 10);
	CHECK(gr.cubic_dist(1).local_size() == 45);
	CHECK(gr.cubic_dist(2).local_size() == 77);

	for(int ix = 0; ix < gr.local_sizes()[0]; ix++){
		for(int iy = 0; iy < gr.local_sizes()[1]; iy++){
			for(int iz = 0; iz < gr.local_sizes()[2]; iz++){

				auto ii = gr.from_symmetric_range(gr.to_symmetric_range(ix, iy, iz));

				CHECK(ii[0] == ix);
				CHECK(ii[1] == iy);
				CHECK(ii[2] == iz);
				
			}
		}
	}

	CHECK(gr.symmetric_range_begin(0) == -60);
	CHECK(gr.symmetric_range_end(0) == 60);

	CHECK(gr.symmetric_range_begin(1) == -22);
	CHECK(gr.symmetric_range_end(1) == 23);
	
	CHECK(gr.symmetric_range_begin(2) == -38);
	CHECK(gr.symmetric_range_end(2) == 39);	

	CHECK(gr.double_grid().enabled() == false);

	auto new_gr = basis::grid(basis::grid(gr), boost::mpi3::environment::get_self_instance());

	CHECK(gr.sizes() == new_gr.sizes());
	CHECK(new_gr.local_sizes() == new_gr.sizes());	
	CHECK(gr.periodic_dimensions() == new_gr.periodic_dimensions());
	
}
#endif

    
#endif
