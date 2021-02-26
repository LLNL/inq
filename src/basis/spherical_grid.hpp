/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__SPHERICAL_GRID
#define INQ__BASIS__SPHERICAL_GRID

/*
 Copyright (C) 2019-2020 Xavier Andrade, Alfredo A. Correa.

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

#include <algorithm> //max, min

#include <math/vector3.hpp>
#include <gpu/atomic.hpp>
#include <gpu/run.hpp>
#include <gpu/reduce.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <cassert>
#include <array>
#include <math/array.hpp>

#include <utils/profiling.hpp>

namespace inq {
namespace basis {

  class spherical_grid {

		//returns the cube that contains the sphere, this makes the initialization O(1) instead of O(N)
		template <class BasisType, typename PosType>
		void static cube(const BasisType & parent_grid, PosType const & pos, double radius, math::vector3<int> & hi, math::vector3<int> & lo){
			for(int idir = 0; idir < 3; idir++){
				lo[idir] = floor((pos[idir] - radius)/parent_grid.rspacing()[idir]) - 1;
				hi[idir] = ceil((pos[idir] + radius)/parent_grid.rspacing()[idir]) + 1;
				
				lo[idir] = std::max<int>(parent_grid.symmetric_range_begin(idir), lo[idir]);
				hi[idir] = std::max<int>(parent_grid.symmetric_range_begin(idir), hi[idir]);
				
				lo[idir] = std::min<int>(parent_grid.symmetric_range_end(idir), lo[idir]);
				hi[idir] = std::min<int>(parent_grid.symmetric_range_end(idir), hi[idir]);
			}
		}

#ifdef ENABLE_CUDA
	public:
#endif
		struct point_data {
			math::vector3<int> coords_;
			float distance_; //I don't think we need additional precision for this, and we get aligned memory
			math::vector3<double> relative_pos_;
		};
		
  public:

		//we need to make an additional public function to make cuda happy
		template <class basis>
		void initialize(const basis & parent_grid, const ions::UnitCell & cell, const math::vector3<double> & center_point, const double radius){
					CALI_CXX_MARK_FUNCTION;
	
      ions::periodic_replicas rep(cell, center_point, parent_grid.diagonal_length());

			long count = 0;

			//FIRST PASS: we count the number of points for allocation
			for(unsigned irep = 0; irep < rep.size(); irep++){

				math::vector3<int> lo, hi;
				cube(parent_grid, rep[irep], radius, hi, lo);
				
				//OPTIMIZATION: this iteration should be done only over the local points
				math::vector3<int> local_sizes = parent_grid.local_sizes();

				count += gpu::run(gpu::reduce(hi[2] - lo[2]), gpu::reduce(hi[1] - lo[1]), gpu::reduce(hi[0] - lo[0]),
													[lo, local_sizes, point_op = parent_grid.point_op(), re = rep[irep], radius] GPU_LAMBDA (auto iz, auto iy, auto ix){
														
														auto ii = point_op.from_symmetric_range({int(lo[0] + ix), int(lo[1] + iy), int(lo[2] + iz)});
														
														utils::global_index ii0(ii[0]);
														utils::global_index ii1(ii[1]);
														utils::global_index ii2(ii[2]);
														
														int ixl = point_op.cubic_dist()[0].global_to_local(ii0);
														int iyl = point_op.cubic_dist()[1].global_to_local(ii1);
														int izl = point_op.cubic_dist()[2].global_to_local(ii2);
														
														if(ixl < 0 or ixl >= local_sizes[0]) return 0;
														if(iyl < 0 or iyl >= local_sizes[1]) return 0;
														if(izl < 0 or izl >= local_sizes[2]) return 0;
														
														auto rpoint = point_op.rvector(ii0, ii1, ii2);
														
														auto n2 = norm(rpoint - re);
														if(n2 > radius*radius) return 0;
														
														return 1;
													});
				
			}

			
			points_.reextent({count});

			if(count == 0) return;

			math::array<long, 1> destination_count(1, long(0));
			assert(destination_count[0] == 0);

			//SECOND PASS: we generate the list of points
			for(unsigned irep = 0; irep < rep.size(); irep++){

				math::vector3<int> lo, hi;
				cube(parent_grid, rep[irep], radius, hi, lo);

				math::vector3<int> local_sizes = parent_grid.local_sizes();
				
				gpu::run(hi[2] - lo[2], hi[1] - lo[1], hi[0] - lo[0],
								 [lo, local_sizes, point_op = parent_grid.point_op(), re = rep[irep], radius,
									des = begin(destination_count), poi = begin(points_)] GPU_LAMBDA (auto iz, auto iy, auto ix){
									 
									 auto ii = point_op.from_symmetric_range({int(lo[0] + ix), int(lo[1] + iy), int(lo[2] + iz)});
									 
									 utils::global_index ii0(ii[0]);
									 utils::global_index ii1(ii[1]);
									 utils::global_index ii2(ii[2]);
									 
									 int ixl = point_op.cubic_dist()[0].global_to_local(ii0);
									 int iyl = point_op.cubic_dist()[1].global_to_local(ii1);
									 int izl = point_op.cubic_dist()[2].global_to_local(ii2);
									 
									 if(ixl < 0 or ixl >= local_sizes[0]) return;
									 if(iyl < 0 or iyl >= local_sizes[1]) return;
									 if(izl < 0 or izl >= local_sizes[2]) return;
									 
									 auto rpoint = point_op.rvector(ii0, ii1, ii2);
									 
									 auto n2 = norm(rpoint - re);
									 if(n2 > radius*radius) return;
									 
									 auto dest = gpu::atomic::add(&des[0], 1);
									 
									 poi[dest].coords_ = {ixl, iyl, izl};
									 poi[dest].distance_ = sqrt(n2);
									 poi[dest].relative_pos_ = rpoint - re;
									 
								 });
			}

			//OPTIMIZATION: order the points for better memory access

		}
		
		
		const static int dimension = 1;
		
		template <class basis>
    spherical_grid(const basis & parent_grid, const ions::UnitCell & cell, const math::vector3<double> & center_point, const double radius):
			volume_element_(parent_grid.volume_element()),
			center_(center_point){

			initialize(parent_grid, cell, center_point, radius);
    }

		auto create_comm(boost::mpi3::communicator & comm) const {
			auto color = MPI_UNDEFINED;
			if(size() != 0) color = 1;
			return comm.split(color, 0);
		}
		
    long size() const {
      return points_.size();
    }
    
    template <class array_3d, class array_1d>
    void gather(const array_3d & grid, array_1d && subgrid) const {

			CALI_CXX_MARK_SCOPE("spherical_grid::gather(3d)");
				
			//DATAOPERATIONS STL TRANSFORM
			std::transform(points_.begin(), points_.end(), subgrid.begin(),
										 [& grid](auto point){
											 return grid[point.coords_[0]][point.coords_[1]][point.coords_[2]];
										 });
		}

		template <class array_4d>
		math::array<typename array_4d::element, 2> gather(const array_4d & grid) const {

			CALI_CXX_MARK_SCOPE("spherical_grid::gather(4d)");
			
			const int nst = std::get<3>(sizes(grid));

			CALI_MARK_BEGIN("spherical_grid::gather(4d)::allocation");
			math::array<typename array_4d::element, 2> subgrid({this->size(), nst});
			CALI_MARK_END("spherical_grid::gather(4d)::allocation");

			math::prefetch(subgrid);
			
			gpu::run(nst, size(),
							 [sgr = begin(subgrid), gr = begin(grid), poi = begin(points_)] GPU_LAMBDA (auto ist, auto ipoint){
								 sgr[ipoint][ist] = gr[poi[ipoint].coords_[0]][poi[ipoint].coords_[1]][poi[ipoint].coords_[2]][ist];
							 });
			
			return subgrid;
    }

    template <class array_2d, class array_4d>
    void scatter_add(const array_2d & subgrid, array_4d && grid) const{
			CALI_CXX_MARK_SCOPE("spherical_grid::scatter_add");

			gpu::run(std::get<1>(sizes(subgrid)), size(),
							 [sgr = begin(subgrid), gr = begin(grid), poi = begin(points_)] GPU_LAMBDA (auto ist, auto ipoint){
								 gr[poi[ipoint].coords_[0]][poi[ipoint].coords_[1]][poi[ipoint].coords_[2]][ist] += sgr[ipoint][ist];
							 });
    }
    
    template <class array_1d, class array_3d>
    void scatter(const array_1d & subgrid, array_3d && grid) const{

			CALI_CXX_MARK_SCOPE("spherical_grid::scatter");
			
      for(int ipoint = 0; ipoint < size(); ipoint++){
				grid[points_[ipoint].coords_[0]][points_[ipoint].coords_[1]][points_[ipoint].coords_[2]] = subgrid[ipoint];
      }
    }

		const double & volume_element() const {
			return volume_element_;
		}

		auto & points(int ii) const {
			return points_[ii].coords_;
		}
		
		auto & distance(int ii) const {
			return points_[ii].distance_;
		}
		
		auto & point_pos(int ii) const {
			return points_[ii].relative_pos_;
		}
		
		friend auto sizes(const spherical_grid & sphere){
			return std::array<long, dimension>{sphere.size()};
		}

		auto & center() const {
			return center_;
		}

		template <typename PointsType>
		struct sphere_ref {

			PointsType points_;

			GPU_FUNCTION auto & points(int ii) const {
				return points_[ii].coords_;
			}

			GPU_FUNCTION auto & distance(int ii) const {
				return points_[ii].distance_;
			}

			GPU_FUNCTION auto & point_pos(int ii) const {
				return points_[ii].relative_pos_;
			}
			
		};

		auto ref() const {
			return sphere_ref<decltype(cbegin(points_))>{cbegin(points_)};
		}		
		
  private:

		math::array<point_data, 1> points_;
		double volume_element_;
		math::vector3<double> center_;
		
  };

}
}

#ifdef INQ_BASIS_SPHERICAL_GRID_UNIT_TEST
#undef INQ_BASIS_SPHERICAL_GRID_UNIT_TEST

#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>
#include <math/array.hpp>
#include <math/complex.hpp>

TEST_CASE("class basis::spherical_grid", "[basis::spherical_grid]") {
	
	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;

	auto comm = boost::mpi3::environment::get_world_instance();
 
  double ll = 10.0;
  
  ions::UnitCell cell(vector3<double>(ll, 0.0, 0.0), vector3<double>(0.0, ll, 0.0), vector3<double>(0.0, 0.0, ll));
  
  double ecut = 20.0;
  
  basis::real_space pw(cell, input::basis::cutoff_energy(ecut), comm);
  
  SECTION("Point 0 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {0.0, 0.0, 0.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);

    math::array<complex, 3> grid({pw.local_sizes()[0], pw.local_sizes()[1], pw.local_sizes()[2]});
    std::vector<complex> subgrid(sphere.size());

    for(long ii = 0; ii < grid.num_elements(); ii++) grid.data_elements()[ii] = 0.0;
    
    sphere.gather(grid, subgrid);

    for(unsigned ii = 0; ii < subgrid.size(); ii++) subgrid[ii] = 1.0; 
    
    sphere.scatter(subgrid, grid);

    double sum = 0.0;
	for(long ii = 0; ii < grid.num_elements(); ii++) sum += real(grid.data_elements()[ii]);
	comm.all_reduce_in_place_n(&sum, 1, std::plus<>{});
	
    CHECK(sum == 257.0_a);
    
  }

  SECTION("Point -l/2 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {-ll/2.0, 0.0, 0.0}, 2.0);
		
		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);
    
    math::array<complex, 4> grid({pw.local_sizes()[0], pw.local_sizes()[1], pw.local_sizes()[2], 20}, 0.0);
    math::array<complex, 2> subgrid({sphere.size(), 20}, 0.0);

    for(long ii = 0; ii < grid.num_elements(); ii++) grid.data_elements()[ii] = 1.0;
    
    sphere.gather(grid, subgrid);

    double sum = 0.0;
	for(long ii = 0; ii < subgrid.num_elements(); ii++) sum += real(subgrid.data_elements()[ii]);
	comm.all_reduce_in_place_n(&sum, 1, std::plus<>{});

    CHECK(sum == Approx(20.0*257.0));
    
    for(long ii = 0; ii < subgrid.num_elements(); ii++) subgrid.data_elements()[ii] = 0.0;
    
	sphere.scatter(subgrid, grid);

    sum = 0.0;
	for(long ii = 0; ii < grid.num_elements(); ii++) sum += real(grid.data_elements()[ii]);
	comm.all_reduce_in_place_n(&sum, 1, std::plus<>{});

	CHECK(sum == Approx(20.0*(pw.size() - size)));

  }

  SECTION("Point l/2 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, 0.0, 0.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);

    math::array<complex, 6> grid({1, pw.local_sizes()[0], pw.local_sizes()[1], pw.local_sizes()[2], 2, 20}, 0.0);
    math::array<complex, 3> subgrid({sphere.size(), 2, 20}, 0.0);

    sphere.gather(grid[0], subgrid);

    sphere.scatter(subgrid, grid[0]);
    
  }

  SECTION("Point -l/2 -l/2 -l/2"){
    
    basis::spherical_grid sphere(pw, cell, {-ll/2.0, -ll/2.0, -ll/2.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);
    
  }

  SECTION("Point l/2 l/2 l/2"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, ll/2.0, ll/2.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);
    
  }

  SECTION("Point l/2 l/2 l/2"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, ll/2.0, ll/2.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);
    
  }

  SECTION("Point l/4 l/4 l/4"){
    
    basis::spherical_grid sphere(pw, cell, {ll/4.0, ll/4.0, ll/4.0}, 2.0);

		auto size = sphere.size();
		comm.all_reduce_in_place_n(&size, 1, std::plus<>{});
    CHECK(size == 257);
    
  }
  
}
#endif

#endif

