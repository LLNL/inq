/* -*- indent-tabs-mode: t -*- */

#ifndef SPHERICAL_GRID_HPP
#define SPHERICAL_GRID_HPP

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

#include <math/d3vector.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <cassert>
#include <array>
#include <math/array.hpp>

namespace basis {

  class spherical_grid {

  public:

		const static int dimension = 1;
		
		template <class basis>
    spherical_grid(const basis & parent_grid, const ions::UnitCell & cell, const math::d3vector & center_point, const double radius):
			volume_element_(parent_grid.volume_element()){

      ions::periodic_replicas rep(cell, center_point, parent_grid.diagonal_length());

			std::vector<std::array<int, 3> > tmp_points;
					
			//DATAOPERATIONS LOOP 4D
      for(int ix = 0; ix < parent_grid.rsize()[0]; ix++){
				for(int iy = 0; iy < parent_grid.rsize()[1]; iy++){
					for(int iz = 0; iz < parent_grid.rsize()[2]; iz++){
						auto rpoint = parent_grid.rvector(ix, iy, iz);
						
						for(unsigned irep = 0; irep < rep.size(); irep++){
							auto n2 = norm(rpoint - rep[irep]);
							if(n2 > radius*radius) continue;
							
							tmp_points.push_back({ix, iy, iz});
							distance_.push_back(sqrt(n2));
							relative_pos_.push_back(rpoint - rep[irep]);
						}
						
					}
				}
      }

			points_.reextent({tmp_points.size()});
			
			for(unsigned ii = 0; ii < tmp_points.size(); ii++){
				for(int jj = 0; jj < 3; jj++){
					points_[ii][jj] = tmp_points[ii][jj];
				}
			}
			
    }

    long size() const {
      return points_.size();
    }
    
    template <class array_3d, class array_1d>
    void gather(const array_3d & grid, array_1d && subgrid) const {
			
			//DATAOPERATIONS STL TRANSFORM
			std::transform(points_.begin(), points_.end(), subgrid.begin(),
										 [& grid](auto point){
											 return grid[point[0]][point[1]][point[2]];
										 });
		}

		template <class array_4d>
    math::array<typename array_4d::element_type, 2> gather(const array_4d & grid) const {
			const int nst = std::get<3>(sizes(grid));
			math::array<typename array_4d::element_type, 2> subgrid({this->size(), nst});

			//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifdef HAVE_CUDA
			gpu::run(nst, size(),
							 [sgr = begin(subgrid), gr = begin(grid), pts = begin(points_)]
							 __device__ (auto ist, auto ipoint){
								 sgr[ipoint][ist] = gr[pts[ipoint][0]][pts[ipoint][1]][pts[ipoint][2]][ist];
							 });
#else
      for(int ipoint = 0; ipoint < size(); ipoint++){
				for(int ist = 0; ist < nst; ist++) subgrid[ipoint][ist] = grid[points_[ipoint][0]][points_[ipoint][1]][points_[ipoint][2]][ist];
      }
#endif
			return subgrid;
    }

    template <class array_2d, class array_4d>
    void scatter_add(const array_2d & subgrid, array_4d && grid) const{
			
			//DATAOPERATIONS LOOP + GPU::RUN 2D
#ifdef HAVE_CUDA
			gpu::run(std::get<1>(sizes(subgrid)), size(),
							 [sgr = begin(subgrid), gr = begin(grid), pts = begin(points_)]
							 __device__ (auto ist, auto ipoint){
								 gr[pts[ipoint][0]][pts[ipoint][1]][pts[ipoint][2]][ist] += sgr[ipoint][ist];
							 });
#else
      for(int ipoint = 0; ipoint < size(); ipoint++){
				for(int i1 = 0; i1 < std::get<1>(sizes(subgrid)); i1++){
					grid[points_[ipoint][0]][points_[ipoint][1]][points_[ipoint][2]][i1] += subgrid[ipoint][i1];
				}
      }
#endif
    }
    
    template <class array_1d, class array_3d>
    void scatter(const array_1d & subgrid, array_3d && grid) const{

			//DATAOPERATIONS LOOP 1D (random access output)
      for(int ipoint = 0; ipoint < size(); ipoint++){
				grid[points_[ipoint][0]][points_[ipoint][1]][points_[ipoint][2]] = subgrid[ipoint];
      }
    }

    auto points() const {
      return points_;
    }

		const double & volume_element() const {
			return volume_element_;
		}

		const auto & distance() const {
			return distance_;
		}

		const auto & point_pos() const {
			return relative_pos_;
		}		
		
		friend auto sizes(const spherical_grid & sphere){
			return std::array<long, dimension>{sphere.size()};
		}

  private:

		math::array<std::array<int, 3>, 1> points_;
		std::vector<float> distance_; //I don't think we need additional precision for this. XA
		std::vector<math::d3vector> relative_pos_;
		double volume_element_;
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>
#include <math/array.hpp>
#include <math/complex.hpp>

TEST_CASE("class basis::spherical_grid", "[basis::spherical_grid]") {
	
  using namespace Catch::literals;
  using math::d3vector;

  double ll = 10.0;
  
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  
  double ecut = 20.0;
  
  basis::real_space pw(cell, input::basis::cutoff_energy(ecut));
  
  SECTION("Point 0 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {0.0, 0.0, 0.0}, 2.0);
						       
    REQUIRE(sphere.size() == 257);

    math::array<complex, 3> grid(pw.rsize());
    std::vector<complex> subgrid(sphere.size());

    for(long ii = 0; ii < grid.num_elements(); ii++) grid.data()[ii] = 0.0;
    
    sphere.gather(grid, subgrid);

    for(unsigned ii = 0; ii < subgrid.size(); ii++) subgrid[ii] = 1.0; 
    
    sphere.scatter(subgrid, grid);

    double sum = 0.0;
    for(long ii = 0; ii < grid.num_elements(); ii++) sum += real(grid.data()[ii]);

    REQUIRE(sum == 257.0_a);
    
  }

  SECTION("Point -l/2 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {-ll/2.0, 0.0, 0.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);
    
    math::array<complex, 4> grid({pw.rsize()[0], pw.rsize()[1], pw.rsize()[2], 20}, 0.0);
    math::array<complex, 2> subgrid({sphere.size(), 20}, 0.0);

    for(long ii = 0; ii < grid.num_elements(); ii++) grid.data()[ii] = 1.0;
    
    sphere.gather(grid, subgrid);

    double sum = 0.0;
    for(long ii = 0; ii < subgrid.num_elements(); ii++) sum += real(subgrid.data()[ii]);

    REQUIRE(sum == Approx(20.0*257.0));
    
    for(long ii = 0; ii < subgrid.num_elements(); ii++) subgrid.data()[ii] = 0.0;
    
    sphere.scatter(subgrid, grid);

    sum = 0.0;
    for(long ii = 0; ii < grid.num_elements(); ii++) sum += real(grid.data()[ii]);

    REQUIRE(sum == Approx(20.0*(pw.size() - sphere.size())));
    
  }

  SECTION("Point l/2 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, 0.0, 0.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);

    math::array<complex, 6> grid({1, pw.rsize()[0], pw.rsize()[1], pw.rsize()[2], 2, 20}, 0.0);
    math::array<complex, 3> subgrid({sphere.size(), 2, 20}, 0.0);

    sphere.gather(grid[0], subgrid);

    sphere.scatter(subgrid, grid[0]);
    
  }

  SECTION("Point -l/2 -l/2 -l/2"){
    
    basis::spherical_grid sphere(pw, cell, {-ll/2.0, -ll/2.0, -ll/2.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);
    
  }

  SECTION("Point l/2 l/2 l/2"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, ll/2.0, ll/2.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);
    
  }

  SECTION("Point l/2 l/2 l/2"){
    
    basis::spherical_grid sphere(pw, cell, {ll/2.0, ll/2.0, ll/2.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);
    
  }

  SECTION("Point l/4 l/4 l/4"){
    
    basis::spherical_grid sphere(pw, cell, {ll/4.0, ll/4.0, ll/4.0}, 2.0);
    
    REQUIRE(sphere.size() == 257);
    
  }
  
}
#endif

#endif

