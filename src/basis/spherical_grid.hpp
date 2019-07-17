/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

namespace basis {

  class spherical_grid {

  public:

		const static int dimension = 1;
		
		template <class basis>
    spherical_grid(const basis & parent_grid, const ions::UnitCell & cell, const math::d3vector & center_point, const double radius){

      ions::periodic_replicas rep(cell, center_point, parent_grid.diagonal_length());
      
      for(int ix = 0; ix < parent_grid.rsize()[0]; ix++){
				for(int iy = 0; iy < parent_grid.rsize()[1]; iy++){
					for(int iz = 0; iz < parent_grid.rsize()[2]; iz++){
						auto rpoint = parent_grid.rvector(ix, iy, iz);
						
						for(int irep = 0; irep < rep.size(); irep++){
							if(norm(rpoint - rep[irep]) <= radius*radius) points_.push_back({ix, iy, iz});
						}
						
					}
				}
      }
			
    }

    long size() const {
      return points_.size();
    }
    
    template <class array_3d, class array_1d>
    void gather(const array_3d & grid, array_1d && subgrid) const {

			//DATAOPERATIONS
			
      for(int ipoint = 0; ipoint < size(); ipoint++){
				subgrid[ipoint] = grid[points_[ipoint][0]][points_[ipoint][1]][points_[ipoint][2]];
      }
    }

    template <class array_2d, class array_4d>
    void scatter_add(const array_2d & subgrid, array_4d && grid) const{

			//DATAOPERATIONS
			
      for(int ipoint = 0; ipoint < size(); ipoint++){
				for(int i1 = 0; i1 < std::get<1>(sizes(subgrid)); i1++){
					grid[points_[ipoint][0]][points_[ipoint][1]][points_[ipoint][2]][i1] += subgrid[ipoint][i1];
				}
      }
    }
    
    template <class array_1d, class array_3d>
    void scatter(const array_1d & subgrid, array_3d && grid) const{

			//DATAOPERATIONS
			
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

		friend auto sizes(const spherical_grid & sphere){
			return std::array<long, dimension>{sphere.size()};
		}
		
  private:
    
    std::vector<std::array<int, 3> > points_;
		double volume_element_;
		
  };
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <ions/unitcell.hpp>
#include <multi/array.hpp>
#include <math/complex.hpp>

TEST_CASE("class basis::spherical_grid", "[spherical_grid]") {
  
  using namespace Catch::literals;
  using math::d3vector;

  double ll = 10.0;
  
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  
  double ecut = 20.0;
  
  basis::real_space pw(cell, ecut);
  
  SECTION("Point 0 0 0"){
    
    basis::spherical_grid sphere(pw, cell, {0.0, 0.0, 0.0}, 2.0);
						       
    REQUIRE(sphere.size() == 257);

    boost::multi::array<complex, 3> grid(pw.rsize());
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
    
    boost::multi::array<complex, 4> grid({pw.rsize()[0], pw.rsize()[1], pw.rsize()[2], 20}, 0.0);
    boost::multi::array<complex, 2> subgrid({sphere.size(), 20}, 0.0);

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

    boost::multi::array<complex, 6> grid({1, pw.rsize()[0], pw.rsize()[1], pw.rsize()[2], 2, 20}, 0.0);
    boost::multi::array<complex, 3> subgrid({sphere.size(), 2, 20}, 0.0);

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

