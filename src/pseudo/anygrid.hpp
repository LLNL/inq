#ifndef PSEUDO_ANYGRID_HPP
#define PSEUDO_ANYGRID_HPP

/*
 Copyright (C) 2018 Xavier Andrade

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

#include <vector>
#include <string>
#include "base.hpp"
#include "spline.h"

namespace pseudopotential {

  class anygrid : public pseudopotential::base {

  public:

    anygrid(bool uniform_grid):
      uniform_grid_(uniform_grid){
    }
    
    double mesh_spacing() const {
      return 0.01;
    }
    
    int mesh_size() const {
      if(uniform_grid_) {
	return mesh_size_;
      } else {
	return grid_.size();
      }
    }
    
    virtual void grid(std::vector<double> & val) const {
      if(uniform_grid_){
	pseudopotential::base::grid(val);
      } else {
	val = grid_;
      }
    }

    virtual void grid_weights(std::vector<double> & val) const {
      if(uniform_grid_){
	pseudopotential::base::grid(val);
      } else {
	val = grid_weights_;
      }
    }

  protected:

    void interpolate(std::vector<double> & function) const {
      if(!uniform_grid_) return;

      std::vector<double> function_in_grid = function;
      
      assert(function.size() == grid_.size());
      
      Spline function_spline;
      function_spline.fit(grid_.data(), function_in_grid.data(), function_in_grid.size(), SPLINE_FLAT_BC, SPLINE_NATURAL_BC);
      
      function.clear();
      for(double rr = 0.0; rr <= grid_[grid_.size() - 1]; rr += mesh_spacing()){
	function.push_back(function_spline.value(rr));
      }
    }

    bool uniform_grid_;
    std::vector<double> grid_;
    std::vector<double> grid_weights_;
    int mesh_size_;
    
  };

}

#endif
