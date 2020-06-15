/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__UTILS__FINITE_DIFFERENCE
#define INQ__UTILS__FINITE_DIFFERENCE

/*
 Copyright (C) 2020 Xavier Andrade

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

namespace inq {
namespace utils {

    template <class FunctionType, class VecType>
    auto finite_difference_gradient1D(FunctionType const & function, VecType const & point){
	double const delta = 0.1;
	VecType gradient = {0.0, 0.0, 0.0};
	for(int idir = 0; idir < 3; idir++){
	    for(auto isign : {-1, 1}){
		auto point_plus_delta = point;
		point_plus_delta[idir] += isign*delta;
		    gradient[idir] += isign*function(point_plus_delta)/(2.0*delta);
	    }
	}
	return gradient;
}

}
}
#endif

