/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef OPERATIONS__PRECONDITION
#define OPERATIONS__PRECONDITION

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

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

#include <basis/coefficients_set.hpp>
#include <cstdlib>

namespace operations {

	class preconditioner {

	public:

		preconditioner(double arg_ecutpr):
			ecutpr_(arg_ecutpr){
		}
		
		template <class type>
		void operator()(basis::coefficients_set<basis::fourier_space, type> & phi) const {

			//DATAOPERATIONS
			for(int ix = 0; ix < phi.basis().gsize()[0]; ix++){
				for(int iy = 0; iy < phi.basis().gsize()[1]; iy++){
					for(int iz = 0; iz < phi.basis().gsize()[2]; iz++){
						
						if(phi.basis().g_is_zero(ix, iy, iz)){
							for(int ist = 0; ist < phi.set_size(); ist++) phi.cubic()[0][0][0][ist] = 0;
							continue;
						}
						double e = 0.5*phi.basis().g2(ix, iy, iz);
						double scal = ( e < ecutpr_ ) ? 0.5/ecutpr_ : 0.5/e;

						for(int ist = 0; ist < phi.set_size(); ist++) phi.cubic()[ix][iy][iz][ist] *= -scal;
					}
				}
			}

		}
		
		template <class type>
		void operator()(basis::coefficients_set<basis::real_space, type> & phi) const {
			
			auto fphi = operations::space::to_fourier(phi);
			operator()(fphi);
			phi = operations::space::to_real(fphi);

		}

	private:
			double ecutpr_;
	};
	
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("function operations::precondition", "[precondition]") {

	using namespace Catch::literals;

}


#endif

#endif
