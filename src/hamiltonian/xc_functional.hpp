/* -*- indent-tabs-mode: t -*- */

#ifndef HAMILTONIAN__XC_FUNCTIONAL
#define HAMILTONIAN__XC_FUNCTIONAL

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

#include <xc.h>

namespace hamiltonian {
	namespace xc_functional {
		
		template <class density_type, class exc_type, class vxc_type>
		void unpolarized(long size, density_type const & density, exc_type & exc, vxc_type & vxc){
			
			xc_func_type func;
			int vmajor, vminor, vmicro, func_id = XC_LDA_X;
			
			xc_version(&vmajor, &vminor, &vmicro);
			printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);
			
			if(xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0){
				fprintf(stderr, "Functional '%d' not found\n", func_id);
				exit(1);
			}
			
			switch(func.info->family)
				{
				case XC_FAMILY_LDA:
					xc_lda_exc_vxc(&func, size, density.data(), exc.data(), vxc.data());
					break;
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					break;
				}
			
			xc_func_end(&func);

		}
	}
}

#endif
