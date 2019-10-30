/* -*- indent-tabs-mode: t -*- */

#ifndef MATH__ARRAY
#define MATH__ARRAY

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

#include <multi/array.hpp>

#ifdef HAVE_CUDA
#include <multi/memory/adaptors/cuda/allocator.hpp>
#endif

namespace math {

  template <class type, size_t dim,
						class allocator = std::allocator<type>
						>
  using array = boost::multi::array<type, dim, allocator>;
	//	using array = boost::multi::array<type, dim, boost::multi::memory::cuda::managed::allocator<type>>;
 
}

#endif
