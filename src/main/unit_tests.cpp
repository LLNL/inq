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

#define UNIT_TEST

#include <catch2/catch.hpp>

#include <config/path.hpp>
#include <math/spline.hpp>
#include <pseudo/element.hpp>
#include <pseudo/pseudopotential.hpp>
#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <ions/interaction.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/plane_wave.hpp>
#include <hamiltonian/atomic_potential.hpp>

