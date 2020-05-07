#    
#    Copyright (C) 2020 Xavier Andrade <xavier@llnl.gov>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation, either
#    version 3 of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public
#    License along with this program. If not, see
#    <https://www.gnu.org/licenses/>.

#Override the defaults to use inc++
CXX=$(top_builddir)/scripts/inc++
CXXLD=$(CXX)
CXXFLAGS=
LDFLAGS=
LIBS=

electron_gas_non_interacting_SOURCES = electron_gas_non_interacting.cpp 
electron_gas_lda_pz_SOURCES = electron_gas_lda_pz.cpp
hydrogen_local_SOURCES = hydrogen_local.cpp
hydrogen_non_local_SOURCES = hydrogen_non_local.cpp
neon_SOURCES = neon.cpp
kripton_SOURCES = kripton.cpp
nitrogen_SOURCES = nitrogen.cpp
silicon_SOURCES = silicon.cpp
h2o_real_time_SOURCES = h2o_real_time.cpp

noinst_PROGRAMS =                                \
	electron_gas_non_interacting                   \
	electron_gas_lda_pz                            \
	hydrogen_local                                 \
	hydrogen_non_local                             \
	neon                                           \
	kripton                                        \
	nitrogen                                       \
	silicon                                        \
	h2o_real_time

