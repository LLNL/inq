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

AM_CXXFLAGS =                             \
	-DNO_CATCH_MAIN                         \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/          \
	-DSHARE_DIR='"$(pkgdatadir)"'           \
	-DPREFIX_DIR='"$(prefix)"'              \
	@CFLAGS_CUDA@

LIBS = @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_CUDA@ @LIBS@

noinst_HEADERS =                             \
	config/path.hpp                            \
	gpu/run.hpp                                \
	utils/match.hpp                            \
	utils/merge_optional.hpp                   \
	utils/partition.hpp                        \
	utils/skeleton_wrapper.hpp                 \
	input/interaction.hpp                      \
	input/atom.hpp                             \
	input/basis.hpp                            \
	input/cell.hpp                             \
	input/config.hpp                           \
	input/scf.hpp                              \
	input/species.hpp                          \
	math/array.hpp                             \
	math/complex.hpp                           \
	math/vec3d.hpp                             \
	math/vector3.hpp                           \
	ions/geometry.hpp                          \
	ions/periodic_replicas.hpp                 \
	ions/interaction.hpp                       \
	ions/unitcell.hpp                          \
	basis/base.hpp                             \
	basis/grid.hpp                             \
	basis/real_space.hpp                       \
	basis/field.hpp                            \
	basis/field_set.hpp                        \
	basis/fourier_space.hpp                    \
	basis/spherical_grid.hpp                   \
	basis/trivial.hpp                          \
	states/ks_states.hpp                       \
	operations/overlap.hpp                     \
	operations/randomize.hpp                   \
	operations/scal.hpp                        \
	operations/space.hpp                       \
	operations/orthogonalize.hpp               \
	operations/shift.hpp                       \
	operations/preconditioner.hpp              \
	operations/diagonalize.hpp                 \
	operations/integral.hpp                    \
	operations/add.hpp                         \
	operations/subspace_diagonalization.hpp    \
	operations/sum.hpp                         \
	operations/transfer.hpp                    \
	operations/matrix_operator.hpp             \
	operations/laplacian.hpp                   \
	density/calculate.hpp                      \
	density/normalize.hpp                      \
	hamiltonian/atomic_potential.hpp           \
	hamiltonian/energy.hpp                     \
	hamiltonian/exchange_operator.hpp          \
	hamiltonian/ks_hamiltonian.hpp             \
	hamiltonian/projector.hpp                  \
	hamiltonian/self_consistency.hpp           \
	hamiltonian/xc_functional.hpp              \
	solvers/poisson.hpp                        \
	solvers/least_squares.hpp                  \
	solvers/linear.hpp                         \
	mixers/linear.hpp                          \
	mixers/pulay.hpp                           \
	eigensolvers/conjugate_gradient.hpp        \
	eigensolvers/steepest_descent.hpp          \
	systems/ions.hpp                           \
	systems/electrons.hpp                      \
	ground_state/calculate.hpp

all_srcs =                                   \
	ions/unitcell.cpp

lib_LIBRARIES = libinq.a

libinq_a_SOURCES = $(all_srcs)

bin_PROGRAMS  =                              \
	@PACKAGE_NAME@_unit_tests

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp main/unit_tests_main.cpp $(all_srcs)

check: @PACKAGE_NAME@_unit_tests
	./@PACKAGE_NAME@_unit_tests
