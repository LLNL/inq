AM_CXXFLAGS =                             \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/          \
	-DSHARE_DIR='"$(pkgdatadir)"'           \
	-DPREFIX_DIR='"$(prefix)"'              \
	@CFLAGS_CUDA@

LIBS = @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_CUDA@ @LIBS@

noinst_HEADERS =                             \
	config/path.hpp                            \
	gpu/run.hpp                                \
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
	math/vec3d.hpp                             \
	math/complex.hpp                           \
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
	operations/calculate_density.hpp           \
	operations/integral.hpp                    \
	operations/add.hpp                         \
	operations/subspace_diagonalization.hpp    \
	operations/sum.hpp                         \
	operations/matrix_operator.hpp             \
	hamiltonian/atomic_potential.hpp           \
	hamiltonian/energy.hpp                     \
	hamiltonian/exchange_operator.hpp          \
	hamiltonian/ks_hamiltonian.hpp             \
	hamiltonian/projector.hpp                  \
	hamiltonian/self_consistency.hpp           \
	hamiltonian/xc_functional.hpp              \
	solvers/poisson.hpp                        \
	solvers/least_squares.hpp                  \
	solvers/linear_mixer.hpp                   \
	solvers/linear.hpp                         \
	solvers/pulay_mixer.hpp                    \
	eigensolvers/conjugate_gradient.hpp        \
	eigensolvers/steepest_descent.hpp          \
	systems/ions.hpp                           \
	systems/electrons.hpp

all_srcs =                                   \
	ions/unitcell.cpp

bin_PROGRAMS  =                              \
	@PACKAGE_NAME@_unit_tests

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp main/unit_tests_main.cpp $(all_srcs)

#TEST PROGRAMS

@PACKAGE_NAME@_electron_gas_non_interacting_SOURCES = tests/electron_gas_non_interacting.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_electron_gas_lda_pz_SOURCES = tests/electron_gas_lda_pz.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_hydrogen_local_SOURCES = tests/hydrogen_local.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_hydrogen_non_local_SOURCES = tests/hydrogen_non_local.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_neon_SOURCES = tests/neon.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_kripton_SOURCES = tests/kripton.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_nitrogen_SOURCES = tests/nitrogen.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_silicon_SOURCES = tests/silicon.cpp main/unit_tests_main.cpp $(all_srcs)

test_programs =                                                 \
	@PACKAGE_NAME@_electron_gas_non_interacting	                  \
	@PACKAGE_NAME@_electron_gas_lda_pz                            \
	@PACKAGE_NAME@_hydrogen_local                                 \
	@PACKAGE_NAME@_hydrogen_non_local                             \
	@PACKAGE_NAME@_neon                                           \
	@PACKAGE_NAME@_kripton                                        \
	@PACKAGE_NAME@_nitrogen                                       \
	@PACKAGE_NAME@_silicon

bin_PROGRAMS += $(test_programs) 

check: @PACKAGE_NAME@_unit_tests $(test_programs)
	./@PACKAGE_NAME@_unit_tests
	./@PACKAGE_NAME@_hydrogen_local
	./@PACKAGE_NAME@_neon
