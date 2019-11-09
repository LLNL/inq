AM_CXXFLAGS =                             \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/          \
	-DSHARE_DIR='"$(pkgdatadir)"'           \
	-DPREFIX_DIR='"$(prefix)"'              \
	$(GSL_CFLAGS)                           \
	@CFLAGS_CUDA@

LIBS += @GSL_LIBS@ @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_CUDA@

noinst_HEADERS =                             \
	gpu/run.hpp                                \
	utils/merge_optional.hpp                   \
	parser/gsl_userdef.h                       \
	parser/input_file.hpp                      \
	parser/liboct_parser.h                     \
	parser/string_f.h                          \
	parser/symbols.h                           \
	input/interaction.hpp                      \
	input/basis.hpp                            \
	input/cell.hpp                             \
	input/config.hpp                           \
	input/species.hpp                          \
	input/atom.hpp                             \
	math/array.hpp                             \
	math/d3vector.hpp                          \
	ions/geometry.hpp                          \
	ions/periodic_replicas.hpp                 \
	ions/interaction.hpp                       \
	basis/grid.hpp                             \
	basis/real_space.hpp                       \
	basis/field.hpp                            \
	basis/field_set.hpp                        \
	basis/spherical_grid.hpp                   \
	basis/trivial.hpp                          \
	states/ks_states.hpp                       \
	operations/overlap.hpp                     \
	operations/randomize.hpp                   \
	operations/scal.hpp                        \
	operations/space.hpp                       \
	operations/orthogonalization.hpp           \
	operations/shift.hpp                       \
	operations/preconditioner.hpp              \
	operations/diagonalize.hpp                 \
	operations/calculate_density.hpp           \
	operations/integral.hpp                    \
	operations/add.hpp                         \
	operations/subspace_diagonalization.hpp    \
	operations/sum.hpp                         \
	hamiltonian/atomic_potential.hpp           \
	hamiltonian/ks_hamiltonian.hpp             \
	hamiltonian/self_consistency.hpp           \
	hamiltonian/energy.hpp                     \
	hamiltonian/xc_functional.hpp              \
	solvers/poisson.hpp                        \
	solvers/steepest_descent.hpp               \
	solvers/linear_mixer.hpp                   \
	solvers/linear.hpp                         \
	solvers/pulay_mixer.hpp                    \
	systems/ions.hpp                           \
	systems/electrons.hpp

all_srcs =                                   \
	parser/symbols.c                           \
	parser/parse.c                             \
	parser/parse_exp.c                         \
	parser/gsl_userdef.c                       \
	ions/unitcell.cpp

bin_PROGRAMS  =                              \
	inq                                        \
	@PACKAGE_NAME@_unit_tests

inq_SOURCES = main/main.cpp $(all_srcs)

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp main/unit_tests_main.cpp $(all_srcs)

#TEST PROGRAMS

@PACKAGE_NAME@_electron_gas_non_interacting_SOURCES = tests/electron_gas_non_interacting.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_electron_gas_lda_pz_SOURCES = tests/electron_gas_lda_pz.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_hydrogen_local_SOURCES = tests/hydrogen_local.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_hydrogen_non_local_SOURCES = tests/hydrogen_non_local.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_neon_SOURCES = tests/neon.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_kripton_SOURCES = tests/kripton.cpp main/unit_tests_main.cpp $(all_srcs)
@PACKAGE_NAME@_nitrogen_SOURCES = tests/nitrogen.cpp main/unit_tests_main.cpp $(all_srcs)

test_programs =                                                 \
	@PACKAGE_NAME@_electron_gas_non_interacting	                  \
	@PACKAGE_NAME@_electron_gas_lda_pz                            \
	@PACKAGE_NAME@_hydrogen_local                                 \
	@PACKAGE_NAME@_hydrogen_non_local                             \
	@PACKAGE_NAME@_neon                                           \
	@PACKAGE_NAME@_kripton                                        \
	@PACKAGE_NAME@_nitrogen

bin_PROGRAMS += $(test_programs) 

check: @PACKAGE_NAME@_unit_tests $(test_programs)
	for test_prog in @PACKAGE_NAME@_unit_tests $(test_programs); do $(top_builddir)/src/$$test_prog; done
