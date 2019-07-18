AM_CXXFLAGS =                                   \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/          \
	-DSHARE_DIR='"$(pkgdatadir)"'           \
	$(GSL_CFLAGS)

LIBS += @GSL_LIBS@ @LIBS_LAPACK@ @LIBS_BLAS@

noinst_HEADERS =                                   \
	parser/gsl_userdef.h                       \
	parser/input_file.hpp                       \
	parser/liboct_parser.h                     \
	parser/string_f.h                          \
	parser/symbols.h                           \
	math/d3vector.hpp                          \
	math/spline.hpp                            \
	math/spherical_harmonic.hpp                \
	pseudo/anygrid.hpp                         \
	pseudo/base.hpp                            \
	pseudo/element.hpp                         \
	pseudo/psml.hpp                            \
	pseudo/qso.hpp                             \
	pseudo/set.hpp                             \
	pseudo/upf.hpp                             \
	pseudo/upf1.hpp                            \
	pseudo/upf2.hpp                            \
	ions/geometry.hpp                          \
	ions/periodic_replicas.hpp                 \
	ions/interaction.hpp                       \
	basis/grid.hpp                             \
	basis/real_space.hpp                       \
	basis/coefficients.hpp                     \
	basis/coefficients_set.hpp                 \
	basis/spherical_grid.hpp                   \
	states/ks_states.hpp                       \
	functionals/LDAFunctional.hpp              \
	functionals/PBESolFunctional.hpp           \
	functionals/BLYPFunctional.hpp             \
	functionals/PBERevFunctional.hpp           \
	functionals/XCFunctional.hpp               \
	functionals/PBEFunctional.hpp              \
	operations/overlap.hpp                     \
	operations/randomize.hpp                   \
	operations/scal.hpp                        \
	operations/space.hpp                       \
	operations/orthogonalization.hpp           \
	operations/shift.hpp                       \
	operations/diagonalize.hpp                 \
	hamiltonian/atomic_potential.hpp           \
	hamiltonian/ks_hamiltonian.hpp             \
	solvers/poisson.hpp                        \
	solvers/steepest_descent.hpp               \
	main/electrons.hpp

all_srcs =                                         \
	parser/symbols.c                           \
	parser/parse.c                             \
	parser/parse_exp.c                         \
	parser/gsl_userdef.c                       \
	math/spline.cpp                            \
	ions/unitcell.cpp                          \
	functionals/PBEFunctional.cpp              \
	functionals/PBERevFunctional.cpp           \
	functionals/BLYPFunctional.cpp             \
	functionals/PBESolFunctional.cpp

bin_PROGRAMS  =                                    \
	inq                                        \
	@PACKAGE_NAME@_unit_tests

inq_SOURCES = main/main.cpp $(all_srcs)

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp main/unit_tests_main.cpp $(all_srcs)

check: $(@PACKAGE_NAME@_unit_tests) install
	@PACKAGE_NAME@_unit_tests -s
