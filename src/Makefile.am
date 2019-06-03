AM_CXXFLAGS =                                   \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/catch2    \
	-I$(top_srcdir)/external_libs/rapidxml  \
	-DSHARE_DIR='"$(pkgdatadir)"'

noinst_HEADERS =                                   \
	math/d3vector.hpp                          \
	math/spline.hpp                            \
	pseudo/anygrid.hpp                         \
	pseudo/base.hpp                            \
	pseudo/element.hpp                         \
	pseudo/psml.hpp                            \
	pseudo/qso.hpp                             \
	pseudo/set.hpp                             \
	pseudo/share_directory.hpp                 \
	pseudo/upf.hpp                             \
	pseudo/upf1.hpp                            \
	pseudo/upf2.hpp                            \
	ions/geometry.hpp                          \
	ions/periodic_replicas.hpp                 \
	ions/interaction.hpp

all_srcs =                                         \
	math/spline.cpp                            \
	ions/unitcell.cpp

bin_PROGRAMS  =                                    \
	inq                                        \
	@PACKAGE_NAME@_unit_tests

inq_SOURCES = main/main.cpp $(all_srcs)

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp main/unit_tests_main.cpp $(all_srcs)

check: $(@PACKAGE_NAME@_unit_tests) install
	@PACKAGE_NAME@_unit_tests -s
