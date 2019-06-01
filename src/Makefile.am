AM_CXXFLAGS =                                   \
	-I$(top_srcdir)/src/                    \
	-I$(top_srcdir)/external_libs/catch2    \
	-I$(top_srcdir)/external_libs/rapidxml  \
	-DSHARE_DIR='"$(pkgdatadir)"'

noinst_HEADERS =                                   \
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
	ions/position.hpp                          \
	ions/geometry.hpp

bin_PROGRAMS  =                                    \
	inq                                        \
	@PACKAGE_NAME@_unit_tests

inq_SOURCES = main/main.cpp

@PACKAGE_NAME@_unit_tests_SOURCES = main/unit_tests.cpp
