#Override the defaults to use inc++
CXX=$(top_builddir)/scripts/inc++
CXXLD=$(CXX)
CXXFLAGS=""
LDFLAGS=""
LIBS=""

electron_gas_non_interacting_SOURCES = electron_gas_non_interacting.cpp 
electron_gas_lda_pz_SOURCES = electron_gas_lda_pz.cpp
hydrogen_local_SOURCES = hydrogen_local.cpp
hydrogen_non_local_SOURCES = hydrogen_non_local.cpp
neon_SOURCES = neon.cpp
kripton_SOURCES = kripton.cpp
nitrogen_SOURCES = nitrogen.cpp
silicon_SOURCES = silicon.cpp

noinst_PROGRAMS =                                \
	electron_gas_non_interacting                   \
	hydrogen_local                                 \
	hydrogen_non_local                             \
	neon                                           \
	kripton                                        \
	nitrogen                                       \
	silicon
