include_directories(.)
include_directories(../external_libs)
include_directories(../src)

link_libraries(inq)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPREFIX_DIR='\"${CMAKE_INSTALL_PREFIX}\"'" )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DSHARE_DIR='\"${CMAKE_INSTALL_PREFIX}/share/inq\"'" )

add_executable(electron_gas_non_interacting electron_gas_non_interacting.cpp) 
add_executable(electron_gas_lda_pz electron_gas_lda_pz.cpp)
add_executable(hydrogen_local hydrogen_local.cpp)
add_executable(hydrogen_non_local hydrogen_non_local.cpp)
add_executable(neon neon.cpp)
add_executable(kripton kripton.cpp)
add_executable(nitrogen nitrogen.cpp)
add_executable(silicon silicon.cpp)
add_executable(h2o_ground_state h2o_ground_state.cpp)
add_executable(h2o_real_time h2o_real_time.cpp)
