include_directories(.)
include_directories(../external_libs)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPREFIX_DIR='\"${CMAKE_INSTALL_PREFIX}\"'" )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DSHARE_DIR='\"${CMAKE_INSTALL_PREFIX}/share/inq\"'" )


file(GLOB INQ_SOURCES "*/*.cpp")
add_library(inq STATIC ${INQ_SOURCES})

add_executable(inq_unit_tests main/unit_tests.cpp main/unit_tests_main.cpp)


