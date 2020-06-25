include_directories(.)
include_directories(../external_libs)

file(GLOB INQ_SOURCES "*/*.${SRC_EXT}")
add_library(inq STATIC ${INQ_SOURCES})

add_executable(inq_unit_tests main/unit_tests.${SRC_EXT} main/unit_tests_main.${SRC_EXT})
target_link_libraries(inq_unit_tests PUBLIC inq)
