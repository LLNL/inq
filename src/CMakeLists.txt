include_directories(.)
include_directories(../external_libs)

add_executable(inq_unit_tests main/unit_tests.${SRC_EXT} main/unit_tests_main.${SRC_EXT})
