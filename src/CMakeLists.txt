include_directories(.)
include_directories(../external_libs)

add_executable(inq_unit_tests main/unit_tests.${SRC_EXT} main/unit_tests_main.${SRC_EXT})

file(GLOB INCLUDE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} density/*.hpp)

foreach(TEST_FILE ${INCLUDE_FILES})
  get_filename_component(dir ${TEST_FILE} DIRECTORY)
  get_filename_component(name ${TEST_FILE} NAME_WE)
  set(TEST_SRC ${dir}/${name}.${SRC_EXT})
  set(TEST_EXE "${dir}_${name}.utx")
  add_custom_command(OUTPUT ${TEST_SRC}
    COMMAND echo ARGS \\\#define INQ_${dir}_${name}_UNIT_TEST > ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#include \\\"${CMAKE_CURRENT_SOURCE_DIR}/${TEST_FILE}\\\" >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#ifndef INQ_TESTED >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#error NO TEST INCLUDED >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#endif >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}    
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_FILE})
  add_executable(${TEST_EXE} ${TEST_SRC} main/unit_tests_main.${SRC_EXT})
endforeach()

