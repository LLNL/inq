# -*- indent-tabs-mode: t -*-
include_directories(.)
include_directories(../external_libs)


add_library(unit_tests_main OBJECT main/unit_tests_main.${SRC_EXT})

file(GLOB INCLUDE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */*.hpp)


foreach(TEST_FILE ${INCLUDE_FILES})
  get_filename_component(dir ${TEST_FILE} DIRECTORY)
  get_filename_component(name ${TEST_FILE} NAME_WE)
  set(TEST_SRC ${dir}/${name}.${SRC_EXT})
  set(TEST_O "${dir}_${name}.o")
  set(TEST_EXE "${dir}_${name}.utx")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  string(TOUPPER INQ_${dir}_${name}_UNIT_TEST test_macro)
  add_custom_command(OUTPUT ${TEST_SRC}
    COMMAND echo ARGS \\\#define ${test_macro} > ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#include \\\"${CMAKE_CURRENT_SOURCE_DIR}/${TEST_FILE}\\\" >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#ifdef ${test_macro} >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#error NO TEST INCLUDED, MACRO ${test_macro} STILL DEFINED >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}
    COMMAND echo ARGS \\\#endif >> ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SRC}    
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_FILE})
  add_library(${TEST_O} OBJECT ${TEST_SRC})
  add_executable(${TEST_EXE})
  target_link_libraries(${TEST_EXE} ${TEST_O})
  target_link_libraries(${TEST_EXE} unit_tests_main)

  add_test(NAME unit_test_${dir}_${name} COMMAND ${PROJECT_SOURCE_DIR}/scripts/run_test.sh ./${TEST_EXE})
endforeach()

