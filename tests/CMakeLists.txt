# Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include_directories(.)
include_directories(../external_libs)
include_directories(../src)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(TEST_FILE ${TEST_SRCS})
  get_filename_component(TEST_EXE ${TEST_FILE} NAME_WE)
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_include_directories(${TEST_EXE} PRIVATE ${catch2_SOURCE_DIR}/src)
  target_link_libraries(${TEST_EXE} PRIVATE Catch2::Catch2 inq)
  add_test(NAME inq::tests::${TEST_EXE} COMMAND ${PROJECT_BINARY_DIR}/scripts/run_test.sh ./${TEST_EXE})
  if(ENABLE_CUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
  endif()
endforeach()

set_tests_properties(inq::tests::h2o_real_time PROPERTIES DEPENDS inq::tests::h2o_ground_state)
set_tests_properties(inq::tests::al4h1_td PROPERTIES DEPENDS inq::tests::al4h1)

#PYTHON TESTS

file(GLOB TEST_SCRIPTS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

foreach(TEST_FILE ${TEST_SCRIPTS})
  add_test(NAME inq::tests::${TEST_FILE} COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_FILE})
  set_property(TEST inq::tests::${TEST_FILE} PROPERTY ENVIRONMENT_MODIFICATION "PYTHONPATH=path_list_append:${PROJECT_BINARY_DIR}")
endforeach()
