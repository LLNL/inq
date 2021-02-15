include_directories(.)
include_directories(../external_libs)
include_directories(../src)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(TEST_FILE ${TEST_SRCS})
  get_filename_component(TEST_EXE ${TEST_FILE} NAME_WE)
  add_executable(${TEST_EXE} ${TEST_FILE})
  add_test(NAME ${TEST_EXE} COMMAND ${PROJECT_SOURCE_DIR}/scripts/run_test.sh ./${TEST_EXE})
  if(ENABLE_CUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
  endif()
endforeach()

set_tests_properties(h2o_real_time PROPERTIES DEPENDS h2o_ground_state)
