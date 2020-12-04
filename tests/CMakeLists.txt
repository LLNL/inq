include_directories(.)
include_directories(../external_libs)
include_directories(../src)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.${SRC_EXT})

foreach(TEST_FILE ${TEST_SRCS})
  get_filename_component(TEST_EXE ${TEST_FILE} NAME_WE)
  add_executable(${TEST_EXE} ${TEST_FILE})
  add_test(NAME ${TEST_EXE} COMMAND ${PROJECT_SOURCE_DIR}/scripts/run_test.sh ./${TEST_EXE})
endforeach()
