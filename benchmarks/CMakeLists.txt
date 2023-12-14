include_directories(.)
include_directories(../external_libs)
include_directories(../src)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(TEST_FILE ${TEST_SRCS})
  get_filename_component(TEST_EXE ${TEST_FILE} NAME_WE)
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_link_libraries(${TEST_EXE} PRIVATE Catch2::Catch2 inq)
  target_include_directories(${TEST_EXE} PRIVATE ${catch2_SOURCE_DIR}/src)
  if(ENABLE_GPU)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE ${GPU_LANGUAGE})
  endif()
endforeach()

