add_subdirectory(pseudopod)

message("\nCONFIGURING LIBXC\n")
set(DISABLE_FXC ON)
set(BUILD_TESTING OFF)
add_subdirectory(libxc)
message("\nDONE CONFIGURING LIBXC\n")

message("\nCONFIGURING CALIPER\n")
string(REPLACE -Werror -Wno-error CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(BUILD_SHARED_LIBS ON)
add_subdirectory(caliper)
message("\nDONE CONFIGURING CALIPER\n")
