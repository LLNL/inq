add_subdirectory(pseudopod)

message("\nCONFIGURING LIBXC\n")
set(DISABLE_FXC ON)
add_subdirectory(libxc)
message("\nDONE CONFIGURING LIBXC\n")

message("\nCONFIGURING CALIPER\n")
set(BUILD_SHARED_LIBS ON)
add_subdirectory(caliper)
message("\nDONE CONFIGURING CALIPER\n")
