add_subdirectory(pseudopod)

message("\nCONFIGURING SPGLIB\n")
add_subdirectory(spglib)

message("\nCONFIGURING LIBXC\n")
set(DISABLE_FXC ON)
set(BUILD_TESTING OFF)
add_subdirectory(libxc)
message("\nDONE CONFIGURING LIBXC\n")

