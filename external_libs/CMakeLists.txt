add_subdirectory(pseudopod)

message("\nCONFIGURING LIBXC\n")
add_subdirectory(libxc)
message("\nDONE CONFIGURING LIBXC\n")

message("\nCONFIGURING CALIPER\n")
add_subdirectory(caliper)
message("\nDONE CONFIGURING CALIPER\n")
