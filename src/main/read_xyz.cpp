#define SHARE_DIR "/home/xavier/files/codes/inq/inq/share/"

#include "ions/geometry.hpp"

#include <fstream>
#include <string>
#include <cmath>


int main(const int argc, const char ** argv){

  std::string xyz_file_name = "file.xyz";

  std::cout << "Opening file " << xyz_file_name << std::endl;

   ions::geometry geo(xyz_file_name);

}
