#source /opt/intel/system_studio_2020/bin/compilervars.sh intel64
rm -rf ./intel.build ./intel.install 
mkdir ./intel.build
mkdir ./intel.install
cd intel.build
rm -f CMakeCache.txt
CC=/opt/intel/system_studio_2020/bin/icc CXX=/opt/intel/system_studio_2020/bin/icpc \
  cmake ../.. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_COMPILER=/opt/intel/system_studio_2020/bin/icc \
    -DCMAKE_CXX_COMPILER=/opt/intel/system_studio_2020/bin/icpc \
    -DCMAKE_LINKER=/opt/intel/system_studio_2020/bin/icpc \
    -DCMAKE_INSTALL_PREFIX=`pwd`/../intel.install \
    -DCMAKE_C_COMPILE_FLAGS="$(mpicxx -showme:compile || mpicxx -cxx= -compile_info)" \
    -DCMAKE_CXX_COMPILE_FLAGS="$(mpicxx -showme:compile || mpicxx -cxx= -compile_info)" \
    -DCMAKE_CXX_FLAGS="-O3 -Wall -Wextra -Wl,-rpath,/opt/intel/system_studio_2020/lib/intel64" \
    -DCMAKE_EXE_LINKER_FLAGS="$(mpicxx -showme:link || mpicxx -cxx= -link_info)" \
    $* && \
    make -j $(($(nproc)/2 + 1)) && \
    make install || exit 1

