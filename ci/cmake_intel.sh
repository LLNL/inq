. /opt/intel/system_studio_2020/bin/compilervars.sh intel64
#export OMPI_CC=icc
#export OMPI_CXX=icpc
#export CXXFLAGS=" "
#export LIBS=$(for x in `mpic++ --showme:libs`; do echo -n -l$x" " ; done)
#export LDFLAGS=$(for x in `mpic++ --showme:libs`; do echo -n -l$x" " ; done)
#export LIBS=`mpicxx.mp
rm -rf ./intel.dir.build ./intel.dir.install 
mkdir ./intel.dir.build
mkdir ./intel.dir.install
cd intel.dir.build
rm -f CMakeCache.txt
cmake ../.. \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_C_COMPILER=/opt/intel/system_studio_2020/bin/icc -DCMAKE_CXX_COMPILER=/opt/intel/system_studio_2020/bin/icpc -DCMAKE_LINKER=/opt/intel/system_studio_2020/bin/icpc \
	-DCMAKE_INSTALL_PREFIX=`pwd`/../intel.dir.install \
	-DCMAKE_C_COMPILE_FLAGS="$(mpicxx -showme:compile || mpicxx -cxx= -compile_info)" \
	-DCMAKE_CXX_COMPILE_FLAGS="$(mpicxx -showme:compile || mpicxx -cxx= -compile_info)" \
	-DCMAKE_CXX_FLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread -O3 -Wall -Wextra -Warray-bounds -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Wno-unused-parameter -Werror -Wl,-rpath,/opt/intel/system_studio_2020/lib/intel64" \
	-DCMAKE_EXE_LINKER_FLAGS="$(mpicxx -showme:link || mpicxx -cxx= -link_info)" \
	$*
make -j $(($(nproc)/2 + 1))  && make install && ctest -j $(($(nproc)/2 + 1)) --output-on-failure

cd src; (INQ_EXEC_ENV="mpirun --oversubscribe -np 2" ctest --output-on-failure --timeout 270 -j8 || exit) ; cd ..
cd src; (INQ_EXEC_ENV="mpirun --oversubscribe -np 3" ctest --output-on-failure --timeout 270 -j5 || exit) ; cd ..
cd src; (INQ_EXEC_ENV="mpirun --oversubscribe -np 4" ctest --output-on-failure --timeout 270 -j4 || exit) ; cd ..

#	-DCMAKE_CXX_COMPILE_FLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread" \
#	-DCMAKE_COMPILE_FLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread" \
#	-DCMAKE_C_FLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread" \

