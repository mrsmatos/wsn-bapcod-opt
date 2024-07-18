cpu=`getconf _NPROCESSORS_ONLN`
cpuX2=$(($cpu*2))
cd Tools && \
tar -xzf boost_1_76_0.tar.gz && \
cd boost_1_76_0 && \
unamestr=`uname` && \
if [ "$unamestr" == "Darwin" ]; then sed -i .old s/.*-fcoalesce-templates.*// tools/build/src/tools/darwin.jam; fi && \
./bootstrap.sh

if [ "$?" -ge 1 ];
then
    if `grep -q "error: implicit declaration of function" bootstrap.log`;
    then
		echo "-------------------------------------------------------------"
		echo "-------------------------------------------------------------"
		echo "  There seems to be a problem about -Wimplicit-function-declaration when compiling Boost."
		echo "  Trying to patch the source code."
		sed -i .old '/timestamp/a\ 
#include "../filesys.h"' tools/build/src/engine/modules/path.c
		echo "  Trying to run bootstrap.sh once again."
		./bootstrap.sh
	fi
fi 
if [ $? -ge 1 ]; 
then
   echo "definitive bootstrap.sh failure"
   exit 1
fi
./b2 cxxflags=" -fPIC " cflags=" -fPIC $CFLAGS " -j$cpuX2 --with-system --with-filesystem --with-program_options --with-thread --with-chrono --with-timer --with-regex install --prefix=./build/ && \
rm -rf libs boost doc tools 
