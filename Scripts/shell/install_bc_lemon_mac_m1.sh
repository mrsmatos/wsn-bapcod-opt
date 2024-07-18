cpu=`getconf _NPROCESSORS_ONLN`
cd Tools &&\
tar -xzf lemon-1.3.1.tar.gz &&\
cd lemon-1.3.1 &&\
ls &&\
mkdir -p build &&\
cd build &&\
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_OSX_ARCHITECTURES=x86_64 &&\
make -j$cpu VERBOSE=1 install &&\
cd ../.. 
