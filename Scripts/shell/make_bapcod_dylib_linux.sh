#!/bin/sh

echo "Compile Bapcod Julia library \n"


BAPCOD_ROOT=/home/eduardo/bapcodframework
BUILD_ROOT="$BAPCOD_ROOT/build"
ARCH=Linux-x64

BAPCOD_INCLUDE_DIR="$BAPCOD_ROOT/Bapcod/include_dev"
LIBBAPCOD="$BUILD_ROOT/Bapcod/lib/Linux-x64/libbapcod_release.a"
OBJECTS_PATH_LOCATION="$BUILD_ROOT/Bapcod/CMakeFiles/bapcod.dir/src"

CPLEXSO_PATH_LOCATION="/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/"
BOOSTSO_PATH_LOCATION="$BAPCOD_ROOT/Tools/boost_1_76_0/build/lib"
#COINORSO_PATH_LOCATION="/Users/ruslansadykov/Programming/Cgl-0.60.3/build/lib"
BCPRCSP_LIB="$BUILD_ROOT/Tools/rcsp/lib/librcsp.a"

if ! [ -e $LIBBAPCOD ]
then
    echo "\033[31m $LIBBAPCOD not found \033[00m"
    echo "Check LIBBAPCOD variable content in run.sh \n"
    exit
fi


if ! [ -e $CPLEXSO_PATH_LOCATION ]
then
    echo "\033[31m $CPLEXSO_PATH_LOCATION not found \033[00m"
    echo "Check CPLEXSO_PATH_LOCATION variable content in run.sh \n"
    exit
fi


if ! [ -e $BOOSTSO_PATH_LOCATION ]
then
    echo "\033[31m $BOOSTSO_PATH_LOCATION not found \033[00m"
    echo "Check BOOSTSO_PATH_LOCATION variable content in run.sh \n"
    exit
fi


if ! [ -e $OBJECTS_PATH_LOCATION ]
then
    echo "\033[31m $OBJECTS_PATH_LOCATION not found \033[00m"
    echo "Check OBJECTS_PATH_LOCATION variable content in run.sh \n"
    exit
fi

if ! [ -e $LIBBAPCOD ]
then
   echo "$LIBBAPCOD does not exist!!"
fi   

export LD_LIBRARY_PATH=$BOOSTSO_PATH_LOCATION:$CPLEXSO_PATH_LOCATION:$COINORSO_PATH_LOCATION:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH

echo "ld paths:"
echo $LD_LIBRARY_PATH

echo "dyld paths:"
echo $DYLD_LIBRARY_PATH

gcc -I$BAPCOD_INCLUDE_DIR -shared -fPIC -o libdevjulia.so -Wl,--whole-archive $LIBBAPCOD $BCPRCSP_LIB $BOOSTSO_PATH_LOCATION/libboost_program_options.a $BOOSTSO_PATH_LOCATION/libboost_system.a $BOOSTSO_PATH_LOCATION/libboost_filesystem.a $BOOSTSO_PATH_LOCATION/libboost_regex.a $BOOSTSO_PATH_LOCATION/libboost_thread.a $BOOSTSO_PATH_LOCATION/libboost_timer.a $BOOSTSO_PATH_LOCATION/libboost_chrono.a -Wl,--no-whole-archive -L$CPLEXSO_PATH_LOCATION -lcplex12100 -lstdc++


#$LIBBAPCOD_DEBUG $CLQSEP_LIB $CVRPSEP_LIB $RHECCSEP_LIB $BOOSTSO_PATH_LOCATION/libboost_program_options.a  $BOOSTSO_PATH_LOCATION/libboost_system.a $BOOSTSO_PATH_LOCATION/libboost_filesystem.a $BOOSTSO_PATH_LOCATION/libboost_regex.a $BOOSTSO_PATH_LOCATION/libboost_thread.a  $BOOSTSO_PATH_LOCATION/libboost_timer.a $BOOSTSO_PATH_LOCATION/libboost_chrono.a -Wl,--no-whole-archive -L$CPLEXSO_PATH_LOCATION -lcplex12100

#gcc -I$BAPCOD_INCLUDE_DIR -shared -fPIC -o libdevjulia.so -Wl,--whole-archive $LIBBAPCOD $BCPRCSP_LIB $BOOSTSO_PATH_LOCATION/libboost_program_options.a $BOOSTSO_PATH_LOCATION/libboost_system.a $BOOSTSO_PATH_LOCATION/libboost_filesystem.a $BOOSTSO_PATH_LOCATION/libboost_regex.a $BOOSTSO_PATH_LOCATION/libboost_thread.a $BOOSTSO_PATH_LOCATION/libboost_timer.a $BOOSTSO_PATH_LOCATION/libboost_chrono.a -Wl,--whole-archive -L$CPLEXSO_PATH_LOCATION -lcplex12100 -lstdc++ -Wl, 
#,-all_load -L$COINORSO_PATH_LOCATION -lClp -lCgl -lOsiClp -lOsi -lCoinUtils
