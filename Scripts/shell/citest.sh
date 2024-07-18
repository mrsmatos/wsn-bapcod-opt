#!/bin/bash
set -e
date ; hostname ; pwd
cpu=`getconf _NPROCESSORS_ONLN`
echo "the PATH variable before source is: "
echo $PATH

#env
BRANCH=$CI_COMMIT_REF_NAME
echo "Current Branch name of base is:"
echo $BRANCH
echo -n "=== Launching boost, lemon compilation in background + all repos cloning ..." ; date
errfile="err-file-main"
rm -f $errfile
(bash -x Scripts/shell/citest_cloneapplications.sh            > cloneapp.log    2>&1 || (echo ; echo "=== ERROR on cloning applications" ; cat cloneapp.log    ; touch $errfile ; exit 1) ) &

#gccnewversion=$(gcc --version | grep ^gcc | sed 's/^.* //g')
#
#if [ -e Tools/gcc-version ]; then
#  gccoldversion=`cat Tools/gcc-version`
#else
#  gccoldversion="no"
#fi

#if [ "$gccoldversion" != "$gccnewversion" ]; then
#  echo "No gcc version saved or different version - $gccnewversion vs $gccoldversion - need to build"
#  rm -rf Tools/lemon-1.3.1 Tools/boost_1_76_0
#fi

if [ -e "$BOOST_ROOT" ]; then
  echo Try to use existing module with BOOST_ROOT=$BOOST_ROOT
  mkdir -p Tools/boost_1_76_0
  ln -s $BOOST_ROOT Tools/boost_1_76_0/build
fi

if [ -e "$LEMON_ROOT" ]; then
  echo Try to use existing module with LEMON_ROOT=$LEMON_ROOT
  mkdir -p Tools/lemon-1.3.1
  ln -s $LEMON_ROOT Tools/lemon-1.3.1/build
fi

if [ ! -e Tools/lemon-1.3.1 ]; then
  echo Need to compile Lemon
(bash    Scripts/shell/install_bc_lemon.sh                    > lemon.log       2>&1 || (echo ; echo "=== ERROR on install Lemon"        ; cat lemon.log       ; touch $errfile ; exit 1) ) &
fi

if [ ! -e Tools/boost_1_76_0 ]; then
  echo Need to compile Boost
(bash    Scripts/shell/install_bc_boost.sh                    > boost.log       2>&1 || (echo ; echo "=== ERROR on install BOOST"        ; cat boost.log       ; touch $errfile ; exit 1) ) &
fi

##DISABLE##(bash -x Scripts/shell/multirepos_cloneandcheckout.sh $BRANCH > clonemaster.log 2>&1 || (echo ; echo "=== ERROR on script multirepos clone checkout" ; cat clonemaster.log ; touch $errfile ; exit 1) ) &
echo -n "=== Waiting ... "
wait

if [ -e $errfile ]; then
    echo "+++++ Error detected in boost/lemon/all repos cloning +++++"
    #cat *.log
    exit 1
fi
echo -n " ... Done "; date

#echo $gccnewversion > Tools/gcc-version
err=0

rm Bapcod/include/bcModelNonPublicCuts.hpp \
   Bapcod/include_dev/bcNonPublicCuts.hpp \
   Bapcod/include_solverInterfaces/bcClpSolverC.hpp \
   Bapcod/include_solverInterfaces/bcGlpkSolverC.hpp \
   Bapcod/include_solverInterfaces/bcGRBSolverC.hpp \
   Bapcod/include_solverInterfaces/bcXpressMPSolverC.hpp \
   Bapcod/src/bcModelNonPublicCuts.cpp \
   Bapcod/src/bcNonPublicCapacityCuts.cpp \
   Bapcod/src/bcNonPublicCglCut.cpp \
   Bapcod/src/bcNonPublicCliqueCuts.cpp \
   Bapcod/src/bcNonPublicHomExtCapCuts.cpp \
   Bapcod/src/bcNonPublicOverlElimCuts.cpp \
   Bapcod/src/bcNonPublicResConsKnapsackCuts.cpp \
   Bapcod/src/bcNonPublicStrongKPathCuts.cpp \
   Bapcod/src_solverInterfaces/bcClpSolverC.cpp \
   Bapcod/src_solverInterfaces/bcGlpkSolverC.cpp \
   Bapcod/src_solverInterfaces/bcGRBSolverC.cpp \
   Bapcod/src_solverInterfaces/bcLpSolveSolverC.cpp \
   Bapcod/src_solverInterfaces/bcXpressMPSolverC.cpp \
   Bapcod/src_solverInterfaces/bcXpressSolverC.cpp


(cd Tools && wget https://gitlab.inria.fr/api/v4/projects/27629/packages/generic/rcsp-bin/0.0.2/rcsp.tgz && tar xzf rcsp.tgz && ln -s librcsp-*-* rcsp)
mkdir -p build_rel_rcsp 
cd build_rel_rcsp
echo "=== cmake ..."
cmake -DCMAKE_BUILD_TYPE=Release ..
echo "=== make -j$cpu bapcod ..."
make -j$cpu bapcod #CuttingStockDemo  GeneralizedAssignmentDemo  VehicleRoutingWithTimeWindowsDemo  VertexColoringDemo \
                   #BinPackingAndCuttingStock  CapacitatedVehicleRouting  GeneralizedAssignment       MultiCommodityFlow         VertexColoring
# TODO get separation-libs ! make -j$cpu bapcod CVRPSEP CliqueSep RheccSep 
echo -n "=== With RCSP - Compiling application and run tests ... " ; date
sh Scripts/shell/build_applications_and_run_tests.sh || err=1
echo -n " ... Done (with rcsp)"; date

cd ..

mv Applications/CapacitatedVehicleRouting         Applications/CapacitatedVehicleRouting-needs-rcsp
mv        Demos/VehicleRoutingWithTimeWindowsDemo        Demos/VehicleRoutingWithTimeWindowsDemo-needs-rcsp
mv Tools/rcsp  Tools/rcsp-disabled

mkdir -p build_rel 
cd build_rel
echo "=== cmake ..."
cmake -DCMAKE_BUILD_TYPE=Release ..
echo "=== make -j$cpu bapcod ..."
make -j$cpu bapcod #CuttingStockDemo  GeneralizedAssignmentDemo  VehicleRoutingWithTimeWindowsDemo  VertexColoringDemo \
                   #BinPackingAndCuttingStock  CapacitatedVehicleRouting  GeneralizedAssignment       MultiCommodityFlow         VertexColoring
# TODO get separation-libs ! make -j$cpu bapcod CVRPSEP CliqueSep RheccSep 
echo -n "=== Without RCSP - Compiling application and run tests ... " ; date
sh Scripts/shell/build_applications_and_run_tests.sh || err=1
echo -n " ... Done (without rcsp)"; date

exit $err
