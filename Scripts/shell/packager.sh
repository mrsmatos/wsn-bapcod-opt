#! /bin/bash

set -e 

rm -rf pa-ck*

git clone git@gitlab.inria.fr:realopt/bapcodframework.git pa-ck-source
cd pa-ck-source

mkdir -p pa-ck

cp -r CMakeLists.txt \
      Demos \
      License.pdf \
      README.md \
      Tools \
      Bapcod \
      Applications \
   pa-ck

rm pa-ck/Bapcod/include/bcModelNonPublicCuts.hpp \
   pa-ck/Bapcod/include_dev/bcNonPublicCuts.hpp \
   pa-ck/Bapcod/include_solverInterfaces/bcCbcSolverC.hpp \
   pa-ck/Bapcod/include_solverInterfaces/bcGlpkSolverC.hpp \
   pa-ck/Bapcod/include_solverInterfaces/bcGRBSolverC.hpp \
   pa-ck/Bapcod/include_solverInterfaces/bcXpressMPSolverC.hpp \
   pa-ck/Bapcod/src/bcModelNonPublicCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicCapacityCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicCglCut.cpp \
   pa-ck/Bapcod/src/bcNonPublicCliqueCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicHomExtCapCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicOverlElimCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicResConsKnapsackCuts.cpp \
   pa-ck/Bapcod/src/bcNonPublicStrongKPathCuts.cpp \
   pa-ck/Bapcod/src_solverInterfaces/bcGlpkSolverC.cpp \
   pa-ck/Bapcod/src_solverInterfaces/bcGRBSolverC.cpp \
   pa-ck/Bapcod/src_solverInterfaces/bcLpSolveSolverC.cpp \
   pa-ck/Bapcod/src_solverInterfaces/bcXpressMPSolverC.cpp \
   pa-ck/Bapcod/src_solverInterfaces/bcXpressSolverC.cpp

mkdir -p pa-ck/Scripts/shell

cp Scripts/shell/{install_bc_lemon.sh,install_bc_boost.sh,install_bc_lemon_mac_m1.sh,install_bc_boost_mac_m1.sh} \
   pa-ck/Scripts/shell

mkdir -p pa-ck/Scripts/python
mkdir -p pa-ck/Scripts/python/TestOutputSummary

cp Scripts/python/TestOutputSummary/testOutputSummary.py \
   pa-ck/Scripts/python/TestOutputSummary

(cd pa-ck; tar -cf - . | gzip -9 > ../../bapcod.tar.gz)

rm -rf pa-ck*

cd ..
rm -rf pa-ck-source
