#!/bin/sh

#This script must be run in BapcodProject folder (the root folder).

# APPLICATIONS="GeneralizedAssignment
# 	BinPackingAndCuttingStock
#               BinPackingWithConflicts
# 	ShiftSchedulingUsingGrammars
#               VertexColoring
#               ParallelMachineScheduling
# 	CapacitatedVehicleRouting
# 	MultiCommodityFlow
# 	UnrelatedParallelMachineTotalWeightedCompletionTime"

APPLICATIONS="
GeneralizedAssignment
VertexColoring
CapacitatedVehicleRouting
MultiCommodityFlow
"

NB_CORES=4 #Nb core to build the applications

# Configure the projects in Release mode
#cmake . -DCMAKE_BUILD_TYPE=Release

if [ ! -e Applications ]; then
    echo " Cloning the repository applications"
    git clone https://gitlab+deploy-token-13:kKNhN3HVX1mfM7y5C8-Y@plmlab.math.cnrs.fr:/realopt/applications.git Applications
    #git clone https://gitlab+deploy-token-8:m6jxG65PhxgHyUe6mPEa@plmlab.math.cnrs.fr:/realopt/applications-old-toremove Applications
fi
    
cd Applications

errfile="err-file-sub"
rm -f $errfile

# for app in ${APPLICATIONS}
# do
#     lc_app=`echo "${app}" | awk '{print tolower($0)}'`
#     if [ ! -e ${app} ]; then
# 	echo " Cloning the repository ${lc_app}"
#         # gitlab+deploy-token-7 bisbA-vJzDnL8AWWcPPR 
# 	( (git clone git@plmlab.math.cnrs.fr:realopt/${lc_app} ${app})  > $app.log    2>&1 || (echo ; echo "=== ERROR on clining $app" ; cat $app.log    ; touch $errfile ; exit 1) ) &
#     fi
# done

if [ ! -e CapacitatedVehicleRouting ]; then
    ( (git clone https://ci-code-read-only:Jv3pdtJzmo3ex_3Kxs6G@gitlab.inria.fr/sadykov/CapacitatedVehicleRouting.git CapacitatedVehicleRouting > CapacitatedVehicleRouting.log    2>&1) || (echo ; echo "=== ERROR on cloning CapacitatedVehicleRouting" ; cat CapacitatedVehicleRouting.log    ; touch $errfile ; exit 1) ) &
fi

if [ ! -e BinPackingAndCuttingStock ]; then
    ( (git clone https://ci-code-read-only:7tfck6nnqENA4-3KKDUf@gitlab.inria.fr/sadykov/BinPackingAndCuttingStock.git BinPackingAndCuttingStock > BinPackingAndCuttingStock.log    2>&1) || (echo ; echo "=== ERROR on cloning BinPackingAndCuttingStock" ; cat BinPackingAndCuttingStock.log    ; touch $errfile ; exit 1) ) &
fi

if [ ! -e GeneralizedAssignement ]; then
    ( (git clone 'https://ci-code-read-only:92uD_BqA8iFjKP6_JhzR@gitlab.inria.fr/sadykov/GeneralizedAssignement.git' GeneralizedAssignement > GeneralizedAssignement.log    2>&1) || (echo ; echo "=== ERROR on cloning GeneralizedAssignement" ; cat GeneralizedAssignement.log    ; touch $errfile ; exit 1) ) &
fi

if [ ! -e VertexColoring ]; then
    ( (git clone 'https://ci-code-read-only:S4UJeDV2ejL-vxh_vjU7@gitlab.inria.fr/sadykov/VertexColoring.git' VertexColoring > VertexColoring.log    2>&1) || (echo ; echo "=== ERROR on cloning VertexColoring" ; cat VertexColoring.log    ; touch $errfile ; exit 1) ) &
fi

if [ ! -e MultiCommodityFlow ]; then
    ( (git clone 'https://ci-code-read-only:8zzN287w6sLmQy86gEJv@gitlab.inria.fr/sadykov/MultiCommodityFlow.git' MultiCommodityFlow > MultiCommodityFlow.log    2>&1) || (echo ; echo "=== ERROR on cloning MultiCommodityFlow" ; cat MultiCommodityFlow.log    ; touch $errfile ; exit 1) ) &
fi

# if [ ! -e VertexColoring ]; then
# ( (git clone https://gitlab+deploy-token-9:wzW61HrSy31_yZ8TLDtj@gitlab.inria.fr:/rsadykov/vertexcoloring.git VertexColoring > VertexColoring.log    2>&1) || (echo ; echo "=== ERROR on cloning VertexColoring" ; cat VertexColoring.log    ; touch $errfile ; exit 1) ) &
# fi

# if [ ! -e MultiCommodityFlow ]; then
# ( (git clone https://gitlab+deploy-token-10:mpWAzxdJP3_MRyvpe-LM@gitlab.inria.fr:/rsadykov/multicommodityflow.git MultiCommodityFlow > MultiCommodityFlow.log    2>&1) || (echo ; echo "=== ERROR on cloning MultiCommodityFlow" ; cat MultiCommodityFlow.log    ; touch $errfile ; exit 1) ) &
# fi

# if [ ! -e BinPackingAndCuttingStock ]; then
# ( (git clone https://gitlab+deploy-token-12:X7bPJzboaxr-2YaYErcM@gitlab.inria.fr:/rsadykov/BinPackingAndCuttingStock.git BinPackingAndCuttingStock > BinPackingAndCuttingStock.log    2>&1) || (echo ; echo "=== ERROR on cloning BinPackingAndCuttingStock" ; cat MultiCommodityFlow.log    ; touch $errfile ; exit 1) ) &
# fi


wait

if [ -e $errfile ]; then
    echo "+++++ Error detected in applications cloning +++++"
    cat *.log
    exit 1
fi


#git clone git@gitlab.inria.fr:realopt/UnrelatedParallelMachineTotalWeightedCompletionTime
