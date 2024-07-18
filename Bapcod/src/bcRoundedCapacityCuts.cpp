/**
 *
 * This file bcRoundedCapacityCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

//
//  bcCapacityCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 11/09/2017.
//
//

#include "bcNetworkBasedCuts.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"

GenericRCSPCapacityCutConstr::GenericRCSPCapacityCutConstr(Model * modelPtr,
                                                           ProbConfig * probConfPtr,
                                                           const std::string & name,
                                                           const Double & nonRootPriorityLevel,
                                                           const Double & rootPriorityLevel,
                                                           const bool & isFacultative,
                                                           const bool & equalityCase,
                                                           const int & maxCapacity,
                                                           const std::vector<int> & demands,
                                                           const int & twoPathCutsResId) :
        GenericCutConstr(modelPtr, probConfPtr, name, isFacultative ? 'F' : 'C',
                         SelectionStrategy::MostViolated, nonRootPriorityLevel, rootPriorityLevel),
        _maxCapacity(maxCapacity), _equalityCase(equalityCase), _demands(demands), _twoPathCutsResId(twoPathCutsResId),
        _interfacePtr(nullptr), _cgSpConfPts()
{

}
GenericRCSPCapacityCutConstr::~GenericRCSPCapacityCutConstr()
{
    delete _interfacePtr;
}

bool GenericRCSPCapacityCutConstr::prepareSeparation()
{
    bcp_rcsp::RoundCapCutsSeparatorParameters cutSepParams;
    cutSepParams.maxNumPerRound = param().RCSPcapCutsMaxNumPerRound();
    cutSepParams.equalityCase = _equalityCase;
    cutSepParams.printLevel = param().DEFAULTPRINTLEVEL();
    cutSepParams.cutViolationTolerance = param().BapCodCutViolationTolerance();
    if (printL(0))
        cutSepParams.printLevel = 1;
    else
        cutSepParams.printLevel = 0;
    std::vector<const bcp_rcsp::GraphData *> graphs;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() != NULL) {
            graphs.push_back((*cgSpConfPtrIt)->rcspGraphPtr());
            int graphId = (*cgSpConfPtrIt)->rcspGraphPtr()->id;
            if (graphId >= (int)_cgSpConfPts.size())
                _cgSpConfPts.resize(graphId + 1);
            _cgSpConfPts[graphId] = *cgSpConfPtrIt;
        }
    }

    std::vector<double> demands(_demands.size());
    for (int index = 0; index < (int)_demands.size(); ++index)
        demands[index] = (double)_demands[index];
    _interfacePtr = bcp_rcsp::createAndPrepareRoundCapCutSeparation(graphs, demands, (double)_maxCapacity,
                                                                    _twoPathCutsResId, cutSepParams);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "BaPCod error : could not prepare rounded capacity cuts separation" << std::endl;
        return false;
    }

    return true;

}

void GenericRCSPCapacityCutConstr
::cutSeparationRoutine(const VarPtrSet & curSol,
                       std::multiset <InstanciatedConstr * , CutSeparationPriorityComp > & generatedCutConstrSet)
{
    bool printCuts = printL(2);

    bcp_rcsp::FractionalMasterSolution rcspFracSolution;
    rcspFracSolution.solPts.reserve(curSol.size());
    rcspFracSolution.values.reserve(curSol.size());
    for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
        if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;

        MastColumn *colPtr = static_cast<MastColumn *>(*varPtrIt);

        auto rcspSolPtr = colPtr->spSol()->rcspSolPtr();
        if (rcspSolPtr != nullptr)
        {
            rcspFracSolution.solPts.push_back(rcspSolPtr);
            rcspFracSolution.values.push_back(colPtr->val());
        }
    }

    std::vector<const bcp_rcsp::RoundedCapacityCut *> cutPts;
    if (!_interfacePtr->separate(rcspFracSolution, cutPts))
        return;

    for (auto * rcspCutPtr : cutPts)
    {
        if (printCuts)
        {
            std::cout << "New capacity cut " << rcspCutPtr->id << " for set ";
            for (const auto setId : rcspCutPtr->setIds)
                std::cout << setId << " ";
            std::cout << ":\n";
        }

        std::stringstream ss;
        ss << "CAP_" << rcspCutPtr->id;
        MultiIndex newCutId(rcspCutPtr->id);
        char inequalitySense = (rcspCutPtr->sense == LESS_OR_EQUAL_SENSE) ? 'L' : 'G';
        InstMasterConstr * newCutPtr = new InstMasterConstr(newCutId, this, probConfPtr(),
                                                            ss.str(), (double)rcspCutPtr->rhs, inequalitySense,
                                                            defaultType(), defaultKind(), defaultFlag());
        for (const auto & triple : rcspCutPtr->varsMembership)
        {
            InstanciatedVar * instVarPtr = _cgSpConfPts[std::get<0>(triple)]->instVarPts()[std::get<1>(triple)];
            newCutPtr->includeMember(instVarPtr, std::get<2>(triple), false);
            if (printCuts)
                std::cout << ((std::get<2>(triple) > 0.0) ? "+" : "") << std::get<2>(triple) << "*"
                          << instVarPtr->name();
        }
        if (printCuts)
            std::cout << ((inequalitySense == 'L') ? " <= " : " >= ") << rcspCutPtr->rhs << std::endl;

        generatedCutConstrSet.insert(newCutPtr);

        delete rcspCutPtr;
    }
    cutPts.clear();
}
#endif

