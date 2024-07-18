/**
 *
 * This file bcModelRCSPSolver.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *r
 **/

//
//  bcModelRCSPSolver.cpp
//  Project
//
//  Created by Ruslan Sadykov on 05/10/2017.
//
//

#include "bcNetworkFlowC.hpp"
#include "bcProbConfigC.hpp"
#include "bcModelC.hpp"
#include "bcModelRCSPSolver.hpp"
#include "bcNetworkBasedCuts.hpp"
#include "bcNetworkBasedBranching.hpp"

#ifdef USE_NON_PUBLIC_CUTS
#include "bcModelNonPublicCuts.hpp"
#include "bcNonPublicCuts.hpp"
#endif

#ifdef BCP_RCSP_IS_FOUND

namespace bcp_rcsp {
  struct GraphData;
}

bool runRCSPoracleAsStandaloneSolver(const BcModel & bcModel, const std::string & fileName)
{
    std::ifstream ifs(fileName.c_str(), std::ios::in);

    if(!ifs)
    {
        std::cerr << "RCSP solver error : cannot find standalone RCSP input file " << fileName <<std::endl;
        return false;
    }
    std::string line;
    for (int lineNum = 0; lineNum < 6; ++lineNum)
    {
        std::getline(ifs, line);
        if (ifs.eof() || ifs.fail())
        {
            std::cerr << "RCSP solver error : standalone RCSP input file " << fileName << " in the wrong format "
                      << std::endl;
            return false;
        }
    }

    int nbMainRes, nbResources;
    ifs >> nbMainRes >> nbResources;
    if (ifs.eof() || ifs.fail())
    {
        std::cerr << "RCSP solver error : standalone RCSP input file " << fileName << " in the wrong format "
                  << std::endl;
        return false;
    }

    ifs.close();

    Model * modelPtr = (Model *)bcModel;
    BcFormulation bcForm;
    BcRCSPFunctor rcspFunctor(bcForm, modelPtr->param());

    return rcspFunctor.runAsStandaloneRCSPsolver(fileName);
}


/*****************************************************************************/
/*****************************************************************************/
/*****                          TODO: RCSP solver                        *****/
/*****************************************************************************/
/*****************************************************************************/

double RCSPLabelExtensionCostFunctor::getResConsDependentCost(const int arcId,
                                                              const std::vector<double> & resConsumption,
                                                              const bool isTailResCons) const
{
    if (_networkPtr == NULL)
        return 0.0;
    const BcArcInfo * arcInfoPtr = _networkPtr->getArcInfoPtr((int)arcId);
    return operator()(arcInfoPtr, &(resConsumption[0]), isTailResCons);
}

RCSPLabelExtensionCostFunctor::RCSPLabelExtensionCostFunctor(const BcFormulation & bcForm, const BcVar & costVar) :
    _networkPtr(NULL), _costInstVarPtr((InstanciatedVar *)costVar)
{
    if (bcForm.isDefined())
      _networkPtr = bcForm.probConfPtr()->networkFlowPtr();
}

RCSPFeasibilityCheckFunctor::RCSPFeasibilityCheckFunctor(const BcFormulation & bcForm) :
        _networkPtr(NULL)
{
    if (bcForm.isDefined())
        _networkPtr = bcForm.probConfPtr()->networkFlowPtr();
}

bool RCSPFeasibilityCheckFunctor::isFeasible(const bcp_rcsp::Solution & rcspSol) const
{
    if ((_networkPtr == NULL) || rcspSol.arcIds.empty())
        return true;

    std::vector<int> orderedVertIds;
    const BcArcInfo * arcInfoPtr = _networkPtr->getArcInfoPtr(rcspSol.arcIds.front());
    orderedVertIds.push_back(arcInfoPtr->tailVertId);
    for (auto arcId : rcspSol.arcIds)
    {
        arcInfoPtr = _networkPtr->getArcInfoPtr(arcId);
        orderedVertIds.push_back(arcInfoPtr->headVertId);
    }
    return operator()(orderedVertIds, rcspSol.waitingTimes, rcspSol.resConsumption);
}

void BcRCSPFunctor::fillRCSPSolverParameters(const ControlParameters & bapcodParams)
{
    _solverParams.checkDominInOtherBuckets = bapcodParams.RCSPcheckDominInOtherBuckets();
    _solverParams.useBidirectSearchInPricing = bapcodParams.RCSPuseBidirectionalSearch();
    if (bapcodParams.RCSPuseBidirectionalSearch())
        _solverParams.useBidirectSearchInEnumeration = RCSP_BIDIR_ENUM_SEARCH;
    else
        _solverParams.useBidirectSearchInEnumeration = RCSP_MONODIR_ENUM_SEARCH;
    _solverParams.useCompletionBoundsInPricing = bapcodParams.RCSPuseCompletionBoundsInPricing();
    _solverParams.domChecksThresholdInPricing = bapcodParams.RCSPdomChecksThresholdInPricing();
    _solverParams.stopCutGenTimeThresholdInPricing = bapcodParams.RCSPstopCutGenTimeThresholdInPricing();
    _solverParams.hardTimeThresholdInPricing = bapcodParams.RCSPhardTimeThresholdInPricing();
    _solverParams.redCostFixingTimeThreshold = bapcodParams.RCSPredCostFixingTimeThreshold();
    _solverParams.maxNumOfColsPerIteration = bapcodParams.RCSPmaxNumOfColsPerIteration();
    _solverParams.maxNumOfColsPerExactIteration = bapcodParams.RCSPmaxNumOfColsPerExactIteration();
    _solverParams.allowRoutesWithSameVerticesSet = bapcodParams.RCSPallowRoutesWithSameVerticesSet();
    _solverParams.heurLabelingStrategy = bapcodParams.RCSPheurLabelingStrategy();
    _solverParams.labelSplitStrategy = bapcodParams.RCSPlabelSplitStrategy();
    _solverParams.numberOfBucketsPerVertex = bapcodParams.RCSPnumberOfBucketsPerVertex();
    _solverParams.dynamicBucketSteps = bapcodParams.RCSPdynamicBucketSteps();
    _solverParams.dynBuckStepAdjustRatioThreshold = bapcodParams.RCSPdynBuckStepAdjustRatioThreshold();
    _solverParams.dynBuckStepAdjustNumBuckArcsThreshold = bapcodParams.RCSPdynBuckStepAdjustNumBuckArcsThreshold();
    _solverParams.applyReducedCostFixing = bapcodParams.RCSPapplyReducedCostFixing();
    _solverParams.redCostFixingFalseGap = bapcodParams.RCSPredCostFixingFalseGap();
    _solverParams.useComplBoundsInRedCostFixing = bapcodParams.RCSPuseComplBoundsInRedCostFixing();
    _solverParams.maxNumOfColumnsInDSSR = bapcodParams.RCSPmaxNumOfColumnsInDSSR();
    _solverParams.maxCycleSizeInDSSR = bapcodParams.RCSPmaxCycleSizeInDSSR();
    _solverParams.dynamicNGmode = std::abs(bapcodParams.RCSPdynamicNGmode());
    _solverParams.maxNumOfColumnsInDynamicNG = bapcodParams.RCSPmaxNumOfColumnsInDynamicNG();
    _solverParams.maxCycleSizeInDynamicNG = bapcodParams.RCSPmaxCycleSizeInDynamicNG();
    _solverParams.initNGneighbourhoodSize = bapcodParams.RCSPinitNGneighbourhoodSize();
    _solverParams.maxNGneighbourhoodSize = bapcodParams.RCSPmaxNGneighbourhoodSize();
    _solverParams.maxNGaverNeighbourhoodSize = bapcodParams.RCSPmaxNGaverNeighbourhoodSize();
    _solverParams.useExactComplBoundsInEnumeration = bapcodParams.RCSPuseExactComplBoundsInEnumeration();
    _solverParams.maxNumOfLabelsInEnumeration = bapcodParams.RCSPmaxNumOfLabelsInEnumeration();
    _solverParams.maxNumOfLabelsInHeurEnumeration = bapcodParams.RCSPmaxNumOfLabelsInHeurEnumeration();
    _solverParams.maxNumOfEnumeratedSolutions = bapcodParams.RCSPmaxNumOfEnumeratedSolutions();
    _solverParams.useVertexIdsInEnumeration = bapcodParams.RCSPuseVertexIdsInEnumeration();
    _solverParams.stopRatioForConcatenationInEnum = bapcodParams.RCSPstopRatioForConcatenationInEnum();
    _solverParams.generatePosRedCostPaths = bapcodParams.InsertNewNonNegColumnsDirectlyInFormRatherThanInPool();
    _solverParams.toleranceInDSSR = bapcodParams.optimalityGapTolerance();
    _solverParams.printLevel = bapcodParams.RCSPprintLevel();
    _solverParams.labelSplitStrategy = bapcodParams.RCSPlabelSplitStrategy();
}

BcRCSPFunctor::BcRCSPFunctor(const BcFormulation & spForm, const ControlParameters & bapcodParams) :
    _probConfPtr(spForm.probConfPtr()), _graphId(spForm.probConfPtr()->id().first()),
    _useMetaSolver(bapcodParams.RCSPuseMetaSolver() > 0),
    _solverParams(), _solverInterface(NULL), _solverColGenPhaseConfig(), _verificationFunctorPtr(NULL),
    _feasibilityCheckFunctorPtr(NULL), _labelExtensionCostFunctorPtr(NULL),
    _messageIdToCutGeneration(PricingSolverCutsMessage::noMessage), _paramFormulationIsBuilt(true),
    _pureCostInstVarPtr(NULL)
{
    fillRCSPSolverParameters(bapcodParams);
}

BcRCSPFunctor::BcRCSPFunctor(const BcFormulation & spForm) :
    BcRCSPFunctor::BcRCSPFunctor(spForm, spForm.probConfPtr()->param())
{
}

BcRCSPFunctor::~BcRCSPFunctor()
{
    clearCutPts();
    delete _solverInterface;
    delete _verificationFunctorPtr;
    delete _labelExtensionCostFunctorPtr;
    delete _feasibilityCheckFunctorPtr;
}

bool BcRCSPFunctor::runAsStandaloneRCSPsolver(const std::string & fileName, const int colGenPhase)
{
    if (!_probConfPtr->fillRCSPGraph())
        return false;

    bcp_rcsp::Data data(*_probConfPtr->rcspGraphPtr(), _solverParams);
    data.colGenPhasesConfig = _solverColGenPhaseConfig;
    /// we do not pass callbacks to the solver, as they are not used when running from file

    return bcp_rcsp::createAndRunFromFile(data, fileName, colGenPhase);
}

bool BcRCSPFunctor::prepareSolver()
{
    if ((_verificationFunctorPtr != NULL) && !_verificationFunctorPtr->prepareSolver())
        return false;

    if (_probConfPtr->rcspGraphPtr() == NULL)
    {
        std::cerr << "RCSP functor error: RCSP graph was not build by BaPCod" << std::endl;
        return false;
    }

    bcp_rcsp::Data data(*_probConfPtr->rcspGraphPtr(), _solverParams);
    if (_verificationFunctorPtr != nullptr)
    {
        data.verificationSolverPtr = _verificationFunctorPtr->_solverInterface;
        _verificationFunctorPtr->_solverInterface = NULL;
    }
    data.labelExtensionCostFunctorPtr = _labelExtensionCostFunctorPtr;
    _labelExtensionCostFunctorPtr = NULL;
    data.feasibilityCheckFunctorPtr = _feasibilityCheckFunctorPtr;
    _feasibilityCheckFunctorPtr = NULL;
    data.colGenPhasesConfig = _solverColGenPhaseConfig;

    if (_useMetaSolver)
        _solverInterface = bcp_rcsp::createAndPrepareMetaSolver(data);
    else
        _solverInterface = bcp_rcsp::createAndPrepareSolver(data);
    return (_solverInterface != NULL);
}

bool BcRCSPFunctor::setupNode(BcFormulation spPtr, const BcSolverOracleInfo * infoPtr)
{
    _messageIdToCutGeneration = PricingSolverCutsMessage::noMessage;

    std::vector<const bcp_rcsp::AccumResConsBranchConstraint *> accumResConsBranchCtrPts;
    getPackSetResConsActiveBranchConstrList(spPtr, accumResConsBranchCtrPts);

    std::vector<const bcp_rcsp::RyanFosterBranchConstraint *> ryanFosterBranchConstrPts;
    getPackSetRyanFosterActiveBranchConstrList(spPtr, ryanFosterBranchConstrPts);

    const bcp_rcsp::SolverRecord * solverStatePtr = NULL;
    const RCSPOracleInfo * rcspInfoPtr = dynamic_cast<const RCSPOracleInfo *>(infoPtr);
    if (rcspInfoPtr != NULL)
        solverStatePtr = rcspInfoPtr->getSolverState();
    if (!_solverInterface->restoreFromRecord(solverStatePtr, spPtr.debugSolutionIsValidAtThisNode(),
                                             accumResConsBranchCtrPts, ryanFosterBranchConstrPts))
        return true; /// true return status means infeasible or "something wrong"
    return false;
}

BcSolverOracleInfo * BcRCSPFunctor::recordSolverOracleInfo(const BcFormulation spPtr)
{
    const bcp_rcsp::SolverRecord * solverRecordPtr = _solverInterface->saveToRecord();
    return new RCSPOracleInfo(solverRecordPtr);
}

void BcRCSPFunctor::addToSolution(const NetworkFlow * networkPtr, const bcp_rcsp::Solution * rcspSolPtr,
                                      const std::unordered_map<int, InstanciatedVar *> & resIdToVarPtrMap,
                                      BcSolution & bcSol) const
{
    bcSol.solutionPtr()->setRcspSolPtr(rcspSolPtr);

    /// variables associated with arcs
    for (const auto & arcId : rcspSolPtr->arcIds)
    {
        const NetworkArc * netArcPtr = networkPtr->netArcPtr(arcId);
        if (netArcPtr->varToCoeffMaps().size() == 1)
        {
            for (auto & varCoeffPair: netArcPtr->varToCoeffMaps().front())
                bcSol.updateVarVal(BcVar(varCoeffPair.first), varCoeffPair.second);
        }
        else
        {
            double minReducedCost = 0.0;
            int bestMappingId = 0;
            auto mappingIt = netArcPtr->varToCoeffMaps().begin();
            for (auto & varCoeffPair: *mappingIt)
                minReducedCost += bcp_rcsp::roundP(varCoeffPair.first->curCost()) * varCoeffPair.second;
            ++mappingIt;
            int mappingId = 1;
            while (mappingIt != netArcPtr->varToCoeffMaps().end())
            {
                double thisReducedCost = 0.0;
                for (auto & varCoeffPair: *mappingIt)
                    thisReducedCost += bcp_rcsp::roundP(varCoeffPair.first->curCost()) * varCoeffPair.second;
                if (minReducedCost > thisReducedCost)
                {
                    minReducedCost = thisReducedCost;
                    bestMappingId = mappingId;
                }
                ++mappingIt;
                mappingId += 1;
            }
            for (auto & varCoeffPair: netArcPtr->varToCoeffMaps()[bestMappingId])
                bcSol.updateVarVal(BcVar(varCoeffPair.first), varCoeffPair.second);
        }
    }

    /// variable associated with pure arc costs
    if (_pureCostInstVarPtr != nullptr)
        bcSol.updateVarVal(BcVar(_pureCostInstVarPtr), rcspSolPtr->pureArcCost);

    /// variable associated with resources
    for (auto & resIsVarIdPair: resIdToVarPtrMap)
    {
        const std::vector<double> & finalResConsumption = rcspSolPtr->resConsumption.back();
        bcSol.updateVarVal(BcVar(resIsVarIdPair.second), finalResConsumption[resIsVarIdPair.first]);
    }
}

bool BcRCSPFunctor::fillRCSPInput(BcFormulation spPtr, const int colGenPhase,
                                  const std::vector<InstanciatedVar *> & instVarPts, bcp_rcsp::SolverInput & input)
{
    input.rollbackStateIsSaved = spPtr.rollbackPointSavedStatus();
    input.zeroRedCostThreshold = spPtr.zeroReducedCostThreshold();
    input.colGenPhase = colGenPhase;
    input.checkDebugSolution = spPtr.debugSolutionIsValidAtThisNode();

    int numVariables = (int)instVarPts.size();
    for (int varId = 0; varId < numVariables; ++varId)
    {
        InstanciatedVar * instVarPtr = instVarPts[varId];
        if ((_paramFormulationIsBuilt && !instVarPtr->inCurProb()) || (instVarPtr->curUb() <= 0))
        {
            input.varRedCosts[varId] = BapcodInfinity;
        }
        else if (instVarPtr->curLb() > Double::precision)
        {
            std::cerr << "RCSP functor error: non-zero lower bounds on variables associated to arcs are not supported"
                      << std::endl;
            return false;
        }
        else
        {
            input.varRedCosts[varId] = instVarPtr->curCost();
        }
    }

    getActiveRankOneCuts(spPtr, input.rankOneCuts);

    getActiveDiscreteCuts(spPtr, _discreteCutPts);
    for (const auto cutPtr : _discreteCutPts)
        input.discreteCuts.push_back(std::make_pair(cutPtr, cutPtr->curDualVal()));

    getActiveRouteLoadKnapsackCuts(spPtr, input.routeLoadKnapCuts);

#ifdef USE_NON_PUBLIC_CUTS
    getActiveStrongKPathCuts(spPtr, input.strongKPathCuts);

    clearCutPts();

#ifdef CLIQUE_SEP_IS_FOUND
    getActiveCliqueCuts(spPtr, _cliqueCutPts);
    for (const auto cutPtr : _cliqueCutPts)
    {
        input.cliqueCuts.push_back(std::make_pair(cutPtr, cutPtr->curDualVal()));
    }
#endif // CLIQUE_SEP_IS_FOUND

#endif // USE_NON_PUBLIC_CUTS

    return true;
}

void BcRCSPFunctor::clearCutPts()
{
#ifdef USE_NON_PUBLIC_CUTS
#ifdef CLIQUE_SEP_IS_FOUND
    for (auto cutPtr : _cliqueCutPts)
        delete cutPtr;
    _cliqueCutPts.clear();
#endif
#endif

    for (auto cutPtr : _discreteCutPts)
        delete cutPtr;
    _discreteCutPts.clear();
}

bool BcRCSPFunctor::operator() (IN_ BcFormulation spPtr,
                                IN_ int colGenPhase,
                                OUT_ double & objVal,
                                OUT_ double & dualBound,
                                OUT_ BcSolution & primalSol)
{
    _messageIdToCutGeneration = PricingSolverCutsMessage::noMessage;

    const std::vector<InstanciatedVar *> & instVarPts = spPtr.probConfPtr()->instVarPts();
    bcp_rcsp::SolverInput input((int)instVarPts.size());
    if (!fillRCSPInput(spPtr, colGenPhase, instVarPts, input))
    {
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        return false;
    }

    if (colGenPhase == 0)
    {
        auto & instVarRedCosts = spPtr.probConfPtr()->instVarRedCosts();
        instVarRedCosts.resize(input.varRedCosts.size());
        for (size_t varIndex = 0; varIndex < input.varRedCosts.size(); ++varIndex)
        {
            instVarRedCosts[varIndex] = input.varRedCosts[varIndex];

//            auto * iVarPtr = spPtr.probConfPtr()->instVarPts()[varIndex];
//            if (iVarPtr->id().first() == 2)
//                std::cout << " RC[" << iVarPtr->name() << "]=" << instVarRedCosts[varIndex];
        }
    }
    
    bcp_rcsp::SolverOutput output;

    if (!_solverInterface->runPricing(input, output))
    {
        if (output.shouldPerformRollback)
        {
            _messageIdToCutGeneration = PricingSolverCutsMessage::doCutsRollback;
        }
        else
        {
            std::cerr << "RCSP functor error: could not properly terminate the pricing" << std::endl;
            _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        }
        objVal = dualBound = BapcodInfinity;
        return false;
    }

    if (output.solPts.empty())
    {
        objVal = dualBound = spPtr.zeroReducedCostThreshold();
        return false;
    }

    NetworkFlow * networkPtr = ((NetworkFlow *)spPtr.network());
    const std::unordered_map<int, InstanciatedVar *> & resIdToVarIdMap = spPtr.probConfPtr()->resIdToVarIdMap();
    auto rcspSolIt = output.solPts.begin();
    addToSolution(networkPtr, *rcspSolIt, resIdToVarIdMap, primalSol);

    /// we now rectify the solution value due to rounding (roundP) used by RCSP solver
    for (auto & varCoeffPair : primalSol.solutionPtr()->solVarValMap())
        output.bestSolutionCost += (varCoeffPair.first->curCost() - bcp_rcsp::roundP(varCoeffPair.first->curCost()))
                                   * varCoeffPair.second;

    for (++rcspSolIt; rcspSolIt != output.solPts.end(); ++rcspSolIt)
    {
        BcSolution additionalSol(spPtr);
        addToSolution(networkPtr, *rcspSolIt, resIdToVarIdMap, additionalSol);
        primalSol.appendSol(additionalSol);
    }

    dualBound = objVal = output.bestSolutionCost;

    return true;
}

void BcRCSPFunctor::reducedCostFixingAndEnumeration(IN_ BcFormulation spPtr,
                                                    IN_ const int & enumerationMode,
                                                    IN_ const double & threshold)
{
    const std::vector<InstanciatedVar *> & instVarPts = spPtr.probConfPtr()->instVarPts();
    bcp_rcsp::SolverInput input(instVarPts.size());
    if (!fillRCSPInput(spPtr, 0, instVarPts, input))
    {
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        return;
    }
    input.fixingThreshold = threshold;

    if (!_solverInterface->runRedCostFixingAndEnumeration(input, enumerationMode))
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
}

int BcRCSPFunctor::getMessageIdToCutGeneration() const
{
    return _messageIdToCutGeneration;
}

bool BcRCSPFunctor::getEnumeratedStatus() const
{
    return _solverInterface->getEnumeratedStatus();
}

int BcRCSPFunctor::getNumberOfEnumeratedSolutions() const
{
    return _solverInterface->getNumberOfEnumeratedSolutions();
}

void BcRCSPFunctor::checkEnumeratedSolutions(BcFormulation spPtr, const std::vector<Solution *> & solPts,
                                             std::vector<bool> & solIsEnumerated)
{
    std::vector<const bcp_rcsp::Solution *> rcspSolsPts;
    rcspSolsPts.reserve(solPts.size());
    for (const auto * solPtr : solPts)
        rcspSolsPts.push_back(solPtr->rcspSolPtr());
    _solverInterface->checkEnumeratedSolutions(rcspSolsPts, solIsEnumerated);
}

bool BcRCSPFunctor::setDebugSolution(const std::vector<std::vector<int> > & ids, bool vertexBased)
{
    return _solverInterface->setDebugPaths(ids, vertexBased);
}

bool BcRCSPFunctor::getDebugSolution(BcFormulation spPtr, BcSolution & primalSol)
{
    std::vector<bcp_rcsp::Solution *> rcspSolsPts;
    if (!_solverInterface->getDebugSolutions(rcspSolsPts))
    {
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        return false;
    }

    NetworkFlow * networkPtr = ((NetworkFlow *)spPtr.network());
    const std::unordered_map<int, InstanciatedVar *> & resIdToVarIdMap = spPtr.probConfPtr()->resIdToVarIdMap();
    for (auto pathPtr : rcspSolsPts)
    {
        BcSolution bcSol(spPtr);
        addToSolution(networkPtr, pathPtr, resIdToVarIdMap, bcSol);
        primalSol += bcSol;
    }

    return true;
}

void BcRCSPFunctor::getEnumeratedSolutions(IN_ BcFormulation spPtr,
                                           IN_ const int & maxNumberOfSolutions,
                                           OUT_ BcSolution & enumeratedSol,
                                           OUT_ std::vector<double> & reducedCosts)
{
    const std::vector<InstanciatedVar *> & instVarPts = spPtr.probConfPtr()->instVarPts();
    bcp_rcsp::SolverInput input((int)instVarPts.size());
    if (!fillRCSPInput(spPtr, 0, instVarPts, input))
    {
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        return;
    }

    std::vector<bcp_rcsp::Solution *> rcspSolPts;
    if (!_solverInterface->getEnumeratedSolutions(input, maxNumberOfSolutions, rcspSolPts, reducedCosts))
    {
        _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
        return;
    }

    if (rcspSolPts.empty())
        return;

    NetworkFlow * networkPtr = ((NetworkFlow *)spPtr.network());
    const std::unordered_map<int, InstanciatedVar *> & resIdToVarIdMap = spPtr.probConfPtr()->resIdToVarIdMap();
    addToSolution(networkPtr, rcspSolPts.front(), resIdToVarIdMap, enumeratedSol);
    /// We add solutions in the reverse order as "enumeratedSol += additionalSol" inserts solution in the
    /// beginning of the solutions chain rooted at enumeratedSol
    auto pathIt = rcspSolPts.end();
    for (--pathIt; pathIt != rcspSolPts.begin(); --pathIt)
    {
        BcSolution additionalSol(spPtr);
        addToSolution(networkPtr, *pathIt, resIdToVarIdMap, additionalSol);
        enumeratedSol += additionalSol;
    }
}

bool BcRCSPFunctor::isProperSolution(IN_ BcFormulation spPtr, IN_ BcSolution & bcSol)
{
    bcp_rcsp::Solution rcspSol(_graphId);
    rcspSol.arcIds = bcSol.orderedIds();
    return _solverInterface->isProper(rcspSol);
}

bool BcRCSPFunctor::solSatisfiesCurrentSpRelaxation(IN_ BcFormulation spPtr, IN_ const BcSolution & bcSol)
{
    bcp_rcsp::Solution rcspSol(_graphId);
    rcspSol.arcIds = bcSol.orderedIds();
    return _solverInterface->satisfiesCurrentRelaxation(rcspSol);
}

bool BcRCSPFunctor
     ::drawPrimalSolutionToDotFile(IN_ BcFormulation spPtr,
                                   IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                   IN_ const std::string & filename)
{
    bcp_rcsp::FractionalMasterSolution rcspFracSolution;
    for (std::vector<std::pair<BcSolution, double> >::const_iterator pairIt = colsInMasterSolution.begin();
         pairIt != colsInMasterSolution.end(); ++pairIt)
    {
        if (pairIt->first.formulation().id().first() != _graphId)
            continue;
        rcspFracSolution.solPts.push_back(pairIt->first.solutionPtr()->rcspSolPtr());
        rcspFracSolution.values.push_back(pairIt->second);
    }

    _solverInterface->drawMasterSolutionInDotFile(rcspFracSolution, filename);

    return true;
}

bool BcRCSPFunctor::improveCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                               IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                               IN_ const bool & masterConverged)
{
    const ControlParameters & params = spPtr.probConfPtr()->param();

    if (params.RCSPdynamicNGmode() == DYNAMIC_NG_MODE_NONE)
        return false;

    if (_solverInterface->getEnumeratedStatus())
        return false;

    /// if params.RCSPdynamicNGmode() is negative then we increase NG only when there are no non-robust cuts
    if (params.RCSPdynamicNGmode() < 0)
    {
        std::vector<std::pair<const bcp_rcsp::RankOneCut *, double> > rankOneCuts;
        getActiveRankOneCuts(spPtr, rankOneCuts);
        if (!rankOneCuts.empty())
            return false;

#ifdef USE_NON_PUBLIC_CUTS
        std::vector<std::pair<const bcp_rcsp::StrongKPathCut *, double> > kPathCuts;
        getActiveStrongKPathCuts(spPtr, kPathCuts);
        if (!kPathCuts.empty())
            return false;
#endif
    }

    lightenCurrentSpRelaxation(spPtr, masterConverged, DSSR_MODE_AFTER_EACH_RELAXATION_IMPROVE);

    bcp_rcsp::FractionalMasterSolution fracSolution;
    for (std::vector<std::pair<BcSolution, double> >::const_iterator pairIt = colsInMasterSolution.begin();
         pairIt != colsInMasterSolution.end(); ++pairIt)
    {
        if (pairIt->first.formulation().id().first() != _graphId)
            continue;
        fracSolution.solPts.push_back(pairIt->first.solutionPtr()->rcspSolPtr());
        fracSolution.values.push_back(pairIt->second);
    }
    return _solverInterface->augmentNGneighbourhoods(fracSolution);
}

bool BcRCSPFunctor::lightenCurrentSpRelaxation(IN_ BcFormulation spPtr, const int & masterConverged,
                                               const int & callMode)
{
    const std::vector<InstanciatedVar *> & instVarPts = spPtr.probConfPtr()->instVarPts();
    const ControlParameters & params = spPtr.probConfPtr()->param();
    int paramDynamicNGmode = abs(params.RCSPdynamicNGmode());
    if (!_solverInterface->getEnumeratedStatus() && masterConverged && (params.RCSPuseDSSRInMode() & callMode))
    {
        bcp_rcsp::SolverInput input((int)instVarPts.size());
        if (!fillRCSPInput(spPtr, 0, instVarPts, input))
        {
            _messageIdToCutGeneration = PricingSolverCutsMessage::interruptSolution;
            return false;
        }

        return _solverInterface->reduceNGneighbourhoods(input);
    }
    return false;
}

void BcRCSPFunctor::setDiscreteCase(const bool & value)
{
    _solverParams.imposeDiscreteCase = value;
}

void BcRCSPFunctor::setVerificationFunctor(BcRCSPFunctor * functorPtr)
{
    _verificationFunctorPtr = functorPtr;
}

void BcRCSPFunctor::setLabelExtensionCostFunctor(RCSPLabelExtensionCostFunctor * functorPtr)
{
    _labelExtensionCostFunctorPtr = functorPtr;
    if (functorPtr != nullptr)
        _pureCostInstVarPtr = functorPtr->costInstVarPtr();
}

void BcRCSPFunctor::setPureCostBcVar(BcVar bcVar)
{
    if ((InstanciatedVar *)bcVar != nullptr)
    {
        _pureCostInstVarPtr = (InstanciatedVar *)bcVar;
    }
}


void BcRCSPFunctor::setFeasibilityCheckFunctor(RCSPFeasibilityCheckFunctor * functorPtr)
{
    _feasibilityCheckFunctorPtr = functorPtr;
}

void BcRCSPFunctor::setParamPhasesConfig(const std::vector<bcp_rcsp::ColGenPhaseConfig> & colGenPhaseConfigVector)
{
    _solverColGenPhaseConfig = colGenPhaseConfigVector;
}

void BcRCSPFunctor::setParamUseMoreTimers(const bool & value)
{
    _solverParams.useMoreTimers = value;
}

void BcRCSPFunctor::imposeSameResConsumptionInBucketCase()
{
    _solverParams.imposeSameResConsumptionInBucketCase = true;
}

void BcRCSPFunctor::setParamFormulationIsBuilt(const bool & value)
{
    _paramFormulationIsBuilt = value;
}

void BcRCSPFunctor::setParamPrintLevel(const int & value)
{
    _solverParams.printLevel = value;
}

void BcRCSPFunctor::saveToStandaloneRCSPfile(const std::string & fileName)
{
    _solverParams.standaloneRCSPfileName = fileName;
}

void BcRCSPFunctor::enumeratedCoreFileName(const std::string & fileName, BcVar outOfEnumCoreVar)
{
    std::cout << "BcRCSPFunctor WARNING : Enumerated core is not supported in this version " << std::endl;
}

/// ATTENTION! This function is only for the CVRP problem,
void BcRCSPFunctor::runWithDuals(BcFormulation spPtr, std::vector<double> & duals)
{
    NetworkFlow * networkPtr = ((NetworkFlow *)spPtr.network());
    if (networkPtr == NULL)
        return;

    const std::vector<InstanciatedVar *> & instVarPts = spPtr.probConfPtr()->instVarPts();
    int numVariables = instVarPts.size();
    bcp_rcsp::SolverInput input(numVariables);
    for (int varId = 0; varId < numVariables; ++varId)
    {
        InstanciatedVar * instVarPtr = instVarPts[varId];
        input.varRedCosts[varId] = instVarPtr->costrhs();
    }

    const std::unordered_map<InstanciatedVar *, int> & instVarToIdMap = spPtr.probConfPtr()->instVarToIdMap();
    for (lemon::ListDigraph::ArcIt lemonArc(networkPtr->digraph()); lemonArc != lemon::INVALID; ++lemonArc)
    {
        NetworkArc * netArcPtr = networkPtr->netArcPtr(lemonArc);
        const BcArcInfo * arcInfo = networkPtr->getArcInfoPtr(netArcPtr->id());

        if ((arcInfo->headVertId >= duals.size()) || (arcInfo->tailVertId >= duals.size()))
            continue;

        for (auto & pair : netArcPtr->varToCoeffMaps().front())
        {
            auto uomapIt = instVarToIdMap.find(pair.first);
            if (uomapIt != instVarToIdMap.end())
            {
                input.varRedCosts[uomapIt->second] -= pair.second * duals[arcInfo->headVertId] / 2.0;
                input.varRedCosts[uomapIt->second] -= pair.second * duals[arcInfo->tailVertId] / 2.0;
            }
        }
    }

    bcp_rcsp::SolverOutput output;

    _solverInterface->runPricing(input, output);
}

void BcRCSPFunctor::columnGenerationTerminated(IN_ BcFormulation spPtr, bool afterRedCostFixing,  int nodeOrder,
                                               int nodeDepth, int cutSepRound, double dualBound, double elapsedTime,
                                               bool masterConverged)
{
    bool softTimeLimitIsReached = false;
    if (!_solverInterface->columnGenerationTerminated(afterRedCostFixing, nodeOrder, nodeDepth, cutSepRound, dualBound,
                                                      elapsedTime, masterConverged, softTimeLimitIsReached))
    {
        std::cerr << "BaPCod RCSP functor error in the column generation termination callback " << std::endl;
        return;
    }

    if (!afterRedCostFixing)
    {
        if ((_messageIdToCutGeneration == PricingSolverCutsMessage::noMessage) && softTimeLimitIsReached)
            _messageIdToCutGeneration = PricingSolverCutsMessage::stopCutGeneration;
        if ((_messageIdToCutGeneration == PricingSolverCutsMessage::stopCutGeneration) && !softTimeLimitIsReached)
            _messageIdToCutGeneration = PricingSolverCutsMessage::noMessage;
    }
}

void BcRCSPFunctor::addToSolution(bcp_rcsp::Solution * rcspSolPtr, BcSolution & bcSol) const
{
    BcFormulation spForm = bcSol.formulation();
    NetworkFlow * networkPtr = ((NetworkFlow *)spForm.network());
    const std::unordered_map<int, InstanciatedVar *> & resIdToVarIdMap = spForm.probConfPtr()->resIdToVarIdMap();

    addToSolution(networkPtr, rcspSolPtr, resIdToVarIdMap, bcSol);
}

#endif /* BCP_RCSP_IS_FOUND */

/*****************************************************************************/
/*****************************************************************************/
/*****                        TODO: Cuts                                 *****/
/*****************************************************************************/
/*****************************************************************************/

BcLimMemRankOneCutConstrArray::BcLimMemRankOneCutConstrArray(const BcFormulation & formulation,
                                                             const double & rootPriorityLevel,
                                                             const double & nonRootPriorityLevel,
                                                             const bool & isFacultative,
                                                             const std::string & spVarName) :
        BcCutConstrArray(), _genCutCostrPtr(NULL)
{
    if (printL(5))
        std::cout << " BcLimMemRankOneCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcLimMemRankOneCutConstrArray =  R1C" << std::endl;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("R1C");
    const ControlParameters & params = formulation.probConfPtr()->param();

    if ((_genericCutConstrPtr == NULL) && (params.RCSPrankOneCutsMaxNumRows() > 0))
    {
        if (printL(5))
            std::cout << "BcLimMemRankOneCutConstrArray() : need to create cut" << std::endl;

#ifdef BCP_RCSP_IS_FOUND
        int memoryType = VERTEX_MEMORY_TYPE;
        if (params.RCSPrankOneCutsMemoryType() == RANK1_ARC_MEMORY_SEPARATION)
            memoryType = ARC_MEMORY_TYPE;
        else if (params.RCSPrankOneCutsMemoryType() == RANK1_FULL_MEMORY_SEPARATION)
            memoryType = FULL_MEMORY_TYPE;
        _genCutCostrPtr
                = new GenericLimMemRankOneCutConstr(formulation.probConfPtr()->modelPtr(), formulation.probConfPtr(),
                                                    "R1C", nonRootPriorityLevel, rootPriorityLevel, spVarName,
                                                    memoryType, isFacultative);
#else
        std::cerr << "BaPCod error : cannot use rank-1 cuts, as BCP_RCSP library is not found."
              << std::endl;
      exit(1);
#endif
        _genericCutConstrPtr = _genCutCostrPtr;
        _genericCutConstrPtr->defaultSense('L');
        _genericCutConstrPtr->defaultCostRhs(1.0);
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcLimMemRankOneCutConstrArray::~BcLimMemRankOneCutConstrArray()
{
}

void BcLimMemRankOneCutConstrArray::runStandAloneSeparation(const std::string & filePath)
{
#ifdef BCP_RCSP_IS_FOUND
    ((GenericLimMemRankOneCutConstr*)_genCutCostrPtr)->prepareSeparation();
    auto start = chrono::steady_clock::now();
    ((GenericLimMemRankOneCutConstr*)_genCutCostrPtr)->cutSeparationRoutine(filePath);
    auto end = chrono::steady_clock::now();
    std::cout << "Standalone R1C separation time : "
              << chrono::duration_cast<chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
#endif
}

BcCapacityCutConstrArray::BcCapacityCutConstrArray(const BcFormulation & formulation,
                                                   const int & maxCapacity,
                                                   const std::vector<int> & demands,
                                                   const bool & isFacultative,
                                                   const bool & equalityCase,
                                                   const int & twoPathCutsResId,
                                                   const double & rootPriorityLevel,
                                                   const double & nonRootPriorityLevel) :
        BcCutConstrArray()
{
    if (printL(5))
        std::cout << " BcCapacityCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcCapacityCutConstrArray = CAP" << std::endl;

    const ControlParameters & params = formulation.probConfPtr()->param();
    if (!params.RCSPuseCapacityCuts())
        return;

    if ((params.RCSPcapacityCutsSeparator() == 0) && !equalityCase)
    {
        if (printL(-1))
            std::cout << "BaPCod warning : RCC separator (CVRPSEP) is not activated "
                      << " as it does not support non-equality case" << std::endl;
        return;
    }

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("CAP");

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << "BcCapacityCutConstrArray() : need to create cut" << std::endl;

        if (params.RCSPcapacityCutsSeparator() == 0) /// CVRPSEP separator
        {
#if defined(CVRPSEP_IS_FOUND) && defined(USE_NON_PUBLIC_CUTS)
            _genericCutConstrPtr = new GenericCapacityCutConstr(formulation.probConfPtr()->modelPtr(),
                                                              formulation.probConfPtr(), "CAP",
                                                              nonRootPriorityLevel, rootPriorityLevel, isFacultative,
                                                              params.RCSPcapCutsMaxNumPerRound(),
                                                              maxCapacity, demands);
#else
            std::cerr << "BaPCod error : cannot use CVRPSEP separator of rounded capacity cuts, "
                      <<  "as CVRPSEP library is not found." << std::endl;
            exit(1);
#endif
        }
        else /// RCSP library separator
        {
#ifdef BCP_RCSP_IS_FOUND
            _genericCutConstrPtr = new GenericRCSPCapacityCutConstr(formulation.probConfPtr()->modelPtr(),
                                                                    formulation.probConfPtr(), "CAP",
                                                                    nonRootPriorityLevel, rootPriorityLevel,
                                                                    isFacultative, equalityCase, maxCapacity, demands,
                                                                    twoPathCutsResId);
#else
            std::cerr << "BaPCod error : cannot use RCSP separator of rounded capacity cuts, "
                    <<  "as BCP_RCSP library is not found." << std::endl;
          exit(1);
#endif
        }
        _genericCutConstrPtr->defaultSense('L');
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    else
    {
        _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("CAP2");
        if (_genericCutConstrPtr == NULL)
        {
            if (printL(5))
                std::cout << "BcCapacityCutConstrArray() : need to create cut" << std::endl;

            if (params.RCSPcapacityCutsSeparator() == 0) /// CVRPSEP separator
            {
#ifdef CVRPSEP_IS_FOUND
                _genericCutConstrPtr = new GenericCapacityCutConstr(formulation.probConfPtr()->modelPtr(),
                                                                formulation.probConfPtr(), "CAP2",
                                                                nonRootPriorityLevel, rootPriorityLevel, isFacultative,
                                                                params.RCSPcapCutsMaxNumPerRound(), maxCapacity,
                                                                demands);
#endif
            }
            else /// RCSP library separator
            {
#ifdef BCP_RCSP_IS_FOUND
                _genericCutConstrPtr = new GenericRCSPCapacityCutConstr(formulation.probConfPtr()->modelPtr(),
                                                                        formulation.probConfPtr(), "CAP2",
                                                                        nonRootPriorityLevel, rootPriorityLevel,
                                                                        isFacultative, equalityCase, maxCapacity, demands,
                                                                        twoPathCutsResId);
#endif
            }
            _genericCutConstrPtr->defaultSense('L');
            _genericCutConstrPtr->defaultFlag('d');
            _genericCutConstrPtr->defaultVal(0);
        }

    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcCapacityCutConstrArray::~BcCapacityCutConstrArray()
{
}


#ifdef BCP_RCSP_IS_FOUND

BcExtendedArcCut::BcExtendedArcCut(const NetworkFlow * networkPtr, ExtendedArcCut * cutPtr):
        _networkPtr(networkPtr), _cutPtr(cutPtr)
{
}

double BcExtendedArcCut::curDualVal() const
{
    long long int scaleFactor = _cutPtr->probConfPtr()->param().SafeDualBoundScaleFactor();
    if (scaleFactor > 0)
        return ceil(_cutPtr->valOrSepPointVal()._val * scaleFactor);
    return _cutPtr->valOrSepPointVal();
}

bool BcExtendedArcCut::inCurProb() const
{
    return _cutPtr->inCurProb();
}

bcp_rcsp::DiscreteCutInterface * BcExtendedArcCut::createCopy() const
{
    return new BcExtendedArcCut(_networkPtr, _cutPtr);
}

int BcExtendedArcCut::id() const
{
    return _cutPtr->id().first();
}

void BcExtendedArcCut::nicePrint() const
{
    _cutPtr->nicePrint();
}

bool BcExtendedArcCut::customCut() const
{
    return (_cutPtr->cutInfoPtr() != NULL);
}

double BcExtendedArcCut::getArcCoefficient(const int & tailVertId, const int & headVertId,
                                           const double * tailResCons) const
{
    return _cutPtr->getArcCoefficient(tailVertId, headVertId, tailResCons);
}

double BcExtendedArcCut::getArcCoefficient(const int & arcId, const double * resCons,
                                           const bool & isTailResCons) const
{
    return _cutPtr->getArcCoefficient(_networkPtr, arcId, resCons, isTailResCons);
}

double BcExtendedArcCut::getVertRouteCoefficient(const std::vector<int> & routeVertIds,
                                                 const std::vector<std::vector<double> > & routeResCons) const
{
    return _cutPtr->getRouteCoefficient(routeVertIds, routeResCons);
}

double BcExtendedArcCut::getArcRouteCoefficient(const std::vector<int> & routeArcIds,
                                                const std::vector<std::vector<double> > & routeResCons) const
{
    return _cutPtr->getRouteCoefficient(_networkPtr, routeArcIds, routeResCons);
}

void BcExtendedArcCut::reserveCut() const
{
    _cutPtr->incrParticipation(2);
}

void BcExtendedArcCut::releaseCut() const
{
    _cutPtr->decrParticipation(2);
}

#endif /* BCP_RCSP_IS_FOUND */

BcCustomExtendedArcCutArray::BcCustomExtendedArcCutArray(const BcFormulation & formulation,
                                                         const std::string & name,
                                                         const char & type,
                                                         const SelectionStrategy & priorityRule,
                                                         const double & rootPriorityLevel,
                                                         const double & nonRootPriorityLevel) :
        _genExtArcCutCostrPtr(NULL)
{
    if (printL(5))
        std::cout << " BcCustomExtendedArcCutArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcCustomExtendedArcCutArray = " << name << std::endl;

    GenericCutConstr * genCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr(name);


    if (genCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcCustomExtendedArcCutArrayFunctor() : need to create cut " << std::endl;

#ifdef BCP_RCSP_IS_FOUND
        _genExtArcCutCostrPtr = new GenericExtendedArcCutConstr(formulation.probConfPtr()->modelPtr(),
                                                                formulation.probConfPtr(), name, nonRootPriorityLevel,
                                                                rootPriorityLevel);
        _genExtArcCutCostrPtr->defaultFlag('d');
        _genExtArcCutCostrPtr->defaultVal(0);
#else
        std::cerr << "BaPCod error : cannot use extended arc cuts, as the BCP_RCSP library is not found."
                  << std::endl;
        exit(1);
#endif
    }
    else
    {
#ifdef BCP_RCSP_IS_FOUND
        _genExtArcCutCostrPtr = dynamic_cast<GenericExtendedArcCutConstr *>(genCutConstrPtr);
#endif
    }
}


BcCustomExtendedArcCutArray::~BcCustomExtendedArcCutArray()
{
}

const BcCustomExtendedArcCutArray & BcCustomExtendedArcCutArray
      ::attach(BcCustomExtendedArcCutSeparationFunctor * separationFunctorPtr)
{
#ifdef BCP_RCSP_IS_FOUND
    _genExtArcCutCostrPtr->setSeparationFunctor(separationFunctorPtr);
    separationFunctorPtr->_genExtArcCutCostrPtr = _genExtArcCutCostrPtr;
#endif
    return *this;
}

BcCustomExtendedArcCutSeparationFunctor::BcCustomExtendedArcCutSeparationFunctor() :
        _genExtArcCutCostrPtr(NULL)
{
}

BcCustomExtendedArcCutSeparationFunctor::~BcCustomExtendedArcCutSeparationFunctor()
{
}

BcConstr BcCustomExtendedArcCutSeparationFunctor::createNewCut(BcCustomExtendedArcCutInfo * cutInfoPtr,
                                                               const char & sense, const double & rhs)
{
#ifdef BCP_RCSP_IS_FOUND
    ExtendedArcCut * newCutPtr = new ExtendedArcCut(_genExtArcCutCostrPtr, _genExtArcCutCostrPtr->probConfPtr(),
                                                    _genExtArcCutCostrPtr->defaultName(), rhs, sense, cutInfoPtr);
    return BcConstr(newCutPtr);
#else
    return BcConstr(nullptr);
#endif
}

int BcCustomExtendedArcCutSeparationFunctor
    ::cutSeparationRoutine(BcFormulation formPtr, BcSolution & projectedSol,
                           std::list<std::pair<double, BcSolution> > & columnsInSol,
                           const double & violationTolerance, std::list<BcConstr> & cutList)
{
    return 0;
}

#ifdef BCP_RCSP_IS_FOUND


void getActiveRankOneCuts(const BcFormulation & spPtr,
                          std::vector<std::pair<const bcp_rcsp::RankOneCut *, double> > & rankOneCuts)
{
    if (spPtr.probConfPtr() == NULL)
    {
        std::cerr << "ERROR Model BcFormulation == NULL in getLimMemRankOneActiveMasterCutsList" << std::endl;
        exit(1);
    }
    rankOneCuts.clear();

    MasterConf * mastConfPtr = spPtr.probConfPtr()->modelPtr()->master();
    GenericCutConstr * genR1CutConstrPtr = mastConfPtr->getGenericCutConstr("R1C");
    if (genR1CutConstrPtr == NULL)
        return;

    long long int scaleFactor = spPtr.probConfPtr()->param().SafeDualBoundScaleFactor();

    const IndexCell2InstancConstrPtrMap & constrPtrMap = genR1CutConstrPtr->indexCell2InstancConstrPtrMap();
    for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
        if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
            && it->second->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
        {
            LimMemRankOneCut * cutPtr = static_cast<LimMemRankOneCut *>(it->second);
            double curDualValue = (scaleFactor > 0) ? (double)ceil(cutPtr->valOrSepPointVal()._val * scaleFactor)
                                                    : (double)cutPtr->valOrSepPointVal();
            rankOneCuts.push_back(std::make_pair(cutPtr->rcspCutPtr(), curDualValue));
        }
}

void getActiveDiscreteCuts(const BcFormulation & spPtr, std::vector<const BcExtendedArcCut *> & cutPts)
{
    if ((spPtr.probConfPtr() == NULL) || (spPtr.probConfPtr()->networkFlowPtr() == NULL))
    {
        std::cerr << "ERROR Model BcFormulation == NULL in getExtendedArcMasterCutsList" << std::endl;
        exit(1);
    }

    MasterConf * mastConfPtr = spPtr.probConfPtr()->modelPtr()->master();

    std::set<GenericCutConstr *, DynamicGenConstrSort> _candidateCutGenericConstr;

    std::set<GenericCutConstr *, DynamicGenConstrSort>::iterator genCutPtrIt;
    for (genCutPtrIt = mastConfPtr->candidateCutGenericConstr().begin();
         genCutPtrIt != mastConfPtr->candidateCutGenericConstr().end(); ++genCutPtrIt)
    {
        GenericExtendedArcCutConstr * genExtArcCutPtr = dynamic_cast<GenericExtendedArcCutConstr *>(*genCutPtrIt);
        if (genExtArcCutPtr != NULL)
        {
            const IndexCell2InstancConstrPtrMap & constrPtrMap = genExtArcCutPtr->indexCell2InstancConstrPtrMap();
            for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin();
                 it != constrPtrMap.end(); ++it)
                if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
                    && it->second->isTypeOf(VcId::ExtendedArcCutConstrMask))
                {
                    ExtendedArcCut * cutPtr = static_cast<ExtendedArcCut *>(it->second);
                    cutPts.push_back(new BcExtendedArcCut(spPtr.probConfPtr()->networkFlowPtr(), cutPtr));
                }
        }
    }
}

#endif /* BCP_RCSP_IS_FOUND */

BcResConsumptionKnapsackCutConstrArray
::BcResConsumptionKnapsackCutConstrArray(const BcFormulation & formulation,
                                         const double & rootPriorityLevel,
                                         const double & nonRootPriorityLevel) :
        BcCutConstrArray(), _genCutCostrPtr(NULL)
{
    if (printL(5))
        std::cout << " BcResConsumptionKnapsackCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcResConsumptionKnapsackCutConstrArray =  RCK" << std::endl;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("RCK");
    const ControlParameters & params = formulation.probConfPtr()->param();

    if ((_genericCutConstrPtr == NULL) && (params.RCSPresConsKnapsackCutsMode() != -1))
    {
        if (printL(5))
            std::cout << "BcResConsumptionKnapsackCutConstrArray() : need to create cut" << std::endl;

#ifdef BCP_RCSP_IS_FOUND
        _genCutCostrPtr = new GenericResConsKnapsackCutConstr(formulation.probConfPtr()->modelPtr(),
                                                              formulation.probConfPtr(), "RCK", nonRootPriorityLevel,
                                                              rootPriorityLevel);
#else
        std::cerr << "BaPCod error : cannot use resource consumption knapsack cuts, as BCP_RCSP library is not found."
              << std::endl;
    exit(1);
#endif
        _genericCutConstrPtr = _genCutCostrPtr;
        _genericCutConstrPtr->defaultSense('L');
        _genericCutConstrPtr->defaultCostRhs(0.0);
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcResConsumptionKnapsackCutConstrArray::~BcResConsumptionKnapsackCutConstrArray()
{
}

#ifdef BCP_RCSP_IS_FOUND
void getActiveRouteLoadKnapsackCuts(const BcFormulation & spPtr,
                                    std::vector<std::pair<const bcp_rcsp::RouteLoadKnapsackCut *, double> > & rlkCuts)
{
    if (!spPtr.isColGenSp())
    {
        std::cerr << "BaPCod error : formulation is not a col.gen.sp. problem in "
                  << "getActiveRouteLoadKnapsackCuts()" << std::endl;
        exit(1);
    }
    rlkCuts.clear();

    MasterConf * mastConfPtr = spPtr.probConfPtr()->modelPtr()->master();
    GenericCutConstr * genConstrPtr = mastConfPtr->getGenericCutConstr("RCK");
    if (genConstrPtr == NULL)
        return;

    auto cgSpConfPtr = static_cast<ColGenSpConf *>(spPtr.probConfPtr());

    long long int scaleFactor = spPtr.probConfPtr()->param().SafeDualBoundScaleFactor();

    const IndexCell2InstancConstrPtrMap & constrPtrMap = genConstrPtr->indexCell2InstancConstrPtrMap();
    for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
        if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
            && it->second->isTypeOf(VcId::ResConsKnapsackCutConstrMask))
        {
            ResConsKnapsackCut * rcKnapCutPtr = static_cast<ResConsKnapsackCut *>(it->second);
            double curDualValue = (scaleFactor > 0) ? (double)ceil(rcKnapCutPtr->valOrSepPointVal()._val * scaleFactor)
                                                    : (double)rcKnapCutPtr->valOrSepPointVal();
            if (rcKnapCutPtr->isRelatedTo(cgSpConfPtr))
                rlkCuts.push_back(std::make_pair(rcKnapCutPtr->rcspCutPtr(), curDualValue));
        }
}
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****                   TODO: Branching constraints                     *****/
/*****************************************************************************/
/*****************************************************************************/

BcPathsPerNetworkBranching::BcPathsPerNetworkBranching(const BcFormulation & formulation,
                                                       const double & priorityLevel,
                                                       const bool & toBeUsedInPreprocessing) :
        _genPathsPerNetworkBranchingConstrPtr(NULL)
{
    std::string name = "PPN";

    if (printL(5))
        std::cout << " BcPathsPerNetworkBranching() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcPathsPerNetworkBranching =  " << name << std::endl;

    GenericBranchingConstr * genBranchConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(name);
    if (genBranchConstrPtr != NULL)
        _genPathsPerNetworkBranchingConstrPtr = dynamic_cast<GenPathsPerNetworkBranchingConstr *>(genBranchConstrPtr);

    if (_genPathsPerNetworkBranchingConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcPathsPerNetworkBranching() : need to create branching  " << std::endl;

        Model * modelPtr = formulation.probConfPtr()->modelPtr();

        _genPathsPerNetworkBranchingConstrPtr =
                new GenPathsPerNetworkBranchingConstr(modelPtr, formulation.probConfPtr(), name,
                                                      SelectionStrategy::MostFractional, priorityLevel, priorityLevel,
                                                      toBeUsedInPreprocessing);
        _genPathsPerNetworkBranchingConstrPtr->defaultFlag('d');
    }
}

BcPathsPerNetworkBranching::~BcPathsPerNetworkBranching()
{
}

BcPackSetAssignBranching::BcPackSetAssignBranching(const BcFormulation & formulation,
                                                   const double & priorityLevel,
                                                   const bool & toBeUsedInPreprocessing) :
        _genPackSetAssignBranchingConstrPtr(NULL)
{
    std::string name = "ESA";

    if (printL(5))
        std::cout << " BcPackSetAssignBranching() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcPackSetAssignBranching =  " << name << std::endl;

    GenericBranchingConstr * genBranchConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(name);
    if (genBranchConstrPtr != NULL)
        _genPackSetAssignBranchingConstrPtr = dynamic_cast<GenPackSetAssignBranchingConstr *>(genBranchConstrPtr);

    if (_genPackSetAssignBranchingConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcPackSetAssignBranching() : need to create branching  " << std::endl;

        Model * modelPtr = formulation.probConfPtr()->modelPtr();

        _genPackSetAssignBranchingConstrPtr =
                new GenPackSetAssignBranchingConstr(modelPtr, formulation.probConfPtr(), name,
                                                    SelectionStrategy::MostFractional, priorityLevel, priorityLevel,
                                                    toBeUsedInPreprocessing);
        _genPackSetAssignBranchingConstrPtr->defaultFlag('d');
    }
}

BcPackSetAssignBranching::~BcPackSetAssignBranching()
{
}

BcPackSetResConsumptionBranching::BcPackSetResConsumptionBranching(const BcFormulation & formulation,
                                                                   const double & priorityLevel)
#ifdef BCP_RCSP_IS_FOUND
        : _PackSetResConsGenBranchConstrPtr(NULL)
#endif
{
    std::string name = "ESRC";

    if (printL(5))
        std::cout << " BcPackSetResConsumptionBranching() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcPackSetResConsumptionBranching =  " << name << std::endl;

    GenericBranchingConstr * genBranchConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(name);
#ifdef BCP_RCSP_IS_FOUND
    if (genBranchConstrPtr != NULL)
        _PackSetResConsGenBranchConstrPtr = dynamic_cast<PackSetResConsGenBranchConstr *>(genBranchConstrPtr);

    if (_PackSetResConsGenBranchConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcPackSetResConsumptionBranching() : need to create branching  " << std::endl;

        Model * modelPtr = formulation.probConfPtr()->modelPtr();

        _PackSetResConsGenBranchConstrPtr =
                new PackSetResConsGenBranchConstr(modelPtr, formulation.probConfPtr(), name,
                                                  SelectionStrategy::MostFractional, priorityLevel);
        _PackSetResConsGenBranchConstrPtr->defaultFlag('d');
    }
#else
    std::cerr << "BaPCod error : cannot use elem. set res. consumption branching, as the BCP_RCSP library is not found."
              << std::endl;
    exit(1);
#endif
}

BcPackSetResConsumptionBranching::~BcPackSetResConsumptionBranching()
{
}

#ifdef BCP_RCSP_IS_FOUND

void getPackSetResConsActiveBranchConstrList(const BcFormulation & spPtr,
                                             std::vector<const bcp_rcsp::AccumResConsBranchConstraint *> & constrPts)
{
    if (spPtr.probConfPtr() == NULL)
    {
        std::cerr << "ERROR Model BcFormulation == NULL in getPackSetResConsActiveBranchConstrList" << std::endl;
        exit(1);
    }
    constrPts.clear();

    ConstrIndexManager & probConstrSet = spPtr.probConfPtr()->modelPtr()->master()->probPtr()->probConstrSet();
    for (ConstrIndexManager::iterator constrPtrIt = probConstrSet.begin(VcIndexStatus::Active, 'd');
         constrPtrIt != probConstrSet.end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    {
        if ((*constrPtrIt)->isTypeOf(VcId::PackSetResConsInstMastBranchConstrMask))
        {
            PackSetResConsInstMastBranchConstr * bcConstrPtr
                    = static_cast<PackSetResConsInstMastBranchConstr *>(*constrPtrIt);
            if (bcConstrPtr != NULL)
                constrPts.push_back(bcConstrPtr->getRcspConstrPtr());
        }
    }
}

#endif

BcPackSetRyanFosterBranching::BcPackSetRyanFosterBranching(const BcFormulation & formulation,
                                                           const double & priorityLevel,
                                                           const bool & usePackingSets)
#ifdef BCP_RCSP_IS_FOUND
        : _packSetRyanFosterGenBranchConstrPtr(NULL)
#endif
{
    std::string name = "PSRF";

    if (printL(5))
        std::cout << " BcPackSetRyanFosterBranching() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcPackSetRyanFosterBranching =  " << name << std::endl;

#ifdef BCP_RCSP_IS_FOUND
    GenericBranchingConstr * genBranchConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(name);
    if (genBranchConstrPtr != NULL)
        _packSetRyanFosterGenBranchConstrPtr = dynamic_cast<PackSetRyanFosterGenBranchConstr *>(genBranchConstrPtr);

    if (_packSetRyanFosterGenBranchConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcPackSetRyanFosterBranching() : need to create branching  " << std::endl;

        Model * modelPtr = formulation.probConfPtr()->modelPtr();

        _packSetRyanFosterGenBranchConstrPtr =
                new PackSetRyanFosterGenBranchConstr(modelPtr, formulation.probConfPtr(), name,
                                                     SelectionStrategy::MostFractional, priorityLevel,
                                                     usePackingSets);
        _packSetRyanFosterGenBranchConstrPtr->defaultFlag('d');
    }
#endif
}

BcPackSetRyanFosterBranching::~BcPackSetRyanFosterBranching()
{
}

#ifdef BCP_RCSP_IS_FOUND

void getPackSetRyanFosterActiveBranchConstrList(const BcFormulation & spPtr,
                                                std::vector<const bcp_rcsp::RyanFosterBranchConstraint *> & constrPts)
{
    if (spPtr.probConfPtr() == NULL)
    {
        std::cerr << "ERROR Model BcFormulation == NULL in getPackSetRyanFosterActiveBranchConstrList" << std::endl;
        exit(1);
    }
    constrPts.clear();

    ConstrIndexManager & probConstrSet = spPtr.probConfPtr()->modelPtr()->master()->probPtr()->probConstrSet();
    for (ConstrIndexManager::iterator constrPtrIt = probConstrSet.begin(VcIndexStatus::Active, 'd');
         constrPtrIt != probConstrSet.end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    {
        if ((*constrPtrIt)->isTypeOf(VcId::PackSetRyanFostInstMastBranchConstrMask))
        {
            PackSetRyanFosterInstMastBranchConstr * bpConstrPtr
                    = static_cast<PackSetRyanFosterInstMastBranchConstr *>(*constrPtrIt);
            constrPts.push_back(bpConstrPtr->getRcspConstrPtr());
        }
    }
}

#endif
