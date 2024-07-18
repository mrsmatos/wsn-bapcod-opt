/**
 *
 * This file bcAlg4EvalByColAndCutGen.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4EvalByColAndCutGen.hpp"
#include "bcStabilizationColgen.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcNetworkBasedCuts.hpp"
#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#ifdef USE_NON_PUBLIC_CUTS
#include "bcNonPublicCuts.hpp"
#endif
#endif

Alg4EvalByColAndCutGen::Alg4EvalByColAndCutGen(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons) :
    Alg4EvalBySimplexBasedColGen(probPtr, masterCommons),
    _saveInfoForAutoRankOneCutsMemory(false)
{
 if (printL(5))
    std::cout << " Alg4EvalByColAndCutGen:: NEW ALG" << std::endl;
}

void Alg4EvalByColAndCutGen::reintroduceArtVarInMast()
{
  if (!_need2reintroduceArtVarInMast)
    return;

  if (printL(5))
    std::cout << " MasterConf::reintroduceArtVarInMast(): ArtificialVar cost update" << std::endl;

  ///  Needed to avoid deleting in the set we are currently scanning
  std::list<Variable *> var2Activate;
  for (std::list<Variable *>::const_iterator it =
       _nonStabArtVarPtrList.begin();
       it != _nonStabArtVarPtrList.end(); ++it)
  {
    if (!_probPtr->probVarSet().count(*it, VcIndexStatus::Active))
      var2Activate.push_back(*it);

    if (printL(5))
      std::cout << " MasterConf::reintroduceArtVarInMast() var " << (*it)->name() << " at cost = " << (*it)->costrhs()
                << std::endl;
  }
  if (!var2Activate.empty())
    _probPtr->addVarSet(var2Activate, 1, 2);

  _need2reintroduceArtVarInMast = false;

  return;
}

bool Alg4EvalByColAndCutGen::setupAlgo(Node * nodePtr)
{
  return Alg4EvalBySimplexBasedColGen::setupAlgo(nodePtr);
}

bool Alg4EvalByColAndCutGen::eval()
{
  Alg4MasterSolAndBounds backupBounds(_masterCommons.objStatus());
  setCutSepRound((_maxNbOfCutRounds > 0) ? 0 : -1);

  std::set<GenericCutConstr *, DynamicGenConstrSort>::const_iterator cutGenPtrIt;
  std::set<double> priorityLevelsSet;
  for (cutGenPtrIt = _masterCommons.candidateCutGenericConstr().begin();
       cutGenPtrIt != _masterCommons.candidateCutGenericConstr().end(); ++cutGenPtrIt)
    priorityLevelsSet.insert((*cutGenPtrIt)->dynamicPriorityLevel(_currentNodePtr->isRoot()));
  if (param().ColGenSpRelaxationImprovementPriority() > 0.0)
      priorityLevelsSet.insert(param().ColGenSpRelaxationImprovementPriority());

  double maxCutSeparatorPriority = (priorityLevelsSet.empty() ? 0.0 : *(--priorityLevelsSet.end()));
  _lastPriorityLevel = maxCutSeparatorPriority;
//  if (_currentNodePtr->isRoot()) {
//    /// we assign the largest priority to the sp. relaxation improvement in the root
//    /// TO DO: allow user to define the priority of the sp. relaxation improvement
//    _lastPriorityLevel += 1.0;
//    priorityLevelsSet.insert(_lastPriorityLevel);
//  }

  for (cutGenPtrIt = _masterCommons.candidateCutGenericConstr().begin();
       cutGenPtrIt != _masterCommons.candidateCutGenericConstr().end(); ++cutGenPtrIt)
    (*cutGenPtrIt)->resetSeparationPhase();

  /// Combined phase I and phase II solution
  if (Alg4EvalBySimplexBasedColGen::eval())
    return true; /// Problem infeasible

  /// we do not cleanup cuts if masterConverged is false, as in this case the dual solution value is not
  /// guaranteed to be close to the primal one
  if (masterConverged())
    cleanupRestrictedMastCuts();

  /// we are saving the rollback point
  if (_saveInfoForAutoRankOneCutsMemory)
    _currentNodePtr->saveAutoRankOneCutsMemoryInfo(false, false);
  backupBoundsExceptIncIpPrimalBound(backupBounds);
  _currentNodePtr->saveProblemAndEvalAlgInfo();
  setRollbackPointSavedStatus(true);
  runColGenSpRelaxationLightening(1); /// first bit (001), called after the first convergence of a node

  if (_pricingSolverCutsMessageId == PricingSolverCutsMessage::stopCutGeneration)
    {
      if (_minNbOfCutRounds == 0)
        {
          if (printL(-1))
            std::cout << "----- Cut generation stopped by the pricing solver" << std::endl;
          return false;
        }
      else
        _pricingSolverCutsMessageId = PricingSolverCutsMessage::noMessage;
    }

  bool addCuts = (_currentNodePtr->depth() <= param().MasterCuttingPlanesDepthLimit());

  if (printL(5))
    std::cout << "ColAndCutGenSolver::solve(): enters cutting plane generation" << std::endl;

  Bound prevNodeLpValue = algCurLpPrimalBound();
  double prevPriorityLevel = BapcodInfinity;
  Double gapImprovement(100);
  int cutTailingOffCounter = 0;
  int nbOfCutRounds = 0;

  while (progStatus().doRun() && !subProbSolutionsEnumeratedToMIP() && (nbOfCutRounds < _maxNbOfCutRounds)
         & !isConquered())
    {
      printAndRecordActiveCutsStats(false);

      bool dualSolutionStored = masterConverged();
      ConstrPtr2DoubleMap storedDualSolution;
      if (dualSolutionStored)
        _probPtr->storeDualSolution(storedDualSolution);

      //      bool relaxationWasImproved = runColGenSpRelaxationImprovement(subprobsWithImprovedRelaxation);
//      if ((_lastPriorityLevel > maxCutSeparatorPriority) && !relaxationWasImproved)
//        runColGenSpRelaxationLightening(2); /// second bit (010), called after the convergence of SP relaxation
//      if (addCuts && ((_lastPriorityLevel <= maxCutSeparatorPriority) || !relaxationWasImproved))
//        relaxationWasImproved = addCutToMaster('F') || relaxationWasImproved;
      if (!addCuts || !addCutToMaster('F', masterConverged()))
          break;

//      if (!relaxationWasImproved)
//        break;

      nbOfCutRounds += 1;
      setCutSepRound(nbOfCutRounds);
      if (printL(5))
        std::cout << "ColAndCutGenSolver::solve()  after cutting plane procedure" << std::endl;

      if (_colGenStabilizationPtr != nullptr)
        {
          StabilizationInfo * stabInfoPtr = _colGenStabilizationPtr->recordStabilizationInfo();
          _colGenStabilizationPtr->setDownStab();
          _colGenStabilizationPtr->setupStab(stabInfoPtr, algIncLpDualBound(),
                                             param().MaxNbOfStagesInColGenProcedure() - 1,
                                             _currentNodePtr->depth() + 1);
          delete stabInfoPtr;
        }

      /// commented by Ruslan : for the moment we do not desactivate artificial variables if pure Phase I succeeded,
      /// as local artificial variables and associated constrains should be synchronized
      /// TO DO : think how to overcome this
      // reintroduceArtVarInMast();

      if (printL(-1))
      {
        _probPtr->printDynamicVarConstrStats(std::cout);
        std::cout << std::endl;
      }

      ///  Combined phase I and phase II solution
      if (Alg4EvalBySimplexBasedColGen::eval())
        return true; /// Problem infeasible

      if ((_pricingSolverCutsMessageId == PricingSolverCutsMessage::doCutsRollback) && dualSolutionStored)
        {
          /// we first restore the dual solutions
          _probPtr->resetDualSolution(storedDualSolution);
          if (runColGenSpRelaxationLightening(4)) /// third bit (100), called on rollbak
            {
              if (printL(-1))
                std::cout << "----- Column generation is repeated with ""lighter"" subproblems" << std::endl;
              if (Alg4EvalBySimplexBasedColGen::eval())
                return true; /// Problem infeasible
            }
        }

      if (_pricingSolverCutsMessageId == PricingSolverCutsMessage::doCutsRollback)
        {
          if (printL(-1))
            std::cout << "----- Column generation is interrupted by the pricing solver, "
                      << "we do rollback to the point before the latest cut generation round" << std::endl;
          restoreBoundsExceptIncIpPrimalBound(backupBounds);
          _currentNodePtr->makeProbSetupInfoObligatoryForFullSetup();
          /// current solution can be integer but there will be duality gap after rollback
          /// so we need to change _solIsInteger
          _solIsInteger = false;

          break;
        }

      if ((algIncIpPrimalBound() - prevNodeLpValue) > std::abs((double)prevNodeLpValue) * 0.1)
      {
        /// we use false gap of 10% if the current gap is larger
        gapImprovement = (algCurLpPrimalBound() - prevNodeLpValue) / (std::abs((double)prevNodeLpValue) * 0.1);
        if (printL(0))
          std::cout << "False gap improvement since the last cut separation : " << gapImprovement
                    << " (" << prevNodeLpValue << ")" << std::endl;
      }
      else
      {
        gapImprovement = (algCurLpPrimalBound() - prevNodeLpValue) / (algIncIpPrimalBound() - prevNodeLpValue);
        if (printL(-1))
          std::cout << "Gap improvement since the last cut separation : " << gapImprovement
                    << " (" << prevNodeLpValue << ")" << std::endl;
      }
      prevNodeLpValue = algCurLpPrimalBound();
      if (_lastPriorityLevel != prevPriorityLevel)
        {
          prevPriorityLevel = _lastPriorityLevel;
          cutTailingOffCounter = 0;
        }

      /// we do not cleanup cuts if masterConverged is false, as in this case the dual solution value is not
      /// guaranteed to be close to the primal one
      if (masterConverged())
        cleanupRestrictedMastCuts();

      /// we are saving the rollback point (rollbackPointSavedStatus is already set to true)
      backupBoundsExceptIncIpPrimalBound(backupBounds);
      _currentNodePtr->saveProblemAndEvalAlgInfo();

      if (_saveInfoForAutoRankOneCutsMemory)
        _currentNodePtr->saveAutoRankOneCutsMemoryInfo(false, _probPtr->rankOneCutsArePresent());

//      std::cout  << "$$$$# "  << _pricingSolverCutsMessageId << " " << nbOfCutRounds << " " << _minNbOfCutRounds
//                 << std::endl;
      if (_pricingSolverCutsMessageId == PricingSolverCutsMessage::stopCutGeneration)
        {
          if (nbOfCutRounds >= _minNbOfCutRounds)
            {
              if (printL(-1))
                std::cout << "----- Cut generation stopped by the pricing solver" << std::endl;
              break;
            }
          else
            _pricingSolverCutsMessageId = PricingSolverCutsMessage::noMessage;
        }
      if (bapcodInit().startTime().getElapsedTime() > param().GlobalTimeLimitInTick())
        {
          if (printL(-1))
            std::cout << "----- Cut generation stopped by global time limit" << std::endl;
          break;
        }
      if (gapImprovement < param().CutTailingOffThreshold())
        {
          bool allCutSeparatorsAreAtMaximumPhase = true;
          for (cutGenPtrIt = _masterCommons.candidateCutGenericConstr().begin();
               cutGenPtrIt != _masterCommons.candidateCutGenericConstr().end(); ++cutGenPtrIt)
            if (((*cutGenPtrIt)->dynamicPriorityLevel(_currentNodePtr->isRoot()) >= _lastPriorityLevel)
                && !(*cutGenPtrIt)->separationPhaseIsMaximum())
              {
                if (printL(0))
                  std::cout << "Separator " << (*cutGenPtrIt)->defaultName() << " not in maximum phase" << std::endl;
                allCutSeparatorsAreAtMaximumPhase = false;
              }

          /// we update tailing off counter only if the separation phase is maximum for all cut separators
          if (allCutSeparatorsAreAtMaximumPhase)
            cutTailingOffCounter += 1;
          else
            for (cutGenPtrIt = _masterCommons.candidateCutGenericConstr().begin();
                 cutGenPtrIt != _masterCommons.candidateCutGenericConstr().end(); ++cutGenPtrIt)
              (*cutGenPtrIt)->increaseSeparationPhase();
          if ((cutTailingOffCounter >= param().CutTailingOffCounterThreshold()) && (nbOfCutRounds >= _minNbOfCutRounds))
            {
              std::set<double>::iterator dblIt = priorityLevelsSet.find(_lastPriorityLevel);
              if (dblIt != priorityLevelsSet.begin())
                {
                  /// tailing off occurs for this priority level, so we decrease it to pass to more expensive separator
                  --dblIt;
                  _lastPriorityLevel = *dblIt;
                  if (printL(0))
                    std::cout << "----- Cut separators priority level decreased to " << *dblIt << "-----" << std::endl;
                  if (_lastPriorityLevel == maxCutSeparatorPriority)
                    runColGenSpRelaxationLightening(2); /// second bit (010), called after the convergence of SP relax.
                }
              else /// current priority level is minumum possible
                {
                  if (_saveInfoForAutoRankOneCutsMemory)
                    _currentNodePtr->saveAutoRankOneCutsMemoryInfo(true);
                  if (printL(-1))
                    std::cout << "----- Cut generation is stopped due to tailing off -----" << std::endl;
                  break;
                }
            }
          else if (printL(0))
            {
                std::cout << "Cut generation tailing off counter increased to " << cutTailingOffCounter << std::endl;
            }
        }
      if ((_maxNbOfCutRounds > 0) && (nbOfCutRounds >= _maxNbOfCutRounds))
        {
          if (printL(-1))
            std::cout << "----- Maximum number of cut rounds is reached " << std::endl;
          break;
        }
      if ((param().CutTailingOffAbsTolerance() > 0)
          && (algIncIpPrimalBound() - algIncLpDualBound() < param().CutTailingOffAbsTolerance()))
        {
          if (printL(-1))
            std::cout << "----- Cut generation is stopped as the absolute tolerance threshold is reached ------  "
                      << std::endl;
          break;
        }
    }

  if (_maxNbOfCutRounds > 0)
    printAndRecordActiveCutsStats(_currentNodePtr->isRoot());

  if (printL(0) && (_maxLevelOfSubProbRestriction == 0) && (param().printMasterPrimalSols() == 3))
    _probPtr->printDetailedPrimalSol();

  if ((_maxLevelOfSubProbRestriction == 0) && (param().fractionalSolutionDotFile() != ""))
      drawPrimalSolutionToDotFile(param().fractionalSolutionDotFile());

  if (printL(0) && (_maxLevelOfSubProbRestriction == 0) && (param().printMasterPrimalSols() == 4)
      && (_currentNodePtr->treatOrder() == 1))
    _probPtr->printDualSol();

  return false;
}

void Alg4EvalByColAndCutGen::printAndRecordActiveCutsStats(bool recordStats)
{
    /// for each type of cuts, the tuple <nb. cuts with zero rhs, nb. cuts with non-zero rhs,
    ///                                   aver. dual value of cuts with non-zero rhs.>
    std::map<std::string, std::tuple<int, int, double> > numActiveCutsPerType;
#ifdef BCP_RCSP_IS_FOUND
    int numDepotCapCutsWithOneY = 0;
    int numDepotCapCutsWithTwoYs = 0;
    std::vector<int> numPackingCuts(RANK1_CUTS_MAX_NUMBER_OF_ROWS + 1, 0);
    std::vector<int> numCoveringCuts(RANK1_CUTS_MAX_NUMBER_OF_ROWS + 1, 0);
    int maxOneKvalue = std::abs(param().RCSPresConsKnapsackCutsMode());
    std::vector<int> numRLKCutsByOneK( maxOneKvalue + 1, 0);
#endif

    for (ConstrIndexManager::iterator constrPtrIt = _probPtr->probConstrSet().begin(VcIndexStatus::Active, 'd');
         constrPtrIt != _probPtr->probConstrSet().end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    {
        if (!(*constrPtrIt)->val().isZero()
            && (*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask)
            && !(*constrPtrIt)->isTypeOf(VcId::BranchingConstrBaseTypeMask))
        {
            double absDualValue = std::abs((double)(*constrPtrIt)->val());
            if (absDualValue < Double::precision)
                continue;
            InstanciatedConstr * instConstrPtr = static_cast<InstanciatedConstr *>(*constrPtrIt);
            std::string cutClass = instConstrPtr->genericName();
            auto mapItPair = numActiveCutsPerType.insert(std::make_pair(cutClass,std::make_tuple(0, 0, 0.0)) );
            if (instConstrPtr->curRhs() == 0)
            {
                std::get<0>(mapItPair.first->second) += 1;  /// number of zero rhs occurrences
            }
            else
            {
                std::get<1>(mapItPair.first->second) += 1;  /// number of non-zero occurrences
                std::get<2>(mapItPair.first->second) += absDualValue * instConstrPtr->curRhs(); /// total dual value
            }
#ifdef BCP_RCSP_IS_FOUND
            if (cutClass == "DCC")
            {
                int numYmembers = 0;
                for (const auto & pair : instConstrPtr->member2coefMap())
                    if (pair.first->genericName() == "Y")
                        numYmembers += 1;
                if (numYmembers == 1)
                    numDepotCapCutsWithOneY += 1;
                else
                    numDepotCapCutsWithTwoYs += 1;
            }
            if ((cutClass == "R1C") && instConstrPtr->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
            {
                    const bcp_rcsp::RankOneCut * rank1CutPtr
                        = static_cast<const LimMemRankOneCut *>(instConstrPtr)->rcspCutPtr();
                    if (rank1CutPtr->cutClass == SET_PACKING_CUT)
                        numPackingCuts[rank1CutPtr->numRows] += 1;
                    else
                        numCoveringCuts[rank1CutPtr->numRows] += 1;
            }
            if ((cutClass == "RCK" && instConstrPtr->isTypeOf(VcId::ResConsKnapsackCutConstrMask)))
            {
                const bcp_rcsp::RouteLoadKnapsackCut * rlkCutPtr
                        = static_cast<const ResConsKnapsackCut *>(instConstrPtr)->rcspCutPtr();
                numRLKCutsByOneK[rlkCutPtr->oneKvalue] += 1;
            }
#endif
        }
    }

    if (numActiveCutsPerType.empty())
        return;

    if (printL(0))
        std::cout << "Current active cuts :";

    for (const auto & pair : numActiveCutsPerType)
    {
        if (recordStats)
        {
            bapcodInit().statistics().incrCounter(std::string("bcCountRootActiveCut") + pair.first,
                                                  std::get<0>(pair.second) + std::get<1>(pair.second));
            bapcodInit().statistics().incrValue(std::string("bcRecRootContribCut") + pair.first,
                                                std::get<2>(pair.second));
        }
        if (printL(0)) {
            std::cout << " " << pair.first << "(";
            if (std::get<0>(pair.second) > 0)
                std::cout << std::get<0>(pair.second) << "+";
            std::cout << std::get<1>(pair.second) << "," <<  std::get<2>(pair.second) << ")";
        }
    }

#ifdef BCP_RCSP_IS_FOUND
    if (numDepotCapCutsWithOneY + numDepotCapCutsWithTwoYs > 0)
    {
        if (recordStats) {
            bapcodInit().statistics().incrCounter(std::string("bcCountRootActiveDCCwithOneY"),
                                                  numDepotCapCutsWithOneY);
            bapcodInit().statistics().incrCounter(std::string("bcCountRootActiveDCCwithTwoYs"),
                                                  numDepotCapCutsWithTwoYs);
        }
        if (printL(0))
            std::cout << " DCCwithOneY(" << numDepotCapCutsWithOneY << ")"
                      << " DCCwithTwoYs(" << numDepotCapCutsWithTwoYs << ")";
    }
    for (int numRows = 1; numRows <= RANK1_CUTS_MAX_NUMBER_OF_ROWS; ++numRows)
    {
        if (numPackingCuts[numRows] > 0)
        {
            if (recordStats)
                bapcodInit().statistics().incrCounter(std::string("bcCountRootActive") + std::to_string(numRows)
                                                      + "rowPackR1C", numPackingCuts[numRows]);
            if (printL(0))
                std::cout << " " <<  numRows << "rowPackR1C(" << numPackingCuts[numRows] << ")";
        }
        if (numCoveringCuts[numRows] > 0)
        {
            if (recordStats)
              bapcodInit().statistics().incrCounter(std::string("bcCountRootActive") + std::to_string(numRows)
                                                    + "rowCovR1C", numCoveringCuts[numRows]);
            if (printL(0))
                std::cout << " " << numRows << "rowCovR1C(" << numCoveringCuts[numRows] << ")";
        }
    }
    for (int oneKvalue = 0; oneKvalue <= maxOneKvalue; ++oneKvalue)
        if (numRLKCutsByOneK[oneKvalue] > 0)
        {
            if (recordStats)
            {
                if (oneKvalue == 0)
                    bapcodInit().statistics().incrCounter(std::string("bcCountRootActiveRLKCsRounding"),
                                                          numRLKCutsByOneK[0]);
                else
                    bapcodInit().statistics().incrCounter(std::string("bcCountRootActiveRLKCsOne")
                                                          + std::to_string(oneKvalue), numRLKCutsByOneK[oneKvalue]);
            }
            if (printL(0))
            {
                if (oneKvalue == 0)
                    std::cout << " RLKCsRounding" << oneKvalue << "(" << numRLKCutsByOneK[oneKvalue] << ")";
                else
                    std::cout << " RLKCs1/" << oneKvalue << "(" << numRLKCutsByOneK[oneKvalue] << ")";
            }
    }
#endif

    if (printL(0))
        std::cout << std::endl;
}

void Alg4EvalByColAndCutGen::drawPrimalSolutionToDotFile(std::string & filename)
{
    std::map<ColGenSpConf *, std::vector<MastColumn *> > colsInMasterSolutionMap;
    for (VarPtrSet::const_iterator varPtrIt = _probPtr->inPrimalLpSol().begin();
         varPtrIt != _probPtr->inPrimalLpSol().end(); ++varPtrIt)
        if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
        {
            MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
            auto mapIt = colsInMasterSolutionMap.find(colPtr->cgSpConfPtr());
            if (mapIt != colsInMasterSolutionMap.end())
            {
                mapIt->second.push_back(colPtr);
            }
            else
            {
                colsInMasterSolutionMap[colPtr->cgSpConfPtr()] = std::vector<MastColumn *>();
                colsInMasterSolutionMap[colPtr->cgSpConfPtr()].push_back(colPtr);
            }
        }

    for (auto mapIt = colsInMasterSolutionMap.begin(); mapIt != colsInMasterSolutionMap.end(); ++mapIt)
        mapIt->first->probPtr()->drawPrimalSolutionToDotFile(mapIt->second, filename);
}

void Alg4EvalByColAndCutGen::setDownAlgo()
{
  Alg4EvalBySimplexBasedColGen::setDownAlgo();
}
