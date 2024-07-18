/**
 *
 * This file bcAlg4DivingHeuristic.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4DivingHeuristic.hpp"
#include "bcModelC.hpp"
#include "bcRestrictedMasterIpHeuristicC.hpp"
#include "bcModelFormulationC.hpp"

#define DivingHeuristicPrintLevel 2

using namespace std;

LocalSearchHeuristic::LocalSearchHeuristic(Problem * probPtr,
                                           MasterCommons4PrimalHeuristic & masterCommons,
                                           MultitokenSelectionStrategyVector & colSelectionCriteria) :
    DivingHeuristic(probPtr, masterCommons, colSelectionCriteria), _maxIterationsNumber(1),
        _fixedVarsRatio(0.5)
{
}

LocalSearchHeuristic::~LocalSearchHeuristic()
{
}

DiveInfo::DiveInfo(const int depth, const int tabuSize, const int performRestrictedMaster_) :
    GenChildNodesInfo(), _tabuVariables(), remainingDepth(depth), maxTabuSize(tabuSize),
    performRestrictedMaster(performRestrictedMaster_)
{
}

DiveInfo::DiveInfo(const VarPtrSet & tabuVars, const int depth, const int tabuSize,
                   const int performRestrictedMaster_) :
    GenChildNodesInfo(), _tabuVariables(), remainingDepth(depth), maxTabuSize(tabuSize),
    performRestrictedMaster(performRestrictedMaster_)
{
    tabuVariables(tabuVars);
}

DiveInfo::DiveInfo(const DiveInfo & that):
    GenChildNodesInfo(), _tabuVariables(), remainingDepth(that.remainingDepth), maxTabuSize(that.maxTabuSize),
    performRestrictedMaster(0)
{
    tabuVariables(that.tabuVariables());
}

DiveInfo::~DiveInfo()
{
    for (VarPtrSet::const_iterator vpsIt = _tabuVariables.begin(); vpsIt != _tabuVariables.end(); ++vpsIt)
    {
        if ((*vpsIt)->isTypeOf(VcId::MastColumnMask))
            (*vpsIt)->decrParticipation(5);
    }
}

const VarPtrSet & DiveInfo::tabuVariables() const
{
    return _tabuVariables;
}

void DiveInfo::tabuVariables(const VarPtrSet & tabuVars)
{
    _tabuVariables = tabuVars;
    for (VarPtrSet::const_iterator vpsIt = _tabuVariables.begin(); vpsIt != _tabuVariables.end(); ++vpsIt)
    {
        if ((*vpsIt)->isTypeOf(VcId::MastColumnMask))
            (*vpsIt)->incrParticipation(10);
    }
}

bool Algorithm4DivingEval::setupAlgo(Node * nodePtr)
{
  if (Alg4EvalByColAndCutGen::setupAlgo(nodePtr))
    return true;

  DivingEvalInfo * divingEvalInfoPtr =
      dynamic_cast<DivingEvalInfo *>(nodePtr->nodeEvalInfoPtr());

  bapcodInit().require(divingEvalInfoPtr != NULL,
                       "BaPCod error: nodeEvalInfo for Algorithm4DivingEval is not of type DivingEvalInfo.");

  _nbNeededProperColumns = divingEvalInfoPtr->nbNeededProperColumns;
  _enumerationWithFalseGapMode = divingEvalInfoPtr->enumerationWithFalseGapMode;
  _runCutGeneration = divingEvalInfoPtr->runCutGeneration;

  return false;
}

NodeEvalInfo * Algorithm4DivingEval::recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr)
{
  DivingEvalInfo * divingEvalInfoPtr = NULL;
  if (nodeEvalInfoPtr != NULL)
    {
      divingEvalInfoPtr = dynamic_cast<DivingEvalInfo *>(nodeEvalInfoPtr);

      bapcodInit().require(divingEvalInfoPtr != NULL,
                           "BaPCod error: nodeEvalInfo passed to Algorithm4DivingEval::recordNodeEvalInfo"
                           " is not of type DivingEvalInfo");

      /// nbNeededProperColumns will be set in the diving heuristic
      divingEvalInfoPtr->nbNeededProperColumns = 0;
    }
  else
     /// nbNeededProperColumns will be set in the diving heuristic
    divingEvalInfoPtr = new DivingEvalInfo();

  return Alg4EvalByColAndCutGen::recordNodeEvalInfo(globalTreatOrder, divingEvalInfoPtr);
}

void Algorithm4DivingEval::addFixedSolutionBasedCutsToMaster()
{
  if ((_currentNodePtr->localFixedSolution() == NULL) || _currentNodePtr->localFixedSolution()->solVarValMap().empty())
    return;

  bool coreCutsArePresent = false;
  std::vector<GenericCutConstr *> coreCutConstrPts;
  for (std::set<GenericCutConstr *, DynamicGenConstrSort>::const_iterator genCutPtrIt
          = _masterCommons.candidateCutGenericConstr().begin();
       genCutPtrIt != _masterCommons.candidateCutGenericConstr().end(); ++genCutPtrIt)
    if ((*genCutPtrIt)->type() == 'C')
      coreCutConstrPts.push_back(*genCutPtrIt);

  if (coreCutConstrPts.empty())
    return;

  VarPtr2DoubleMap oldPartialSolution(_probPtr->partialSolution());
  VarPtr2DoubleMap fixedSolution(_currentNodePtr->localFixedSolution()->solVarValMap());
  for (VarPtr2DoubleMap::iterator mapIt = fixedSolution.begin(); mapIt != fixedSolution.end(); ++mapIt)
  {
    oldPartialSolution.erase(mapIt->first);
  }

  Time startCp;

  std::map<std::string, int> numCutGeneratedPerType; /// for printing purposes
  std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> generatedCutConstrSet;

  for (std::vector<GenericCutConstr *>::iterator genCutPtrIt = coreCutConstrPts.begin();
       genCutPtrIt != coreCutConstrPts.end(); ++genCutPtrIt)
  {
    int cutsBefore = generatedCutConstrSet.size();
    (*genCutPtrIt)->cutSeparationBasedOnFixedSol(oldPartialSolution, fixedSolution, generatedCutConstrSet);
    if (generatedCutConstrSet.size() - cutsBefore > 0)
      numCutGeneratedPerType[(*genCutPtrIt)->defaultName()] = generatedCutConstrSet.size() - cutsBefore;
  }

  double cutSepTime = startCp.getElapsedTime();
  bapcodInit().statistics().incrTimer("bcTimeCutSeparation", cutSepTime);

  addCutsToProblem('C', cutSepTime, numCutGeneratedPerType, generatedCutConstrSet);
}

bool Algorithm4DivingEval::eval()
{
  /// we first try to complete the current partial solution with a heuristic solution for the residual problem
  if (_probPtr->runMasterHeuristicFunctor() && checkIfCurSolIsMasterLpFeasible())
    {
      _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());
        updateAlgPrimalLpBounds();
      if (checkIfCurSolIsInteger() && !addCutToMaster('C'))
        {
          if (printL(0))
            std::cout << "Master heuristic function returned solution of value " << _algCurLpPrimalBound << std::endl;
          updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
        }
    }

  addFixedSolutionBasedCutsToMaster();

  if ((_enumerationWithFalseGapMode == 1)
      || ((_enumerationWithFalseGapMode == 2) && (_masterCommons.totalNumberOfEnumeratedSubprobSolutions() < 0)))
  {
      double paramMaxRelativeGap = 0.1; /// make a parameter??

      int saveDoRedCostFixingAndEnumeration = _doRedCostFixingAndEnumeration;
      _doRedCostFixingAndEnumeration = 0;
      if (Alg4EvalBySimplexBasedColGen::eval())
          return true; /// Problem infeasible

      _doRedCostFixingAndEnumeration = 2;
      int numReductions = 0;
      bool repeatColumnGeneration = false;
      double falseGap = algIncIpPrimalBound() - algCurLpDualBound();
      if ((algCurLpDualBound().positive())
          && ((double)algIncIpPrimalBound() > (double)(algCurLpDualBound() * (1 + paramMaxRelativeGap))))
          falseGap = paramMaxRelativeGap * algCurLpDualBound();
      do
      {
          if (printL(0))
              std::cout << "False gap = " << falseGap << std::endl;

          repeatColumnGeneration = runReducedCostFixingAndEnumeration(_enumerationWithFalseGapMode, falseGap);

          falseGap *= param().EnumHeuristicFalseGapMultiplier();
          numReductions += 1;
      }
      while ((numReductions < param().EnumHeuristicNumberOfTries())
             && (_masterCommons.totalNumberOfEnumeratedSubprobSolutions() < 0));

      if (!repeatColumnGeneration)
          return false;

      //_doRedCostFixingAndEnumeration = saveDoRedCostFixingAndEnumeration;
  }


  if (_runCutGeneration && Alg4EvalByColAndCutGen::eval())
      return true;
  if (!_runCutGeneration && Alg4EvalBySimplexBasedColGen::eval())
      return true;

    /// if only proper columns are generated, we stop
  if (param().GenerateProperColumns() || _solIsInteger || isConquered()
      || (_masterCommons.totalNumberOfEnumeratedSubprobSolutions() >= 0))
    return false;

  /// otherwise, as the solution may consist of only non-proper columns, we might want to generate proper
  /// columns heuristically
  _currentNodePtr->recordPrimalSol(_probPtr->inPrimalLpSol());

  /// the user might want to add additional columns in the solution so that there is a larger choice of proper columns
  /// this is also used in the case column contain a partial information: in this case the user generate here
  /// the columns with complete information based on columns with partial information in the solution
  VarPtrSet tabuVariables;
  std::vector<MastColumn *> columnsToAdd;
  DiveInfo * diveInfoPtr = dynamic_cast<DiveInfo *>(_currentNodePtr->genChildNodesInfoPtr());
  if (diveInfoPtr != NULL)
    tabuVariables = diveInfoPtr->tabuVariables();
  _probPtr->runDivingColCutGenTerminatedFunctor(_currentNodePtr->primalSol(), tabuVariables, columnsToAdd);
  for (std::vector<MastColumn *>::iterator colPtrIt = columnsToAdd.begin(); colPtrIt != columnsToAdd.end();
       ++colPtrIt)
    {
       (*colPtrIt)->val(2 * Double::precision);
       _currentNodePtr->addToPrimalSol(*colPtrIt);
    }

  /// we find the maximum priority level of variables present in the solution
  /// the heuristic pricing solver may be called only for subproblems of this priority level or higher
  SelectionStrategy::PriorityEnum firstSelectedRule = param().RoundingColSelectionCriteria().front().selectedRule();
  Double maxPriorityLevelInSolution(0);
  if (firstSelectedRule == SelectionStrategy::HighestPriority)
    {
      for (SolutionVarInfoPtrList::const_iterator infoIt = _currentNodePtr->primalSol().begin();
           infoIt != _currentNodePtr->primalSol().end(); infoIt++)
        {
          if (maxPriorityLevelInSolution < (*infoIt)->priorityLevel)
            maxPriorityLevelInSolution = (*infoIt)->priorityLevel;
        }
    }

  std::list<ColGenSpConf *> cgSpConfWithSufficientPriorityList;
  for (std::vector<ColGenSpConf *>::const_iterator cgspIt = _masterCommons.colGenSubProbConfPts().begin();
       cgspIt != _masterCommons.colGenSubProbConfPts().end(); ++cgspIt)
    {
      if ( (*((*cgspIt)->upperBoundPtr()) > 0) && ((*cgspIt)->priorityLevel() >= maxPriorityLevelInSolution) )
        cgSpConfWithSufficientPriorityList.push_back(*cgspIt);
    }

  int nbProperColumnsNeeded = _nbNeededProperColumns;
  for (SolutionVarInfoPtrList::const_iterator infoIt = _currentNodePtr->primalSol().begin();
       (infoIt != _currentNodePtr->primalSol().end()) && (nbProperColumnsNeeded > 0); infoIt++)
    {
      if (param().FixIntValBeforeRoundingHeur() && !(*infoIt)->value.fractional())
        {
          nbProperColumnsNeeded = 0;
        }
      else if (((*infoIt)->canRoundDown || (*infoIt)->canRoundUp)
               && ((*infoIt)->priorityLevel >= maxPriorityLevelInSolution)
               && (*infoIt)->value.fractional() )
        {
          nbProperColumnsNeeded -= 1;
        }
    }

  if (nbProperColumnsNeeded == 0)
    return false;

  /// there is no variable to fix or round, so we need to generate a proper column

  Time start;
  std::set<MastColumn *, SortMastColumnPerNonDecreasingRedCost> colPtsSortedByRedCost;
  for (std::list<ColGenSpConf *>::iterator cgspIt = cgSpConfWithSufficientPriorityList.begin();
       cgspIt != cgSpConfWithSufficientPriorityList.end(); ++cgspIt)
    {
      int maxLevelOfSubProbRestriction(-1);
      Problem * spProbPtr = (*cgspIt)->probPtr();
      spProbPtr->setPrimalLpBound(Bound::infPrimalBound(_masterCommons.objStatus()));
      spProbPtr->setDualBound(Bound::infDualBound(_masterCommons.objStatus()));
      spProbPtr->resetSolution();
      spProbPtr->clearRecordedSol();
      if (spProbPtr->solveProb(maxLevelOfSubProbRestriction) <= 0)
        continue;
      MastColumn * colPtr = (*cgspIt)->recordSubproblemSolution(spProbPtr->retrieveCurPrimalLpSol(), false, 3);
      colPtr->computeReducedCost();
      colPtsSortedByRedCost.insert(colPtr);

      const std::list<Solution *> & spSolList = spProbPtr->recordedSolList();
      for (std::list <Solution * >::const_iterator spSolPt = spSolList.begin(); spSolPt != spSolList.end(); spSolPt++)
        {
          colPtr = (*cgspIt)->recordSubproblemSolution(*spSolPt, false,3);
          colPtr->computeReducedCost();
          colPtsSortedByRedCost.insert(colPtr);
        }
      /// we just clear _tempColPtrList4Insertion of (*cgspIt) instead of insertColumnsInMaster
      /// otherwise the segmentation fault at the end
      (*cgspIt)->clearColPtrList4Insertion();
    }
  bapcodInit().statistics().incrTimer("bcTimeHeurSpOracle", start.getElapsedTime_dbl());
  int colCounter = 0;
  Double maxReducedCost;
  std::set<MastColumn *, SortMastColumnPerNonDecreasingRedCost>::iterator colPtrIt = colPtsSortedByRedCost.begin();
  while ((colPtrIt != colPtsSortedByRedCost.end()) && (colCounter < nbProperColumnsNeeded))
    {
      maxReducedCost = (*colPtrIt)->reducedCost();
      colCounter += 1;
      ++colPtrIt;
    }
  colPtsSortedByRedCost.erase(colPtrIt, colPtsSortedByRedCost.end());
  if (colPtsSortedByRedCost.size() < nbProperColumnsNeeded)
    maxReducedCost = BapcodInfinity;

  /// we look over all active columns in the master to check whether there exists a proper column
  /// with small enough reduced cost
  solveRestrictedMastLP();
  _probPtr->retrieveRedCosts();   /// fast reduced cost retrieval
  /// reduced cost of columns generated by the heuristic were nullified, we need to generate them again
  for (colPtrIt = colPtsSortedByRedCost.begin(); colPtrIt != colPtsSortedByRedCost.end(); ++colPtrIt)
    {
      (*colPtrIt)->computeReducedCost();
      if (printL(0))
        std::cout << "Column " << (*colPtrIt)->name() << " is generated by heuristic with reduced cost "
                  << (*colPtrIt)->reducedCost() << ", status is " << (*colPtrIt)->vcIndexStatus() << std::endl;
    }
  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); varPt++)
    if ((*varPt)->isTypeOf(VcId::MastColumnMask) && (*varPt)->val().isZero())   /// columns with positive value
                                                                                /// are in the solution and are
                                                                                /// already checked
      {
        MastColumn * colPtr = static_cast<MastColumn *>(*varPt);
        /// TO DO: we need also to check here whether this variable is in the tabu list or not
        if ((colPtr->cgSpConfPtr()->priorityLevel() >= maxPriorityLevelInSolution)
            && colPtr->suitableToFixValue(1.0) && (colPtr->reducedCost() < maxReducedCost))
          {
            if (printL(0))
              std::cout << "Column " << colPtr->name() << " is active in the master with reduced cost "
                        << colPtr->reducedCost() << std::endl;
            colPtsSortedByRedCost.insert(colPtr);
            if (colPtsSortedByRedCost.size() > nbProperColumnsNeeded)
              {
                colPtsSortedByRedCost.erase(--colPtsSortedByRedCost.end());
                colPtrIt = --colPtsSortedByRedCost.end();
                maxReducedCost = (*colPtrIt)->reducedCost();
              }
          }
      }

  if (!colPtsSortedByRedCost.empty())
    {
      maxReducedCost = (*(--colPtsSortedByRedCost.end()))->reducedCost();
      Double minReducedCost = (*(colPtsSortedByRedCost.begin()))->reducedCost();
      if (printL(0))
        std::cout << colPtsSortedByRedCost.size() << " (from "<< nbProperColumnsNeeded
                  << " needed) proper columns generated by heuristic "
                  << "or found in the master problem, their reduced cost is from "
                  << minReducedCost << " to " << maxReducedCost << std::endl;
      for (colPtrIt = colPtsSortedByRedCost.begin(); colPtrIt != colPtsSortedByRedCost.end(); ++colPtrIt)
        {
          (*colPtrIt)->val(2 * Double::precision);
          _currentNodePtr->addToPrimalSol(*colPtrIt);
        }
    }
  else if (printL(0))
    {
      std::cout << "No proper columns generated by heuristic or found in the master problem, we needed "
                << nbProperColumnsNeeded << std::endl;
    }

  return false;
}


bool DiveAlgorithm::setupAlgo(Node * nodePtr)
{
  if (Alg4GenChildrenOfNode::setupAlgo(nodePtr))
    return true;

  bapcodInit().require(_currentNodePtr->genChildNodesInfoPtr() != NULL,
                       "BaPCod error: genChildNodesInfoPtr for DiveAlgorithm is null.");

  DiveInfo * diveInfoPtr = dynamic_cast<DiveInfo *>(_currentNodePtr->genChildNodesInfoPtr());

  bapcodInit().require(diveInfoPtr != NULL,
                       "BaPCod error: genChildNodesInfoPtr for DiveAlgorithm is not of type DiveInfo.");

  _tabuVariables = diveInfoPtr->tabuVariables();
  _remainingDepth = diveInfoPtr->remainingDepth;
  _maxTabuSize = diveInfoPtr->maxTabuSize;

  return false;
}

bool DiveAlgorithm::fixVariables()
{
  if (!param().FixIntValBeforeRoundingHeur)
    return false;

  bool thereAreFixedVariables = false;
  Solution *solPtr = new Solution();

  SelectionStrategy::PriorityEnum firstSelectedRule = _colSelectionCriteria.front().selectedRule();
  Double maxPriorityLevel(0);
  if (firstSelectedRule == SelectionStrategy::HighestPriority)
    {
      /// we find the maximum priorityLevel of variables which can be fixed or rounded
      for (SolutionVarInfoPtrList::const_iterator infoIt = _currentNodePtr->primalSol().begin();
           infoIt != _currentNodePtr->primalSol().end(); infoIt++)
        if (((*infoIt)->canRoundDown || (*infoIt)->canRoundUp)
            && ((*infoIt)->priorityLevel > maxPriorityLevel)
            && !_tabuVariables.count((*infoIt)->varPtr))
          {
            maxPriorityLevel = (*infoIt)->priorityLevel;
          }
    }

  for (SolutionVarInfoPtrList::const_iterator infoIt = _currentNodePtr->primalSol().begin();
       infoIt != _currentNodePtr->primalSol().end(); infoIt++)
    {
      if (!param().GenerateProperColumns() && thereAreFixedVariables)
          break;

      if ((*infoIt)->varPtr->type() == 'C')
        continue;

      if (param().UseDivingHeurOnMasterColOnly() || param().UseDivingHeurOnPureMastVarOnly() )
      {
        if ((*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
        {
          if (param().UseDivingHeurOnPureMastVarOnly() )
            continue;
        }
        else if (param().UseDivingHeurOnMasterColOnly())
        {
          continue;
        }
      }

      if (!(*infoIt)->value.fractional() && !(*infoIt)->value.isZero() && (*infoIt)->canRoundDown
          && (*infoIt)->canRoundUp)
        {
          if (printL(0))
            std::cout << "Fixed variable " << (*infoIt)->varPtr->name() << " at value " << (*infoIt)->value << std::endl;
          solPtr->includeVar((*infoIt)->varPtr, (*infoIt)->value, true);
          thereAreFixedVariables = true;
        }
        else if (((*infoIt)->value > 1.0) && (*infoIt)->canRoundDown)
        {
            if (printL(0))
              std::cout << "Fixed variable " << (*infoIt)->varPtr->name() << " at value "
                        << Dfloor((*infoIt)->value) << std::endl;
            solPtr->includeVar((*infoIt)->varPtr, Dfloor((*infoIt)->value), true);
            thereAreFixedVariables = true;
        }
    }

  if (thereAreFixedVariables)
    {
      std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
      Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                        tmpLocalNodeBrConstrList, solPtr);
      newChildNodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_tabuVariables, _remainingDepth, _maxTabuSize));

      if (!param().GenerateProperColumns())
        {
          DivingEvalInfo * divingEvalInfoPtr = dynamic_cast<DivingEvalInfo *>(newChildNodePtr->nodeEvalInfoPtr());
          bapcodInit().require(divingEvalInfoPtr != NULL,
                               "BaPCod error: nodeEvalInfo of a diving node is not of type divingEvalInfo");
          divingEvalInfoPtr->nbNeededProperColumns = ((_remainingDepth > 0) ? _maxTabuSize + 1
                                                                            : _tabuVariables.size() + 1);
        }

      _currentNodePtr->sons().push_back(newChildNodePtr);
      return true;
    }

  delete solPtr;
  return false;
}

bool DiveAlgorithm::solutionCausesInfeasibility(int & globalTreatOrder, Solution * solPtr)
{
 if (!param().DivingHeurPreprocessBeforeChoosingVar())
    return false;

  bool infeasibleAfterPreprocessing = false;

  std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
  Node * tempNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                tmpLocalNodeBrConstrList, solPtr->clone());
  tempNodePtr->setPreprocessor(new Algorithm4PreprocessingInDive(_masterCommons.problemList()));
  tempNodePtr->setEvalAlg(new Alg4EvalByNone(_masterCommons.masterCommons4EvalAlg()));
  if (tempNodePtr->probSetupInfoPtr()->treatOrderId == globalTreatOrder)
    tempNodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupOfNode(_masterCommons.masterCommons4ProblemSetup()));
  else
    tempNodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup()));
  tempNodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(_masterCommons.masterCommons4ProblemSetup()));

  tempNodePtr->evaluation(globalTreatOrder, _currentNodePtr->nodeIncIpPrimalBound());

  infeasibleAfterPreprocessing = tempNodePtr->infeasible();
  delete tempNodePtr;

  return infeasibleAfterPreprocessing;
}

/// the implementation of this function is not the most efficient
/// in a more efficient version, we need to sort rounding candidates
/// according to the selection criteria, and then take the first candidate
/// which does not cause infeasibility (by preprocessing)
Solution * DiveAlgorithm::roundVariable(int & globalTreatOrder)
{
  VarPtrSet varsCausingInfeasibility;
  Double incumbentRoundedValue(0);
  SolutionVarInfo * incumbentInfoPtr = NULL;
  bool roundingCausesInfeasibility = false;
  do
    {
      for (std::vector<SelectionStrategy>::iterator ssIt = _colSelectionCriteria.begin();
           ssIt != _colSelectionCriteria.end(); ++ssIt)
        ssIt->initializeIncumbent();
      incumbentInfoPtr = NULL;
      for (SolutionVarInfoPtrList::const_iterator infoIt = _currentNodePtr->primalSol().begin();
           infoIt != _currentNodePtr->primalSol().end(); infoIt++)
        {
          if ((*infoIt)->varPtr->type() == 'C')
            continue;

          if (param().UseDivingHeurOnMasterColOnly() || param().UseDivingHeurOnPureMastVarOnly())
            {
              if ((*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
                {
                  if (param().UseDivingHeurOnPureMastVarOnly)
                    continue;
                }
              //TODO: add param RhMinRoundingPriorityForPurMastVars
              else if (param().UseDivingHeurOnMasterColOnly || (*infoIt)->varPtr->priority() < 0 )
                {
                  continue;
                }
            }

          if (_tabuVariables.count((*infoIt)->varPtr) || varsCausingInfeasibility.count((*infoIt)->varPtr))
            continue;

          if (printL(DivingHeuristicPrintLevel))
            std::cout << "DiveAlgorithm::roundVariable() considers var "
                      << (*infoIt)->varPtr->name() << " with value = " << (*infoIt)->value
                      << " floor = " << Dfloor((*infoIt)->value)
                      << " ceil = " << Dceil((*infoIt)->value)
                      << " canRoundDown = " << (*infoIt)->canRoundDown
                      << " canRoundUp = " << (*infoIt)->canRoundUp
                      << " priority = "<< (*infoIt)->varPtr->priority()<< std::endl;

          if (param().IgnoreIntValWhenSelectingRoundedVarInRH() && !(*infoIt)->value.fractional())
            {
              if (printL(DivingHeuristicPrintLevel))
                std::cout << "DiveAlgorithm::roundVariable() IGNORES col  "
                          << (*infoIt)->varPtr->name() << " whose cur val is " << (*infoIt)->value << std::endl;
              continue;
            }

          int status = 0;
          Double challengerRoundedValue;
          std::vector<SelectionStrategy>::iterator ssIt;
          for (ssIt = _colSelectionCriteria.begin(); ssIt != _colSelectionCriteria.end(); ++ssIt)
            {
              status = ssIt->computeCriteriaAndValToCompareToIncumbent(*infoIt, challengerRoundedValue,
                                                                       DivingHeuristicPrintLevel);
              /// status == -1 =>  incumbent doninates challenger
              /// status == 1 =>   challenger doninates incumbent, therefore it defines the new incumbent
              /// status == 0 =>  no dominance

              if (status == 1)
                {
                  incumbentInfoPtr = *infoIt;
                  incumbentRoundedValue = challengerRoundedValue;
                  if (printL(DivingHeuristicPrintLevel))
                    std::cout << "selectNextCol4DivingHeur: new incumbent is " << incumbentInfoPtr->varPtr->name()
                              << " incumbentRoundedValue = " << incumbentRoundedValue << std::endl;
                }
              if (status != 0)
                break;
            }

          /**
           * If incumbent has been replaced (status == 1),
           * update criteria value of incumbent for next criteria
           */
          if ((status == 1) && (ssIt != _colSelectionCriteria.end()))
            {
              for (++ssIt; ssIt != _colSelectionCriteria.end(); ++ssIt)
                {
                  ssIt->initializeIncumbent();
                  ssIt->computeCriteriaToCompareToIncumbent(*infoIt, challengerRoundedValue,
                                                            DivingHeuristicPrintLevel);
                }
            }
        } /// end of for(VarPtrSet::const_iterator mastVarPt =

      if (incumbentInfoPtr != NULL)
        {
          Solution * localSolPtr = new Solution();
          localSolPtr->includeVar(incumbentInfoPtr->varPtr, incumbentRoundedValue, true);
          /// this function runs preprocessing algorithm to verify whether an infeasibility can be determined
          roundingCausesInfeasibility = solutionCausesInfeasibility(globalTreatOrder, localSolPtr);
          if (roundingCausesInfeasibility)
            varsCausingInfeasibility.insert(incumbentInfoPtr->varPtr);
          delete localSolPtr;
        }
    }
  while ((incumbentInfoPtr != NULL) && roundingCausesInfeasibility);

  if (incumbentInfoPtr == NULL)
    {
      if (printL(0))
        std::cout << "DiveAlgorithm::roundVariable() SELECTS NO var  " << std::endl;
      return NULL;
    }

  if (printL(0))
  {
    std::cout << "Selected variable " << incumbentInfoPtr->varPtr->name()
              << " with priority " << incumbentInfoPtr->priorityLevel
              << " with value " << incumbentInfoPtr->value
              << " and reduced cost " << incumbentInfoPtr->reducedCost << std::endl;
    if (incumbentInfoPtr->varPtr->isTypeOf(VcId::MastColumnMask))
      static_cast<MastColumn *>(incumbentInfoPtr->varPtr)->spSol()->printOrderedSolution(std::cout);
  }

  Solution * solPtr = new Solution();
  solPtr->includeVar(incumbentInfoPtr->varPtr, incumbentRoundedValue, true);
  return solPtr;
}

void DiveAlgorithm::prepareCandidateNodeInStrongDive(Node * nodePtr, int & globalTreatOrder)
{
    nodePtr->setPreprocessor(new Algorithm4PreprocessingInDive(_masterCommons.problemList()));

  Algorithm4DivingEval * evalAlgPtr = new Algorithm4DivingEval(_masterCommons.problemList().front(),
                                                               _masterCommons.masterCommons4EvalAlg());
  /// we need to set this option to false, otherwise there the algorithm will be stopped
  /// as there may be a gap and an integer solution
  evalAlgPtr->setOptionNeedBeConqueredIfSolIsInteger(false);
  evalAlgPtr->setOptionMaxNbOfPenaltyUpdates(param().MaxNbOfPenaltUpdatesDuringRH());
  evalAlgPtr->setOptionNonExactEvaluation(true);

  const StrongBranchingPhaseParameter & evalAlgParams = param().EvalAlgParamsInDiving();
  if (evalAlgParams.active())
  {
    evalAlgPtr->setOptionMaxNbOfCgIterations(evalAlgParams.maxNumOfColGenIterations());
    evalAlgPtr->setOptionMaxNbOfCutRounds(evalAlgParams.maxNumCutRounds());
    if (evalAlgParams.exact())
      evalAlgPtr->setOptionMinNbOfCutRounds(param().MinNumOfCutRoundsBeforeStopBySp());
    else
      evalAlgPtr->setOptionMinNbOfCutRounds(evalAlgParams.minNumCutRounds());
    if (evalAlgParams.minLevelOfSpRestriction() < param().MaxNbOfStagesInColGenProcedure())
      evalAlgPtr->setOptionMinLevelOfSpRestriction(evalAlgParams.minLevelOfSpRestriction());
    else
      evalAlgPtr->setOptionMinLevelOfSpRestriction(param().MaxNbOfStagesInColGenProcedure() - 1);
    if (evalAlgParams.doRedCostFixingAndEnumeration())
      evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(2);
    else
      evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(0);
    evalAlgPtr->setOptionLogPrintFrequency(0);
  }
  else /// evalAlgParams is not active
  {
    evalAlgPtr->setOptionMaxNbOfCgIterations(param().MaxNbOfCgIteDuringRh());
    evalAlgPtr->setOptionLogPrintFrequency(0);
  }

   nodePtr->setEvalAlg(evalAlgPtr);
    if (nodePtr->probSetupInfoPtr()->treatOrderId == globalTreatOrder)
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupOfNode(_masterCommons.masterCommons4ProblemSetup()));
    else
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup()));
    nodePtr->setProblemSetDownAlgorithm(new ProblemFullSetDownAlgorithm(_masterCommons.masterCommons4ProblemSetup()));
}

void DiveAlgorithm::runStrongDive(int & globalTreatOrder, const int maxNumberOfCandidates,
                                  const int maxNumberOfDives, const bool inheritParentDualBound)
{
  std::vector<Node *> candidateNodes;
  VarPtrSet savedTabuVariables(_tabuVariables);
    for (VarPtrSet::iterator vpsIt = savedTabuVariables.begin(); vpsIt != savedTabuVariables.end(); ++vpsIt)
        if ((*vpsIt)->isTypeOf(VcId::MastColumnMask))
            (*vpsIt)->incrParticipation(11);
  int numNodesWithSameLpValue = 0;
  int candidateNumber = 0;

  Problem * mastProbPtr = _masterCommons.problemList().front();
  do
    {
      Solution * solPtr = NULL;
      if (mastProbPtr->divingFixingFunctorDefined())
        solPtr = mastProbPtr->runDivingFixingFunctor(_currentNodePtr->probSetupInfoPtr()->masterPartialSolutionInfo,
                                                     _currentNodePtr->primalSol(), _tabuVariables);
      else
        solPtr = roundVariable(globalTreatOrder);

      if (solPtr == NULL)
        break;

      Variable * firstRoundedVarPtr = NULL;
      if (!solPtr->solVarValMap().empty())
        firstRoundedVarPtr = solPtr->solVarValMap().begin()->first;
      for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
           mapIt != solPtr->solVarValMap().end(); ++mapIt)
        _tabuVariables.insert(mapIt->first);

      std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
      Node * newCandidateNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                            tmpLocalNodeBrConstrList, solPtr, inheritParentDualBound);
      prepareCandidateNodeInStrongDive(newCandidateNodePtr, globalTreatOrder);
      candidateNumber += 1;

      newCandidateNodePtr->evaluation(globalTreatOrder, _currentNodePtr->nodeIncIpPrimalBound());
      if (newCandidateNodePtr->infeasible())
        {
          if (printL(0))
            std::cout << "Strong diving candidate " << candidateNumber
                      << " ( " << firstRoundedVarPtr->name() << " ) is infeasible" << std::endl;
          delete newCandidateNodePtr;
        }
      else
        {
          if (newCandidateNodePtr->primalBoundIsUpdated())
            _currentNodePtr->updateNodeIncPrimalSolution(newCandidateNodePtr->nodeIncIpPrimalSolPtr());

          if (newCandidateNodePtr->isConquered())
            {
              if (printL(0))
                std::cout << "Strong diving candidate " << candidateNumber
                          << " ( " << firstRoundedVarPtr->name() << " ) is conquered" << std::endl;
              delete newCandidateNodePtr;
            }
          else
            {
              if (newCandidateNodePtr->nodeIncLpDualBound() <= _currentNodePtr->nodeIncLpDualBound())
                numNodesWithSameLpValue += 1;
              candidateNodes.push_back(newCandidateNodePtr);

              if (printL(0))
                std::cout << "Strong diving candidate " << candidateNumber
                          << " ( " << firstRoundedVarPtr->name() << " ) lp value is "
	  		              << std::setprecision(10) << newCandidateNodePtr->nodeIncLpPrimalBound() << std::endl;
            }
        }
    }
  while ((candidateNumber < maxNumberOfCandidates) && (numNodesWithSameLpValue < maxNumberOfDives));

  if ((candidateNumber < maxNumberOfCandidates) && (numNodesWithSameLpValue == maxNumberOfDives) && printL(0))
    std::cout << "Strong diving is stopped prematurely as other candidats could not improve found ones" << std::endl;

  std::stable_sort(candidateNodes.begin(), candidateNodes.end(), SmallestNodeLpValue());

  _tabuVariables = savedTabuVariables;
  for (VarPtrSet::iterator vpsIt = savedTabuVariables.begin(); vpsIt != savedTabuVariables.end(); ++vpsIt)
    if ((*vpsIt)->isTypeOf(VcId::MastColumnMask))
      (*vpsIt)->decrParticipation(6);
  for (int it = 0; it < candidateNodes.size(); it++)
    {
      Node * curNodePtr = candidateNodes[it];
      if ((it < maxNumberOfDives) && !curNodePtr->treated())
        {
          curNodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_tabuVariables, _remainingDepth, _maxTabuSize));
          _currentNodePtr->sons().push_back(curNodePtr);
          const VarPtr2DoubleMap & nodeSolVarValMap = curNodePtr->localFixedSolution()->solVarValMap();
          for (VarPtr2DoubleMap::const_iterator mapIt = nodeSolVarValMap.begin();
               mapIt != nodeSolVarValMap.end(); ++mapIt)
            _tabuVariables.insert(mapIt->first);
        }
      else
        delete curNodePtr;
    }
}

void DiveAlgorithm::run(int & globalTreatOrder)
{
  if (_currentNodePtr->nodeIncLpPrimalBound() >= _currentNodePtr->nodeIncIpPrimalBound())
  {
    if (printL(0))
      std::cout << "Diving node is prunned by primal bound (" << _currentNodePtr->nodeIncLpPrimalBound()
                << " >= " <<  _currentNodePtr->nodeIncIpPrimalBound() << ")" << std::endl;
    return;
  }

  Problem * mastProbPtr = _masterCommons.problemList().front();
  if (!mastProbPtr->divingFixingFunctorDefined() && fixVariables())
    return;

  _remainingDepth -= 1;

  std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;

  if (param().StrongDivingCandidatesNumber() > 1)
    {
      int maxNumberOfDives = _maxTabuSize + 1 - _tabuVariables.size();
      if (_remainingDepth < 0)
        maxNumberOfDives = 1;
      runStrongDive(globalTreatOrder, param().StrongDivingCandidatesNumber(), maxNumberOfDives, true);
      return;
    }
  VarPtrSet varsCausingInfesibility;
  do
    {
      Solution * solPtr = NULL;
      if (mastProbPtr->divingFixingFunctorDefined())
        solPtr =  mastProbPtr->runDivingFixingFunctor(_currentNodePtr->probSetupInfoPtr()->masterPartialSolutionInfo,
                                                      _currentNodePtr->primalSol(), _tabuVariables);
      else
        solPtr = roundVariable(globalTreatOrder);
      if (solPtr == NULL)
        {
          if (param().UseDivingHeurOnMasterColOnly())
          {
            /// we create the node which will only run the restricted master
            Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                              tmpLocalNodeBrConstrList);
            newChildNodePtr->associateGenChildNodesInfoPtr(new DiveInfo(0, 0, 2));
            _currentNodePtr->sons().push_back(newChildNodePtr);
          }
          break;
        }

      Variable * firstRoundedVarPtr = (solPtr->solVarValMap().empty() ? NULL
                                                                      : solPtr->solVarValMap().begin()->first);

      bool requireReoptimizationWithColGen = false;
      if ( !param().runColGenAfterFixingPureMastVarInDiving()
           && firstRoundedVarPtr->isTypeOf(VcId::MastColumnMask)
           && (_currentNodePtr->localFixedSolution() != NULL)
           && (_currentNodePtr->localFixedSolution()->solVarValMap().begin())->first->isTypeOf(VcId::InstMasterVarMask))
         /// previous optimization was without colGen (since dive was on pure mast var).
         /// we need additional node without any solution to fix , just to reoptimize with colGen
        {
          if (printL(0))
            std::cout << "DiveAlgorithm::run() COLGEN REOPT REQUIRED. UNSELECTS var "
                      << firstRoundedVarPtr->name() << std::endl;
          delete solPtr;
          solPtr = NULL;
          firstRoundedVarPtr = NULL;
          requireReoptimizationWithColGen = true;
        }

      Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                        tmpLocalNodeBrConstrList, solPtr);

      if (!param().GenerateProperColumns())
        {
          DivingEvalInfo * divingEvalInfoPtr = dynamic_cast<DivingEvalInfo *>(newChildNodePtr->nodeEvalInfoPtr());
          bapcodInit().require(divingEvalInfoPtr != NULL,
                               "BaPCod error: nodeEvalInfo of a diving node is not of type divingEvalInfo");
          divingEvalInfoPtr->nbNeededProperColumns = ((_remainingDepth > 0) ? _maxTabuSize + 1
                                                                            : _tabuVariables.size() + 1);
        }

      newChildNodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_tabuVariables, _remainingDepth, _maxTabuSize, 0));
      _currentNodePtr->sons().push_back(newChildNodePtr);

      if (firstRoundedVarPtr == NULL) //the last child added for the restrMip is unique.
        {
          break;
        }
      else
        {
          for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
               mapIt != solPtr->solVarValMap().end(); ++mapIt)
            _tabuVariables.insert(mapIt->first);
        }
    }
  while ((_remainingDepth >= 0) && (_tabuVariables.size() <= _maxTabuSize));
}

void DiveAlgorithm::setDownAlgo()
{
  _tabuVariables.clear();
  Alg4GenChildrenOfNode::setDownAlgo();
}

void DivingHeuristic::setOptionMaxDepth(const int value)
{
  _maxDepth = value;
}

void DivingHeuristic::setOptionMaxDiscrepancy(const int value)
{
  _maxDiscrepancy = value;
}

void DivingHeuristic::setOptionStopAfterFirstSolutionFound(const bool value)
{
    _stopAfterFirstSolutionFound = value;
}

void DivingHeuristic::setOptionRunRestrMasterAfterFalseGapEnumeration(const bool value)
{
    _runRestrMasterAfterFalseGapEnumeration = value;
}

DivingHeuristic::DivingHeuristic(Problem * probPtr,
                                 MasterCommons4PrimalHeuristic & masterCommons,
                                 MultitokenSelectionStrategyVector & colSelectionCriteria) :
    Alg4PrimalHeuristicOfNode(probPtr, masterCommons), _stopAfterFirstSolutionFound(false), _maxDepth(0),
    _maxDiscrepancy(0), _runRestrMasterAfterFalseGapEnumeration(false), _colSelectionCriteria(colSelectionCriteria)
{
}

DivingHeuristic::~DivingHeuristic()
{
}

Alg4EvalOfNode * DivingHeuristic::createEvaluationAlgorithm(const bool runColGen)
{
      Algorithm4DivingEval * divingEvalAlgPtr = new Algorithm4DivingEval(_problemPtr,
                                                                         _masterCommons.masterCommons4EvalAlg());
      /// we need to set this option to false, otherwise there the algorithm will be stopped
      /// as there may be a gap and an integer solution
      divingEvalAlgPtr->setOptionNeedBeConqueredIfSolIsInteger(false);
      divingEvalAlgPtr->setOptionMaxNbOfPenaltyUpdates(param().MaxNbOfPenaltUpdatesDuringRH());
      divingEvalAlgPtr->setOptionNonExactEvaluation(true);

      const StrongBranchingPhaseParameter & evalAlgParams = param().EvalAlgParamsInDiving();
      if (evalAlgParams.active())
      {
        if (runColGen)
          divingEvalAlgPtr->setOptionMaxNbOfCgIterations(evalAlgParams.maxNumOfColGenIterations());
        else
          divingEvalAlgPtr->setOptionMaxNbOfCgIterations(0);
        divingEvalAlgPtr->setOptionMaxNbOfCutRounds(evalAlgParams.maxNumCutRounds());
        if (evalAlgParams.exact())
          divingEvalAlgPtr->setOptionMinNbOfCutRounds(param().MinNumOfCutRoundsBeforeStopBySp());
        else
          divingEvalAlgPtr->setOptionMinNbOfCutRounds(evalAlgParams.minNumCutRounds());
        if (evalAlgParams.minLevelOfSpRestriction() < param().MaxNbOfStagesInColGenProcedure())
          divingEvalAlgPtr->setOptionMinLevelOfSpRestriction(evalAlgParams.minLevelOfSpRestriction());
        else
          divingEvalAlgPtr->setOptionMinLevelOfSpRestriction(param().MaxNbOfStagesInColGenProcedure() - 1);
        if (evalAlgParams.doRedCostFixingAndEnumeration())
          divingEvalAlgPtr->setOptionDoRedCostFixingAndEnumeration(1);
        else
          divingEvalAlgPtr->setOptionDoRedCostFixingAndEnumeration(0);
        if (evalAlgParams.exact())
        {
          int logFrequency = param().ColGenLogFrequency();
          if (!printL(0) && (logFrequency < 10))
            logFrequency = 10;
          divingEvalAlgPtr->setOptionLogPrintFrequency(logFrequency);
        }
        else
        {
          int logFrequency = evalAlgParams.logPrintFrequency();
          if (!printL(0) && (logFrequency < 10))
            logFrequency = 10;
          divingEvalAlgPtr->setOptionLogPrintFrequency(logFrequency);
        }
      }
      else /// evalAlgParams is not active
      {
        if (runColGen)
          divingEvalAlgPtr->setOptionMaxNbOfCgIterations(param().MaxNbOfCgIteDuringRh());
        else
          divingEvalAlgPtr->setOptionMaxNbOfCgIterations(0);
        int logFrequency = param().ColGenLogFrequency();
        if (!printL(0) && (logFrequency < 10))
          logFrequency = 10;
        divingEvalAlgPtr->setOptionLogPrintFrequency(logFrequency);
      }

      return divingEvalAlgPtr;
}

void DivingHeuristic::replaceNodeInNonProperStrongDiving(Node **nodePtr)
{
  std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
  int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
  Node * newNodePtr = new Node(newNodeRef, *nodePtr, tmpLocalNodeBrConstrList, NULL);

  /// we add information about the number of needed proper columns to the evalAlgInfo
  ColGenEvalInfo * colGenEvalInfoPtr = dynamic_cast<ColGenEvalInfo *>(newNodePtr->nodeEvalInfoPtr());
  bapcodInit().require(colGenEvalInfoPtr != NULL,
                       "BaPCod error: nodeEvalInfo in diving heuristic is not of type ColGenEvalInfo.");
  DivingEvalInfo * divingEvalInfoPtr = new DivingEvalInfo(*colGenEvalInfoPtr, param().StrongDivingCandidatesNumber);
  newNodePtr->removeNodeEvalInfoAssociation();
  newNodePtr->associateNodeEvalInfoPtr(divingEvalInfoPtr);

  /// we copy dive algorithm information
  newNodePtr->associateGenChildNodesInfoPtr((*nodePtr)->genChildNodesInfoPtr());

  delete (*nodePtr);
  *nodePtr = newNodePtr;
}

bool DivingHeuristic::prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder)
{
  bool nodeShouldBeTreated = true;

  DiveInfo * diveInfoPtr = dynamic_cast<DiveInfo *>(nodePtr->genChildNodesInfoPtr());
  bapcodInit().require(diveInfoPtr != NULL,
                       "BaPCod error: genChildNodesInfoPtr for DiveAlgorithm is not of type DiveInfo.");

  if (nodePtr->isToBePruned(_currentBaPNodePtr->nodeIncIpPrimalBound()))
    {
      nodePtr->prunedAtBeginningOfTreatNode(true);
      return nodeShouldBeTreated = false;
    }

  /// we just create the node which runs the restricted master without any evaluation algorithm
  if (diveInfoPtr->performRestrictedMaster == 2)
  {
    nodePtr->setSolved(true);
  }

  if(!nodePtr->solved())
    {
      nodePtr->setPreprocessor(new Algorithm4PreprocessingInDive(_masterCommons.problemList()));

      /// set problem setup algorithm
      if ((nodePtr->probSetupInfoPtr()->treatOrderId == globalNodesTreatOrder)
          && !nodePtr->probSetupInfoPtr()->fullSetupIsObligatory)
        /// the problem was not altered after solving the parent node, so we do not do problem setup
        /// (partial solution fixing will taken into account in the preprocessing)
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupOfNode(_masterCommons.masterCommons4ProblemSetup()));
      else
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup()));

       /// set problem set down algorithm and evalAlg
        bool runColGen = true;
        if ( !param().runColGenAfterFixingPureMastVarInDiving() && (nodePtr->localFixedSolution() != NULL)
            && nodePtr->localFixedSolution()->solVarValMap().begin()->first->isTypeOf(VcId::InstMasterVarMask))
        {
           runColGen = false;
        }
        nodePtr->setEvalAlg(createEvaluationAlgorithm(runColGen));

      if ((diveInfoPtr->performRestrictedMaster == 0)
          && (((diveInfoPtr->remainingDepth > 0) && (diveInfoPtr->tabuVariables().size() < diveInfoPtr->maxTabuSize))
              || (param().StrongDivingCandidatesNumber() > 1)
              || param().DivingHeurPreprocessBeforeChoosingVar()
              || _problemPtr->divingFixingFunctorDefined()))
        /// we do full set down, if we can expect that there will be more than one child node,
        /// or if the candidates are testes for infeasibility by the preprocessing
        /// or if diving fixing functor is defined
        /// (for it, we need the current partial solution stored in the full problem set down)
        nodePtr->setProblemSetDownAlgorithm(new ProblemFullSetDownAlgorithm(
                                            _masterCommons.masterCommons4ProblemSetup()));
      else
        nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(
                                            _masterCommons.masterCommons4ProblemSetup()));
    }

  /// set the child nodes generation algorithm
  if (diveInfoPtr->performRestrictedMaster > 0)
  {
    RestrictedMasterIpHeuristic * restrictedMasterIpHeuristicPtr =
                                            new RestrictedMasterIpHeuristic(_problemPtr, _masterCommons);
    restrictedMasterIpHeuristicPtr->setOptionActivateAllColumns(param().ActivateAllColumnsForRestrictedMasterIpHeur());
    nodePtr->addPrimalHeuristic(restrictedMasterIpHeuristicPtr);

    nodePtr->setGenChildNodesAlgorithm(NULL);
  }
  else
  {
    nodePtr->setGenChildNodesAlgorithm(new DiveAlgorithm(_masterCommons.masterCommons4GenChildNodes(),
                                                         _colSelectionCriteria));
  }

  return nodeShouldBeTreated;
}

void DivingHeuristic::printDivingNodeInformation(const Node * nodePtr, const int diveNumber)
{
  if (nodePtr->localFixedSolution() == NULL)
      return;

  DiveInfo * diveInfoPtr = dynamic_cast<DiveInfo *>(nodePtr->genChildNodesInfoPtr());
  bapcodInit().require(diveInfoPtr != NULL,
                       "BaPCod error: genChildNodesInfoPtr for DivingHeuristic is not of type DiveInfo.");

  if (printL(-1))
  {
    std::cout << "---- Diving heuristic node with dive number = " << diveNumber
              << ", fix depth = " << (_maxDepth - diveInfoPtr->remainingDepth)
              << " and tabu list size = " << diveInfoPtr->tabuVariables().size() << std::endl;
    nodePtr->printFixedSolution(std::cout, printL(0));
  }
}


/// Implementation of the diving heuristic is similar to the implementation of
/// the branch-and-price algorithm (which is in MasterConf for now),
/// so we might think of creation a generic SearchTreeAlgorithm class
bool DivingHeuristic::runDiving(int & globalTreatOrder, Node * rootNodePtr)
{
  bool solutionImproved = false;
  int diveNumber = 0;
  int prevNodeDepth = 0;

  std::priority_queue<Node *, std::vector<Node *>, LargestNodeDepth> searchTree;

  Node * curNodePtr = rootNodePtr;
  searchTree.push(curNodePtr);
  _generatedNodes.push_back(curNodePtr);

  while (!searchTree.empty())
    {
      curNodePtr = searchTree.top();
      searchTree.pop();

      if (curNodePtr->depth() <= prevNodeDepth)
      {
        if (solutionImproved && _stopAfterFirstSolutionFound)
          break;
        diveNumber += 1;
      }
      prevNodeDepth = curNodePtr->depth();


      if (curNodePtr->solved() && !param().GenerateProperColumns())
        {
           /// "non-proper strong diving" case
           printDivingNodeInformation(curNodePtr, diveNumber);
           replaceNodeInNonProperStrongDiving(&curNodePtr);
        }

      if (prepareNodeForTreatment(curNodePtr, globalTreatOrder))
        {
          if (curNodePtr->localFixedSolution() != NULL)
          {
            Variable * firstVarInFixedSol = curNodePtr->localFixedSolution()->solVarValMap().begin()->first;

            bool previousDiveWasOnPureVar = false;
            BcVar pureVar(NULL);
            BcSolution spSolution(NULL);

            if (firstVarInFixedSol->isTypeOf(VcId::MastColumnMask))
            {
              MastColumn * colPtr = static_cast<MastColumn *>(firstVarInFixedSol);
              spSolution = BcSolution(colPtr->spSol());
            }
            else
            {
              previousDiveWasOnPureVar = true;
              pureVar = BcVar(static_cast<InstanciatedVar*>(firstVarInFixedSol));
            }
          }

          printDivingNodeInformation(curNodePtr, diveNumber);

          Bound primalBound(_currentBaPNodePtr->nodeIncIpPrimalBound());
          if (!curNodePtr->treat(globalTreatOrder, primalBound))
            {
              if (printL(0))
                std::cout << "ERROR: Branch-and-Price is interrupted" << std::endl;
              break;
            }

          /// the output of the treated node are the generated child nodes
          /// and possibly the updated bounds and the updated solution

          for (std::list<Node *>::const_iterator nodeIt =
              curNodePtr->sons().begin(); nodeIt != curNodePtr->sons().end();
              ++nodeIt)
            {
              searchTree.push(*nodeIt);
              _generatedNodes.push_back(*nodeIt);
            }

          /// nodePtr->nodeIncIpDualBound() is not used here,
          /// as the dual bound of a heuristic node is not a valid global bound

          if (curNodePtr->primalBoundIsUpdated())
            {
              solutionImproved = true;
              _currentBaPNodePtr->updateNodeIncPrimalSolution(curNodePtr->nodeIncIpPrimalSolPtr());
              if ((diveNumber > 0) && _stopAfterFirstSolutionFound)
                break;
            }
        }
    }
    return solutionImproved;
}

void DivingHeuristic::runBody(int & globalTreatOrder)
{
  std::list<BranchingConstrBaseType *> localNodeBrConstrList; /// this should be empty

  /// we make a "clone" of the node which launched the heuristic
  int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
  Node * rootNodePtr = new Node(newNodeRef, _currentBaPNodePtr, localNodeBrConstrList, NULL);

  /// We need to re-run the evaluation algorithm, as Algorithm4DivingEval::eval() should be called
  /// now Algorithm4DivingEval is the only evaluation algorithm used in the col. gen. based diving
  ColGenEvalInfo * colGenEvalInfoPtr = dynamic_cast<ColGenEvalInfo *>(rootNodePtr->nodeEvalInfoPtr());
  bapcodInit().require(colGenEvalInfoPtr != NULL,
                       "BaPCod error: nodeEvalInfo in DivingHeuristic is not of type ColGenEvalInfo.");

  int numberOfNeededProperColumns = _maxDiscrepancy + 1;
  if (param().StrongDivingCandidatesNumber() > numberOfNeededProperColumns)
    numberOfNeededProperColumns = param().StrongDivingCandidatesNumber();
  int enumerationWithFalseGapMode = 0;
  if (param().RCSPmaxNumOfLabelsInHeurEnumeration() > 0)
      enumerationWithFalseGapMode = (_runRestrMasterAfterFalseGapEnumeration ? 1 : 2);
      //enumerationWithFalseGapMode = 2;
  bool runCutGeneration = false; /// cut generation is skipped, we just need to re-run column generation
  DivingEvalInfo * divingEvalInfoPtr = new DivingEvalInfo(*colGenEvalInfoPtr, numberOfNeededProperColumns,
                                                          runCutGeneration, enumerationWithFalseGapMode);
  rootNodePtr->removeNodeEvalInfoAssociation();
  rootNodePtr->associateNodeEvalInfoPtr(divingEvalInfoPtr);

  int performRestrictedMaster = (_runRestrMasterAfterFalseGapEnumeration ? 1 : 0);
  rootNodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_maxDepth, _maxDiscrepancy, performRestrictedMaster));

  if (printL(0))
    std::cout << "------------------------------------------------" << std::endl;
  if (printL(-1))
  {
    if (performRestrictedMaster)
      std::cout << "-- Enumeration for restr. mast. heur. started --" << std::endl;
    else if (enumerationWithFalseGapMode > 0)
      std::cout << "---- Diving heur. with enumeraiton started -----" << std::endl;
    else
      std::cout << "----------- Diving heuristic started -----------" << std::endl;
  }
  if (printL(0))
    std::cout << "------------------------------------------------" << std::endl;

  runDiving(globalTreatOrder, rootNodePtr);
}

Solution * LocalSearchHeuristic::fixPartialSolution()
{
  Solution * incSolutionPtr = _currentBaPNodePtr->nodeIncIpPrimalSolPtr();
  if (incSolutionPtr == NULL)
    incSolutionPtr = _currentBaPNodePtr->probConfPtr()->primalSolutionPtr();
  if (incSolutionPtr == NULL)
    return NULL;
  const VarPtr2DoubleMap & varValMap = incSolutionPtr->solVarValMap();

  std::vector<std::pair<Variable *, Double> > referenceSolutionVector;
  for (VarPtr2DoubleMap::const_iterator it = varValMap.begin(); it != varValMap.end(); ++it)
    referenceSolutionVector.push_back(std::make_pair(it->first, it->second));
  int referenceSolutionSize = referenceSolutionVector.size();

  std::vector<bool> fixedVariables(referenceSolutionSize, false);
  int numberOfVariablesToFix(Iceil(_fixedVarsRatio * referenceSolutionSize));
  if (numberOfVariablesToFix >= referenceSolutionSize)
    return NULL;

  int numberOfFixedVariables = 0;
  Solution * solPtr = new Solution();
  while (numberOfFixedVariables < numberOfVariablesToFix)
    {
      int randomIndex = rand() % referenceSolutionSize;
      if (!fixedVariables[randomIndex])
        {
          fixedVariables[randomIndex] = true;
          numberOfFixedVariables += 1;
          solPtr->includeVar(referenceSolutionVector[randomIndex].first,
              referenceSolutionVector[randomIndex].second, false);
        }
    }

  return solPtr;
}

void LocalSearchHeuristic::runBody(int & globalTreatOrder)
{

  if (printL(0))
  {
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "-------Local search heuristic started-------" << std::endl;
  }

  /// first dive to generation the first reference solution
  std::list<BranchingConstrBaseType *> localNodeBrConstrList; /// this should be empty
  int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
  Node * nodePtr = new Node(newNodeRef, _currentBaPNodePtr, localNodeBrConstrList, NULL);

  ColGenEvalInfo * colGenEvalInfoPtr = dynamic_cast<ColGenEvalInfo *>(nodePtr->nodeEvalInfoPtr());
  bool runCutGeneration = false; /// cut generation is skipped, we just need to re-run column generation
  DivingEvalInfo * divingEvalInfoPtr = new DivingEvalInfo(*colGenEvalInfoPtr, 0, runCutGeneration, false);
  nodePtr->associateNodeEvalInfoPtr(divingEvalInfoPtr);

  nodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_maxDepth = 10000, _maxDiscrepancy = 1));

  _stopAfterFirstSolutionFound = true;
  runDiving(globalTreatOrder, nodePtr);
  _stopAfterFirstSolutionFound = false;

  int iterationNumber = 0;
  while (!_currentBaPNodePtr->isConquered() && (iterationNumber < _maxIterationsNumber))
    {
      iterationNumber += 1;
      if (printL(0))
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "Local search heuristic iteration " << iterationNumber << std::endl;
      }

      Solution * fixedSolPtr = fixPartialSolution();
      if (fixedSolPtr == NULL)
        break;
      int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
      nodePtr = new Node(newNodeRef, _currentBaPNodePtr, localNodeBrConstrList, fixedSolPtr);

      ColGenEvalInfo * colGenEvalInfoPtr = dynamic_cast<ColGenEvalInfo *>(nodePtr->nodeEvalInfoPtr());
      bool runCutGeneration = false; /// cut generation is skipped, we just need to re-run column generation
      DivingEvalInfo * divingEvalInfoPtr = new DivingEvalInfo(*colGenEvalInfoPtr, 0, runCutGeneration, false);
      nodePtr->associateNodeEvalInfoPtr(divingEvalInfoPtr);

      nodePtr->associateGenChildNodesInfoPtr(new DiveInfo(_maxDepth = 0, _maxDiscrepancy = 0));
      runDiving(globalTreatOrder, nodePtr);
    }
}

void LocalSearchHeuristic::setOptionMaxIterationsNumber(const int value)
{
  _maxIterationsNumber = value;
}

void LocalSearchHeuristic::setOptionFixedVarsRatio(const double value)
{
  _fixedVarsRatio = value;
}