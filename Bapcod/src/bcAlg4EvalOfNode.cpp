/**
 *
 * This file bcAlg4EvalOfNode.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4EvalByMip.hpp"
#include "bcAlg4Master.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcNetworkBasedCuts.hpp"
#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#ifdef USE_NON_PUBLIC_CUTS
#include "bcNonPublicCuts.hpp"
#endif
#endif

//#define DEBUG_MISSING_SOLUTION

Alg4EvalOfNode::Alg4EvalOfNode(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons) :
  Alg4Master(probPtr), _masterCommons(masterCommons), _lastPriorityLevel(BapcodInfinity),
  _pricingSolverCutsMessageId(PricingSolverCutsMessage::noMessage), _addCutToMasterFirstCall(false),
  _solIsInteger(false), _needBeConqueredIfSolIsInteger(true), _nonStabArtVarPtrList()
{
}

void Alg4EvalOfNode::addCutsToProblem(char const & flag, const double & cutSepTime,
                                      std::map<std::string, int> & numCutGeneratedPerType,
                                      std::multiset<InstanciatedConstr *,
                                                    CutSeparationPriorityComp> & generatedCutConstrSet)
{
  if (generatedCutConstrSet.empty())
    return;

  Time startAc;
  bapcodInit().statistics().incrCounter("bcCountCutInMaster", (long int) generatedCutConstrSet.size());

  int nbOfCutsKept = 0;
  if (printL(5))
    std::cout << "MasterConf::addCutToMaster(): generated cut nb = " << generatedCutConstrSet.size() << std::endl;

  std::multiset<InstanciatedConstr *, CutSeparationPriorityComp>::iterator itCpt;
  for (itCpt = generatedCutConstrSet.begin(); (itCpt != generatedCutConstrSet.end())
                                              && (nbOfCutsKept < param().MaxNbOfCutGeneratedAtEachIter()); ++itCpt)
  {
    /// cast the cut without adding it to the problem
    InstanciatedConstr * instConstrPtr = _masterCommons.castAndAddConstraint(*itCpt, false);
    /// cuts do not participate in preprocessing, as for the moment there is no mechanism to keep their slacks
    /// up to date
    instConstrPtr->toBeUsedInPreprocessing(false);
    /// add the cut to the master problem, its local artificial variables will be also added
    _probPtr->addConstraintToProblem(instConstrPtr);

    InstMasterConstr * instMastConstrPtr = static_cast<InstMasterConstr *>(instConstrPtr);
    if (_probPtr->curNodePtr() != NULL)
      instMastConstrPtr->treatOrderId(_probPtr->curNodePtr()->treatOrder());

    /// activate the cut and set it to the formulation
    /// its local artificial variables will be also activated and set to the formulation
    instConstrPtr->activateConstraint(true);

    if (instConstrPtr->posLocalArtVarPtr() != NULL)
      _nonStabArtVarPtrList.push_back(instConstrPtr->posLocalArtVarPtr());
    if (instConstrPtr->negLocalArtVarPtr() != NULL)
      _nonStabArtVarPtrList.push_back(instConstrPtr->negLocalArtVarPtr());

    if ((_masterCommons.debugSolution() != NULL) && _currentNodePtr->debugSolutionAtThisNode())
    {
        Double lhsValue = instConstrPtr->computeLhs(_masterCommons.debugSolution()->solVarValMap());
        if (((instConstrPtr->sense() == 'G') && (lhsValue + Double::precision < instConstrPtr->costrhs()))
            || ((instConstrPtr->sense() == 'L') && (lhsValue - Double::precision > instConstrPtr->costrhs())))
        {
            std::cerr << "BaPCod WARNING : new cut " << instConstrPtr->name() << " is violated by the debug solution "
                      << std::endl;
            if (printL(-1))
            {
                std::cout << "BaPCod WARNING : new cut " << instConstrPtr->name()
                          << " is violated by the debug solution "
                          << std::endl;
                InstMasterConstr *imConstrPtr = static_cast<InstMasterConstr *>(instConstrPtr);
                std::cout << "New cut " << instConstrPtr->name() << " : ";
                for (MapSubProbVariablePtr2Double::const_iterator mapIt = imConstrPtr->subProbVarMember2coefMap().begin();
                     mapIt != imConstrPtr->subProbVarMember2coefMap().end(); ++mapIt) {
                    std::cout << "+" << mapIt->second << "*" << mapIt->first->name();
                }
                for (ConstVarConstrPtr2Double::const_iterator mapIt = imConstrPtr->member2coefMap().begin();
                     mapIt != imConstrPtr->member2coefMap().end(); ++mapIt) {
                    if (!mapIt->first->isTypeOf(VcId::MastColumnMask)) {
                        std::cout << "+" << mapIt->second << "*" << mapIt->first->name();
                    }
                }
                std::cout << (instConstrPtr->sense() == 'G' ? " >= " : " <= ") << instConstrPtr->costrhs()
                          << "(debug sol. LHS is " << lhsValue << ")" << std::endl;
            }
        }
    }

    ++nbOfCutsKept;
  }

  /// add the cuts and their local artificial variables to the formulation
  _probPtr->addConstrInForm();
  _probPtr->addVarInForm();

  /// we delete cuts which were not kept
  for (std::multiset<InstanciatedConstr *, CutSeparationPriorityComp>::iterator itCpt2 = itCpt;
       itCpt2 != generatedCutConstrSet.end(); ++itCpt2)
    delete *itCpt2;
  generatedCutConstrSet.erase(itCpt, generatedCutConstrSet.end());

  /// print information about the added cuts
  if (printL(-1))
  {
    std::cout << "----- Add " << (flag == 'F' ? "fac." : "core") << " cuts :";
    for (std::map<std::string, int>::iterator mapIt = numCutGeneratedPerType.begin();
         mapIt != numCutGeneratedPerType.end(); ++mapIt)
      std::cout << " " << mapIt->first << "(" << mapIt->second << ")";
  }
  Double maxViolation(0.0);
  Double totalViolation(0.0);
  int numWithZeroViolation = 0;
  for (itCpt = generatedCutConstrSet.begin(); itCpt != generatedCutConstrSet.end(); ++itCpt)
  {
    (*itCpt)->computeViolation((*itCpt)->computeLhs(_probPtr->inPrimalLpSol()));
    if (maxViolation < (*itCpt)->violation())
      maxViolation = (*itCpt)->violation();
    totalViolation += (*itCpt)->violation();
    if ((*itCpt)->violation().isZero())
      numWithZeroViolation += 1;
  }
  if (printL(-1))
  {
    std::cout << ", max.viol = " << maxViolation << ", aver.viol = "
              << (double) totalViolation / generatedCutConstrSet.size();
    if (numWithZeroViolation > 0)
      std::cout << ", zero viol = " << numWithZeroViolation;
  }

  double cutAddTime = startAc.getElapsedTime_dbl();
  bapcodInit().statistics().incrTimer("bcTimeAddCutToMaster", cutAddTime);

  if (printL(-1))
    std::cout << ", sep/add took " << cutSepTime/100 << "/" << cutAddTime/100 << " sec." << " -----" << std::endl;
}

void Alg4EvalOfNode::recordCutSeparationStats(const std::string & cutClass,
                                              const std::multiset<InstanciatedConstr *,
                                                                  CutSeparationPriorityComp> & generatedCutConstrSet,
                                              int numCutsBefore, double sepTime)
{
    int numNewCuts = generatedCutConstrSet.size() - numCutsBefore;

    bapcodInit().statistics().incrCounter(std::string("bcCountCut") + cutClass, numNewCuts);
    bapcodInit().statistics().incrTimer(std::string("bcTimeCutSep") + cutClass, sepTime);
#ifdef BCP_RCSP_IS_FOUND
    if (cutClass == "R1C")
    {
        std::vector<int> numPackingCuts(RANK1_CUTS_MAX_NUMBER_OF_ROWS + 1, 0);
        std::vector<int> numCoveringCuts(RANK1_CUTS_MAX_NUMBER_OF_ROWS + 1, 0);
        for (const auto * cutPtr : generatedCutConstrSet)
        {
            if (cutPtr->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
            {
                const bcp_rcsp::RankOneCut * rank1CutPtr = static_cast<const LimMemRankOneCut *>(cutPtr)->rcspCutPtr();
                if (rank1CutPtr->cutClass == SET_PACKING_CUT)
                    numPackingCuts[rank1CutPtr->numRows] += 1;
                else
                    numCoveringCuts[rank1CutPtr->numRows] += 1;
            }
        }
        for (int numRows = 1; numRows <= RANK1_CUTS_MAX_NUMBER_OF_ROWS; ++numRows)
        {
            if (numPackingCuts[numRows] > 0)
                bapcodInit().statistics().incrCounter(std::string("bcCountCut") + std::to_string(numRows)
                                                      + "rowPackR1C", numPackingCuts[numRows]);
            if (numCoveringCuts[numRows] > 0)
                bapcodInit().statistics().incrCounter(std::string("bcCountCut") + std::to_string(numRows)
                                                      + "rowCovR1C", numCoveringCuts[numRows]);
        }
    }
    if (cutClass == "DCC")
    {
        int numDepotCapCutsWithOneY = 0;
        int numDepotCapCutsWithTwoYs = 0;
        for (const auto * cutPtr : generatedCutConstrSet)
        {
            if (cutPtr->genConstrPtr()->defaultName() != "DCC")
                continue;

            int numYmembers = 0;
            for (const auto &pair: cutPtr->member2coefMap())
                if (pair.first->genericName() == "Y")
                    numYmembers += 1;
            if (numYmembers == 1)
                numDepotCapCutsWithOneY += 1;
            else
                numDepotCapCutsWithTwoYs += 1;
        }
        if (numDepotCapCutsWithTwoYs > 0)
            bapcodInit().statistics().incrCounter(std::string("bcCountCutDCCwithTwoYs"), numDepotCapCutsWithTwoYs);
        if (numDepotCapCutsWithOneY > 0)
            bapcodInit().statistics().incrCounter(std::string("bcCountCutDCCwithOneY"), numDepotCapCutsWithOneY);
    }
    if (cutClass == "RCK")
    {
        int maxOneKvalue = std::abs(param().RCSPresConsKnapsackCutsMode());
        std::vector<int> numRLKCutsByOneK(maxOneKvalue + 1, 0);
        for (const auto * cutPtr : generatedCutConstrSet)
        {
            if (cutPtr->isTypeOf(VcId::ResConsKnapsackCutConstrMask))
            {
                const bcp_rcsp::RouteLoadKnapsackCut * rlkCutPtr
                        = static_cast<const ResConsKnapsackCut *>(cutPtr)->rcspCutPtr();
                numRLKCutsByOneK[rlkCutPtr->oneKvalue] += 1;
            }
        }
        for (int oneKvalue = 0; oneKvalue <= maxOneKvalue; ++oneKvalue)
            if (numRLKCutsByOneK[oneKvalue] > 0)
            {
                if (oneKvalue == 0)
                    bapcodInit().statistics().incrCounter(std::string("bcCountCutRLKCRounding"),
                                                          numRLKCutsByOneK[0]);
                else
                    bapcodInit().statistics().incrCounter(std::string("bcCountCutRLKCsOne")
                                                          + std::to_string(oneKvalue), numRLKCutsByOneK[oneKvalue]);
            }
    }
#endif
}

void Alg4EvalOfNode
     ::removeColumnsNotSatisfyingColGenSpRelaxation(const std::set<ColGenSpConf *> & subprobsWithImprovedRelaxation)
{
    if (subprobsWithImprovedRelaxation.empty())
        return;

    std::list<Variable *> varsToRemoveFromForm;
    for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
         varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); )
    {
        if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
        {
            MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
            ++varPtrIt;             /// this should be done before making the column unsuitable

            ColGenSpConf * cgSpConfPtr = colPtr->cgSpConfPtr();
            if (subprobsWithImprovedRelaxation.count(cgSpConfPtr)
                && !cgSpConfPtr->probPtr()->solSatisfiesCurrentSpRelaxation(colPtr->spSol()))
            {
                _probPtr->probVarSet().insert(colPtr, VcIndexStatus::Unsuitable);
                colPtr->desactivate();
                varsToRemoveFromForm.push_back(colPtr);
            }

        }
        else
            ++varPtrIt;
    }
    /// now we do not need the solution anymore, as cuts have been generated already
    /// we can now clean up columns not satisfying improved subproblems relaxation
    if (printL(0))
        std::cout << "Removed " << varsToRemoveFromForm.size() << " columns not satisfying improved "
                  << "subproblem relaxation" << std::endl;
    _probPtr->resetSolution('d');         /// sol. may involve columns which will now be removed from the formulation
    _probPtr->delVarsSimplyInForm(varsToRemoveFromForm);
    _probPtr->removeUnusedDynamicVarsFromMemory();
}


/// returns true if column generation should be repeated
bool Alg4EvalOfNode
     ::runColGenSpRelaxationImprovement(bool masterConverged, std::set<ColGenSpConf *> & subprobsWithImprovedRelaxation)
{
    subprobsWithImprovedRelaxation.clear();

    /// in the current implementation subprob. restriction level is embedded into subprob. relaxation level
    if (!progStatus().doRun() || !_solIsMasterLpFeasible)
        return false;

    std::map<ColGenSpConf *, std::vector<MastColumn *> > colsInSolutionMap;
    std::map<ColGenSpConf *, std::vector<MastColumn *> >::iterator mapIt;
    for (VarPtrSet::const_iterator varPtrIt = _probPtr->inPrimalLpSol().begin();
         varPtrIt != _probPtr->inPrimalLpSol().end(); ++varPtrIt)
        if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
        {
            MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
            mapIt = colsInSolutionMap.find(colPtr->cgSpConfPtr());
            if (mapIt != colsInSolutionMap.end())
            {
                mapIt->second.push_back(colPtr);
            }
            else
            {
                colsInSolutionMap[colPtr->cgSpConfPtr()] = std::vector<MastColumn *>();
                colsInSolutionMap[colPtr->cgSpConfPtr()].push_back(colPtr);
            }
        }

    for (mapIt = colsInSolutionMap.begin(); mapIt != colsInSolutionMap.end(); ++mapIt)
    {
        if (!mapIt->first->enumeratedStatus()
            && mapIt->first->probPtr()->improveCurrentSpRelaxation(mapIt->second, masterConverged))
            subprobsWithImprovedRelaxation.insert(mapIt->first);
    }

    return (!subprobsWithImprovedRelaxation.empty());
}

bool Alg4EvalOfNode::addCutToMaster(char flag, bool masterConverged)
{
    if (printL(5))
        std::cout << "MasterConf::addCutToMaster() flag = " << flag << "  candidateCutGenericConstr().size() = "
                  << _masterCommons.candidateCutGenericConstr().size() << std::endl;

    if (flag == 'F')
      {
        if (_addCutToMasterFirstCall)
          {
            _addCutToMasterFirstCall = false;

            bool status(_masterCommons.pcDelayedConstrPtrList().size() > 0);

            for (std::list<Constraint *>::iterator it = _masterCommons.pcDelayedConstrPtrList().begin();
                 it != _masterCommons.pcDelayedConstrPtrList().end(); ++it)
              {
                (*it)->kind('E');
                if (printL(5))
                    std::cout << "MasterConf::addCutToMaster() add delayed constraint " << (*it);

                _masterCommons.pcConstrPtrList().push_back(*it);
              }

            _probPtr->addConstrSet(_masterCommons.pcDelayedConstrPtrList(), 1, 2);
            _masterCommons.pcDelayedConstrPtrList().clear();
            if (status)
                return (true);
          }
      }

    if (_masterCommons.candidateCutGenericConstr().empty() && param().ColGenSpRelaxationImprovementPriority() == 0.0)
      return false;

    Time startCp;

    std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> generatedCutConstrSet;

    std::map<std::string, int> numCutGeneratedPerType; /// for printing purposes

    std::vector<std::tuple<double, std::string, GenericCutConstr *> > separators;
    for (auto * cutGenConstrPtr : _masterCommons.candidateCutGenericConstr())
        separators.emplace_back(cutGenConstrPtr->dynamicPriorityLevel(_currentNodePtr->isRoot()),
                                cutGenConstrPtr->defaultName(), cutGenConstrPtr);
    if (param().ColGenSpRelaxationImprovementPriority() > 0.0)
        separators.emplace_back(param().ColGenSpRelaxationImprovementPriority(), "", nullptr);

    std::sort(separators.rbegin(), separators.rend());

    /// generic cut constraints
    auto lastPriorityLevel(std::get<0>(separators.front()));
    if (!param().CutDynamicPriorityLevel() && (lastPriorityLevel > _lastPriorityLevel))
      lastPriorityLevel = _lastPriorityLevel;

    bool relaxationImproved = false;
    std::set<ColGenSpConf *> subprobsWithImprovedRelaxation;

    /// we check if we can branch on candidates with larger priority than the current priority level of cut separation
    /// if yes, we prefer to do branching and not generate cuts
    /// deactivated for the moment
//    BranchGeneratorsSet generatedBrConstrGeneratorSet;
//    for (std::set<DynamicGenericConstr *, DynamicGenConstrSort>::const_iterator dgcIt =
//         _masterCommons.candidateBranchingGenericConstr().begin();
//         dgcIt != _masterCommons.candidateBranchingGenericConstr().end(); ++dgcIt)
//      if ((*dgcIt)->priorityLevel() > lastPriorityLevel)
//        {
//          (*dgcIt)->branchingSeparationFindCandidates(_currentNodePtr->listOfFractMastCol(),
//                                                      generatedBrConstrGeneratorSet);
//        }
//    bool branchCandidatesFound = !generatedBrConstrGeneratorSet.empty();
//    for (BranchGeneratorsSet::iterator setIt = generatedBrConstrGeneratorSet.begin();
//         setIt != generatedBrConstrGeneratorSet.end(); ++setIt)
//      delete *setIt;
//    generatedBrConstrGeneratorSet.clear();

//    if (!branchCandidatesFound)
    if (true)
    {
        for (auto & sepTuple : separators)
        {
            auto currPriorityLevel = std::get<0>(sepTuple);
            auto cutGenConstrPtr = std::get<2>(sepTuple);

            if (flag == 'C')
            {
                if (cutGenConstrPtr != nullptr && cutGenConstrPtr->type() == 'F')
                    continue;
                /// we still need to launch the sp. relaxation improvement, as
                /// the columns in the integer solution may not be proper
//                if (cutGenConstrPtr == nullptr)
//                    continue;
            }

            if (currPriorityLevel < lastPriorityLevel)
            {
                if (!relaxationImproved)
                    lastPriorityLevel = currPriorityLevel;
                else
                    break;
            }

            if (cutGenConstrPtr != nullptr)
            {
                double thisSepTime = startCp.getElapsedTime_dbl();
                int cutsBefore = (int)generatedCutConstrSet.size();
                cutGenConstrPtr->cutSeparationRoutine(_probPtr->inPrimalLpSol(), generatedCutConstrSet);
                int numNewCuts = (int)generatedCutConstrSet.size() - cutsBefore;
                thisSepTime = startCp.getElapsedTime_dbl() - thisSepTime;
                if (numNewCuts > 0)
                {
                    relaxationImproved = true;
                    numCutGeneratedPerType[cutGenConstrPtr->defaultName()] = numNewCuts;
                    recordCutSeparationStats(cutGenConstrPtr->defaultName(), generatedCutConstrSet, cutsBefore,
                                             thisSepTime);
                }
            }
            else
            {
                /// sp. relaxation improvement (strengthening the subproblems' relaxation)
                if (runColGenSpRelaxationImprovement(masterConverged, subprobsWithImprovedRelaxation))
                    relaxationImproved = true;
            }

            if (printL(5))
                std::cout << "MasterConf::addCutToMaster() new  candidateCutGenericConstr().size() = "
                          << _masterCommons.candidateCutGenericConstr().size() << std::endl;
        }
    }

    if (flag == 'F')
        _lastPriorityLevel = lastPriorityLevel;

    double cutSepTime = startCp.getElapsedTime_dbl();
    bapcodInit().statistics().incrTimer("bcTimeCutSeparation", cutSepTime);

    if (relaxationImproved)
    {
        if (!generatedCutConstrSet.empty())
            addCutsToProblem(flag, cutSepTime, numCutGeneratedPerType, generatedCutConstrSet);
        if (printL(-1) && !subprobsWithImprovedRelaxation.empty())
            std::cout << "----- Relaxation was improved for " << subprobsWithImprovedRelaxation.size()
                      << " subproblems -----" << std::endl;
    }
    else if (printL(-1) && flag == 'F')
    {
        std::cout << "----- Add fac. cuts : sep. took " << cutSepTime/100 << " sec.";
//        if (branchCandidatesFound)
//          std::cout << "----- and stopped, larger priority branching candidate found ";
//        else
        std::cout << "----- no cuts found ";
        std::cout << " -----" << std::endl;
    }

    /// this should be done after the separation of cuts
    /// (deactivation of columns participating in the solution affects cut separation, need to check why)
    removeColumnsNotSatisfyingColGenSpRelaxation(subprobsWithImprovedRelaxation);

    if ((flag == 'C') && relaxationImproved)
      _solIsInteger = false;

    return relaxationImproved;
}

bool Alg4EvalOfNode::checkIfCurSolIsInteger()
{
  int printlevel = 5;

  _solIsInteger = false;

  /// Test Integrality of aggregate Solution
  if (param().TestAggregateMasterSol4Integrality())
  {
    VarPtr2DoubleMap curAggregateMastSol;
    _solIsInteger = true; /// unless proved otherwise

    if (printL(printlevel))
      std::cout << "ColGenEvalAlg::checkIfCurSolIsInteger? checking inPrimalSol " << std::endl;

    for (VarPtrSet::const_iterator sPtr = _probPtr->inPrimalLpSol().begin();
         sPtr != _probPtr->inPrimalLpSol().end(); sPtr++)
    {
      if ((*sPtr)->isTypeOf(VcId::ArtificialVarMask) && !(*sPtr)->val().isZero())
      {
        if (printL(-1))
          std::cout << "BaPCod WARNING: ColGenEvalAlg::checkIfCurSolIsInteger() mast. sol. should not include art var "
                    << (*sPtr)->name() << " with value " << (*sPtr)->val() << " at this stage" << std::endl;
        std::cerr << "BaPCod WARNING: master solution should not include art var at this stage" << std::endl;
      }

      ValueRecord rec((*sPtr)->val(), param().BapCodIntegralityTolerance);

      if ( (((*sPtr)->type() == 'B') || ((*sPtr)->type() == 'I')) && rec._isFractional
          && !(*sPtr)->isTypeOf(VcId::MastColumnMask))
        _solIsInteger = false;

      /// Show current master LP solution
      if (printL(3))
        std::cout << "ColGenEvalAlg::checkIfCurSolIsInteger? MastCol[ " << (*sPtr)->name() << " ] = "
                  << (*sPtr)->val() << std::endl;
      (*sPtr)->fillAggregateSol(curAggregateMastSol, (*sPtr)->val());
    }

    if (printL(printlevel))
      std::cout << "ColGenEvalAlg::checkIfCurSolIsInteger? checking curAggregateMastSol " << std::endl;

    for (VarPtr2DoubleMap::const_iterator it = curAggregateMastSol.begin(); it != curAggregateMastSol.end(); ++it)
    {
      if (it->first->genVarConstrPtr()->priorityRule() == SelectionStrategy::NotConsideredForIntegralityCheck)
        continue;

      if ((it->first->type() == 'B') || (it->first->type() == 'I'))
      {
        Double intPart = Dfloor(it->second);
        if (printL(4))
          std::cout << "sol[" << it->first->name() << "] = " << it->second << ", round = " << intPart << std::endl;

        if (it->second != intPart)
        {
          _solIsInteger = false;
          break;
        }
      }
    }
  } // end if (param().TestProjectedMasterSol4Integrality())

  /// if observing the projected solution this has not been disproved
  if (_solIsInteger)
  {
    if (printL(1))
      std::cout << "ColGenSolver::checkIfCurSolIsInteger() current solution to the master is integer " << std::endl;
    return true;
  }

  /// Now we test integrality of the projection of the master LP solution
  if (printL(printlevel))
    std::cout << "ColGenSolver::checkIfCurSolIsInteger? checking inPrimalSol " << std::endl;

  bool thereAreFractionalMastColumns = false;
  bool thereAreFractionalPureMastVars = false;
  for (VarPtrSet::const_iterator colPt = _probPtr->inPrimalLpSol().begin();
       colPt != _probPtr->inPrimalLpSol().end(); colPt++)
  {
    /**
     * Consider nonzero solution: need High precision here, otherwise rounding error do accumulate
     */
    double veryHighPrecision = Double::precision * 1e-4; /// used to be VERYHIGHPRECISION parameter with value 1e-10
    if (((*colPt)->val() > veryHighPrecision) || ((*colPt)->val() < -veryHighPrecision))
    {
      if ((*colPt)->isTypeOf(VcId::ArtificialVarMask) && !(*colPt)->val().isZero())
      {
        if (printL(0))
          std::cout << "BaPCod info : ColGenEvalAlg::checkIfCurSolIsInteger() mast. sol. should not include art var "
                    << (*colPt)->name() << " with value " << (*colPt)->val() << " at this stage(1)" << std::endl;
        std::cerr << "BaPCod WARNING: master solution should not include art var at this stage(1)" << std::endl;
      }

      ValueRecord rec((*colPt)->val(), param().BapCodIntegralityTolerance);

      if ((((*colPt)->type() == 'B') || ((*colPt)->type() == 'I')) && rec._isFractional)
      {
        if ((*colPt)->isTypeOf(VcId::MastColumnMask))
          thereAreFractionalMastColumns = true;
        else
          thereAreFractionalPureMastVars = true;
      }

      if (printL(printlevel))
        std::cout << "ColGenSolver::checkIfCurSolIsInteger() founds that lexicographically sorted master var[ "
                  << (*colPt)->name() << " ] is fractional = " << rec._isFractional
                  << " and has value " << (*colPt)->val() << " and tmp value " << rec._lfracValue << std::endl
                  << " thereAreFractionalMastColumns " << thereAreFractionalMastColumns
                  << " thereAreFractionalPureMastVars " << thereAreFractionalPureMastVars << std::endl;
    }
  }

  if ((param().VerifyColsIntegralityInTestSolForIntegrality() && thereAreFractionalMastColumns)
      || thereAreFractionalPureMastVars)
  {
    _solIsInteger = false;
    return false;
  }

  /// Compute the projection of the master solution in the original variable space
  _solIsInteger = true;  /// unless shown otherwise

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    (*spcPt)->listOfFractMastColInColGenSp().clear();

  for (VarPtrSet::const_iterator sPtr = _probPtr->inPrimalLpSol().begin();
       sPtr != _probPtr->inPrimalLpSol().end(); sPtr++)
  {
    ValueRecord rec((*sPtr)->val(), param().BapCodIntegralityTolerance);
    if ((*sPtr)->isTypeOf(VcId::MastColumnMask) && rec._isFractional)
      (*sPtr)->cgSpConfPtr()->listOfFractMastColInColGenSp().push_back(*sPtr, rec);
  }

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
  {
    /**
     * Sort fract col in ILO (lexicographic order dictated by the set of branching constraints)
     */
    ComponentSequence curClassCompBoundSet(*spcPt);
    MasterColSolution nonSortedListOfMastCol((*spcPt)->listOfFractMastColInColGenSp());
    (*spcPt)->listOfFractMastColInColGenSp().clear();
    CompBoundSetGenBranchConstr::ILOsortMastColumn(_currentNodePtr->treeOfColClasses(),
                                                   nonSortedListOfMastCol, curClassCompBoundSet,
                                                   (*spcPt)->listOfFractMastColInColGenSp());

    if (printL(printlevel))
      std::cout << "ColGenSpConf " << (*spcPt)->name() << " has sorted list of frac col has size "
                << (*spcPt)->listOfFractMastColInColGenSp().size() << std::endl;
    /**
     * Test dissagregated solution (OVF solution)
     */
    std::map<int, VarPtr2DoubleMap> curDisaggregateSol;
    if (printL(printlevel))
      std::cout << "ColGenSolver::checkIfCurSolIsInteger(): fract columns in ILO order " << std::endl;

    Double cumVal(0);
    int tindex(0);
    MastColumn * mastColPtr(NULL);

    /**
     * _listOfFractMastCol includes integer variables up to value 1.999; columns are sorted in ILO
     */
    for (MasterColSolution::const_iterator colPt = (*spcPt)->listOfFractMastColInColGenSp().begin();
         colPt != (*spcPt)->listOfFractMastColInColGenSp().end(); colPt++)
    /**
     * Try to take convex combination of frac col to get an nteger solution
     */
    {
      mastColPtr = colPt->first;
      if (printL(printlevel))
      {
        std::cout << "col[ " << mastColPtr->name() << " ] has fractional part " << colPt->second._lfracValue
                  << " cumVal = " << cumVal << std::endl;
        mastColPtr->printColVector();
      }

      Double tmpColVal(colPt->second._lfracValue);
      Double lambdaR(0);
      while (tmpColVal > 0)
      {
        /// Slice is filled
        if (cumVal == (tindex + 1))
        {
          /// Check integrality in mapped solution
          for (VarPtr2DoubleMap::const_iterator it = curDisaggregateSol[tindex].begin();
               it != curDisaggregateSol[tindex].end(); ++it)
          {
            if (printL(printlevel))
              std::cout << "curDisaggregateSol[" << tindex << "][" << it->first->name() << "] = "
                        << it->second << " " << it->first->type() << std::endl;

            InstanciatedVar * instanciatedVarPtr = dynamic_cast<InstanciatedVar *> (it->first);
            if ((instanciatedVarPtr != NULL) && (instanciatedVarPtr->genVarConstrPtr()->priorityRule()
                                                 == SelectionStrategy::NotConsideredForIntegralityCheck))
              continue;

            if (((it->first->type() == 'B') || (it->first->type() == 'I'))
                && (it->second > Dfloor(it->second)))
            {
              if (printL(3))
                std::cout << "ColGenSolver::checkIfCurSolIsInteger():  detected fractional component "
                << "curDisaggregateSol[" << tindex << "][" << it->first->name()
                << "] = " << it->second << "  solIsInteger = false " << std::endl;
              _solIsInteger = false;
              return false;
            }
          }

          /// Goto next slice
          tindex++;
        }

        lambdaR = tmpColVal;
        if (cumVal + lambdaR > tindex + 1)
          lambdaR = (tindex + 1 - cumVal);

        mastColPtr->fillMapOfIntSpVar(curDisaggregateSol[tindex], lambdaR);
        tmpColVal -= lambdaR;
        cumVal += lambdaR;

        if (printL(5))
          std::cout << "ColGenSolver::checkIfCurSolIsInteger(): tindex = " << tindex << " map size "
                    << curDisaggregateSol[tindex].size() << " cumVal = " << cumVal << std::endl;
      }
    }

    /**
     * Check integrality in mapped solution, do not accept a fractional number of slice,
     * or else one could only pay a fraction of the setup cost associated with a column
     */
    if (cumVal > Dfloor(cumVal) + param().BapCodIntegralityTolerance)  /*not an integer number of slices */
    {
      /// Check integrality in mapped solution in last slice
      for (VarPtr2DoubleMap::const_iterator it = curDisaggregateSol[tindex].begin();
           it != curDisaggregateSol[tindex].end(); ++it)
      {
        if (printL(5))
          std::cout << "curDisaggregateSol[" << tindex << "][" << it->first->name() << "] = " << it->second
                    << std::endl;

        InstanciatedVar * instanciatedVarPtr = dynamic_cast<InstanciatedVar *> (it->first);
        if ((instanciatedVarPtr != NULL) && (instanciatedVarPtr->genVarConstrPtr()->priorityRule()
                                             == SelectionStrategy::NotConsideredForIntegralityCheck))
          continue;

        if (((it->first->type() == 'B') || (it->first->type() == 'I')) && (it->second > Dfloor(it->second)))
        {
          if (printL(3))
            std::cout << "ColGenSolver::checkIfCurSolIsInteger(): detected fractional slice and component "
                      << "curDisaggregateSol[" << tindex << "][" << it->first->name() << "] = " << it->second
                      << "  solIsInteger = false " << std::endl;
          _solIsInteger = false;;
          return false;
        }
      }
    }

    /// Else  check last slice (the others have been check on creation)
    for (VarPtr2DoubleMap::const_iterator it = curDisaggregateSol[tindex].begin();
         it != curDisaggregateSol[tindex].end(); ++it)
    {
      if (printL(5))
        std::cout << "last curDisaggregateSol[" << tindex << "][" << it->first->name() << "] = " << it->second
                  << std::endl;

      InstanciatedVar * instanciatedVarPtr = dynamic_cast<InstanciatedVar *> (it->first);
      if ((instanciatedVarPtr != NULL) && (instanciatedVarPtr->genVarConstrPtr()->priorityRule()
                                           == SelectionStrategy::NotConsideredForIntegralityCheck))
        continue;

      if (((it->first->type() == 'B') || (it->first->type() == 'I')) && (it->second > Dfloor(it->second)))
      {
        if (printL(3))
          std::cout << "ColGenSolver::checkIfCurSolIsInteger(): detected fractional slice and component "
                    << "curDisaggregateSol[" << tindex << "][" << it->first->name() << "] = " << it->second
                    << "  solIsInteger = false " << std::endl;
        _solIsInteger = false;
        return false;
      }
    }
  }

  if (printL(1) && _solIsInteger)
    std::cout << "ColGenSolver::checkIfCurSolIsInteger() current solution to the master is integer " << std::endl;

  return _solIsInteger;
}

void Alg4EvalOfNode::setRollbackPointSavedStatus(const bool & status)
{
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    (*spcPt)->setRollbackPointSavedStatus(status);
}

bool Alg4EvalOfNode::checkIfCurSolIsMasterLpFeasible()
{
  _solIsMasterLpFeasible = _probPtr->Problem::primalSolIsFeasible();
  
  return _solIsMasterLpFeasible;
}

void Alg4EvalOfNode::setOptionNeedBeConqueredIfSolIsInteger(const bool value)
{
  _needBeConqueredIfSolIsInteger = value;
}

bool Alg4EvalOfNode::setupAlgo(Node * nodePtr)
{

  bool ret = Alg4Master::setupAlgo(nodePtr);

  _pricingSolverCutsMessageId = PricingSolverCutsMessage::noMessage;
  _solIsMasterLpFeasible = false;
  _solIsInteger = false;

  setRollbackPointSavedStatus(false);

  return ret;
}

NodeEvalInfo * Alg4EvalOfNode::recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr)
{
  if (nodeEvalInfoPtr != NULL)
  {
    nodeEvalInfoPtr->treatOrderId = _currentNodePtr->treatOrder();
  }
  else
    nodeEvalInfoPtr = new NodeEvalInfo(_currentNodePtr->treatOrder());

  return nodeEvalInfoPtr;
}

void Alg4EvalOfNode::setDownAlgo()
{
  Alg4Master::setDownAlgo();
  _nonStabArtVarPtrList.clear();
  
  if (_probPtr == NULL)
      return;
  
  bool timeLimitIsReached = bapcodInit().startTime().getElapsedTime() > param().GlobalTimeLimitInTick();
  bool needBeConqueredIfSolIsInteger = _needBeConqueredIfSolIsInteger && param().StopWhenBranchingFails()
                                       && !timeLimitIsReached;

  bapcodInit().check(!isConquered() && _solIsInteger && _solIsMasterLpFeasible & needBeConqueredIfSolIsInteger,
                     "BaPCod error in Alg4EvalOfNode::setDownAlgo() : primal solution is integer"
                     " after node evaluation but the node duality gap is non-zero.\n "
                     " May be you should increase optimalityGapTolerance or relOptimalityGapTolerance parameter.");
}

