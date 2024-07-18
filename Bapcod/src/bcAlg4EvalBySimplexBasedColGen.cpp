/**
 *
 * This file bcAlg4EvalBySimplexBasedColGen.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <vector>
#include <unordered_map>

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalBySimplexBasedColGen.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcStabilizationColgen.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcFormC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif /* BCP_RCSP_IS_FOUND */

using namespace std;

ColGenEvalInfo::ColGenEvalInfo(const ColGenEvalInfo & that):
        NodeEvalInfo(), stabilizationInfoPtr(NULL), masterLpBasisPtr(new LpBasisRecord(*(that.masterLpBasisPtr))),
        latestReducedCostFixingGap(that.latestReducedCostFixingGap)
{
    if (that.stabilizationInfoPtr != NULL)
        stabilizationInfoPtr = new StabilizationInfo(*(that.stabilizationInfoPtr));
}

ColGenEvalInfo::~ColGenEvalInfo()
{
    if (printL(5))
        std::cout << "ColGenEvalInfo with " << *masterLpBasisPtr << " is deleted " << std::endl;
    if (stabilizationInfoPtr != NULL)
        delete stabilizationInfoPtr;
    if (masterLpBasisPtr != NULL)
        delete masterLpBasisPtr;
}

std::ostream & ColGenEvalInfo::print(std::ostream & os) const
{
    os << "ColGenEvalInfo with number of nodes = " << numberOfNodes
       << ", latestReducedCostFixingGap = " << latestReducedCostFixingGap << std::endl;
    if (masterLpBasisPtr != NULL)
        os << *masterLpBasisPtr;
    if (stabilizationInfoPtr != NULL)
        os << *stabilizationInfoPtr;
    return os;
}

Alg4EvalBySimplexBasedColGen::Alg4EvalBySimplexBasedColGen(Problem* const probPtr,
                                                           MasterCommons4EvalAlg& masterCommons) :
  Alg4EvalByLagrangianDuality(probPtr, masterCommons),
  _subProbSolutionsEnumeratedToMIP(false), _masterConverged(false), _canRunReducedCostFixing(false), _cutSepRound(0)
{

  if ((param().colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
      || (param().colGenDualPriceSmoothingAlphaFactor() > 0))
    _colGenStabilizationPtr = new ColGenStabilization(probPtr, _masterCommons.colGenSubProbConfPts(),
						                              bapcodInit().param());
}

Alg4EvalBySimplexBasedColGen::~Alg4EvalBySimplexBasedColGen()
{
}


bool Alg4EvalBySimplexBasedColGen::setPurePhaseI()
{
  bool update(false);
  _currentlyPerformingPhase1 = true;

  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 's');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); varPt++)
    {
      (*varPt)->resetCost(true);

      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->resetObjCoef(*varPt);

      if (printL(2))
        std::cout << "setPurePhaseI set cost to zero for Var " << (*varPt)->name() << std::endl;

      update = true;
    }

  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); varPt++)
    {
      (*varPt)->resetCost(true);

      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->resetObjCoef(*varPt);

      if (printL(2))
        std::cout << "setPurePhaseI set cost to zero for Var " << (*varPt)->name() << std::endl;

      update = true;
    }

  if (update)
    {
      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->updateObjectiveInFormulation();
    }

  return (update);
}

bool Alg4EvalBySimplexBasedColGen::unsetPurePhaseI()
{
  if (printL(5))
    std::cout << "Alg4EvalBySimplexBasedColGen::unsetPurePhaseI()" << std::endl;

  _need2reintroduceArtVarInMast = true;
  _currentlyPerformingPhase1 = false;

  bool update(false);
  /// Needed to avoid deleting in the set we are currently scanning
  VarPtrSet var2Desactivate;

  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 's');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); varPt++)
    {
      (*varPt)->resetCost(false);
      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->resetObjCoef(*varPt);

      if (printL(5))
        std::cout << "unsetPurePhaseI reset cost  for Var " << (*varPt)->name() << std::endl;

      update = true;
    }

  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); varPt++)
    {
      (*varPt)->resetCost(false);
      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->resetObjCoef(*varPt);

      if (printL(5))
        std::cout << "unsetPurePhaseI reset cost  for Var " << (*varPt)->name() << std::endl;

      update = true;
    }

  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'a');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'a'); varPt++)
    {
      (*varPt)->resetCost(false);
      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->resetObjCoef(*varPt);

      if (printL(5))
        std::cout << "unsetPurePhaseI reset cost  for Var " << (*varPt)->name() << std::endl;

      update = true;
    }

  if (update)
    {
      if (_probPtr->primalFormulationPtr() != NULL)
        _probPtr->primalFormulationPtr()->updateObjectiveInFormulation();
    }
  if (!var2Desactivate.empty())
    _probPtr->delVarSet(var2Desactivate, 1, 2);

  return (update);
}

/// returns true if we need to terminate the column generation algorithm
bool Alg4EvalBySimplexBasedColGen::solveMastLpPrimPh1A2(int & nbOfPenaltiesUpdates, int & nbCgIterations)
{
  /**
   * Implement column generation solution of master LP relaxation.
   * Combined phase I and phase II solution.
   */
  do
  /**
   * Iterate while soft constraints are not met
   * (while artificial variables still in master solution implying phase I is incomplete)
   */
    {
      bool interruptAlgorithm = solveMastLpPrimPh2(nbCgIterations);

      if (_colGenStabilizationPtr != NULL)
        _colGenStabilizationPtr->resetOnColGenTermination();

      if (interruptAlgorithm)
        return true;
      
      if (printL(2))
        std::cout << "solveMastLpPrimPh1A2() nbOfPenaltiesUpdates = " << nbOfPenaltiesUpdates << std::endl;
    }
  while ((nbCgIterations < _maxNbOfCgIterations) && (++nbOfPenaltiesUpdates <= _maxNbOfPenaltyUpdates)
         && updatePenalties(param().ArtVarPenaltyUpdateFactor()));

  if (isConquered())
    {
      if (printL(2))
        std::cout << "solveMastLpPrimPh1A2() : phase 1 cut by bound or infeasibility " << std::endl;
      return true;
    }

  if (nbCgIterations >= _maxNbOfCgIterations)
  {
    if (printL(0))
      std::cout << " Max. number of col. gen. iterations is reached " << std::endl;
      updateAlgPrimalLpBounds();
    return false;
  }

  if (phase1isCompleted())
    return false;

  if (_maxLevelOfSubProbRestriction > 0)
    {
      if (printL(0))
        std::cout << " Artificial variables in the solution : decreasing phase of stage " << std::endl;
      return false;
    }

  if (printL(0))
    std::cout << "  MaxNbOfPenaltiesUpdates is reached : do switch to pure phase 1 " << std::endl;

  if (_colGenStabilizationPtr != NULL)
    _colGenStabilizationPtr->deactivate();     /// stabilization is not implemented for the pure Phase I
  setPurePhaseI();

  /**
   * Phase II loop for artificial problem:
   * Iterate while can generate new columns and termination by bound does not apply
   */
  do
    {
      if (solveRestrictedMastLP() <= 0)     /// Problem infeasible
        {
          if (printL(2))
            std::cout << "ColGenSolver::solveMastLpPrimPh1A2() could not solveRestrictedMastLP()): problem infeasible"
                      << std::endl;
          markInfeasible();
          unsetPurePhaseI();
          return true;
        }
      nbCgIterations += 1;

      /// We manage to get all artificial column out
      if (phase1isCompleted())
        {
          /// Reset master variable costs
          unsetPurePhaseI();

          if (printL(2))
            std::cout << "solveMast: after Resolve master without artificial column " << std::endl;

          /// Resolve master without artificial column
          return solveMastLpPrimPh2(nbCgIterations);
        }

      /**
       * Search pool (master inactive columns) for MastCol with neg red Cost and add them to master
       */
      int nbAddedNegRedCostCol = param().Search4NegRedCostColInInactivePool() ? searchNegRedCostInactiveCol() : 0;
      int nbAllAddedCol = nbAddedNegRedCostCol;

      /**
        * If none were found in the column pool, go for a round of Col generation
        */
      if (nbAddedNegRedCostCol <= 0)
        {
          if (printL(2))
            std::cout << "solveMast: need to generate new MastCol" << std::endl;
          cleanupRestrictedMastColumns(nbCgIterations);
          int status = genNewColumns(_maxLevelOfSubProbRestriction, nbAllAddedCol, nbAddedNegRedCostCol);
          if (status < 0)
            {
              markInfeasible();
              unsetPurePhaseI();
              return true;
            }
        }

      long elapsedTime(0);
      if (printL(-1))
        printIntermediateStatistics(std::cout, _maxLevelOfSubProbRestriction, nbAllAddedCol, nbCgIterations,
                                    elapsedTime);

      if (printL(2))
        std::cout << "solveMast at CG iteration " << nbCgIterations
                  << " found " << nbAddedNegRedCostCol << " MastCol with neg red cost" << std::endl;

      if (nbAddedNegRedCostCol == 0)
        break;

      if (elapsedTime > param().GlobalTimeLimitInTick())
        {
          if (printL(0))
            std::cout << "Global time limit is reached during phase 1 of column generation" << std::endl;
          unsetPurePhaseI();
          return true;
        }
    }
  while (nbCgIterations < _maxNbOfCgIterations);

  if ((nbCgIterations >= _maxNbOfCgIterations) && !_nonExactEvaluation)
  {
    if (printL(-1))
      std::cout << "BaPCod WARNING : maximum number of col.gen. iterations reached, no guarantee for optimality"
                << std::endl;
    std::cerr << "BaPCod WARNING : maximum number of col.gen. iterations reached, no guarantee for optimality"
              << std::endl;
  }

  /// Reset master variable costs
  unsetPurePhaseI();

  if ((_maxLevelOfSubProbRestriction == 0) && !phase1isCompleted())
    {
      if (printL(0))
        std::cout << "Pure phase I determined infeasibility" << std::endl;
      markInfeasible();
    }

  return !_solIsMasterLpFeasible;
}

/// returns true if we need to terminate the column generation algorithm
bool Alg4EvalBySimplexBasedColGen::solveMastLpPrimPh2(int & nbCgIterations)
{
  if (!progStatus().doRun())
    return true;

  int curMaxLevelOfSubProbRestriction(_maxLevelOfSubProbRestriction);

  /**
   *  Phase II loop: Iterate while can generate new columns and termination by bound does not apply
   */
  do /// column generation iterations
    {
      bool phaseOneIsCompleted = false;
      int iterationOfCuttingPlaneIfIntegerSol(0);
      do /// cutting plane algorithm iterations for core constraints needed for the validity of the formulation
        {
          int solverReturnStatus = solveRestrictedMastLP();
          phaseOneIsCompleted = phase1isCompleted();

          if (_colGenStabilizationPtr != NULL)
            _colGenStabilizationPtr->initializationAfterSolvingRestrictedMaster(computeOptimGap(), nbCgIterations,
                                                                                curMaxLevelOfSubProbRestriction);

          if (solverReturnStatus <= 0) /// Problem infeasible
            {
              if (printL(2))
                std::cout << "solveMastLpPrimPh2() could not solveRestrictedMastLP()): problem infeasible"
                          << std::endl;
              markInfeasible();
              return true;
            }

          if (!phaseOneIsCompleted)
            {
              if (printL(0))
                std::cout << "#";
              if (printL(2))
                std::cout << " ColGenSolver::solveMastLpPrimPh2:  phase 1 is not yet completed:"
                          << " doing phase 2 of a given phase 1 stage" << std::endl;
              break; /// exit cutting plane procedure for core constraints
            }

          /// no artificial variables in the restricted master LP solution, we can update incumbent LP primal bound
            updateAlgPrimalLpBounds();

          bool potentiallyImrovingSolution = (algCurLpPrimalBound() < algIncIpPrimalBound());

          if (!potentiallyImrovingSolution || !checkIfCurSolIsInteger())
            break; /// exit cutting plane procedure for core constraints

          if (addCutToMaster('C'))
            {
              /// we need to reset the best Lp value if cuts are added
              resetAlgIncLpPrimalBound(_masterCommons.objStatus());
            }
          else
            {
              updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
              break; /// exit cutting plane procedure for core constraints
            }
        }
      while (bapcodInit().require(++iterationOfCuttingPlaneIfIntegerSol < 100,
                                  "MaxNb Of core cut generation (set to 100) is exceeded"));
      nbCgIterations += 1;

      long elapsedTime = bapcodInit().startTime().getElapsedTime();

      /**
       * Early Termination Test of Phase I and II:
       * terminate column generation if node can be pruned (due to new incumbent solution)
       */
      if (isConquered())
        {
          if (printL(-1))
            printIntermediateStatistics(std::cout, curMaxLevelOfSubProbRestriction, 0, nbCgIterations,
                                        elapsedTime, true, true);
          if (printL(1))
              std::cout << "ColGen early termination due to IP bounds : [" << std::setprecision(12)
                        << algIncIpDualBound() << "," << algIncIpPrimalBound() << "]" << std::setprecision(6)
                        << std::endl;
          return true;
        }

      /// Early Termination Test of Phase II: since master LP value may have changed
      if (phaseOneIsCompleted && earlyCGtermType1())
        {
            if (printL(-1))
                printIntermediateStatistics(std::cout, curMaxLevelOfSubProbRestriction, 0, nbCgIterations,
                                            elapsedTime, true,
                                            (_maxLevelOfSubProbRestriction == 0));
          /// we do not set _masterConverged to true as the best dual bound does not necessarily comes
          /// from this run of the evaluation algorithm
          return false;
        }

      cleanupRestrictedMastColumns(nbCgIterations);

      int nbAddedNegRedCostCol = param().Search4NegRedCostColInInactivePool() ? searchNegRedCostInactiveCol() : 0;
      int nbAllAddedCol = nbAddedNegRedCostCol;

      if ((nbAddedNegRedCostCol > 0) && (param().colgeninfo_file() != ""))
         std::cout << "p";

      ///  If none were found in the column pool, go for a round of Col generation
      if (nbAddedNegRedCostCol <= 0)
        {
          do     /// iterate on smoothing phase with column generation and lagrangian bound computation
            {
              if (printL(2))
                std::cout << "solveMast: need to generate new MastCol" << std::endl;

              if ((_colGenStabilizationPtr != NULL) && _colGenStabilizationPtr->solValueSmoothingIsActive())
                  /// if smoothing is active, the dual solutions will be different at the moment of reduced cost
                  /// fixing, therefore reduced cost fixing should not be run
                  _canRunReducedCostFixing = false;
              else
                  _canRunReducedCostFixing = true;

              int status = genNewColumns(curMaxLevelOfSubProbRestriction, nbAllAddedCol, nbAddedNegRedCostCol);

              ///  In case subproblem infeasibility results in master infeasibility
              if (status < 0)
                {
                  markInfeasible();
                  return true;
                }

              if (_pricingSolverCutsMessageId == PricingSolverCutsMessage::doCutsRollback)
                return true;

              if (_pricingSolverCutsMessageId == PricingSolverCutsMessage::interruptSolution)
                {
                  progStatus().setStat(ProgStatus::terminate);
                  if (printL(-1))
                    std::cout << "SEARCH IS INTERRUPTED due to the pricing problem return status. " << std::endl;
                  /// current solution can be integer but there will be duality gap at the end of the node
                  /// which will trigger an error, so we need to change _solIsInteger
                  _solIsInteger = false;
                  return true;
                }

                updateLagrangianDualBound(true);
            }
          while ((_colGenStabilizationPtr != NULL)
                 && (_colGenStabilizationPtr->updateAfterPricingProblemSolution(nbAddedNegRedCostCol)));

          if (curMaxLevelOfSubProbRestriction != _maxLevelOfSubProbRestriction)
            _maxLevelOfSubProbRestriction = curMaxLevelOfSubProbRestriction;

          if (printL(-1))
            printIntermediateStatistics(std::cout, curMaxLevelOfSubProbRestriction, nbAllAddedCol, nbCgIterations,
                                        elapsedTime);

          if (_colGenStabilizationPtr != NULL)
            _colGenStabilizationPtr->updateAfterColGenIteration();

          /**
           * Early Termination Test of Phase I and II: due to new dual bound
           */
          if (isConquered())
            {
              if (printL(2))
                std::cout << "ColGenSolver early termination : node is pruned by Bound " << std::endl;
              return true;
            }
          /**
           * Early Termination Test of Phase II: since master dual bound may have changed
           */
          if (phaseOneIsCompleted && earlyCGtermType1())
            {
              _masterConverged = (_maxLevelOfSubProbRestriction == 0);
              return false;
            }
        }
      else if (printL(-1))
        {
          printIntermediateStatistics(std::cout, curMaxLevelOfSubProbRestriction, nbAllAddedCol, nbCgIterations,
                                      elapsedTime);
        }

      if (printL(2))
        std::cout << "CG iteration " << nbCgIterations << " : inserted " << nbAllAddedCol << " columns" << std::endl;

      if (nbAddedNegRedCostCol == 0)
        {
          _masterConverged = (_maxLevelOfSubProbRestriction == 0) && phaseOneIsCompleted;
          return false;
        }

      if (elapsedTime > param().GlobalTimeLimitInTick())
        {
          std::cout << "Global time limit is reached during column generation" << std::endl;
          return true;
        }
    }
  while (nbCgIterations < _maxNbOfCgIterations);

  if ((nbCgIterations >= _maxNbOfCgIterations) && !_nonExactEvaluation)
  {
    if (printL(-1))
      std::cout << "BaPCod WARNING : maximum number of col.gen. iterations reached, no guarantee for optimality"
                << std::endl;
    std::cerr << "BaPCod WARNING : maximum number of col.gen. iterations reached, no guarantee for optimality"
              << std::endl;
  }

  if (printL(2))
    std::cout << "solveMastLpPrimPh2() is finished" << std::endl;

  return false;
}

int Alg4EvalBySimplexBasedColGen::solveRestrictedMastLP()
{
  int formPrintLevel = 2;

  Time startbcTimeMastMPsol;
  int maxLevelOfRestriction(0);

  int solverReturnStatus = 0;
  if (_maxNbOfCgIterations == 0)
  {
    /// If there is no pricing subproblem then we call solveProb() directly
    /// because it is overloaded for MipProblem. Therefore, if masterSolMode==MipSolver
    /// (in Classical Benders for example) then MipProblem::solveProb() will be
    /// called instead of Problem::solveProb() (which always solves the LP).
    solverReturnStatus = _probPtr->solveProb(maxLevelOfRestriction, ' ', printL(formPrintLevel));
  }
  else
  {
    /// 'd' means do record the dual solution
    solverReturnStatus = _probPtr->Problem::solveProb(maxLevelOfRestriction, 'd', printL(formPrintLevel));
  }
  
  // Added by Guillaume. Retrieve reduced costs used for stabilization in
  // computation of the master contribution to the Lagrangian Bound.
    if ((param().masterSolMode().status() == SolutionMethod::lpSolver || 
            param().masterSolMode().status() ==  SolutionMethod::mipSolver) && 
        ((param().colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
        || (param().colGenDualPriceSmoothingAlphaFactor() > 0) || param().SplitColIntoDissagregateSpVar )
        && !_probPtr->probVarSet().empty(VcIndexStatus::Active, 's') ) {
      _probPtr->retrieveRedCosts();
    }

  bapcodInit().statistics().incrTimer("bcTimeMastMPsol", startbcTimeMastMPsol.getElapsedTime_dbl());

  _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());

  if (printL(0) && (param().printMasterPrimalSols() == 1))
    _probPtr->printDetailedPrimalSol();

  if (printL(1))
      std::cout << " Restricted master LP is solved in " << startbcTimeMastMPsol.getElapsedTime()/100.0
                << " seconds" << std::endl;

  bapcodInit().statistics().incrRecord("bcAverageDualSolSize", int(_probPtr->inDualSol().size()));
  bapcodInit().statistics().incrCounter("bcCountMastSol");
  bapcodInit().statistics().setCounter("bcCountPrimalSolSize", int(_probPtr->inPrimalLpSol().size()));
  bapcodInit().statistics().incrRecord("bcAveragePrimalSolSize", int(_probPtr->inPrimalLpSol().size()));

  return (solverReturnStatus);
}

bool Alg4EvalBySimplexBasedColGen::updatePenalties(const Double& factor)
{
  if (printL(2))
    std::cout << "Alg4EvalBySimplexBasedColGen::updatePenalties(); factor =  " << factor << std::endl;

  if ((_maxLevelOfSubProbRestriction > 0) && gapSmallerThanTol(algCurLpDualBound(), algIncIpPrimalBound(), param()))
    return false;

  if (_colGenStabilizationPtr != NULL)
    {
      if (printL(0) && (param().colgeninfo_file() != "") && (_logPrintFrequency > 0))
        std::cout << "# ";
      if (_colGenStabilizationPtr->updateOnArtVarsInFinalSolution())
        return (true);
    }

  if (_nonStabArtVarPtrList.empty())
    return (false);

  /// we verify where there are feasibility related artificial variables in the solution
  /// if yes, we increase their cost by factor

  std::list<Variable *>::const_iterator varIt;
  bool atLeastOneArtVarInSolution(false);

    for (varIt = _nonStabArtVarPtrList.begin(); varIt != _nonStabArtVarPtrList.end(); ++varIt)
    {
        if (_probPtr->inPrimalLpSol().count(*varIt))
        {
            atLeastOneArtVarInSolution = true;
            break;
        }
    }
    if (atLeastOneArtVarInSolution)
    {
        std::list<Variable *> varsToUpdateList;
        for (varIt = _nonStabArtVarPtrList.begin(); varIt != _nonStabArtVarPtrList.end(); ++varIt)
            if ((*varIt)->vcIndexStatus() == VcIndexStatus::Active)
            {
                (*varIt)->resetCurCostByValue((*varIt)->curCost() * factor);
                varsToUpdateList.push_back(*varIt);
            }
        _probPtr->updateObjCoeffsInForm(varsToUpdateList);
    }

  return (atLeastOneArtVarInSolution);
}

bool Alg4EvalBySimplexBasedColGen::setupAlgo(Node * nodePtr)
{
  if (Alg4EvalOfNode::setupAlgo(nodePtr))
    return true;

  ColGenEvalInfo * colGenEvalInfoPtr = static_cast<ColGenEvalInfo *>(nodePtr->nodeEvalInfoPtr());

  bapcodInit().require(colGenEvalInfoPtr != NULL,
                       "BaPCod error: NodeEvalInfo for ColGenEvalAlg is not of type colGenSolverInfo.");

  for (std::vector<ColGenSpConf *>::const_iterator spIt = _masterCommons.colGenSubProbConfPts().begin();
       spIt != _masterCommons.colGenSubProbConfPts().end(); ++spIt)
    if (_currentNodePtr->cgSpConfTreeOfColClassesMap().count(*spIt))
      {
        ColClassesVector & cgSpTreeOfColClasses = _currentNodePtr->cgSpConfTreeOfColClassesMap()[*spIt];
        if (!cgSpTreeOfColClasses.empty())
          {
            (*spIt)->spConfHasClassInducedSpVarBounds(true);
            (*spIt)->treeOfColClasses() = cgSpTreeOfColClasses;
          }
      }

  _isLastSpcPtInitialized = false;

  /// fill _nonStabArtVarPtrList
  for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'a');
       varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'a'); ++varPtrIt)
    {
      if ((*varPtrIt)->isTypeOf(VcId::LocalArtificialVarMask))
        {
          LocalArtificialVar * artVarPtr = static_cast<LocalArtificialVar *>(*varPtrIt);
          LocalArtificialVar::LocalArtClassId locClassId = artVarPtr->localClassId();
          if ((locClassId != LocalArtificialVar::NegLocalId) && (locClassId != LocalArtificialVar::PosLocalId))
            continue;
        }
      _nonStabArtVarPtrList.push_back(*varPtrIt);
    }

  _savedMasterLPTime = bapcodInit().statistics().getTime("bcTimeMastMPsol");
  _savedPricingTime = bapcodInit().statistics().getTime("bcTimeCgSpOracle");
  _savedNbColumns = 0;

  /// we do not update basis if this node is treated immediately after its parent
  /// (so the basis in the LP solver remains unchanged)
  if ((_currentNodePtr->treatOrder() != colGenEvalInfoPtr->treatOrderId + 1)
      && (colGenEvalInfoPtr->masterLpBasisPtr != NULL))
    {
      /// we need to add to the basis the local branching constraints, as there were not in the problem
      /// when the basis was retrieved, this is except the case when a node was obtained from previous phase node
      /// in strong branching, in this case the updateBasisWithLocBrConstrsOnSetup will be set to false
      /// (desactivated for the moment), TO DO : move basis reloading to the problem setup
      LpBasisRecord * basisPtrToReload = new LpBasisRecord(*(colGenEvalInfoPtr->masterLpBasisPtr));
      for (std::list<BranchingConstrBaseType *>::const_iterator
           constrIt = _currentNodePtr->localNodeBrConstrList().begin();
           constrIt != _currentNodePtr->localNodeBrConstrList().end(); ++constrIt)
        {
          Constraint * constrPtr = dynamic_cast<Constraint *>(*constrIt);
          if ((constrPtr != NULL) && (!constrPtr->isTypeOf(VcId::InstSubProbBranchingConstrMask)))
            basisPtrToReload->_constrInBasis.push_back
            (ConstrPtr_MpFormIndexStatus(constrPtr, MathProgSolverInterface::BasicVarConstr));
        }
      _probPtr->reloadMemorizedBasis(basisPtrToReload);
      delete basisPtrToReload;
    }

  _probPtr->updateInDualSol(); /// residue from ColGenSolver::setupProbConf, to check if we still need this
  _need2reintroduceArtVarInMast = false; /// residue from ColGenSolver::setupProbConf, to check if we still need this

  if (_colGenStabilizationPtr != NULL)
    _colGenStabilizationPtr->setupStab(colGenEvalInfoPtr->stabilizationInfoPtr, algIncLpDualBound(),
                                       param().MaxNbOfStagesInColGenProcedure() - 1, _currentNodePtr->depth());

  _latestReducedCostFixingGap = colGenEvalInfoPtr->latestReducedCostFixingGap;

  _addCutToMasterFirstCall = nodePtr->isRoot();

  if (_masterCommons.colGenSubProbConfPts().empty())
    _maxNbOfCgIterations = 0;

  _subProbSolutionsEnumeratedToMIP = false;

  return false;
}

bool Alg4EvalBySimplexBasedColGen::eval()
{
  bapcodInit().require(_currentNodePtr != NULL, "_currentNodePtr of Alg4EvalBySimplexBasedColGen should not be NULL");
  bapcodInit().require(_probPtr != NULL, "_probPtr of Alg4EvalBySimplexBasedColGen should not be NULL");

  _solIsMasterLpFeasible = true; /// in column generation master solution is always LP feasible

  int nbCgIterations(0);
  /// we should repeat column generation if the master problem was changed
  /// after some subproblem passed to the enumerated state or after subrob. relaxation level increase //
  do
    {
      _canRunReducedCostFixing = true;
      _masterConverged = false;
      /// we need to reset the best Lp value after a round of cuts or after the reduced cost fixing
      resetAlgIncLpPrimalBound(_masterCommons.objStatus());

      _maxLevelOfSubProbRestriction = param().MaxNbOfStagesInColGenProcedure;

      bapcodInit().check(_maxLevelOfSubProbRestriction < 1,
                         "Alg4EvalBySimplexBasedColGen::eval() ZERO PHASE in stage column generation");

      /// iterate on different stages of the column generation procedure
      while (_maxLevelOfSubProbRestriction > _minLevelOfSpRestriction)
        {
          _maxLevelOfSubProbRestriction -= 1;

          if (printL(2))
            std::cout << "Alg4EvalBySimplexBasedColGen::eval() STAGE " << _maxLevelOfSubProbRestriction << std::endl;

          /// reset stage dual bound
          resetAlgIncStageLpDualBound();
          _algCurLpDualBound = algIncStageLpDualBound();

          if (_maxNbOfCgIterations == 0)
            {
              if (solveRestrictedMastLP() <= 0)
                markInfeasible();
              _algCurLpDualBound = _algCurLpPrimalBound; /// we use here strong duality
              updateAlgPrimalLpBounds();
              updateAlgDualBounds();

              /// deactivated as core cuts are not separated, TO DO : implement it
              // if (checkIfCurSolIsInteger())
              //  updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
              break; /// no reduced cost fixing if no column generation
            }
          else
            {
              Time start;
              int nbOfPenaltiesUpdates = 0;
              bool interruptAlgorithm = solveMastLpPrimPh1A2(nbOfPenaltiesUpdates, nbCgIterations);

              if (printL(2))
                std::cout << "solveMastLpPrimPh1A2() is finished" << std::endl;

              bapcodInit().statistics().incrTimer("bcTimeColGen", start.getElapsedTime_dbl());

              if (nbOfPenaltiesUpdates > bapcodInit().statistics().getCounter("bcMaxPUpd"))
                bapcodInit().statistics().setCounter("bcMaxPUpd", nbOfPenaltiesUpdates);

              if (printL(-1) && (_logPrintFrequency > 1) && (nbCgIterations % _logPrintFrequency != 0))
              {
                long elapsedTime(0);
                printIntermediateStatistics(std::cout, _maxLevelOfSubProbRestriction, 0, nbCgIterations, elapsedTime,
                                            true, (_maxLevelOfSubProbRestriction == 0));
              }

              if (_maxLevelOfSubProbRestriction == _minLevelOfSpRestriction)
                for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
                     spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
                  (*spcPt)->probPtr()->callColGenTerminationCallBack(false, _currentNodePtr->treatOrder(),
                                                                     _currentNodePtr->depth(), _cutSepRound,
                                                                     algIncLpDualBound(),
                                                                     bapcodInit().startTime().getElapsedTime(),
                                                                     _masterConverged);

              if (interruptAlgorithm)
                break;

              if (_minLevelOfSpRestriction > 0)
              {
                /// if all subproblems are enumerated, we go to the exact pricing in any case,
                /// otherwise heuristics may not be efficient
                bool allSubproblemsAreEnumerated = true;
                for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
                     allSubproblemsAreEnumerated && (spcPt != _masterCommons.colGenSubProbConfPts().end()); ++spcPt)
                  allSubproblemsAreEnumerated = allSubproblemsAreEnumerated && (*spcPt)->enumeratedStatus();
                if (allSubproblemsAreEnumerated)
                  _minLevelOfSpRestriction = 0;
              }

            }
        }

      if (printL(0) && (_maxLevelOfSubProbRestriction == 0) && (param().printMasterPrimalSols() == 2))
        _probPtr->printDetailedPrimalSol();
    }
  while (_canRunReducedCostFixing && runReducedCostFixingAndEnumeration());

    for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
        (*spcPt)->probPtr()->callColGenTerminationCallBack(true, _currentNodePtr->isRoot(), _currentNodePtr->depth(),
                                                           _cutSepRound, algIncLpDualBound(),
                                                           bapcodInit().startTime().getElapsedTime(),
                                                           _masterConverged);

    updatePricingSolverCutsMessageId();

  if (printL(-1) && (_maxNbOfCgIterations != 0) && (_logPrintFrequency == 0))
    {
      long elapsedTime(0);
      printIntermediateStatistics(std::cout, _minLevelOfSpRestriction, 0, nbCgIterations,
                                  elapsedTime, true);
    }

  /// added by Ruslan : dual bound can be slightly more than master solution value
  /// because of floating point calculation, we make the dual bound equal to lp value in this case
  if (phase1isCompleted())
    {
      Bound buBound = algIncLpDualBound();
      rectifyIncumbentLpValue();
    }

  int nbOfColumnInMaster = _probPtr->probVarSet().size(VcIndexStatus::Active,'s')
                           + _probPtr->probVarSet().size(VcIndexStatus::Active, 'd');
  if (nbOfColumnInMaster > bapcodInit().statistics().getCounter("bcMaxNbColInMastLp"))
    bapcodInit().statistics().setCounter("bcMaxNbColInMastLp", nbOfColumnInMaster);

  int nbOfRowsInMaster = _probPtr->probConstrSet().size(VcIndexStatus::Active, 's')
                         + _probPtr->probConstrSet().size(VcIndexStatus::Active, 'd');
  if (nbOfRowsInMaster > bapcodInit().statistics().getCounter("bcMaxNbRowInMastLp"))
    bapcodInit().statistics().setCounter("bcMaxNbRowInMastLp", nbOfRowsInMaster);

  if (nbCgIterations > bapcodInit().statistics().getCounter("bcMaxIterCg"))
    bapcodInit().statistics().setCounter("bcMaxIterCg", nbCgIterations);

  if (printL(-1) && (_maxNbOfCgIterations != 0) && (_minLevelOfSpRestriction == 0) && (_logPrintFrequency > 0))
    {
      if (param().SafeDualBoundScaleFactor() > 0)
        std::cout << std::setprecision(20) << "ColGenEvalAlg final dual bound: " << algIncLpDualBound()
                  << " (rounded: " << algIncIpDualBound() << ")" << std::setprecision(6) << std::endl;
      else
        std::cout << "ColGenEvalAlg final dual bound: " << algIncLpDualBound()
                  << " (rounded: " << algIncIpDualBound() << ")" << std::endl;
    }

  if (printL(5))
    std::cout << " Alg4EvalBySimplexBasedColGen::solve() EXITING " << std::endl;

  /// should return false problem if solution feasible
  return (!_solIsMasterLpFeasible);
}

/// returns true if lightening was successeful for at least one subproblem
bool Alg4EvalBySimplexBasedColGen::runColGenSpRelaxationLightening(const int callMode)
{
    bool atLeastOneSubProbWasLightened = false;
    for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
        atLeastOneSubProbWasLightened = (*spcPt)->performSpRelaxationLigntening(masterConverged(), callMode)
                                        || atLeastOneSubProbWasLightened;

    updatePricingSolverCutsMessageId();

    return atLeastOneSubProbWasLightened;
}


struct ColGenSpConfCmpById
{
    bool operator()(const ColGenSpConf * const aPtr, const ColGenSpConf * const bPtr) const
    {
      return (aPtr->id() < bPtr->id());
    }
};


bool Alg4EvalBySimplexBasedColGen::addEnumColumnsToMaster(const int maxNumOfColumns)
{
  long int numEnumeratedSolutions = _masterCommons.totalNumberOfEnumeratedSubprobSolutions();
  if ((numEnumeratedSolutions >= 0) && (numEnumeratedSolutions <= maxNumOfColumns))
  {
    for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
      Solution * enumSolutionPtr = new Solution(*spcPt);
      std::vector<double> dummyVector;
      (*spcPt)->probPtr()->getEnumeratedSolutions(-1, enumSolutionPtr, dummyVector);
      Solution * solPtr = enumSolutionPtr;
      while (solPtr != NULL)
      {
        /// we do not accept empty solutions if the number of enumerated solutions is zero,
        /// the first solution will be empty
        if (!solPtr->solVarValMap().empty())
          (*spcPt)->recordSubproblemSolution(solPtr, false, 1);
        solPtr = solPtr->nextSolPtr();
      }
      (*spcPt)->insertAllColumnsInMaster();
      recordColInForm();
      enumSolutionPtr->deleteSolutionsChain();
      delete enumSolutionPtr;
    }
    return true;
  }
  return false;
}

bool Alg4EvalBySimplexBasedColGen::runEnumSolBasedUserHeuristicFunctor()
{
    if ((_masterCommons.totalNumberOfEnumeratedSubprobSolutions() < 0)
        || (param().MaxNumEnumSolsInUserHeuristicFunctor() <= 1)
        || !_probPtr->enumSolBasedHeuristicFunctorDefined())
        return false;

    int maxNumOfColumns = param().MaxNumEnumSolsInUserHeuristicFunctor();
    std::vector<ColGenSpConf *>::const_iterator spcPt;
    std::map<ColGenSpConf *, Solution *, ColGenSpConfCmpById> enumeratedSolPtrMap;
    std::vector<std::pair<double, Solution *> > solsVector;
    for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
        Solution * solPtr = enumeratedSolPtrMap[*spcPt] = new Solution(*spcPt);
        std::vector<double> redCosts;
        (*spcPt)->probPtr()->getEnumeratedSolutions(maxNumOfColumns, solPtr, redCosts);
        int solId = 0;
        while (solPtr != NULL)
        {
            if (!solPtr->solVarValMap().empty())
                solsVector.push_back(std::make_pair(redCosts[solId] + (*spcPt)->fixedDualCost(), solPtr));
            solPtr = solPtr->nextSolPtr();
            solId += 1;
        }
    }
    std::stable_sort(solsVector.begin(), solsVector.end());
    long int nbOfEnumSolsToPass = (std::min)(maxNumOfColumns, (int)solsVector.size());
    if (printL(0))
        std::cout << "Passing " << nbOfEnumSolsToPass << " enum. solutions to the user heuristic functor "
                  << std::endl;

    std::vector<Solution *> solPts;
    solPts.reserve(nbOfEnumSolsToPass);
    for (int solId = 0; solId < nbOfEnumSolsToPass; ++solId)
    {
        solPts.push_back(solsVector[solId].second);
    }

    if (_probPtr->runEnumSolBasedHeuristicFunctor(solPts, algIncIpPrimalBound()))
    {
        _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());
        updateAlgPrimalLpBounds();
        if (checkIfCurSolIsInteger() && !addCutToMaster('C'))
        {
            if (printL(-1))
                std::cout << "Enum. sol. based user heur. functor returned solution of value "
                          << _algCurLpPrimalBound << std::endl;
            updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
        }
        _masterConverged = false; /// as solution is changed
    }

    for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
        (*spcPt)->insertAllColumnsInMaster();
        recordColInForm();
        enumeratedSolPtrMap[*spcPt]->deleteSolutionsChain();
        delete enumeratedSolPtrMap[*spcPt];
    }
    return true;
}

bool Alg4EvalBySimplexBasedColGen::addEnumColumnsWithSmallestRedCostInMaster(const int maxNumOfColumns)
{
    if ((_masterCommons.totalNumberOfEnumeratedSubprobSolutions() < 0) || (maxNumOfColumns <= 0))
        return false;

    std::vector<ColGenSpConf *>::const_iterator spcPt;
    std::map<ColGenSpConf *, Solution *, ColGenSpConfCmpById> enumeratedSolPtrMap;
    std::vector<std::pair<double, Solution *> > solsVector;
    for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
        Solution * solPtr = enumeratedSolPtrMap[*spcPt] = new Solution(*spcPt);
        std::vector<double> redCosts;
        (*spcPt)->probPtr()->getEnumeratedSolutions(maxNumOfColumns, solPtr, redCosts);
        int solId = 0;
        while (solPtr != NULL)
        {
            if (!solPtr->solVarValMap().empty())
                solsVector.push_back(std::make_pair(redCosts[solId] + (*spcPt)->fixedDualCost(), solPtr));
            solPtr = solPtr->nextSolPtr();
            solId += 1;
        }
    }
    std::stable_sort(solsVector.begin(), solsVector.end());
    long int nbOfColumnsToGenerate = (std::min)(maxNumOfColumns, (int)solsVector.size());
    if (printL(0))
        std::cout << "Added " << nbOfColumnsToGenerate << " enum. columns to the heuristic restricted master "
                  << std::endl;
    for (int solId = 0; solId < nbOfColumnsToGenerate; ++solId)
    {
        Solution * solPtr = solsVector[solId].second;
        ColGenSpConf * cgSpConfPtr = static_cast<ColGenSpConf *>(solPtr->probConfPtr());
        MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false,1);
    }
    for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
         spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
        (*spcPt)->insertAllColumnsInMaster();
        recordColInForm();
        enumeratedSolPtrMap[*spcPt]->deleteSolutionsChain();
        delete enumeratedSolPtrMap[*spcPt];
    }

    _masterConverged = false; /// as solution is changed
    return true;
}


/// the red.cost fixing and enumeration will be done with the current gap if falseGap < 0
/// enumerationMode == 0 if standard enumeration
/// enumerationMode == 1 if enumeraiton with false gap for heuristics
///                      + master is completed with ColumnCleanupThreshold columns with the smallest reduced cost
/// enumerationMode == 2 if enumeraiton with false gap for heuristics
///                      + red. cost. and enum. are not performed for already enumerated subproblems
///                      + no enumeration to the MIP
///                      + also columns with reduced cost larger than the falseGap are not removed from the master
/// enumeration with false gap if falseGap > 0, otherwise standard enumeration with the current gap
/// returns true if the column generation should be repeated
/// TO DO : for the moment, valid only for the minimization objective
bool Alg4EvalBySimplexBasedColGen::runReducedCostFixingAndEnumeration(const int enumerationMode, const double falseGap)
{
  /// in order to run the standard reduced cost fixing (enumerationMode == 0),
  /// we need full convergence of column generation (_masterConverged = true)
  if (!progStatus().doRun() || !_solIsMasterLpFeasible || !_doRedCostFixingAndEnumeration
      || (!_masterConverged && (enumerationMode == 0) ) )
    return false;

  Time startRedCostTime;

  bool repeatColumnGeneration = false;

  std::vector<ColGenSpConf *>::const_iterator spcPt;
  bool atLeastOneOfSubproblemsIsEnumerated = false;
  for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    if ((*spcPt)->enumeratedStatus())
      atLeastOneOfSubproblemsIsEnumerated = true;

  Bound incumbentValue = algIncIpPrimalBound();
  long long int scaleFactor = param().SafeDualBoundScaleFactor();
  BcObjStatus::MinMaxIntFloat objStatus = _masterCommons.objStatus();
  if (objStatus == BcObjStatus::minInt)
    incumbentValue = Bound(floor((double)(incumbentValue - Double::precision))
                           + (scaleFactor > 0 ? 0 : Double::precision), objStatus);
  Double currentGap = falseGap > 0 ? Double(falseGap) : incumbentValue - algCurLpDualBound();
  bool currentGapIsBelowThreshold = (currentGap / _latestReducedCostFixingGap < param().ReducedCostFixingThreshold());
  if (atLeastOneOfSubproblemsIsEnumerated || currentGapIsBelowThreshold || (_doRedCostFixingAndEnumeration == 2))
    {
      if (currentGapIsBelowThreshold)
        _latestReducedCostFixingGap = currentGap;
      std::set<ColGenSpConf *, ColGenSpConfCmpById>  justEnumeratedSubproblems;
      for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
           spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
        {
            if ((*spcPt)->upperBoundPtr() != nullptr && (*spcPt)->upperBoundMastConstrPtr()->curRhs() <= 0)
                continue;
          bool wasEnumerated = (*spcPt)->enumeratedStatus();
          if (!wasEnumerated && !currentGapIsBelowThreshold && (_doRedCostFixingAndEnumeration == 1))
              continue;
          if (wasEnumerated && (enumerationMode == 2))
              continue;
          int enumerationModeToPass = enumerationMode;
          if (enumerationModeToPass > 1)
            enumerationModeToPass = 1;
          if ((enumerationMode == 0)
              && (currentGap / std::abs(incumbentValue._val) > param().RCSPmaxGapToRunEnumeration()))
            enumerationModeToPass = -1;
          Double threshold(0.0);
          if (falseGap > 0)
          {
            threshold = Double(falseGap) - (*spcPt)->fixedDualCost(); /// previous version, not always correct
          }
          else {
              if (scaleFactor > 0)
                  threshold = incumbentValue - algCurLpDualBound()
                              + (*spcPt)->probPtr()->dualBound()._val / scaleFactor;
              else
                  threshold = incumbentValue - algCurLpDualBound() + (*spcPt)->probPtr()->dualBound();
          }
          if (scaleFactor > 0)
            threshold._val = ceil(threshold._val * scaleFactor) + 1;
          (*spcPt)->probPtr()->reducedCostFixingAndEnumeration(enumerationModeToPass, threshold);
          if (!wasEnumerated && (*spcPt)->enumeratedStatus())
            justEnumeratedSubproblems.insert(*spcPt);

        }
      if (atLeastOneOfSubproblemsIsEnumerated || !justEnumeratedSubproblems.empty())
        {
          /// this may (rarely) cause Cplex error 1017 (CPXERR_NOT_FOR_MIP), I do not understand why
          /// so for the moment below we use colPtr->computeReducedCost() instead of colPtr->reducedCost()
          //_probPtr->retrieveRedCosts(); /// fast reduced cost retrieval
          int numEnumColsRemoved = 0;
          int numColsRegeneratedWithEnumFlag = 0;

          //std::vector<Variable *> activeProblemColumns;
          std::map<ColGenSpConf *, std::vector<MastColumn *> > activeProblemColumns;
          for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
               spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
              activeProblemColumns[*spcPt] = std::vector<MastColumn *>();

          for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
               varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); ++varPtrIt)
            if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            {
                MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
                activeProblemColumns[colPtr->cgSpConfPtr()].push_back(colPtr);
            }

          std::list<Variable *> varsToRemoveFromForm;
          int numEnumSols =  _masterCommons.totalNumberOfEnumeratedSubprobSolutions();
          for (spcPt = _masterCommons.colGenSubProbConfPts().begin();
               spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
            {
                if (justEnumeratedSubproblems.count(*spcPt))
                {
                    /// this subproblems has just passed to the enumerated status,
                    /// so we need to regenerate some columns with the enumerated flag
                    int numCols = (int)activeProblemColumns[*spcPt].size();
                    std::vector<bool> toCreateEnumCopy(numCols, false);
                    std::vector<MastColumn *> & colPts = activeProblemColumns[*spcPt];

                    /// here we check which columns are not in the set of enumerated solutions, this might be because
                    /// 1) the column is non-proper
                    /// 2) the column has a reduced cost larger than the threshold,
                    /// 3) the column has an alternative column with a better cost and with the same vector of
                    ///    coefficients in the essential constraints
                    /// for other columns, we create the enumerated copy
                    /// (as the coefficients may change due to non-robust cuts)

                    /// deactivated, as condition 3) just above is too restrictive, and impacts heavily the
                    /// convergence of column generation just after successful enumeration
//                    std::vector<Solution *> solPts;
//                    solPts.reserve(numCols);
//                    for (auto * colPtr : activeProblemColumns[*spcPt])
//                        solPts.push_back(colPtr->spSol());
//                    (*spcPt)->probPtr()->checkEliminatedEnumeratedPaths(solPts, toCreateEnumCopy);

                    /// instead, we just verify conditions 1) and 2) to decide whether a column should be regenerated
                    for (int colId = 0; colId < numCols; ++colId)
                    {
                        if ((*spcPt)->probPtr()->isProperSolution(colPts[colId]->spSol())
                            && ((enumerationMode == 2) || (colPts[colId]->computeReducedCost() < currentGap)))
                            toCreateEnumCopy[colId] = true;
                    }

                    for (int colId = 0; colId < numCols; ++colId)
                    {
                        if (toCreateEnumCopy[colId])
                        {
                            (*spcPt)->recordSubproblemSolution(colPts[colId]->spSol(), false,
                                                               1,NULL, true);

                            numColsRegeneratedWithEnumFlag += 1;
                        }
                        else
                        {
                            numEnumColsRemoved += 1;
                        }
                        MastColumn * colPtr = activeProblemColumns[*spcPt][colId];
                        _probPtr->probVarSet().insert(colPtr, VcIndexStatus::Unsuitable);
                        colPtr->desactivate();
                        varsToRemoveFromForm.push_back(colPtr);
                    }
                }
                else if ((*spcPt)->enumeratedStatus() && (enumerationMode != 2))
                {
                    for (auto colPtr : activeProblemColumns[*spcPt])
                        if (colPtr->computeReducedCost() > currentGap)
                        {
                            numEnumColsRemoved += 1;
                            _probPtr->probVarSet().insert(colPtr, VcIndexStatus::Unsuitable);
                            colPtr->desactivate();
                            varsToRemoveFromForm.push_back(colPtr);
                        }
                }
            }

          _probPtr->removeVarsNotInProblemFromNonZeroRedCostVars(); /// otherwise removed columns may remain in
                                                                    /// problem's_nonZeroRedCostVars

          /// we need to do it due to tolerance problems of the CLP solver
          if (param().solverName() == "CLP_SOLVER")
            _probPtr->removeVarsNotInProblemFromPrimalSolution();

          /// remove columns without enumerated status and with large reduced cost
          if (!varsToRemoveFromForm.empty())
            {
              if (!justEnumeratedSubproblems.empty())
                {
                  /// sol. may involve cols which will now be removed from the formulation
                  /// we should not reset dual solution (it may be used in the enumeration with the false gap)
                  _probPtr->resetSolution();
                  repeatColumnGeneration = true;
                }
              _probPtr->delVarsSimplyInForm(varsToRemoveFromForm);
              _probPtr->removeUnusedDynamicVarsFromMemory();
              varsToRemoveFromForm.clear();
            }

          /// insert newly generated columns with enumerated status
          for (std::set<ColGenSpConf *,ColGenSpConfCmpById>::iterator cgSpConfPtrIt = justEnumeratedSubproblems.begin();
               cgSpConfPtrIt != justEnumeratedSubproblems.end(); ++cgSpConfPtrIt)
            (*cgSpConfPtrIt)->insertAllColumnsInMaster();
          if (!justEnumeratedSubproblems.empty())
            recordColInForm();

          if (printL(0))
          {
            if (numEnumColsRemoved > 0)
              std::cout << "Removed " << numEnumColsRemoved
                        << " columns (not in the enumerated set) from the formulation" << std::endl;
            if (numColsRegeneratedWithEnumFlag > 0)
              std::cout << "Regenerated " << numColsRegeneratedWithEnumFlag
                        << " columns with the 'enumerated' flag" << std::endl;
          }
        }
        if ((enumerationMode == 0) && addEnumColumnsToMaster(param().RCSPmaxNumOfEnumSolutionsForMIP()))
        {
          _pricingSolverCutsMessageId = PricingSolverCutsMessage::stopCutGeneration;
          _masterConverged = false; /// as solution is changed
          _subProbSolutionsEnumeratedToMIP = true;
          repeatColumnGeneration = false;
        }
        if ((enumerationMode == 1) && runEnumSolBasedUserHeuristicFunctor())
        {
            repeatColumnGeneration = false;
        }
        if ((enumerationMode == 1)
            && addEnumColumnsWithSmallestRedCostInMaster(param().MaxNumEnumSolsInRestrictedMasterIpHeur()))
        {
            repeatColumnGeneration = false;
        }
    }
  if (printL(0) && (param().ReducedCostFixingThreshold() > 0) && !currentGapIsBelowThreshold && (enumerationMode == 0))
    std::cout << "    Full reduced cost fixing is not called (gap ratio is "
              << currentGap / _latestReducedCostFixingGap << ")" << std::endl;

  bapcodInit().statistics().incrTimer("bcTimeRedCostFixAndEnum", startRedCostTime.getElapsedTime_dbl());
  if (repeatColumnGeneration && (_colGenStabilizationPtr != NULL))
    {
      StabilizationInfo * stabInfoPtr = _colGenStabilizationPtr->recordStabilizationInfo();
      _colGenStabilizationPtr->setDownStab();
      _colGenStabilizationPtr->setupStab(stabInfoPtr, algIncLpDualBound(),
                                         param().MaxNbOfStagesInColGenProcedure() - 1,
                                         _currentNodePtr->depth() + 1);
      delete stabInfoPtr;
    }

    
  return repeatColumnGeneration;
}

NodeEvalInfo * Alg4EvalBySimplexBasedColGen::recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr)
{
  /// copy fractional solution, the copy will be used in branching and heuristic algorithms
  _currentNodePtr->recordPrimalSol(_probPtr->inPrimalLpSol());

  LpBasisRecord * lpBasisPtr = new LpBasisRecord(std::string("BasisN") + _currentNodePtr->ref());
  _probPtr->retrieveBasis(lpBasisPtr);

  StabilizationInfo * stabInfo = NULL;
  if (_colGenStabilizationPtr != NULL)
    stabInfo = _colGenStabilizationPtr->recordStabilizationInfo();

  ColGenEvalInfo * colGenSolverInfoPtr = NULL;
  if (nodeEvalInfoPtr != NULL)
    {
      colGenSolverInfoPtr = dynamic_cast<ColGenEvalInfo *>(nodeEvalInfoPtr);

      bapcodInit().require(colGenSolverInfoPtr != NULL,
                           "BaPCod error: nodeEvalInfoPtr passed to ColGenEvalAlg::recordNodeEvalInfo"
                           " is not of type ColGenEvalInfo");

      colGenSolverInfoPtr->masterLpBasisPtr = lpBasisPtr;
      colGenSolverInfoPtr->stabilizationInfoPtr = stabInfo;
      colGenSolverInfoPtr->latestReducedCostFixingGap = _latestReducedCostFixingGap;
    }
  else
    colGenSolverInfoPtr = new ColGenEvalInfo(stabInfo, lpBasisPtr, _latestReducedCostFixingGap);

  return Alg4EvalOfNode::recordNodeEvalInfo(globalTreatOrder, colGenSolverInfoPtr);
}

void Alg4EvalBySimplexBasedColGen::setDownAlgo()
{
  for (std::vector<ColGenSpConf*>::const_iterator spIt = _masterCommons.colGenSubProbConfPts().begin();
      spIt != _masterCommons.colGenSubProbConfPts().end(); ++spIt)
    {
      (*spIt)->spConfHasClassInducedSpVarBounds(false);
      if (!(*spIt)->treeOfColClasses().empty())
        {
          /// Ruslan : do we need to save col gen sp trees of col classes?
          /// if we do not save, order of classes is different from the previous implementation
          /// Is this order important for branching?
          _currentNodePtr->cgSpConfTreeOfColClassesMap()[*spIt] = (*spIt)->treeOfColClasses();
          (*spIt)->treeOfColClasses().clear();
        }
    }

  _isLastSpcPtInitialized = false;

  if (_colGenStabilizationPtr != NULL)
    _colGenStabilizationPtr->setDownStab();

  Alg4EvalOfNode::setDownAlgo();
}


bool Alg4EvalBySimplexBasedColGen::checkIfSubProbSolutionsEnumeratedToMIP()
{
  if (!_subProbSolutionsEnumeratedToMIP && (_doRedCostFixingAndEnumeration == 1)
      && addEnumColumnsToMaster(param().RCSPmaxNumOfEnumSolsForEndOfNodeMIP()))
  {
    _subProbSolutionsEnumeratedToMIP = true;
    _pricingSolverCutsMessageId = PricingSolverCutsMessage::stopCutGeneration;
  }
  return _subProbSolutionsEnumeratedToMIP;
}
