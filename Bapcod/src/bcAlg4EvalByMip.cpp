/**
 *
 * This file bcAlg4EvalByMip.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcFormC.hpp"
#include "bcMasterConfC.hpp"

int Alg4EvalByMip::solveRestrictedMastIP() //This name is probably not the most adequate at it is not only for Restricted Mast.
{
  Time start;

  if (!progStatus().doRun())
    return (false);
  int maxLevelOfRestriction(0);
  int solverReturnStatus(_probPtr->solveProb(maxLevelOfRestriction, ' ', printL(2)));

  bapcodInit().statistics().incrTimer("bcTimeSolveRM", start.getElapsedTime_dbl());

  bapcodInit().statistics().incrCounter("bcCountMastIpSol");

  int status = *(_probPtr->probStatus().begin());
  if (printL(0))
    std::cout << "Solution status = " << _probPtr->probStatus().stat2string(status) << std::endl;

  if (_probPtr->probStatus().count(SolutionStatus::Optimum))
    _algCurLpDualBound = _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());
  else if (_probPtr->probStatus().count(SolutionStatus::Infeasible))
    /// if MIP is infeasible, this does not mean that the node is infeasible, as cutoff value was set
    _algCurLpDualBound = _algCurLpPrimalBound = algIncIpPrimalBound();
  else if (_probPtr->probStatus().count(SolutionStatus::PrimalFeasSolFound))
    _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());

  if (!printL(0) && printL(-1)
      && (_probPtr->probStatus().count(SolutionStatus::Optimum)
          || _probPtr->probStatus().count(SolutionStatus::PrimalFeasSolFound)))
  {
    std::cout << "MIP found solution with value " << _algCurLpPrimalBound << std::endl;
  }

    updateAlgDualBounds();
    updateAlgPrimalLpBounds();

  return solverReturnStatus;
}

bool Alg4EvalByMip::setupAlgo(Node * nodePtr)
{
  if (Alg4EvalOfNode::setupAlgo(nodePtr))
    return true;

  /// added by Issam for restricted MIP heuristic. In fact even when masterSol is MIP the cplex formulation turn into
  /// an LP formulation so we need to change it back to be an MIP formulation.
  if (param().masterSolMode().status() == SolutionMethod::mipSolver)
    {
      MIPform * mipFormulation = dynamic_cast<MIPform*>(_probPtr->formulationPtr());
      Bound incumbentValue(_currentNodePtr->nodeIncIpPrimalBound() - _probPtr->partialSolutionValue(),
                           _masterCommons.objStatus());
      if (_masterCommons.objStatus() == BcObjStatus::minInt)
        incumbentValue = Bound(floor((double)(incumbentValue - Double::precision)) + Double::precision,
                               _masterCommons.objStatus());

      mipFormulation->resetMIPpartOfFormulation(incumbentValue, _paramExactSolution);

      /// we calculate the remaining time
      long remainingTimeInTick = param().GlobalTimeLimitInTick() - bapcodInit().startTime().getElapsedTime();
      double timeToGiveToMIP = (((double)remainingTimeInTick) / 100 ) * 0.98;

      if (timeToGiveToMIP > _paramMaxTime)
          timeToGiveToMIP = _paramMaxTime;

      mipFormulation->setTimeLimit(timeToGiveToMIP);
    }
  else
    {
      std::cerr << " BaPCod ERROR when solving master as MIP: masterSolMode must be of type SolutionMethod::mipSolver"
                << std::endl;
      exit(1);
    }

  std::list<Variable*> varsToDeleteFromForm;

  for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'a');
       varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'a');
       ++varPtrIt)
    {
      if (printL(5))
        std::cout << "var to be deleted from form : " << (*varPtrIt)->name() << std::endl;
      varsToDeleteFromForm.push_back(*varPtrIt);
    }

  _probPtr->delVarSet(varsToDeleteFromForm, 1, 2);

  if (printL(3))
    {
      std::cout << "Current MIP : " << std::endl;
      _probPtr->printForm();
    }

  /// added by Boris : change requiredStatus of Problem to prevent infeasible MIPs from causing application termination
  ((MipProblem*) _probPtr)->getMIPRequiredStatus(_problemSolutionRequiredStatus);
  int tmp[] = { 0, 1, 2, 3, 4, 5 };
  std::vector<int> statusVect(tmp, tmp + 5);
  SolutionStatus status(statusVect);
  if (printL(3))
    std::cout << "BaPCod info: Changing required solution status to " << status << std::endl;
  ((MipProblem*) _probPtr)->setMIPRequiredStatus(status);

  static_cast<LPform*>(_probPtr->formulationPtr())->interfacePtr()->setScreenOutput(printL(0));

    MasterConf * masterProbConf = NULL;
    for (auto cutGenConstrPtr : _masterCommons.candidateCutGenericConstr())
        if (cutGenConstrPtr->type() == 'C')
        {
            /// if there are core cuts, the lazy constraints callback should be defined
            /// this should be done after deleting artificial variables, as we need the exact number of
            /// variables in the formulation
            MIPform * mipFormulation = dynamic_cast<MIPform*>(_probPtr->formulationPtr());
            MasterConf * masterConfPtr = dynamic_cast<MasterConf *>(nodePtr->probConfPtr());
            mipFormulation->setLazyConstraintsCallback(masterConfPtr);
            break;
        }

    return false;
}

void Alg4EvalByMip::setDownAlgo()
{
  static_cast<LPform*>(_probPtr->formulationPtr())->interfacePtr()->setScreenOutput(printL(1));

  std::list<Variable*> varsToAddToForm;

  for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Inactive, 'a');
       varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Inactive, 'a');
       ++varPtrIt)
    {
      if (printL(5))
        std::cout << "var to be added to form : " << (*varPtrIt)->name() << std::endl;
      varsToAddToForm.push_back(*varPtrIt);
    }

  _probPtr->addVarSet(varsToAddToForm, 1, 2);

  if (printL(3))
    std::cout << "Changing required solution status back to " << _problemSolutionRequiredStatus << std::endl;
  ((MipProblem*) _probPtr)->setMIPRequiredStatus(_problemSolutionRequiredStatus);
    
  MIPform * mipFormulation = dynamic_cast<MIPform*>(_probPtr->formulationPtr());
  mipFormulation->resetAfterMIP();

  /// we restore the maximum time
  mipFormulation->setTimeLimit(1e+75);

  Alg4EvalOfNode::setDownAlgo();
}

bool Alg4EvalByMip::eval()
{
#ifdef BC_MORE_TIMERS
  Time start;
#endif

  if (!progStatus().doRun())
    return false;

  bool coreCutAdded = false;
  do {
    coreCutAdded = false;
    if (solveRestrictedMastIP() > 0)
    {
      /// Test if LP solution is integer and,
      if (checkIfCurSolIsInteger())
      {
        /**
         * If so and no more delayed constraint need to be added,
         * memorise it if it beats the incumbent
         */
        if (addCutToMaster('C'))
        {
          coreCutAdded = true;
          resetAlgIncLpPrimalBound(_masterCommons.objStatus());
        }
        else
        {
          if (printL(0) && (param().printMasterPrimalSols() == 3))
            _probPtr->printDetailedPrimalSol();
          updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
          return false;
        }
      }
      else
      {
        if (printL(5))
          std::cout << "restrictedMasterIpHeuristic :  Problem sol non integer " << std::endl;
        return false;
      }
      _solIsMasterLpFeasible = true;
    }
    else
    {
      /// if MIP is infeasible, this does not mean that the node is infeasible, as cutoff value was set
      /// thus, we do nothing
    }
  } while (_paramRepeatIfCoreCutAddedAfterMIPHeur && coreCutAdded);

#ifdef BC_MORE_TIMERS
  bapcodInit().statistics().incrTimer("bcTimeRestictMastIpHeur", start.getElapsedTime_dbl());
#endif

  return !_solIsMasterLpFeasible;
}

