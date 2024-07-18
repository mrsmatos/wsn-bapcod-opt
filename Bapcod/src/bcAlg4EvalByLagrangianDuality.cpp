/**
 *
 * This file bcAlg4EvalByLagrangianDuality.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalByLagrangianDuality.hpp"
#include "bcStabilizationColgen.hpp"
#include "bcColGenSpConfC.hpp"

#define PRICING_STARTEGY_PRINTL 2

using namespace std;

Alg4EvalByLagrangianDuality::Alg4EvalByLagrangianDuality(Problem* const probPtr,
                                                         MasterCommons4EvalAlg & masterCommons) :
  Alg4EvalOfNode(probPtr, masterCommons),
  _savedMasterLPTime(0),
  _savedPricingTime(0),
  _savedNbColumns(0),
  _maxNbOfPenaltyUpdates(param().ArtVarMaxNbOfPenaltyUpdates()),
  _maxNbOfCgIterations(param().MaxNbOfCgIterations()),
  _minLevelOfSpRestriction(0),
  _minNbOfCutRounds(0),
  _maxNbOfCutRounds(10000000),
  _doRedCostFixingAndEnumeration(1),
  _logPrintFrequency(10),
  _nonExactEvaluation(false),
  _lastSpcPt(std::vector<ColGenSpConf *>::iterator()),
  _isLastSpcPtInitialized(false),
  _colGenStabilizationPtr(NULL),
  _need2reintroduceArtVarInMast(false),
  _currentlyPerformingPhase1(false),
  _promisingOrOpenSpIndices(),
  _unpromisingSpIndices(),
  _pricingStrat(param().PricingStrategy)
{
    for (int i = 0; i< masterCommons.colGenSubProbConfPts().size(); i ++)
    {
        _promisingOrOpenSpIndices.push_back(i);
    }
}

Alg4EvalByLagrangianDuality::~Alg4EvalByLagrangianDuality()
{
  if (_colGenStabilizationPtr != NULL)
    delete _colGenStabilizationPtr;
}

void Alg4EvalByLagrangianDuality::printIntermediateStatistics(std::ostream & os,
                                                              const int & curMaxLevelOfSubProbRestriction,
                                                              const int & nbNumColumns,
                                                              const int & nbCgIterations, long & elapsedTime,
                                                              const bool & doNotCheckLogFrequency,
                                                              const bool & printIncumbentDual)
{
  _savedNbColumns += nbNumColumns;
  if (doNotCheckLogFrequency || ((_logPrintFrequency > 0) && (nbCgIterations % _logPrintFrequency == 0)))
    {
      long masterLPTime = bapcodInit().statistics().getTime("bcTimeMastMPsol");
      long pricingTime = bapcodInit().statistics().getTime("bcTimeCgSpOracle");
      elapsedTime = bapcodInit().startTime().getElapsedTime();
      if (param().colgeninfo_file() == "")
        {
          double stabAlphaValue = -1.0;
          if (_colGenStabilizationPtr != NULL)
            stabAlphaValue = _colGenStabilizationPtr->curAlphaValue();
          if (param().MaxNbOfStagesInColGenProcedure() > 1)
            os << "<DWph=" << curMaxLevelOfSubProbRestriction << "> ";
          os << "<it=" << std::setfill(' ') << std::setw(3) << nbCgIterations << "> "
             << "<et=" << std::fixed << std::setprecision(2) << elapsedTime / (double) 100 << "> "
             << "<Mt=" << std::setfill(' ') << std::setw(5) << (masterLPTime - _savedMasterLPTime) / 100.0 << "> "
             << "<Spt=" << std::setfill(' ') << std::setw(5) << (pricingTime - _savedPricingTime) / 100.0 << "> "
             << "<nCl=" << std::setw(3) << std::setprecision(0) << _savedNbColumns << "> ";
          if (printL(0) && (stabAlphaValue >= 0))
            os << std::setprecision(2) << "<al=" << stabAlphaValue << "> ";
          if (printIncumbentDual)
            os << std::setprecision(4) << "<DB=" << std::setw(10) << algIncLpDualBound() << "> ";
          else
            os << std::setprecision(4) << "<DB=" << std::setw(10) << algCurLpDualBound() << "> ";
          os  << "<Mlp=" << std::setw(10) << algCurLpPrimalBound() << "> ";
          os.unsetf(std::ios_base::floatfield);
          os << "<PB=" << std::setprecision(8) << algIncIpPrimalBound() << std::setprecision(6) <<  "> "
             << std::endl;
        }
      else if (_colGenStabilizationPtr != NULL)
        _colGenStabilizationPtr->printDetailedStabilizationInformation(os, nbCgIterations + 1,
                                                                       elapsedTime,algCurLpPrimalBound());
      _savedMasterLPTime = masterLPTime;
      _savedPricingTime = pricingTime;
      _savedNbColumns = 0;
    }
}

bool Alg4EvalByLagrangianDuality::earlyCGtermType1()
{
  if (param().runColGenUntilFullConvergence())
    return false;

  bool rvalue(false);

  if (_maxLevelOfSubProbRestriction < 1)
    {
      if (param().TerminateCgWhenRoundedDbCannotImprove)
        rvalue = gapSmallerThanTol(algIncIpDualBound(), algIncLpPrimalBound(), param());
      else
        rvalue = gapSmallerThanTol(algIncLpDualBound(), algIncLpPrimalBound(), param());
    }
  else
    {
      rvalue = gapSmallerThanTol(algIncStageLpDualBound(), algIncLpPrimalBound(), param());
    }

  if (rvalue)
    {
      if (printL(2))
        std::cout << "Alg4EvalByLagrangianDuality: early termination of type 1" << std::endl;
      bapcodInit().statistics().incrCounter("bcCountCgT1");
    }


  return (rvalue);
}

int Alg4EvalByLagrangianDuality::searchNegRedCostInactiveCol()
{
  int nbNewCol(0);
#ifdef BC_MORE_TIMERS
  Time start;
#endif

  ///  Needed to avoid adding in the set we are currently scanning
  std::list<Variable *> var2addInProb;

  Double bestRedCost(BapcodInfinity);
  Double currentRedCost(0);
  MastColumn * bestColPtr(NULL);
  MastColumn * colPtr(NULL);

  if (!_candidatColWithNegRedCost.empty())
    {
      std::list< MastColumn *> tmpCandidatColWithNegRedCost;

      for (std::list< MastColumn *>::iterator colPt = _candidatColWithNegRedCost.begin();
	       colPt != _candidatColWithNegRedCost.end(); ++colPt)
      {
          currentRedCost = (*colPt)->computeReducedCost();
          if (currentRedCost.negative(param().BapCodReducedCostTolerance))
          {
              if (printL(2))
                  std::cout << "found in candidatList an inactive  MastCol[" << (*colPt)->name()
                            << "] with neg red cost rc = " << (*colPt)->reducedCost() << std::endl;

              tmpCandidatColWithNegRedCost.push_back(*colPt);

              if (currentRedCost < bestRedCost)
              {
                  bestRedCost = currentRedCost;
                  bestColPtr = *colPt;
                  if (printL(2))
                      std::cout << "new incumbent  MastCol[" << bestColPtr->name()
                                << "] with neg red cost rc = " << bestColPtr->reducedCost() << std::endl;
              }
          }
      }
        _candidatColWithNegRedCost.clear();  // delete previous list, and reset

        _candidatColWithNegRedCost = tmpCandidatColWithNegRedCost;
    }


  if (bestColPtr == NULL) // no neg red vost foundin previous list
    {
      _candidatColWithNegRedCost.clear();  // delete previous list, and reset

      for (VarIndexManager::iterator it = _probPtr->probVarSet().begin(VcIndexStatus::Inactive, 'd');
	       it != _probPtr->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
          if ((*it)->isTypeOf(VcId::MastColumnMask))
          {
              colPtr = static_cast<MastColumn *>(*it);
              currentRedCost = colPtr->computeReducedCost();
              if (currentRedCost.negative(param().BapCodReducedCostTolerance))
              {
                  if (printL(2))
                      std::cout << "found in probVarSet an inactive  MastCol[" << colPtr->name()
                                << "] with neg red cost rc = " << colPtr->reducedCost()
                                << std::endl;

                  _candidatColWithNegRedCost.push_back(colPtr);

                  if (currentRedCost < bestRedCost)
                  {
                      bestRedCost = currentRedCost;
                      bestColPtr = colPtr;
                      if (printL(2))
                          std::cout << "new incumbent  MastCol[" << bestColPtr->name()
                                    << "] with neg red cost rc = " << bestColPtr->reducedCost() << std::endl;
                  }
              }

              if (printL(6))
                  std::cout << "MastCol[" << colPtr->name() << "]  rc = " << colPtr->reducedCost() << std::endl;
          }
    }

  if (bestColPtr != NULL)
    {
      if (printL(2))
	    std::cout << "ADDING  MastCol[" << bestColPtr->name() << "] with neg red cost rc = "
	              << bestColPtr->reducedCost() << std::endl;

      var2addInProb.push_back(bestColPtr);
      nbNewCol++;
    }

  _probPtr->addVarSet(var2addInProb, 1, 2);
#ifdef BC_MORE_TIMERS
  bapcodInit().statistics().incrTimer("bcTimeSearchColPool", start.getElapsedTime());
#endif

  return (nbNewCol);
}

void Alg4EvalByLagrangianDuality::compMastDualBoundContrib(Bound & mastDualBoundContrib)
{
  /// Check that there are no artificial col in the solution when updating dual bound
  int printLevel = 5;
  ConstrPtrSet::const_iterator constrPt;
  long long int scaleFactor = param().SafeDualBoundScaleFactor();

  if ((scaleFactor <= 0) && !param().SplitColIntoDissagregateSpVar()
      && ((_colGenStabilizationPtr == NULL) || !_colGenStabilizationPtr->isActive()))
  {
      /// if there is no stabilization, master dual bound contribution is equal to the current master LP value
      mastDualBoundContrib = Bound(totalObjVal(), _masterCommons.objStatus());
      return;
  }

  if (scaleFactor > 0)
    mastDualBoundContrib = Bound(ceil(_probPtr->partialSolutionValue()._val * scaleFactor), _masterCommons.objStatus());
  else
    mastDualBoundContrib = Bound(_probPtr->partialSolutionValue(), _masterCommons.objStatus());
  for (constrPt = _probPtr->inDualSol().begin(); constrPt != _probPtr->inDualSol().end(); constrPt++)
    if ((*constrPt)->inCurProb() && ((*constrPt)->type() != 'S')
        && (((*constrPt)->type() != 'X') || !param().SplitColIntoDissagregateSpVar))
      {
        if (scaleFactor > 0)
          mastDualBoundContrib -= ceil((*constrPt)->valOrSepPointVal()._val * (*constrPt)->curRhs()._val * scaleFactor);

        else
          mastDualBoundContrib -= (*constrPt)->valOrSepPointVal() * (*constrPt)->curRhs();

        if (printL(printLevel))
          std::cout << " Alg4EvalByLagrangianDuality::compDualBoundContrib() explicitly: constr "
                    << (*constrPt)->name() << " valOrSepPointVal " << (*constrPt)->valOrSepPointVal() << " rhs "
                    << (*constrPt)->curRhs() << " mastDualBoundContrib " << mastDualBoundContrib << std::endl;
      }

  /// in principle, pure master variables can be dynamic, TO DO : implement the support for this
  if (!_probPtr->probVarSet().empty(VcIndexStatus::Active, 's'))
    {
      if (scaleFactor > 0)
      {
        std::cerr << "BaPCod error : safe dual bound cannot be computed as the stabilization is active and "
                  << "pure master variables are present" << std::endl;
        exit(1);
      }
      /// then we retrieve reduced costs from LP solver
      /// Guillaume : if stabilization is active, reduced costs are retrieved
      /// after restricted master solving.
      
      _probPtr->retrieveRedCosts(); /// fast reduced cost retrieval for \pi_out

      if ((_colGenStabilizationPtr != NULL) && _colGenStabilizationPtr->solValueSmoothingIsActive())
        {
          std::map<VarConstr *, double> mastVarPtr2modRedCostMap;
          for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 's');
               varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); ++varPtrIt)
          {
            mastVarPtr2modRedCostMap.insert(std::make_pair(*varPtrIt, (*varPtrIt)->reducedCost()));
          }

          _colGenStabilizationPtr->changePureMasterVarsReducedCostUsingSepValues(mastVarPtr2modRedCostMap);

           for (std::map<VarConstr *, double>::iterator mapIt = mastVarPtr2modRedCostMap.begin();
                mapIt != mastVarPtr2modRedCostMap.end();  ++mapIt)
            {
              Double modRedCost(mapIt->second);
              /// mult which is calculated here is also needed by the automatic smoothing
              if (modRedCost.negative(param().BapCodReducedCostTolerance()))
                mapIt->first->mult(mapIt->first->curUb());
              else
                mapIt->first->mult(mapIt->first->curLb());

              mastDualBoundContrib += modRedCost * mapIt->first->mult();
              if (printL(printLevel))
                std::cout << " Alg4EvalByLagrangianDuality::compDualBoundContrib() pure master variable"
                          << " under dual price smoothing: var " << mapIt->first->name()
                          << "- rc = " << modRedCost << ", val = " << mapIt->first->val()
                          << ", bounds = [" << mapIt->first->curLb() << " , " << mapIt->first->curUb()
                          << ", contrib = " << modRedCost * mapIt->first->mult()
                          << "], mastDualBoundContrib = " << mastDualBoundContrib << std::endl;
            }
        }
      else
        {
          /// here dual bound is adjusted by pure master variables with non-zero reduced cost
          for (VarPtrSet::iterator varPt = _probPtr->nonZeroRedCostVars().begin();
               varPt != _probPtr->nonZeroRedCostVars().end(); ++varPt)
            if ((*varPt)->isTypeOf(VcId::InstMasterVarMask)) /// if it is a pure master variable
              {
                if ((*varPt)->reducedCost().negative(param().BapCodReducedCostTolerance()))
                  mastDualBoundContrib += (*varPt)->reducedCost() * (*varPt)->curUb();
                else
                  mastDualBoundContrib += (*varPt)->reducedCost() * (*varPt)->curLb();
              }
        }
    }

  if (scaleFactor > 0)
    mastDualBoundContrib._val /= scaleFactor;
}

struct SortConstraintsByViolation
{
    bool operator()(const Constraint * a, const Constraint * b) const
    {
        if (a->violation() < b->violation())
            return true;
        if (a->violation() > b->violation())
            return false;
        return (a->ref() < b->ref());
    }
};

void Alg4EvalByLagrangianDuality::cleanupRestrictedMastCuts()
{
    int numDynamicConstrs = _probPtr->probConstrSet().size(VcIndexStatus::Active, 'd');
    if (param().CutCleanupThreshold() <= 0 || numDynamicConstrs <= param().CutCleanupThreshold())
        return;

    int nbOfCutsToKeep = (std::min)((int)(numDynamicConstrs * param().CutCleanupRatio()), (int)param().CutCleanupThreshold());
    bool keepZeroSlackCuts = (param().CutCleanupRatio() > Double::precision);

    Solution solution(_probPtr->probConfPtr());
    solution.includeVarSet(_probPtr->inPrimalLpSol());

    LpBasisRecord * lpBasisPtr = new LpBasisRecord(std::string("tmp"));
    _probPtr->retrieveBasis(lpBasisPtr, false, true);
    std::vector<Constraint *> vectOfCandConstrPtr;
    int nbOfCutsCannotDelete = 0;
    int nbOfBasisCutsWithZeroSlack = 0;

    for (ConstrIndexManager::iterator constrPt = _probPtr->probConstrSet().begin(VcIndexStatus::Active, 'd');
         constrPt != _probPtr->probConstrSet().end(VcIndexStatus::Active, 'd'); constrPt++)
        if (!(*constrPt)->isTypeOf(VcId::BranchingConstrBaseTypeMask)
            && (*constrPt)->isTypeOf(VcId::InstMasterConstrMask))
        {
            /// all master instantiated dynamic constraints which are not branching ones are considered to be cuts
            if ((*constrPt)->isInBasis())
            {
                Double lhs = (*constrPt)->computeLhs(_probPtr->inPrimalLpSol());
                Double weightedSlack = 0.0;
                if ((*constrPt)->sense() == 'G')
                    weightedSlack = (lhs - (*constrPt)->curRhs());
                else if ((*constrPt)->sense() == 'L')
                    weightedSlack = (*constrPt)->curRhs() - lhs;
                if (std::abs((double)(*constrPt)->curRhs()) > 1.0)
                    weightedSlack /= (*constrPt)->curRhs();

                (*constrPt)->violation(weightedSlack); /// we use _violation field to store the slack of the constraint

//                std::cout << "Basis cut " << (*constrPt)->name() << " with violation = " << (*constrPt)->violation()
//                          << ", val = " << (*constrPt)->val() << ", smoothVal = " << (*constrPt)->valOrSepPointVal()
//                          << std::endl;

                if (keepZeroSlackCuts && ((*constrPt)->violation() < Double::precision))
                {
                    nbOfBasisCutsWithZeroSlack += 1;
                    nbOfCutsCannotDelete += 1;
                }
                else
                {
                    vectOfCandConstrPtr.push_back(*constrPt);
                }

                if (!(*constrPt)->val().isZero())
                    std::cerr << "BaPCod WARNING : cut " << (*constrPt)->name()
                              << " is being cleaned but its dual value is " << (*constrPt)->val() << std::endl;
            }
            else
            {
                nbOfCutsCannotDelete += 1;
//                std::cout << "Non-basis cut " << (*constrPt)->name() << " with val = " << (*constrPt)->val()
//                          << ", smoothVal = " << (*constrPt)->valOrSepPointVal() << std::endl;
            }
        }

    lpBasisPtr->clear(false, true);
    delete lpBasisPtr;

//    std::cout << "NumDynamicConstrs = " << numDynamicConstrs << ", nb of cut candidates = " << vectOfCandConstrPtr.size()
//              << ", nbOfCutsCannotDelete = " << nbOfCutsCannotDelete << ", nbOfCutsToKeep = " << nbOfCutsToKeep
//              << ", nbOfBasisCutsWithZeroSlack = " << nbOfBasisCutsWithZeroSlack << std::endl;

    if (nbOfCutsToKeep - nbOfCutsCannotDelete > 0)
    {
        if (vectOfCandConstrPtr.size() <= nbOfCutsToKeep - nbOfCutsCannotDelete)
            return;

        stable_sort(vectOfCandConstrPtr.begin(), vectOfCandConstrPtr.end(), SortConstraintsByViolation());

        for (auto constrPtrIt = vectOfCandConstrPtr.begin() + nbOfCutsToKeep - nbOfCutsCannotDelete;
             constrPtrIt != vectOfCandConstrPtr.end(); ++constrPtrIt)
            (*constrPtrIt)->desactivateConstraint(VcIndexStatus::Unsuitable);
    }
    else
    {
        for (auto * constrPtr : vectOfCandConstrPtr)
            constrPtr->desactivateConstraint(VcIndexStatus::Unsuitable);
    }

    _probPtr->delConstrInForm();
    _probPtr->delVarInForm();
    _probPtr->removeVarsNotInProblemFromNonZeroRedCostVars(); /// otherwise removed vars may remain in _nonZeroRedCostVars
    _probPtr->removeUnusedDynamicConstrsFromMemory();

    if (printL(1))
        std::cout << "Clean-up is done, remains "
                  << _probPtr->probConstrSet().size(VcIndexStatus::Unsuitable, 'd')
                  << " unsuitable dyn. constrs and "
                  << _probPtr->probConstrSet().size(VcIndexStatus::Active, 'd')
                  << " active dyn. constrs" << std::endl;
}

void Alg4EvalByLagrangianDuality::cleanupRestrictedMastColumns(const int& nbCgIterations)
{
  if (param().ColumnCleanupThreshold() <= 0)
      return;

  _candidatColWithNegRedCost.clear();

  int maxNumberOfColumns = param().ColumnCleanupThreshold();

  int numDynamicVars = _probPtr->probVarSet().size(VcIndexStatus::Active, 'd');

  if (numDynamicVars <= maxNumberOfColumns)
    return;

  if (printL(1))
      std::cout << "Master columns clean-up started with " << numDynamicVars << " dynamic variables in the master "
                << std::endl;

  std::vector<MastColumn *> vectOfCandColumnPtr;

  /// we retrieve reduced costs only when it is needed, unlike primal and dual solutions
  _probPtr->retrieveRedCosts(); /// fast reduced cost retrieval

  LpBasisRecord * lpBasisPtr = new LpBasisRecord(std::string("tmp"));
  _probPtr->retrieveBasis(lpBasisPtr, true, false);

  if (param().solverName() == "CLP_SOLVER")
  {
      /// we artificially enter variables in the primal solution into basis
      /// (so that they are not removed by the cleanup procedure)
      /// non-basic columns may have a non-zero value due to the fact that the solution of CLP solver
      /// does not always respect its tolerance
      for (auto * varPtr : _probPtr->inPrimalLpSol())
      {
          if (!varPtr->isInBasis())
          {
              if (printL(1))
                std::cout << "WARNING: variable " << varPtr->name() << " in not in basis but in the primal sol with  "
                          << " value " << varPtr->val() << std::endl;
              varPtr->isInBasis(true);
          }
      }
  }

  int numberOfZeroReducedCostColumns = 0;
  for (VarIndexManager::iterator varPt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
       varPt != _probPtr->probVarSet().end(VcIndexStatus::Active, 'd'); varPt++)
    if ((*varPt)->isTypeOf(VcId::MastColumnMask))
      {
        if (!(*varPt)->isInBasis())
            vectOfCandColumnPtr.push_back(static_cast<MastColumn *>(*varPt));
        else
            numberOfZeroReducedCostColumns += 1;
      }

  lpBasisPtr->clear(true, false);
  delete lpBasisPtr;

    if (param().solverName() == "CLP_SOLVER")
    {
        /// clear variables artificially in the basis
        for (auto * varPtr : _probPtr->inPrimalLpSol())
            if (varPtr->isInBasis())
                varPtr->isInBasis(false);
    }


    if (static_cast<int>(vectOfCandColumnPtr.size()) + numberOfZeroReducedCostColumns <= maxNumberOfColumns)
    return;

  int numberOfColumnsToKeep = numDynamicVars * param().ColumnCleanupRatio();

  if (numberOfZeroReducedCostColumns >= numberOfColumnsToKeep)
    return;

  int nbOfColumnInMaster = _probPtr->probVarSet().size(VcIndexStatus::Active,'s')
                           + _probPtr->probVarSet().size(VcIndexStatus::Active, 'd');

  if (nbOfColumnInMaster > bapcodInit().statistics().getCounter("bcMaxNbColInMastLp"))
    bapcodInit().statistics().setCounter("bcMaxNbColInMastLp", nbOfColumnInMaster);

  stable_sort(vectOfCandColumnPtr.begin(), vectOfCandColumnPtr.end(), SortMastColumnPerNonDecreasingRedCost());

  std::vector<MastColumn *>::iterator colPtrIt = vectOfCandColumnPtr.begin();
  colPtrIt += (numberOfColumnsToKeep - numberOfZeroReducedCostColumns);

  std::list<Variable *> varsToRemoveFromForm;
  for (; colPtrIt != vectOfCandColumnPtr.end(); ++colPtrIt)
    {
      if (param().Search4NegRedCostColInInactivePool())
        _probPtr->probVarSet().insert(*colPtrIt, VcIndexStatus::Inactive);
      else
        _probPtr->probVarSet().insert(*colPtrIt, VcIndexStatus::Unsuitable);
      (*colPtrIt)->desactivate();
      varsToRemoveFromForm.push_back(*colPtrIt);
    }
  _probPtr->delVarsSimplyInForm(varsToRemoveFromForm);
  _probPtr->removeVarsNotInProblemFromNonZeroRedCostVars(); /// otherwise removed cols may remain in _nonZeroRedCostVars
  _probPtr->removeUnusedDynamicVarsFromMemory();

  if (printL(1))
    std::cout << "Clean-up is done, remains "
              << _probPtr->probVarSet().size(VcIndexStatus::Unsuitable, 'd')
              << " unsuitable columns and "
              << _probPtr->probVarSet().size(VcIndexStatus::Active, 'd')
              << " active columns" << std::endl;

  if ((param().colgeninfo_file() != "") && (_logPrintFrequency > 0))
        std::cout << "c";

  /// we need also to keep the list of inactive columns of a reasonable size
  /// TO DO : implement here definitive removal of some inactive columns

  return;
}

void Alg4EvalByLagrangianDuality::recordColInForm()
{
  if (printL(6))
    std::cout << "MasterConf::recordColInForm()" << std::endl;

  _probPtr->addVarInForm();

  return;
}

void Alg4EvalByLagrangianDuality::updatePricingSolverCutsMessageId()
{
  _pricingSolverCutsMessageId = PricingSolverCutsMessage::noMessage;
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
      /// TO DO : if a pricing solver asks for a rollback, we do not need to run pricing solvers for other subproblems
      int thisPricingSolverMessageId = (*spcPt)->probPtr()->getMessageIdToCutGeneration();
      if (thisPricingSolverMessageId == PricingSolverCutsMessage::interruptSolution)
      {
        _pricingSolverCutsMessageId = PricingSolverCutsMessage::interruptSolution;
      }
      else if (thisPricingSolverMessageId == PricingSolverCutsMessage::doCutsRollback)
        {
          if ((*spcPt)->getRollbackPointSavedStatus())
            {
              _pricingSolverCutsMessageId = PricingSolverCutsMessage::doCutsRollback;
            }
          else if (printL(-1))
            {
              std::cerr << "BaPCod WARNING : pricing problem sent 'doCutsRollback' message, "
                        << "but 'rollback point' is not saved, so we interrupt the solution " << std::endl;
              _pricingSolverCutsMessageId = PricingSolverCutsMessage::interruptSolution;
            }
        }
      else if ((thisPricingSolverMessageId == PricingSolverCutsMessage::stopCutGeneration)
               && (_pricingSolverCutsMessageId != PricingSolverCutsMessage::doCutsRollback))
        {
          _pricingSolverCutsMessageId = PricingSolverCutsMessage::stopCutGeneration;
        }
    }
}

int Alg4EvalByLagrangianDuality::genNewColumns(int & maxLevelOfSubProbRestriction, int & allAddedColumns,
                                               int & addedNegRedCostColumns)
{
  if (!_isLastSpcPtInitialized)
    {
      _lastSpcPt = _masterCommons.colGenSubProbConfPts().begin();
      _isLastSpcPtInitialized = true;
    }

  bapcodInit().statistics().incrCounter("bcCountCg");
  Time start;

  if (param().CyclicSpScanning())
    {
      if (printL(6))
        std::cout << "ColGenSolver::genNewColumns(): lastSpcPt " << (*_lastSpcPt)->id() << std::endl;

      bool loop(false);
      int nbNewNegRedCostCol = 0;
      std::vector<ColGenSpConf *>::const_iterator spcPt = _lastSpcPt;
      do
        {
          nbNewNegRedCostCol -= addedNegRedCostColumns;
          int status = (*spcPt)->genNewColumns(_currentlyPerformingPhase1, maxLevelOfSubProbRestriction,
                                               allAddedColumns, addedNegRedCostColumns);
          nbNewNegRedCostCol += addedNegRedCostColumns;

          ///  In case number of subproblem solutions is already used up, goto the next SP
          if (status == -2)
            continue;

          ///  In case subproblem infeasibility results in master infeasibility
          if (status < 0)
            {
              recordColInForm();
              return status;
            }

          if (printL(5))
            std::cout << "ColGenSolver::genNewColumns(): in solving  " << (*spcPt)->id()
                      << " found " << nbNewNegRedCostCol << " neg. red. cost. columns " << std::endl;

          if (++spcPt == _masterCommons.colGenSubProbConfPts().end())
            {
              spcPt = _masterCommons.colGenSubProbConfPts().begin();
              loop = true;
            }
        }
      while ( (spcPt != _lastSpcPt) && (nbNewNegRedCostCol <= 0) );
      
      ///  Found no column
      if ((spcPt == _lastSpcPt) && (nbNewNegRedCostCol <= 0))
        {
          recordColInForm();
          return 0;
        }

      /// Record empty solutions and dual bound contrib for Non Visited Sp
      if (loop)
        {
          for (std::vector<ColGenSpConf *>::const_iterator spcPtNV = spcPt; spcPtNV != _lastSpcPt; ++spcPtNV)
            (*spcPtNV)->computeSpDualBoundContrib();
        }
      else // (_lastSpcPt <= spcPt)
        {
          std::vector<ColGenSpConf *>::const_iterator spcPtNV = _masterCommons.colGenSubProbConfPts().begin();
          for (; spcPtNV != _lastSpcPt; ++spcPtNV)
            (*spcPtNV)->computeSpDualBoundContrib();

          for (spcPtNV = spcPt; spcPtNV != _masterCommons.colGenSubProbConfPts().end(); ++spcPtNV)
            (*spcPtNV)->computeSpDualBoundContrib();
        }
      _lastSpcPt = spcPt;
      if (printL(5))
        std::cout << "ColGenSolver::genNewColumns(): new lastSpcPt " << (*_lastSpcPt)->id() << std::endl;
    }
  else if (param().PricingStrategy() == 0)
    {
      for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
  	       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
        {
          int status = (*spcPt)->genNewColumns(_currentlyPerformingPhase1, maxLevelOfSubProbRestriction,
                                               allAddedColumns, addedNegRedCostColumns);

          ///  In case number of subproblem solution is already sued up, goto the next SP
          if (status == -2)
            continue;

          ///  In case subproblem infeasibility results in  master infeasibility
          if (status < 0)
            {
              recordColInForm();
              return status;
            }
        }
    }
  else
    {
      int pricingStrategyPrintL = 5;

      int maxNbPromisingSpFound = param().MaxNbPromisingSpFound();
      int maxNbUnpromisingSpFound = param().MaxNbUnpromisingSpFound();
      int nbNewNegRedCostCol = 0;
      int count = 0;

      std::list<int> nextPromisingSpIndices;
      std::list<int> nextUnpromisingSpIndices;

      std::list<int> nonPricedPreviouslyPromisingOrOpenSpIndices;
      std::list<int> nonPricedPreviouslyUnpromisingSpIndices;

      for (; count < _masterCommons.colGenSubProbConfPts().size(); count++)
        {
          int idx = -1;

          if (!_promisingOrOpenSpIndices.empty())
	        {
	          idx = _promisingOrOpenSpIndices.front();
	          _promisingOrOpenSpIndices.pop_front();
            }
          else
	        {
	          idx = _unpromisingSpIndices.front();
	          _unpromisingSpIndices.pop_front();
	        }
          if (printL(pricingStrategyPrintL))
            std::cout << "ColGenSolver::genNewColumns(): --------Considering sp " << idx << std::endl;

          nbNewNegRedCostCol -= addedNegRedCostColumns;
          int status = _masterCommons.colGenSubProbConfPts()[idx]->genNewColumns(_currentlyPerformingPhase1,
                                                                                 maxLevelOfSubProbRestriction,
                                                                                 allAddedColumns,
                                                                                 addedNegRedCostColumns);
          nbNewNegRedCostCol += addedNegRedCostColumns;

          ///  In case number of subproblem solution is already sued up, goto the next SP
          if (status == -2)
	        {
	          nextUnpromisingSpIndices.push_back(idx);
	          continue;
	        }

          ///  In case subproblem infeasibility results in  master infeasibility
          if (status < 0)
	        {
	          recordColInForm();
	          nextUnpromisingSpIndices.push_back(idx);
	          return status;
	        }

          if (status > 0)
	        {
   	          nextPromisingSpIndices.push_back(idx);
  	          if (--maxNbPromisingSpFound <= 0)
		        {
		          count++;
		          break;
		        }
	        }
          else
	        {
	          if (printL(pricingStrategyPrintL))
                std::cout << "ColGenSolver::genNewColumns(): No new colum ADDED by sp " << idx << std::endl;
	          nextUnpromisingSpIndices.push_back(idx);

	          if ((nbNewNegRedCostCol > 0) && (--maxNbUnpromisingSpFound <= 0))
		        {
		          count++;
		          break;
                }
	        }
        }

      /// Note: count is not really needed here. but we keep it since this double checking doesn't cost much.
      for (; count < _masterCommons.colGenSubProbConfPts().size(); ++count)
        {
          if (!_promisingOrOpenSpIndices.empty())
	        {
	          int idx = _promisingOrOpenSpIndices.front();
              _promisingOrOpenSpIndices.pop_front();
	          if (printL(pricingStrategyPrintL))
                std::cout << "ColGenSolver::genNewColumns(): contrib. open or prom sp " << idx << std::endl;
	          _masterCommons.colGenSubProbConfPts()[idx]->updateConf(_currentlyPerformingPhase1);
              _masterCommons.colGenSubProbConfPts()[idx]->computeSpDualBoundContrib();
	          nonPricedPreviouslyPromisingOrOpenSpIndices.push_back(idx);
	        }
          else
	        {
	          int idx = _unpromisingSpIndices.front();
	          _unpromisingSpIndices.pop_front();
	          if (printL(pricingStrategyPrintL))
                std::cout << "ColGenSolver::genNewColumns(): contrib. nonPromising sp " << idx << std::endl;
	          _masterCommons.colGenSubProbConfPts()[idx]->updateConf(_currentlyPerformingPhase1);
	          _masterCommons.colGenSubProbConfPts()[idx]->computeSpDualBoundContrib();
	          nonPricedPreviouslyUnpromisingSpIndices.push_back(idx);
	        }
        }

      if (!_promisingOrOpenSpIndices.empty() && !_unpromisingSpIndices.empty())
        {
	      std::cerr << "Alg4EvalByLagrangianDuality::genNewColumns "
		            << " Error : some SP was not evaluated" << std::endl;
          exit(1);
        }

      if (_pricingStrat.selectedRule() == PricingStrategy::Diversification)
        {
          _promisingOrOpenSpIndices.splice(_promisingOrOpenSpIndices.end(),
                                           nonPricedPreviouslyPromisingOrOpenSpIndices);
          _promisingOrOpenSpIndices.splice(_promisingOrOpenSpIndices.end(), nextPromisingSpIndices);
          _unpromisingSpIndices.splice(_unpromisingSpIndices.end(), nonPricedPreviouslyUnpromisingSpIndices);
          _unpromisingSpIndices.splice(_unpromisingSpIndices.end(), nextUnpromisingSpIndices);
        }
      else if (_pricingStrat.selectedRule() == PricingStrategy::Intensification)
        {
          _promisingOrOpenSpIndices.splice(_promisingOrOpenSpIndices.end(), nextPromisingSpIndices);
          _promisingOrOpenSpIndices.splice(_promisingOrOpenSpIndices.end(),
                                           nonPricedPreviouslyPromisingOrOpenSpIndices);
          _unpromisingSpIndices.splice(_unpromisingSpIndices.end(), nonPricedPreviouslyUnpromisingSpIndices);
          _unpromisingSpIndices.splice(_unpromisingSpIndices.end(), nextUnpromisingSpIndices);
        }
    }
  
  recordColInForm();

  updatePricingSolverCutsMessageId();

  double delta = start.getElapsedTime_dbl();
  bapcodInit().statistics().incrTimer("bcTimeCgSpOracle",  delta);
  if (maxLevelOfSubProbRestriction > 0)
    bapcodInit().statistics().incrTimer("bcTimeHeurSpOracle",  delta);

  return 0;
}

void Alg4EvalByLagrangianDuality::updateLagrangianDualBound(bool updateDualBnd)
{
  int localPrintLevel = 1;

  Bound mastCurLagrangianBoundVal(_masterCommons.objStatus());
  compMastDualBoundContrib(mastCurLagrangianBoundVal);

  /// Subproblem contributions
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    {
      mastCurLagrangianBoundVal += (*spcPt)->dualBoundContrib();
      if (printL(localPrintLevel))
        std::cout << " master dual bound: contrib of SP[" << (*spcPt)->name()
		          << "]  = " << (*spcPt)->dualBoundContrib()
                  << " mastCurLagrangianBoundVal = " << mastCurLagrangianBoundVal << std::endl;
    }

  /// realization of the "do not split while in the head-in" option
  if (param().SplitColIntoDissagregateSpVar() && !param().SplitColIntoDissagregateSpVarInHeadIn
      && (mastCurLagrangianBoundVal > 0.0))
    for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
   	     spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
       (*spcPt)->toSplit();

  if (printL(2) && (_maxLevelOfSubProbRestriction < 1))
    {
      std::cout << "UPDATED CURRENT DUAL BOUND : " << "   objVal() = " << totalObjVal()
                << "   mastCurLagrangianBoundVal = " << mastCurLagrangianBoundVal << std::endl;
    }

  _algCurLpDualBound = mastCurLagrangianBoundVal;

  if (updateDualBnd)
      updateAlgDualBounds();

  if (_colGenStabilizationPtr != NULL)
    _colGenStabilizationPtr->updateOnLagrBoundChange(mastCurLagrangianBoundVal, _maxLevelOfSubProbRestriction);

  return;
}

bool Alg4EvalByLagrangianDuality::phase1isCompleted()
{

  if ((_colGenStabilizationPtr != NULL) && (_colGenStabilizationPtr->stabVarsInSolution()))
    return (false);

  /// Check for the presence of artificial variables in the solution
  if (_nonStabArtVarPtrList.empty())
    return (true);

  for (VarPtrSet::const_iterator sPtr = _probPtr->inPrimalLpSol().begin(); sPtr != _probPtr->inPrimalLpSol().end();
       sPtr++)
  {
      if ((*sPtr)->isTypeOf(VcId::ArtificialVarMask) && !(*sPtr)->val().isZero())
          return false;
  }

//  /// function count of VarPtrSet does not work properly (when containing both mast columns and pure mast variables)
//  std::list<Variable *>::const_iterator it;
//  for (it = _nonStabArtVarPtrList.begin(); it != _nonStabArtVarPtrList.end(); ++it)
//    {
//      if (_probPtr->inPrimalLpSol().count(*it) && !(*it)->val().isZero())
//        return (false);
//    }

  if (printL(5))
    std::cout << "Master::phaseI is completed" << std::endl;

  return (true);
}


void Alg4EvalByLagrangianDuality::setOptionLogPrintFrequency(const int value)
{
  _logPrintFrequency = value;
}

void Alg4EvalByLagrangianDuality::setOptionMaxNbOfCgIterations(const int value)
{
  _maxNbOfCgIterations = value;
}

void Alg4EvalByLagrangianDuality::setOptionMaxNbOfPenaltyUpdates(const int value)
{
  _maxNbOfPenaltyUpdates = value;
}

void Alg4EvalByLagrangianDuality::setOptionMinLevelOfSpRestriction(const int value)
{
  _minLevelOfSpRestriction = value;
}

void Alg4EvalByLagrangianDuality::setOptionMinNbOfCutRounds(const int value)
{
  _minNbOfCutRounds = value;
}

void Alg4EvalByLagrangianDuality::setOptionMaxNbOfCutRounds(const int value)
{
  _maxNbOfCutRounds = value;
}

void Alg4EvalByLagrangianDuality::setOptionDoRedCostFixingAndEnumeration(const int value)
{
  _doRedCostFixingAndEnumeration = value;
}

void Alg4EvalByLagrangianDuality::setOptionNonExactEvaluation(const bool value)
{
  _nonExactEvaluation = value;
}


