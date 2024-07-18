/**
 *
 * This file bcProblemC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcErrorC.hpp"
#include "bcFormC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "bcNetworkBasedCuts.hpp"
#endif //BCP_RCSP_IS_FOUND

#include "bcMasterConfC.hpp"
#include "bcModelRCSPSolver.hpp"
#include "bcProblemC.hpp"
#include "bcProbConfigC.hpp"
#include "bcPrintC.hpp"
#include "bcVarConstrC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcStabilizationInfo.hpp"
#include "bcModelBranchingConstrC.hpp"
#include "bcStabilizationColgen.hpp"

using namespace std;
using namespace VcIndexStatus;

void LpBasisRecord::clear(const bool removeMarksInVars, const bool removeMarksInConstrs)
{
    if (removeMarksInVars)
    {
        for (std::vector< VarPtr_MpFormIndexStatus >::iterator it = _varInBasis.begin();
             it != _varInBasis.end(); ++it)
        {
            it->_varPtr->isInBasis(false);
        }
    }
    _varInBasis.clear();
    if (removeMarksInConstrs)
    {
        for (std::vector< ConstrPtr_MpFormIndexStatus >::iterator it = _constrInBasis.begin();
             it != _constrInBasis.end(); ++it)
        {
            it->_constrPtr->isInBasis(false);
        }
    }
    _constrInBasis.clear();
    return;
}

void VariableSolInfo::applyVarInfo() const
{
  varPtr->problemPtr()->updatePartialSolution(varPtr, value);
}

Problem::Problem(const int & ref,
                 const double & rightHandSideZeroTol,
                 const double & reducedCostTolerance,
                 const BcObjStatus::MinMaxIntFloat & minmaxStatus,
                 const SolutionMethod & solMode,
                 const std::string & name,
                 const SolutionStatus & requiredStatus,
                 const bool & LPpreprocessorOn,
                 const bool & LPprobingOn,
                 const char & LPsolverSelect) :
        _ref(ref),
        _name(name),
        _probConfPtr(NULL),
        _solverOracleFunctorIsDefined(false),
        _solverOracleFunctorPtr(NULL),
        _fracSolBasedHeuristicFunctorPtr(NULL),
        _masterHeuristicFunctorPtr(NULL),
        _divingFixingFunctorPtr(NULL),
        _enumSolBasedHeuristicFunctor(NULL),
        _objStatus(minmaxStatus),
        _probInfeasibleFlag(false),
        _solMode(solMode),
        _minCost(- minmaxStatus * BapcodInfinity),
        _maxCost(minmaxStatus * BapcodInfinity),
        _probIsBuilt(false),
        _formulationPtr(NULL),
        _primalFormulationPtr(NULL),
        _posGlobalArtVarPtr(NULL),
        _negGlobalArtVarPtr(NULL),
        _objVal(0.0),
        _primalBound(minmaxStatus * BapcodInfinity),
        _dualBound(- minmaxStatus * BapcodInfinity),
        _probConstrManager(),
        _probVarManager(),
        _primalSolPtr(NULL),
        _dualSolPtr(NULL),
        _partialSolutionValue(0),
        _partialSolution(),
        _requiredStatus(requiredStatus),
        _probStatus(SolutionStatus::UnSolved),
        _LPpreprocessorOn(LPpreprocessorOn),
        _LPprobingOn(LPprobingOn),
        _LPsolverSelect(LPsolverSelect),
        _curNodePtr(NULL),
        _isRetrievedRedCosts(false),
        _rightHandSideZeroTol(rightHandSideZeroTol),
        _reducedCostTolerance(reducedCostTolerance)
{
}

const ProgStatus& Problem::progStatus() const
{
  return _probConfPtr->progStatus();
}

const ControlParameters& Problem::param() const
{
  return _probConfPtr->param();
}

void Problem::setMIPRequiredStatus(const SolutionStatus & newStatus)
{
    std::cout << "Problem::setMipRequiredStatus : Formulation will be solved as an LP. "
                 "Cannot set the MIP required status." << std::endl;
}

Problem::~Problem()
{
   ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Active, 'd');
   Constraint * constrPtr = *constrPtrIt;

  if (printL(6))
    std::cout << "~Problem() " << name() << std::endl;

  if (_primalSolPtr != NULL)
    {
      delete _primalSolPtr;
      _primalSolPtr = NULL;
    }

    if (_dualSolPtr != NULL)
    {
        delete _dualSolPtr;
        _dualSolPtr = NULL;
    }

  if (_solverOracleFunctorPtr != NULL)
    {
      delete _solverOracleFunctorPtr;
      _solverOracleFunctorPtr = NULL;
    }

  if (_fracSolBasedHeuristicFunctorPtr != NULL)
    {
      delete _fracSolBasedHeuristicFunctorPtr;
      _fracSolBasedHeuristicFunctorPtr = NULL;
    }

  if (_masterHeuristicFunctorPtr != NULL)
    {
      delete _masterHeuristicFunctorPtr;
      _masterHeuristicFunctorPtr= NULL;
    }

  if (_divingFixingFunctorPtr != NULL)
    {
      delete _divingFixingFunctorPtr;
      _divingFixingFunctorPtr= NULL;
    }

  if (_enumSolBasedHeuristicFunctor != NULL)
  {
      delete _enumSolBasedHeuristicFunctor;
      _enumSolBasedHeuristicFunctor = NULL;
  }

  clearRecordedSol();
  deleteForm();

  return;
}

void Problem::probConfPtr(ProbConfig * ptr)
{
  _probConfPtr = ptr;
}

void Problem::defineFormulation()
{
  if (printL(6))
    std::cout << "Prob name = " << name() << " _solMode.status() = " << _solMode.status() << std::endl;

  switch (_solMode.status()) {
  case SolutionMethod::lpSolver:
    {
      _formulationPtr = _primalFormulationPtr = new LPform(this);
      break;
    }
  case SolutionMethod::mipSolver:
  case SolutionMethod::custom2mipSolver:
    {
      _formulationPtr = _primalFormulationPtr = new MIPform(this);
      break;
    }
  case SolutionMethod::customSolver:
  {
    bool testOracles = param().CheckOracleOptimality() || param().CheckSpOracleFeasibility();
    if(testOracles)
      {
	     _formulationPtr = _primalFormulationPtr = new MIPform(this);
      }
    break;
  }
  case SolutionMethod::none:
    break;
  case SolutionMethod::undefined:
    bapcodInit().check(true, "Problem::~defineFormulation(): ERROR undefined solution method");
    break;
  }
}

Solution * Problem::retrieveCurPrimalLpSol(const bool & recordPartialSolOnly) //const
{
    int printlevel = 6;

    if (printL(printlevel))
        std::cout << "Problem::retrieveCurPrimalLpSol() Problem = " << name()
                  << " recordPartialSolOnly = " << recordPartialSolOnly << std::endl;

    Solution * solPtr(NULL);

    if (!recordPartialSolOnly)
    {
        if (_primalSolPtr != NULL)
        {
            solPtr = _primalSolPtr->clone();
        }
        else
        {
            solPtr = new Solution(_probConfPtr);

            if (inPrimalLpSol().empty() && param().GenerateProperColumns())
            {
                Solution * tmpSolPtr = extractIncumbent();
                if (tmpSolPtr != NULL)
                {
                    if (printL(printlevel))
                        std::cout << "Problem::retrieveCurPrimalLpSol() inPrimalSol is empty recordPartialSolOnly: Incumbent extracted" << std::endl;

                    solPtr->includeVars(tmpSolPtr->solVarValMap(), true);
                    delete tmpSolPtr;
                }
            }
            else
            {
                if (printL(printlevel))
                    std::cout << "Problem::retrieveCurPrimalLpSol() inPrimalSol extracted" << std::endl;
                solPtr->includeVarSet(inPrimalLpSol());
            }
        }
        solPtr->cost(partialSolutionValue() + primalBound());

    }
    else /// recordPartialSolOnly
    {
        solPtr = new Solution(_probConfPtr);
        solPtr->cost(partialSolutionValue());
    }

    solPtr->includeVars(_partialSolution, true);

    return (solPtr);
}

Solution * Problem::recordSolution(Solution * solPtr)
{
  bapcodInit().require(solPtr != NULL, "Problem::recordSolution() solution is not defined");

  if (printL(6))
    std::cout << "Problem::recordSolution(): sol  to insert has ref " << solPtr->ref() << std::endl
	          << " with cost = " << solPtr-> cost() << std::endl;

  /// insert in good position
  std::list <Solution * >::iterator it = _recordedSolPtr.begin();

  if (printL(6) && (it != _recordedSolPtr.end()))
    std::cout << "Problem::recordSolution():  first sol in record is  " << (**it) << std::endl;

  while ((it != _recordedSolPtr.end()) && ((**it) < (*solPtr)))
    {
      ++it;
      if (printL(6) && (it != _recordedSolPtr.end()))
	    std::cout << "Problem::recordSolution():  cur sol in record is  " << (**it) << std::endl;
    }
  if (it == _recordedSolPtr.end())
    {
      if (printL(6))
	    std::cout << "Problem::recordSolution(): add sol  at the end " << std::endl;

      _recordedSolPtr.push_back(solPtr);
    }
  else
    {
      if (printL(6)) std::cout << "Problem::recordSolution(): record  sol before " << (**it) <<  std::endl;

      _recordedSolPtr.insert(it, solPtr);
    }

  return (solPtr);
}

void Problem::recordPosGlobalArtVar(GlobalArtificialVar * globalArtVarPtr)
{
  _posGlobalArtVarPtr = globalArtVarPtr;
}

void Problem::recordNegGlobalArtVar(GlobalArtificialVar * globalArtVarPtr)
{
  _negGlobalArtVarPtr = globalArtVarPtr;
}

Solution * Problem::incumbentSol() const
{
  if (_recordedSolPtr.empty())
    return NULL;
  else
    return *(_recordedSolPtr.begin());
}

Solution * Problem::extractIncumbent()
{
  if (_recordedSolPtr.empty())
    return (NULL);
  else
  {
    Solution * ptr =  *(_recordedSolPtr.begin());
    _recordedSolPtr.erase(_recordedSolPtr.begin());
    return (ptr);
  }
}

void Problem::addConstrInForm()
{
  if (printL(5))
    std::cout << "Problem::addContrInForm()  "  << std::endl;

  if (_primalFormulationPtr != NULL)
    _primalFormulationPtr->addConstr2Formulation();

  return;
}

void Problem::addVarInForm()
{
  if (_primalFormulationPtr != NULL)
    _primalFormulationPtr->addVar2Formulation();

  return;
}

void Problem::delConstrInForm()
{
  if (_primalFormulationPtr != NULL)
    _primalFormulationPtr->delConstrFromFormulation();

  return;
}

void Problem::delVarInForm()
{
  if (_primalFormulationPtr != NULL)
    _primalFormulationPtr->delVarFromFormulation();

  return;
}

void Problem::addConstraintToProblem(Constraint * constrPtr)
{
    if (constrPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
    {
        _probNonLinearConstrSet.insert(constrPtr);
    }
    constrPtr->addToProb(this);

    constrPtr->curRhs(constrPtr->costrhs());

    Variable * varPtr = constrPtr->posLocalArtVarPtr();

    if (varPtr != NULL)
        varPtr->addToProb(this);
    varPtr = constrPtr->negLocalArtVarPtr();

    if (varPtr != NULL)
        varPtr->addToProb(this);

    if (constrPtr->stabInfoPtr() != NULL)
    {
        varPtr = constrPtr->stabInfoPtr()->negInnerArtVarPtr();
        if (varPtr != NULL)
            varPtr->addToProb(this);
        varPtr = constrPtr->stabInfoPtr()->negOuterArtVarPtr();
        if (varPtr != NULL)
            varPtr->addToProb(this);
        varPtr = constrPtr->stabInfoPtr()->posInnerArtVarPtr();
        if (varPtr != NULL)
            varPtr->addToProb(this);
        varPtr = constrPtr->stabInfoPtr()->posOuterArtVarPtr();
        if (varPtr != NULL)
            varPtr->addToProb(this);
    }

}

void Problem::insertActiveConstr(Constraint * constrPtr, const int & updateForm)
{

  insertConstr(constrPtr, Active);

  constrPtr->activate();

  if (constrPtr->kind() == 'E')
  {
    if (updateForm >= 1)
    {
      if (printL(6))
        std::cout << "Problem::insertActiveConstr() setConstr2Form "
              << constrPtr->name() << std::endl;

      setConstr2Form(constrPtr);
    }
    if (updateForm >= 2)
    {
      if (printL(6))
        std::cout << "Problem::insertActiveConstr() addConstrInForm() "
              << constrPtr->name() << std::endl;

      addConstrInForm();
    }
  }
}

void Problem::insertConstr(Constraint * constrPtr,
                           const VcIndexStatus::VcStatus & status)
{
  if (status != Undefined)
    _probConstrManager.insert(constrPtr, status);

  if (constrPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
  {
    _probNonLinearConstrSet.insert(constrPtr);
  }

}

void Problem::addConstr(Constraint * constrPtr, const int & flag, const int & updateForm)
{
  if (printL(6))
    std::cout << "Problem " << name() << " addConstr() " << constrPtr->name()
    << " constrPtr->vcIndexStatus() " << constrPtr->vcIndexStatus()
    << " flag " << flag
    << " updateForm " << updateForm << std::endl ;

  constrPtr->resetRhs();

  switch (flag) {
  case 1:
  {
    //if (_probConstrManager.count(constrPtr, Active))
    if (constrPtr->vcIndexStatus() == Active)
      return;
    bapcodInit().check(constrPtr->inCurProb(), "Problem::addConstr(): constr not in _probConstrManager should not be active");

    if (printL(6))
      std::cout << "Problem::addConstr() insert "
            << constrPtr->name() << std::endl;

    if (constrPtr->vcIndexStatus() == Undefined)
    {

      addConstr2Prob(constrPtr);

      /// need to di

      if (printL(6))
        std::cout << "Problem::addConstr() addConstr2Prob "
              << constrPtr->name() << std::endl;
    }
    // else
    //   {
    //     insertActiveConstr(constrPtr, updateForm);
    //   }

    insertActiveConstr(constrPtr, updateForm);

    return;
  }
  case 2:
  {
    // if (count(constrPtr, 3))
    //   return;
    // _inactiveConstrSet.insert(constrPtr);
    // if (_unsuitableConstrSet.count(constrPtr))
    //   _unsuitableConstrSet.erase(constrPtr);
    // else
    //   addConstr2Prob(constrPtr);

    // return;

    if (constrPtr->vcIndexStatus() == Inactive)
      //if (count(constrPtr, 3))
      return;

    if (constrPtr->vcIndexStatus() == Undefined)
    {
      addConstr2Prob(constrPtr); //_probConstrManager.insert(constrPtr);
    }

    insertConstr(constrPtr, Inactive);

    //        if (_unsuitableConstrSet.count(constrPtr))
    //	  _unsuitableConstrSet.erase(constrPtr);
    //        else


    return;
  }
  case 3:
  {
    if (constrPtr->vcIndexStatus() == Unsuitable)
      return;

    if (constrPtr->vcIndexStatus() == Undefined)
    {
      addConstr2Prob(constrPtr);
    }
    insertConstr(constrPtr, Unsuitable);
    return;
  }
  default:
    bapcodInit().check(1, "Problem::addConstr(): flag is not valid");
    return;
  }
}

void Problem::addConstr2Prob(Constraint * constrPtr)
{
  if (printL(6))
    std::cout << "Problem:addConstr2Prob() " << constrPtr->name()  << std::endl;
  constrPtr->addToProb(this);

  return;
}

void Problem::addVarsSimplyInForm(std::list<Variable *> & varsToAddToForm)
{
  std::list<Variable *>::iterator varPtrIt;
  if (!varsToAddToForm.empty())
  {
    for (varPtrIt = varsToAddToForm.begin(); varPtrIt != varsToAddToForm.end(); ++varPtrIt)
      if ((*varPtrIt)->kind() == 'E')
        setVar2Form(*varPtrIt);
    addVarInForm();
  }
}

void Problem::addVarsSimplyInForm(std::list<InstanciatedVar *> & varsToAddToForm)
{
  std::list<InstanciatedVar *>::iterator varPtrIt;
  if (!varsToAddToForm.empty())
  {
    for (varPtrIt = varsToAddToForm.begin(); varPtrIt != varsToAddToForm.end(); ++varPtrIt)
      if ((*varPtrIt)->kind() == 'E')
        setVar2Form(*varPtrIt);
    addVarInForm();
  }
}

void Problem::addConstrsSimplyInForm(std::list<Constraint *> & constrsToAddToForm)
{
  std::list<Constraint *>::iterator constrPtrIt;
  if (!constrsToAddToForm.empty())
  {
    for (constrPtrIt = constrsToAddToForm.begin(); constrPtrIt != constrsToAddToForm.end(); ++constrPtrIt)
      if ((*constrPtrIt)->kind() == 'E')
        setConstr2Form(*constrPtrIt);
    addConstrInForm();
  }
}

void Problem::delVarsSimplyInForm(std::list<Variable *> & varsToRemoveFromForm)
{
  std::list<Variable *>::iterator varPtrIt;
  if (!varsToRemoveFromForm.empty())
    {
      for (varPtrIt = varsToRemoveFromForm.begin(); varPtrIt != varsToRemoveFromForm.end();
          ++varPtrIt)
        {
          if ((*varPtrIt)->kind() == 'E')
            unsetVar2Form(*varPtrIt);
        }
      delVarInForm();
    }
}

void Problem::delConstrsSimplyInForm(std::list<Constraint *> & constrsToRemoveFromForm)
{
  std::list<Constraint *>::iterator constrPtrIt;
  if (!constrsToRemoveFromForm.empty())
    {
      for (constrPtrIt = constrsToRemoveFromForm.begin();
          constrPtrIt != constrsToRemoveFromForm.end(); ++constrPtrIt)
        if ((*constrPtrIt)->kind() == 'E')
          {
            unsetConstr2Form(*constrPtrIt);
          }
      delConstrInForm();
    }
}

int Problem::addVar(Variable * varPtr, const int & flag, const int & updateForm)
{
  int printlevel = 5;
  if (printL(printlevel))
    std::cout << "flag = " << flag << std::endl;

  varPtr->resetCost(false);

  if (printL(printlevel))
    std::cout << varPtr << std::endl;

  ///@todo: Implement record artificial variables
  switch (flag)
  {
    case 1:
    {
      switch (_solMode.status())
      {
        case SolutionMethod::lpSolver:
        case SolutionMethod::mipSolver:
        case SolutionMethod::custom2mipSolver:
        case SolutionMethod::customSolver:
        {
          if (varPtr->vcIndexStatus() == Active)
            return 0;  /// zero variable added to active formulation

          if (printL(printlevel))
            cout << " Adding  var " << varPtr->name() << " in Active Set " << " varPtr->vcIndexStatus() = "
                 << varPtr->vcIndexStatus() << std::endl;


          bapcodInit().check(varPtr->inCurProb(),
                             "Problem::addVar(): var not in Active Set, hence it should not be inCurProb");

          if (varPtr->vcIndexStatus() == VcIndexStatus::Undefined)
          {
            addVar2Prob(varPtr);
          }

          _probVarManager.insert(varPtr, VcIndexStatus::Active);

          varPtr->activate();

          if (varPtr->kind() == 'E')
          {
            if (updateForm >= 1)
              setVar2Form(varPtr);

            if (updateForm >= 2)
              addVarInForm();
          }

          return 1; /// one variable added to active formulation
        }
        case SolutionMethod::none:
          break;
        case SolutionMethod::undefined:
          bapcodInit().check(true, "Problem::~Problem(): ERROR undefined solution method");
          break;
      }
      break;
    }

    case 2:
    {
      if (printL(printlevel))
        std::cout << "flag = 2"
        << " varPtr->vcIndexStatus() = "
        << varPtr->vcIndexStatus()
        << std::endl;

      if (varPtr->vcIndexStatus() == VcIndexStatus::Inactive)
        return 0;

      if (printL(printlevel))
        std::cout << "inactiveVarSet.insert(varPtr)" << std::endl;

      if (varPtr->vcIndexStatus() == VcIndexStatus::Undefined)
        addVar2Prob(varPtr);

      _probVarManager.insert(varPtr, VcIndexStatus::Inactive);
      return 0;
    }
    case 3:
    {
      if (printL(printlevel))
        std::cout << "flag = 3" << " varPtr->vcIndexStatus() = " << varPtr->vcIndexStatus() << std::endl;

      if (varPtr->vcIndexStatus() == VcIndexStatus::Unsuitable)
        return 0;

      if (varPtr->vcIndexStatus() == VcIndexStatus::Undefined)
        addVar2Prob(varPtr);

      _probVarManager.insert(varPtr, VcIndexStatus::Unsuitable);

      return (true);
    }
    default:
      bapcodInit().check(1, "Problem::addVar(): flag is not valid");
      break;
    }
  return 0;
}

void Problem::addVar2Prob(Variable * varPtr)
{
  varPtr->addToProb(this);

  return;
}

void Problem::delVarFromProb(Variable * varPtr)
{
  _probVarManager.erase(varPtr);
  varPtr->deleteFromProb();
  return;
}

void Problem::delConstrFromProb(Constraint * constrPtr)
{
  _probConstrManager.erase(constrPtr);

  constrPtr->deleteFromProb();

  return;
}

void Problem::removeActiveConstr(Constraint * constrPtr, const int & updateForm)
{
  if (constrPtr->vcIndexStatus() != Active)
    return;

  constrPtr->desactivate();

  if (constrPtr->kind() == 'E')
    {
        if (updateForm >= 1)
            unsetConstr2Form(constrPtr);

        if (updateForm >= 2)
            delConstrInForm();
    }
}

void Problem::delConstr(Constraint * constrPtr, const int & flag, const int & updateForm)
{
    switch (flag) {
        case 1:
        {
            /// Delete from Active probConstrSet only, i.e. desactivate
            if (constrPtr->vcIndexStatus() == Active)
            {
                removeActiveConstr(constrPtr, updateForm);
                _probConstrManager.insert(constrPtr, Inactive);
            }
            break;
        }
        case 2:
        {
            /// Delete from probConstrSet and setOfInactiveConstr only, i.e setAside

            if (constrPtr->vcIndexStatus() == Active)
            {
                removeActiveConstr(constrPtr, updateForm);
                _probConstrManager.insert(constrPtr, Unsuitable);
            } else if (constrPtr->vcIndexStatus() == Inactive)
            {
                _probConstrManager.insert(constrPtr, Unsuitable);
            }
            break;
        }
        case 3:
        {

            if (constrPtr->vcIndexStatus() == Active)
            {
                removeActiveConstr(constrPtr, updateForm);
                delConstrFromProb(constrPtr);
            }
            else if ((constrPtr->vcIndexStatus() == Inactive) || (constrPtr->vcIndexStatus() == Unsuitable))
            {
                delConstrFromProb(constrPtr);
            }

            break;
        }
        default:
            bapcodInit().check(1, "Problem::delConstr(): flag is not valid");
            break;
    }
}

template<typename Container>
void Problem::delConstrSet(const Container & oldConstrSet, const int & flag, const int & updateForm)
{
    if (oldConstrSet.empty())
        return;

    typename Container::const_iterator it;
    std::list<Variable *> varPtrList;
    for (it = oldConstrSet.begin(); it != oldConstrSet.end(); ++it)
    {
        updateLocArtVarList(*it, varPtrList);
    }
    delVarSet(varPtrList, flag, updateForm);

    for (it = oldConstrSet.begin(); it != oldConstrSet.end(); ++it)
    {
        delConstr(*it, flag, (updateForm > 0 ? 1 : 0));
    }
    if (updateForm >= 2)
        delConstrInForm();
}

void Problem::updateLocArtVarList(Constraint * constrPtr, std::list<Variable *> & varPtrList)
{
  VarConstrStabInfo * stabPtr = constrPtr->stabInfoPtr();
  if (stabPtr == NULL)
    return;

  /// TO DO: first check whether this is a master constraint, otherwise not artificial variables linked
  if (stabPtr->negOuterArtVarPtr() != NULL)
    varPtrList.push_back(stabPtr->negOuterArtVarPtr());
  if (stabPtr->posOuterArtVarPtr() != NULL)
    varPtrList.push_back(stabPtr->posOuterArtVarPtr());
  if (stabPtr->negInnerArtVarPtr() != NULL)
    varPtrList.push_back(stabPtr->negInnerArtVarPtr());
  if (stabPtr->posInnerArtVarPtr() != NULL)
    varPtrList.push_back(stabPtr->posInnerArtVarPtr());
}

template<typename Container>
void Problem::addConstrSet(const Container & newConstrSet, const int & flag, const int & updateForm)
{
    if (newConstrSet.empty())
        return;

    typename Container::const_iterator it;
    std::list<Variable *> varPtrList;
    for (it = newConstrSet.begin(); it != newConstrSet.end(); ++it)
    {
        updateLocArtVarList(*it, varPtrList);
    }
    addVarSet(varPtrList, flag, updateForm);

    for (it = newConstrSet.begin(); it != newConstrSet.end(); ++it)
    {
        if ((*it)->name().find("vub") != std::string::npos)
        {
            continue;
        }
        if ((*it)->name().find("su") != std::string::npos)
        {
            continue;
        }
        addConstr(*it, flag, (updateForm > 0 ? 1 : 0));

        if (printL(7))
            std::cout << "Problem::addConstrSet(): added constr " << *it << std::endl;
    }
    if (updateForm >= 2)
        addConstrInForm();

    return;
}

template void Problem::addConstrSet<std::list<Constraint *>>(const std::list<Constraint *> & newConstrSet,
                                                             const int & flag, const int & updateForm);


/**
 * _probVarSet if flag=1,
 * _inactiveVarSet if flag=2,
 * _unsuitableVarSet if flag=3
 * updateForm == 0 => don't set it in formulation yet,
 * updateForm >= 1 => prepare it for insertion formulation,
 * updateForm >= 2 => set it in formulation .
 * return true if variable is really added
 */
template<typename Container>
void Problem::addVarSet(const Container & newVarSet, const int & flag, const int & updateForm)
{
    if (newVarSet.empty())
        return;

    typename Container::const_iterator it;
    for (it = newVarSet.begin(); it != newVarSet.end(); ++it)
    {
        addVar(*it, flag, (updateForm > 0 ? 1 : 0));
        if (printL(7))
            std::cout << "Problem::addVarSet(): added var " << *it << std::endl;
    }
    if (updateForm >= 2)
        addVarInForm();
}

template void Problem::addVarSet<std::list<SubProbVariable *>>(const std::list<SubProbVariable *> & newVarSet,
                                                               const int & flag, const int & updateForm);
template void Problem::addVarSet<std::list<Variable *>>(const std::list<Variable *> & newVarSet,
                                                        const int & flag, const int & updateForm);

/**
 * _probVarSet if flag=1,
 * up to _inactiveVarSet if flag=2,
 * up to _unsuitableVarSet  if flag=3
 */
template<typename Container>
void Problem::delVarSet(const Container & oldVarSet, const int & flag, const int & updateForm)
{
    if (oldVarSet.empty())
        return;

    typename Container::const_iterator it;
    for (it = oldVarSet.begin(); it != oldVarSet.end(); ++it)
    {
        delVar(*it, flag, (updateForm > 0 ? 1 : 0));
    }

    if (updateForm >= 2)
        delVarInForm();
}

template void Problem::delVarSet<std::list<Variable *>>(const std::list<Variable *> & oldVarSet, const int & flag,
                                                        const int & updateForm);

template void Problem::delVarSet<std::set<Variable *, VarConstrSort>>
              (const std::set<Variable *, VarConstrSort> & oldVarSet, const int & flag, const int & updateForm);

void Problem::removeVar(Variable * varPtr, const int & updateForm)
{
  if (printL(6))
    std::cout << "Problem::removeVar(" << varPtr->name() << "," << updateForm << ")" << std::endl;

  if (varPtr->vcIndexStatus() != Active)
    return;

  _probVarManager.erase(varPtr);

  varPtr->desactivate();

  if (varPtr->kind() == 'E')
  {
    if (updateForm >= 1)
      unsetVar2Form(varPtr);

    if (updateForm >= 2)
      delVarInForm();
  }
}

void Problem::delVar(Variable * varPtr, const int & flag, const int & updateForm)
{
  if (printL(6))
    std::cout << "Problem::delVar(" << varPtr->name() << "," << flag << "," << updateForm << ")" << std::endl;

  switch (flag) {
  case 1: // delete from probVarSet only, i.e. deactivate
  {
    if (varPtr->vcIndexStatus() == Active)
    {
      removeVar(varPtr, updateForm);
      _probVarManager.insert(varPtr, Inactive);
    }
    return;
  }
  case 2: // delete from probVarSet and setOfInactiveVar only, i.e setAside
  {
    if (varPtr->vcIndexStatus() == Active)
    {
      removeVar(varPtr, updateForm);
      _probVarManager.insert(varPtr, Unsuitable);
    } else if (varPtr->vcIndexStatus() == Inactive)
    {
      _probVarManager.insert(varPtr, Unsuitable);
    }

    return;
  }
  case 3: // delete from probVarSet and setOfInactiveVar and _unsuitableVarSet
  {
    if (varPtr->vcIndexStatus() == Active)
    {
      removeVar(varPtr, updateForm);
      delVarFromProb(varPtr);
    }
    else if ((varPtr->vcIndexStatus() == Inactive) || (varPtr->vcIndexStatus() == Unsuitable))
    {

      delVarFromProb(varPtr);
    }
    return;
  }
  default:
    bapcodInit().check(1, "Problem::delVar(): flag is not valid");
    return;
  }
}

void Problem::updateObjCoeffsInForm(const std::list<Variable *> & varPtrList)
{
  for (std::list<Variable *>::const_iterator it = varPtrList.begin();
       it != varPtrList.end(); ++it)
  {
    Variable * varPtr = *it;

    if (varPtr->index() >= 0)
    {
      if (printL(6))
        std::cout << "Problem::updateObjCoeffsInForm, var " << varPtr->name() << ", cost = " << varPtr->curCost()
                  << std::endl;

      if (primalFormulationPtr() != NULL)
        primalFormulationPtr()->resetObjCoef(varPtr);
    }
  }

  if ( !varPtrList.empty() )
  {
    if (primalFormulationPtr() != NULL)
      primalFormulationPtr()->updateObjectiveInFormulation();
  }
}

void Problem::updateBoundsInForm(const std::list<Variable *> & varPtrList)
{
  for (std::list<Variable *>::const_iterator it = varPtrList.begin();
       it != varPtrList.end(); ++it)
  {
    Variable * varPtr = *it;
    if (varPtr->index() >= 0)
    {
      if (primalFormulationPtr() != NULL)
        primalFormulationPtr()->resetBounds(varPtr);

      if (printL(6))
        std::cout << "Problem::resetBoundsInForm, var " << varPtr->name() << std::endl;
    }
  }

  if ( !varPtrList.empty() )
  {
    if (primalFormulationPtr() != NULL)
      primalFormulationPtr()->updateBoundsInFormulation();
  }
}

void Problem::updateConstrRhsInForm(const std::list<Constraint *> & constrPtrList)
{
  for (std::list<Constraint *>::const_iterator it = constrPtrList.begin();
       it != constrPtrList.end(); ++it)
  {
    Constraint * constrPtr = *it;
    if (constrPtr->index() >= 0)
    {
      if (primalFormulationPtr() != NULL)
        primalFormulationPtr()->resetRhs(constrPtr);

      if (printL(6))
        std::cout << "Problem::resetRhsInForm, constr " << constrPtr->name() << std::endl;
    }
  }

  if ( !constrPtrList.empty() )
  {
    if (primalFormulationPtr() != NULL)
      primalFormulationPtr()->updateConstrRhsInFormulation();
  }
}

void Problem::clearPreprocessingLists()
{
  VarPtrList::iterator varPtrIt;
  for (varPtrIt = _preprocessedVarsList.begin();
      varPtrIt != _preprocessedVarsList.end(); ++varPtrIt)
    (*varPtrIt)->removeFromPreprocessedList();
  _preprocessedVarsList.clear();

  ConstrPtrList::iterator constrPtrIt;
  for (constrPtrIt = _preprocessedConstrsList.begin();
      constrPtrIt != _preprocessedConstrsList.end(); ++constrPtrIt)
    (*constrPtrIt)->removeFromPreprocessedList();
  _preprocessedConstrsList.clear();
}

bool Problem::updateProbConstr(char flag)
{
    for (ConstrIndexManager::iterator constrPt = _probConstrManager.begin(Active, flag);
         constrPt != _probConstrManager.end(Active, flag); ++constrPt)
    {
        (*constrPt)->resetRhs();
        if ((*constrPt)->kind() == 'E' && (_primalFormulationPtr != NULL))
            _primalFormulationPtr->resetRhs(*constrPt);
    }
    return false;
}

bool Problem::updateProbVar(bool inPurePhaseOne, int printlevel, char flag)
{
  for (VarIndexManager::iterator varPt = _probVarManager.begin(Active, flag);
       varPt != _probVarManager.end(Active, flag); ++varPt)
  {
    if (printL( printlevel))
      std::cout << "Problem::updateProbVar():   consider var " << (*varPt)->name()
                << " inCurProb ?" << (*varPt)->inCurProb() << std::endl;
    bapcodInit().require((*varPt)->inCurProb(),
                         "Problem::updateProb(): var in _probVarManager should be marked as inCurProb",
                         ProgStatus::run  );
    {
      if (printL( printlevel))
        std::cout << "Problem::updateProbVar():  var " << (*varPt)->name() << " in [" << (*varPt)->curLb()
                  << ", " << (*varPt)->curUb() << " ] " << std::endl;

      if ((*varPt)->infeasible())
      {
        if (printL(3))
          std::cout << "Problem::updateProbVar(): infeasibility detected, due to variable " << (*varPt)->name()
                    << std::endl;

        _probInfeasibleFlag = true;
        /// Problem infeasible
        /// early exit
        {
          if (_primalFormulationPtr != NULL)
          {
            _primalFormulationPtr->clearColFormulationDataStruct();
            _primalFormulationPtr->clearRowFormulationDataStruct();
          }
          return true;
        }
      }

      (*varPt)->resetCost(inPurePhaseOne);

      if (printL( printlevel))
        std::cout << "   var " << (*varPt)->name() << " has cost " <<  (*varPt)->curCost() << std::endl;

      if ((*varPt)->kind() == 'E')
      {
        if (_primalFormulationPtr != NULL)
          _primalFormulationPtr->resetObjCoef(*varPt);
      }
      /// Set marks to update bounds
      {
        if (_primalFormulationPtr != NULL)
          _primalFormulationPtr->resetBounds(*varPt);
      }


      if ((*varPt)->curCost().negative())
      {
        _minCost += (*varPt)->curCost() * (*varPt)->curUb();
        _maxCost += (*varPt)->curCost() * (*varPt)->curLb();
      } else
      {
        _minCost += (*varPt)->curCost() * (*varPt)->curLb();
        _maxCost += (*varPt)->curCost() * (*varPt)->curUb();
      }
    }
  }
  return false;
}

bool Problem::rankOneCutsArePresent()
{
    GenericCutConstr *genR1CutConstrPtr = probConfPtr()->getGenericCutConstr("R1C");
    if (genR1CutConstrPtr != NULL) {
        const IndexCell2InstancConstrPtrMap &constrPtrMap = genR1CutConstrPtr->indexCell2InstancConstrPtrMap();
        for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
        {
            if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
                && it->second->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
                return true;
        }
    }
    return false;
}

void Problem::updateProbVarCostWithSoftConflictCuts()
{
    MasterConf * mastConfPtr = _probConfPtr->mastConfPtr();
    std::set<GenericCutConstr *, DynamicGenConstrSort>::iterator genCutPtrIt;
    std::set<Variable *> varsToResetObjCoeff;

    for (genCutPtrIt = mastConfPtr->candidateCutGenericConstr().begin();
         genCutPtrIt != mastConfPtr->candidateCutGenericConstr().end(); ++genCutPtrIt)
    {
        GenericSoftConflictsCutConstr * genSoftConfCutConstrPtr
                = dynamic_cast<GenericSoftConflictsCutConstr *>(*genCutPtrIt);
        if (genSoftConfCutConstrPtr != NULL)
        {
            GenericVar * genIndicVarPtr = _probConfPtr->getGenericVar(genSoftConfCutConstrPtr->defaultName() + "V");
            const IndexCell2InstancConstrPtrMap & constrPtrMap
                    = genSoftConfCutConstrPtr->indexCell2InstancConstrPtrMap();
            IndexCell2InstancConstrPtrMap::const_iterator mapIt;
            for (mapIt = constrPtrMap.begin(); mapIt != constrPtrMap.end(); ++mapIt)
                if (mapIt->second->vcIndexStatus() == VcIndexStatus::Active)
                {
                    SoftConflictsCut * cutPtr = static_cast<SoftConflictsCut *>(mapIt->second);
                    if (cutPtr->cutType() == 0)
                    {
                        std::vector<std::pair<SubProbVariable *, SubProbVariable *> >::const_iterator pairIt;
                        for (pairIt = cutPtr->conflicts().begin(); pairIt != cutPtr->conflicts().end(); ++pairIt)
                        {
                            if (pairIt->first->cgSpConfPtr() == _probConfPtr)
                            {
                                MultiIndex id(pairIt->first->id().first(), pairIt->second->id().first());
                                InstanciatedVar *ivPtr = genIndicVarPtr->checkIfInstanciationAlreadyExist(id);
                                if (ivPtr == NULL)
                                {
                                    std::cerr << "BaPCod error : no indicator variable for conflict "
                                              << id.first() << "<->" << id.second() << " in problem "
                                              << _probConfPtr->name() << " for cut " << cutPtr->name() << std::endl;
                                    exit(1);
                                }
                                SubProbVariable *spVarPtr = static_cast<SubProbVariable *>(ivPtr);
                                spVarPtr->updateCurCostWithConstraint(cutPtr);
                                varsToResetObjCoeff.insert(ivPtr);
                            }
                        }
                    }
                    else /// cutType == 1
                    {
                        MultiIndex id(mapIt->second->id().first());
                        InstanciatedVar *ivPtr = genIndicVarPtr->checkIfInstanciationAlreadyExist(id);
                        if (ivPtr != NULL)
                        {
                            SubProbVariable *spVarPtr = static_cast<SubProbVariable *>(ivPtr);
                            spVarPtr->updateCurCostWithConstraint(cutPtr);
                            varsToResetObjCoeff.insert(ivPtr);
                        }
                    }

                }
        }
#ifdef BCP_RCSP_IS_FOUND
        GenericLimMemRankOneCutConstr * genRankOneCutConstrPtr
                = dynamic_cast<GenericLimMemRankOneCutConstr *>(*genCutPtrIt);
        if ((genRankOneCutConstrPtr != NULL) && (genRankOneCutConstrPtr->spVarName() != ""))
        {
            GenericVar * genIndicVarPtr = _probConfPtr->getGenericVar(genRankOneCutConstrPtr->defaultName() + "V");
            const IndexCell2InstancConstrPtrMap & constrPtrMap
                    = genRankOneCutConstrPtr->indexCell2InstancConstrPtrMap();
            IndexCell2InstancConstrPtrMap::const_iterator mapIt;
            for (mapIt = constrPtrMap.begin(); mapIt != constrPtrMap.end(); ++mapIt)
                if (mapIt->second->vcIndexStatus() == VcIndexStatus::Active)
                {
                    MultiIndex id(mapIt->second->id().first());
                    InstanciatedVar * ivPtr = genIndicVarPtr->checkIfInstanciationAlreadyExist(id);
                    if (ivPtr != NULL)
                    {
                        SubProbVariable *spVarPtr = static_cast<SubProbVariable *>(ivPtr);
                        spVarPtr->updateCurCostWithConstraint(mapIt->second);
                        varsToResetObjCoeff.insert(ivPtr);
                    }
                }
        }
#endif //BCP_RCSP_IS_FOUND
    }
    for (std::set<Variable *>::iterator varPtrIt = varsToResetObjCoeff.begin();
         varPtrIt != varsToResetObjCoeff.end(); ++varPtrIt)
    {
        _primalFormulationPtr->replaceObjCoef(*varPtrIt);
        if ((*varPtrIt)->curCost().negative())
        {
            _minCost += (*varPtrIt)->curCost() * (*varPtrIt)->curUb();
            _maxCost += (*varPtrIt)->curCost() * (*varPtrIt)->curLb();
        } else
        {
            _minCost += (*varPtrIt)->curCost() * (*varPtrIt)->curLb();
            _maxCost += (*varPtrIt)->curCost() * (*varPtrIt)->curUb();
        }
    }
};

/**
 *
 * @return True if  problem infeasible
 */
bool Problem::updateProbForColGen(bool inPurePhaseOne)
{
  _minCost = 0;
  _maxCost = 0;

  int printlevel = 6;

  if (printL(printlevel))
    std::cout << "Problem::updateProbForColGen(): _probVarManager.size() = " << _probVarManager.size() << std::endl;

  if (updateProbVar(inPurePhaseOne, printlevel, 's'))
    return true;

  /// updating reduced cost of indicator variables of the soft conflict cuts
  if (param().colGenSubProbSolMode().getStatusAsInteger() != SolutionMethod::customSolver)
    updateProbVarCostWithSoftConflictCuts();

  if (printL( printlevel))
    std::cout << "   _minCost = " << _minCost << "   _maxCost = " << _maxCost << std::endl;

  if (_primalFormulationPtr != NULL)
    {
      _primalFormulationPtr->updateObjectiveInFormulation();
      _primalFormulationPtr->updateBoundsInFormulation();
    }
  return (false);
}

void Problem::hardResetObjective(char flag)
{
    for (VarIndexManager::iterator varPt = _probVarManager.begin(Active, flag);
         varPt != _probVarManager.end(Active, flag); ++varPt)
    {
        (*varPt)->costrhs(0.0);
    }
}

void Problem::hardResetConstraintRHS(char flag)
{
    for (ConstrIndexManager::iterator constrPt = _probConstrManager.begin(Active, flag);
         constrPt != _probConstrManager.end(Active, flag); ++constrPt)
    {
        (*constrPt)->costrhs(0.0);
    }
}

bool Problem::updateProblem()
{
    resetSolution();
    clearRecordedSol();

    _minCost = 0;
    _maxCost = 0;
    int printlevel = 6;

    if (updateProbVar(false, printlevel, 's'))
        return true;

    updateProbConstr('s');

    if (_primalFormulationPtr != NULL)
    {
        _primalFormulationPtr->updateObjectiveInFormulation();
        _primalFormulationPtr->updateBoundsInFormulation();
        _primalFormulationPtr->updateConstrRhsInFormulation();
    }

    setDualBound(Bound::infDualBound(objStatus()));
    setPrimalLpBound(Bound::infPrimalBound(objStatus()));

    return (false);
}

void Problem::setVar2Form(Variable * varPtr)
{
  bapcodInit().require(varPtr->kind() == 'E', "Problem::setVar2Form(): implicit var should not be set in formulation");
  bapcodInit().require(varPtr->inCurProb(), "Problem::setVar2Form():  var should have been activated");

  varPtr->setInForm();

  if (_primalFormulationPtr  != NULL)
  {
    _primalFormulationPtr->setVar2Form(varPtr);
  }
}

void Problem::unsetVar2Form(Variable * varPtr)
{
  bapcodInit().require(varPtr->kind() == 'E', "Problem::unsetVar2Form(): implicit var should not be set in formulation");
  bapcodInit().require(!varPtr->inCurProb(), "Problem::setVar2Form():  var should have been desactivated");

  varPtr->unsetInForm();

  if (_primalFormulationPtr  != NULL)
  {
    _primalFormulationPtr->unsetVar2Form(varPtr);
  }
}

void Problem::setConstr2Form(Constraint * constrPtr)
{
  if (printL(6))
    std::cout << "Problem::setConstr2Form() constr " << constrPtr->name() << std::endl;

  bapcodInit().require(constrPtr->kind() == 'E',
                       std::string("Problem::setConstr2Form():  constraint" + constrPtr->name()
                                   + " is implicit and therefore should not added to the formulation").c_str());

  bapcodInit().require(constrPtr->inCurProb(), std::string("Problem::setConstr2Form():  constraint "
                                                           + constrPtr->name()
                                                           + " should have been activated").c_str());

  bapcodInit().require((bool)_probConstrManager.count(constrPtr, Active),
                       "Problem::setConstr2Form() active Constr should be in probConstrSet");

  constrPtr->setInForm();

  if (_primalFormulationPtr  != NULL)
  {
    _primalFormulationPtr->setConstr2Form(constrPtr);
  }
}

void Problem::unsetConstr2Form(Constraint * constrPtr)
{
  if (printL(6))
    std::cout << "unset constr " << constrPtr->name() << std::endl;

  bapcodInit().require(constrPtr->kind() == 'E',
                       "Problem::unsetConstr2Form():  "
                       "constraint  is implicit and therefore should not added to the formulation");

  if (printL(7) && constrPtr->inCurProb())
  {
      cout << "constr shouldave been deactivated: " << constrPtr->name() << " at "
           << hex << (long)constrPtr << dec << endl;
  }

  bapcodInit().require(!constrPtr->inCurProb(),
                       "Problem::unsetConstr2Form():  constraint should have been desactivated:");

  constrPtr->unsetInForm();

  if (_primalFormulationPtr  != NULL)
  {
    _primalFormulationPtr->unsetConstr2Form(constrPtr);
  }
}

void Problem::buildProblem()
{
  // set variables and constraints for formulation , then build formulation

  if (_probIsBuilt)
    return;
  _probIsBuilt = true;

  for (ConstrIndexManager::iterator itc = _probConstrManager.begin(Active, 's');
       itc != _probConstrManager.end(Active, 's'); ++itc)
  {
    (*itc)->activate();
    if ((*itc)->kind() == 'E') setConstr2Form(*itc);
  }

  for (ConstrIndexManager::iterator itc = _probConstrManager.begin(Active, 'd');
       itc != _probConstrManager.end(Active, 'd'); ++itc)
  {
    (*itc)->activate();
    if ((*itc)->kind() == 'E')
      setConstr2Form(*itc);
  }

  for (VarIndexManager::iterator itv =  _probVarManager.begin(Active, 's');
       itv != _probVarManager.end(Active, 's'); ++itv)
  {
    (*itv)->activate();
    if ((*itv)->kind() == 'E')
      setVar2Form(*itv);
  }

  for (VarIndexManager::iterator itv =  _probVarManager.begin(Active, 'd');
       itv != _probVarManager.end(Active, 'd'); ++itv)
  {
    (*itv)->activate();
    if ((*itv)->kind() == 'E')
    {
      setVar2Form(*itv);
    }
  }

  for (VarIndexManager::iterator itv =  _probVarManager.begin(Active, 'a');
       itv != _probVarManager.end(Active, 'a'); ++itv)
  {
    (*itv)->activate();
    if ((*itv)->kind() == 'E')
    {
      setVar2Form(*itv);
    }
  }
  //end

  if (printL(6)) printProb();

  if (_primalFormulationPtr != NULL)
    _primalFormulationPtr->buildFormulation();

  return;
}

void Problem::storeDualSolution(ConstrPtr2DoubleMap & dualSolution) const
{
  dualSolution.clear();
  for (ConstrPtrSet::const_iterator constrPtrIt = _inDualSol.begin(); constrPtrIt != _inDualSol.end(); ++constrPtrIt)
    {
      dualSolution[*constrPtrIt] = (*constrPtrIt)->val();
    }
}

void Problem::resetDualSolution(const ConstrPtr2DoubleMap & dualSolution)
{
  resetSolution('d');
  for (ConstrPtr2DoubleMap::const_iterator mapIt = dualSolution.begin(); mapIt != dualSolution.end(); ++mapIt)
    {
      mapIt->first->val(mapIt->second);
      _inDualSol.insert(mapIt->first);
    }
}

void Problem::resetSolution(const char & flag)
{
  _probStatus = SolutionStatus::UnSolved;

    if (printL(1))
        std::cout << "resetSolution started " << std::endl;

    if (_primalSolPtr != NULL)
  {
    delete _primalSolPtr;
    _primalSolPtr = NULL;
  }

    for (VarPtrSet::const_iterator vit = _inPrimalSol.begin(); vit != _inPrimalSol.end(); ++vit)
      (*vit)->val(0);
    _inPrimalSol.clear();

    if (printL(1))
        std::cout << "primal solution cleared " << std::endl;


    for (VarPtrSet::const_iterator vit = _nonZeroRedCostVars.begin();
       vit != _nonZeroRedCostVars.end(); ++vit)
    (*vit)->reducedCost(0);

  _nonZeroRedCostVars.clear();

    if (printL(1))
        std::cout << "reduced costs are cleared " << std::endl;

    if (flag == 'd')
    {
      for (ConstrPtrSet::const_iterator cit = _inDualSol.begin(); cit != _inDualSol.end(); ++cit)
        {
          (*cit)->val(0);
        }
      _inDualSol.clear();
    }

    if (printL(1))
        std::cout << "dual solution cleared " << std::endl;

    _isRetrievedRedCosts = false;

  return;
}

void Problem::resetPartialSolution()
{
    _partialSolutionValue = 0;
    for (VarPtr2DoubleMap::iterator mapIt = _partialSolution.begin(); mapIt != _partialSolution.end(); ++mapIt)
        if (mapIt->first->isTypeOf(VcId::MastColumnMask))
            mapIt->first->decrParticipation(9);
    _partialSolution.clear();
}

/**
 * Record current  constraint val in inDualSol
 */
void Problem::updateInDualSol()
{
  int printlevel = 6;

  _inDualSol.clear();

  for (ConstrIndexManager::iterator constrPt = _probConstrManager.begin(Active, 's');
       constrPt != _probConstrManager.end(Active, 's'); ++constrPt)
  {
    if (!zeroTest((*constrPt)->valOrSepPointVal()))
    {
      _inDualSol.insert(*constrPt);
      if (printL(printlevel))
        std::cout << "Problem::updateInDualSol() DualSol[" << (*constrPt)->name() << "] = "
                  << (*constrPt)->valOrSepPointVal() << std::endl;
    }
  }

  for (ConstrIndexManager::iterator constrPt = _probConstrManager.begin(Active, 'd');
       constrPt != _probConstrManager.end(Active, 'd'); ++constrPt)
  {
    if (!zeroTest((*constrPt)->valOrSepPointVal()))
    {
      _inDualSol.insert(*constrPt);
      if (printL(printlevel))
        std::cout << "Problem::updateInDualSol() DualSol[" << (*constrPt)->name() << "] = "
                  << (*constrPt)->valOrSepPointVal() << std::endl;
    }
  }
}


void Problem::deleteForm()
{
  if (_primalFormulationPtr != NULL)
  {
    delete _primalFormulationPtr;
    _primalFormulationPtr = NULL;
  }

  _formulationPtr = NULL;
  return;
}

void Problem::deactivateAndRemoveAllVarsAndConstrsFromMemory()
{
   resetSolution('d');

   /// we deactivate all dynamic variables and constraints
   ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Active, 'd');
   while (constrPtrIt != _probConstrManager.end(VcIndexStatus::Active, 'd'))
    {
       Constraint * constrPtr = *constrPtrIt;
       ++constrPtrIt;
       constrPtr->desactivateConstraint(VcIndexStatus::Unsuitable);
    }
    VarIndexManager::iterator varPtrIt = _probVarManager.begin(VcIndexStatus::Active, 'd');
    while (varPtrIt != _probVarManager.end(VcIndexStatus::Active, 'd'))
    {
        Variable * varPtr = *varPtrIt;
        varPtrIt++;
        varPtr->desactivateVariable(VcIndexStatus::Unsuitable);
    }
    varPtrIt = probVarSet().begin(VcIndexStatus::Inactive, 'd');
    while (varPtrIt != _probVarManager.end(VcIndexStatus::Inactive, 'd'))
    {
        Variable * varPtr = *varPtrIt;
        varPtrIt++;
        varPtr->desactivateVariable(VcIndexStatus::Unsuitable, false);
    }
    delConstrInForm();
    delVarInForm();

    /// we remove all dynamic variables and constraints from the memory
    /// (as well as artificial variables where were associated with the dynamic constraints)

    /// constraints should be removed first (their associated artificial variables are checked during deletion)
    constrPtrIt = _probConstrManager.begin(Unsuitable, 'd');
    int counter = 0;
    while (constrPtrIt != _probConstrManager.end(Unsuitable, 'd'))
    {
        Constraint * constrPtr = *constrPtrIt;
        if (printL(1) && (constrPtr->participation() != 0))
            std::cout << "BaPCod warning : participation of constraint " << constrPtr->name() << " is "
                      << constrPtr->participation() << std::endl;
        constrPtrIt++;
        _probConstrManager.erase(constrPtr);
        delete constrPtr;
        counter++;
    }
    if (printL(3))
        std::cout << "Removed " << counter << " dynamic constraints from memory " << std::endl;

    varPtrIt = _probVarManager.begin(Unsuitable, 'd');
    counter = 0;
    while (varPtrIt != _probVarManager.end(Unsuitable, 'd'))
    {
        Variable * varPtr = *varPtrIt;
        if (printL(1) && (varPtr->participation() != 0))
            std::cout << "BaPCod warning : participation of variable " << varPtr->name() << " is "
                      << varPtr->participation() << std::endl;
        varPtrIt++;
        _probVarManager.erase(varPtr);
        delete varPtr;
        counter++;
    }
    if (printL(3))
      std::cout << "Removed " << counter << " dynamic variables from memory " << std::endl;

    varPtrIt = _probVarManager.begin(Unsuitable, 'a');
    counter = 0;
    while (varPtrIt != _probVarManager.end(Unsuitable, 'a'))
    {
        Variable * varPtr = *varPtrIt;
        varPtrIt++;
        _probVarManager.erase(varPtr);
        delete varPtr;
        counter++;
    }
    if (printL(3))
        std::cout << "Removed " << counter << " artificial variables from memory " << std::endl;
}

void Problem::removeUnusedDynamicVarsFromMemory(const bool & checkArtVars)
{
  if ((param().MaxTimeForRestrictedMasterIpHeur() > 0) && param().ActivateAllColumnsForRestrictedMasterIpHeur())
    return;

  if (param().Search4NegRedCostColInInactivePool())
    return;

  if(printL(7))
    std::cout << "probVarPts.size(Unsuitable, 'd') = "  <<
      _probVarManager.size(Unsuitable, 'd') << std::endl;

  VarIndexManager::iterator varPtrIt = _probVarManager.begin(Unsuitable, 'd');
  while (varPtrIt != _probVarManager.end(Unsuitable, 'd'))
    {
      Variable * varPtr = *varPtrIt;
      varPtrIt++;
      if (printL(7))
        {
          std::cout << "col " << std::endl;
          std::cout << "in 0x" << hex  << (long) varPtr << dec << std::endl;
          std::cout << "with indexInProb " << varPtr->vcIndexInProb() << std::endl;
          std::cout << varPtr->name() << std::endl;
          std::cout << "is being tested for removal from problem" << std::endl;
        }

      if (varPtr->isTypeOf(VcId::MastColumnMask) && (varPtr->participation() == 0))
        {
          _probVarManager.erase(varPtr);
          if (printL(7))
            {
              std::cout << "col " << std::endl;
              std::cout << "in 0x" << hex  << (long) varPtr << dec << std::endl;
              std::cout << "with indexInProb " << varPtr->vcIndexInProb() << std::endl;
              std::cout <<varPtr->name() << std::endl;
              std::cout << "has been removed from problem" << std::endl;
            }
          delete varPtr;
        }
    }
  /// checkArtVars should be false if the procedure is called inside an evaluation algorithm
  if (checkArtVars)
    {
      VarIndexManager::iterator varPtrIt = _probVarManager.begin(Unsuitable, 'a');
      while (varPtrIt != _probVarManager.end(Unsuitable, 'a'))
      {
          if ((*varPtrIt)->isTypeOf(VcId::LocalArtificialVarMask))
          {
              LocalArtificialVar * locArtVarPtr = static_cast<LocalArtificialVar *>(*varPtrIt);
              ++varPtrIt;
              if (locArtVarPtr->canBeDeleted())
              {
                  _probVarManager.erase(locArtVarPtr);
                  delete locArtVarPtr;
              }
          }
          else
          {
              ++varPtrIt;
          }
      }
    }
  if (printL(7))
    std::cout << "probVarPts.size(Unsuitable, 'd') = "  << _probVarManager.size(Unsuitable, 'd') << std::endl;
}

void Problem::removeUnusedDynamicConstrsFromMemory()
{
    ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(Unsuitable, 'd');
    while (constrPtrIt != _probConstrManager.end(Unsuitable, 'd'))
    {
        Constraint * constrPtr = (*constrPtrIt);
        ++constrPtrIt;

        if (printL(7))
          std::cout << "constraint " << (constrPtr)->name() << " is being tested for removal from problem" << std::endl;

        if (constrPtr->isTypeOf(VcId::InstMasterConstrMask) && (constrPtr->participation() == 0))
          {
          /// we cannot now delete associated local art. var as they can participate
          /// in the _nonStabArtVarPtrList of the eval. alg.
          /// they will be deleted in removeUnusedDynamicVarsFromMemory(true)
          /// when it is called from Alg4ProblemSetDownOfNode::run()
          Variable * posLocalArtVarPtr = constrPtr->posLocalArtVarPtr();
          if ((posLocalArtVarPtr != NULL) && (posLocalArtVarPtr->isTypeOf(VcId::LocalArtificialVarMask)))
            {
              LocalArtificialVar * locArtVarPtr = static_cast<LocalArtificialVar *>(posLocalArtVarPtr);
              /// by doing this, we allow deletion of the var. in removeUnusedDynamicVarsFromMemory
              locArtVarPtr->setConstraintPtr(NULL);
            }
          Variable * negLocalArtVarPtr = constrPtr->negLocalArtVarPtr();
          if ((negLocalArtVarPtr != NULL) && (negLocalArtVarPtr->isTypeOf(VcId::LocalArtificialVarMask)))
            {
              LocalArtificialVar * locArtVarPtr = static_cast<LocalArtificialVar *>(negLocalArtVarPtr);
              /// by doing this, we allow deletion of the var. in removeUnusedDynamicVarsFromMemory
              locArtVarPtr->setConstraintPtr(NULL);
            }

          if (printL(7))
              std::cout << "constraint " << (constrPtr)->name() << " is being removed from problem" << std::endl;
          _probConstrManager.erase(constrPtr);
          if (constrPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
            _probNonLinearConstrSet.erase(constrPtr);
          delete constrPtr;
        }
    }
}

void Problem::retrieveRedCosts()
{
  if (!_isRetrievedRedCosts)
  {
     _primalFormulationPtr->LPform::retrieveRedCosts(false, _nonZeroRedCostVars);
     _isRetrievedRedCosts = true;
  }
}

void Problem::setProbStatus(const SolutionStatus & stat)
{
  if (printL(5))
    std::cout << "Problem::setProbStatus() for " << name() << ",  stat = " << stat << std::endl;

  _probStatus = stat;
}

void MipProblem::setProbStatus(const SolutionStatus & stat)
{
  if (printL(5))
    std::cout << "MipProblem::setProbStatus() for " << name() << ",  stat = " << stat << std::endl;

  _probStatus = stat;
  _mipProbStatus = stat;
}

void Problem::setProbStatus(const int & stat)
{
  if (printL(5))
    std::cout << "Problem::setProbStatus() for " << name() << ",  stat = " << stat << std::endl;

  _probStatus = stat;
}

void MipProblem::setProbStatus(const int & stat)
{
  if (printL(5))
    std::cout << "Problem::setProbStatus() for " << name() << ",  stat = " << stat << std::endl;

  _probStatus = stat;
  _mipProbStatus = stat;
}

void Problem::printDynamicVarConstrStats(std::ostream & os, const bool & completePrint)
{
  int numActiveDynVars = _probVarManager.size(VcIndexStatus::Active,'d');
  int numInactiveDynVars = _probVarManager.size(VcIndexStatus::Inactive,'d');
  int numDynVars = numActiveDynVars + _probVarManager.size(VcIndexStatus::Inactive,'d')
                   + _probVarManager.size(VcIndexStatus::Unsuitable,'d');
  int numActiveDynConstrs = _probConstrManager.size(VcIndexStatus::Active,'d');
  int numDynConstrs = numActiveDynConstrs + _probConstrManager.size(VcIndexStatus::Unsuitable,'d');
  int numActiveArtVars = _probVarManager.size(VcIndexStatus::Active, 'a');
  int numArtVars = numActiveArtVars + _probVarManager.size(VcIndexStatus::Unsuitable, 'a');
  std::cout << numDynVars << " columns (" << numActiveDynVars << " active";
  if (numInactiveDynVars > 0)
    std::cout << ", " << numInactiveDynVars << " inactive";
  std::cout << "), " << numDynConstrs << " dyn. constrs. (" << numActiveDynConstrs << " active), "
            << numArtVars << " art. vars. (" << numActiveArtVars << " active)";

  if (!completePrint)
    return;

  std::map<int, int> treatOrderMap;
  for (VarIndexManager::iterator varPtrIt = _probVarManager.begin(VcIndexStatus::Unsuitable, 'd');
       varPtrIt != _probVarManager.end(VcIndexStatus::Unsuitable, 'd'); ++varPtrIt)
    {
      MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
      int treatOrderId = colPtr->treatOrderId();
      if (treatOrderMap.count(treatOrderId))
        treatOrderMap[treatOrderId] += 1;
      else
        treatOrderMap[treatOrderId] = 1;
    }
  if (!treatOrderMap.empty())
    {
      std::cout << std::endl << "Unsuitable columns : ";
      for (std::map<int, int>::iterator mapIt = treatOrderMap.begin(); mapIt != treatOrderMap.end(); ++mapIt)
        std::cout << mapIt->first << "(" << mapIt->second << ") ";
    }
  treatOrderMap.clear();
  for (ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Unsuitable, 'd');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Unsuitable, 'd'); ++constrPtrIt)
    if ((*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask))
      {
        InstMasterConstr * instMastConstrPtr = static_cast<InstMasterConstr *>(*constrPtrIt);
        if (!instMastConstrPtr->isTypeOf(VcId::BranchingConstrBaseTypeMask))
          continue;
        int treatOrderId = instMastConstrPtr->treatOrderId();
        if (treatOrderMap.count(treatOrderId))
          treatOrderMap[treatOrderId] += 1;
        else
          treatOrderMap[treatOrderId] = 1;
      }
  if (!treatOrderMap.empty())
    {
      std::cout << std::endl << "Unsuitable branching constraints : ";
      for (std::map<int, int>::iterator mapIt = treatOrderMap.begin(); mapIt != treatOrderMap.end(); ++mapIt)
        std::cout << mapIt->first << "(" << mapIt->second << ") ";
    }
  treatOrderMap.clear();
  for (ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Unsuitable, 'd');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Unsuitable, 'd'); ++constrPtrIt)
    if ((*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask))
      {
        InstMasterConstr * instMastConstrPtr = static_cast<InstMasterConstr *>(*constrPtrIt);
        if (instMastConstrPtr->isTypeOf(VcId::BranchingConstrBaseTypeMask))
          continue;
        int treatOrderId = instMastConstrPtr->treatOrderId();
        if (treatOrderMap.count(treatOrderId))
          treatOrderMap[treatOrderId] += 1;
        else
          treatOrderMap[treatOrderId] = 1;
      }
  if (!treatOrderMap.empty())
    {
      std::cout << std::endl << "Unsuitable cuts : ";
      for (std::map<int, int>::iterator mapIt = treatOrderMap.begin(); mapIt != treatOrderMap.end(); ++mapIt)
        std::cout << mapIt->first << "(" << mapIt->second << ") ";
    }

  std::cout << std::endl;
  /// we compute now estimated memory size of the memberships of dynamic variables and constraints
  unsigned long numberOfMembers = 0;
  for (VarIndexManager::iterator varPtrIt = _probVarManager.begin(VcIndexStatus::Active, 'd');
       varPtrIt != _probVarManager.end(VcIndexStatus::Active, 'd'); ++varPtrIt)
    numberOfMembers += (*varPtrIt)->member2coefMap().size();
  for (VarIndexManager::iterator varPtrIt = _probVarManager.begin(VcIndexStatus::Unsuitable, 'd');
       varPtrIt != _probVarManager.end(VcIndexStatus::Unsuitable, 'd'); ++varPtrIt)
    numberOfMembers += (*varPtrIt)->member2coefMap().size();
  for (ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Active, 'a');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Active, 'a'); ++constrPtrIt)
    numberOfMembers += (*constrPtrIt)->member2coefMap().size();
  for (ConstrIndexManager::iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Unsuitable, 'd');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Unsuitable, 'd'); ++constrPtrIt)
    numberOfMembers += (*constrPtrIt)->member2coefMap().size();
  std::cout << "Estimated memory usage of dynamic var/constr membership is "
            << (numberOfMembers * 32.0) / 1000000 << " Mb.";
}

int Problem::enumerateAllColumns(Solution * solPtr)
{
    int nbEnumSols = 0;
    if (!_solverOracleFunctorIsDefined)
        return -1;

    _solverOracleFunctorPtr->reducedCostFixingAndEnumeration(BcFormulation(_probConfPtr),
                                                             1 /* enumeration for heuristic mode */,
                                                             BapcodInfinity);
    if (!_probConfPtr->enumeratedStatus())
        return -1;

    nbEnumSols += getNumberOfEnumeratedSolutions();
    std::vector<double> dummyVector;
    Solution * localSolPtr = new Solution(_probConfPtr);
    getEnumeratedSolutions(-1, localSolPtr, dummyVector);
    localSolPtr->previousSolPtr(solPtr);

    return nbEnumSols;
}


void Problem::reducedCostFixingAndEnumeration(const int & enumerationMode, const Double & threshold)
{
  if (_solverOracleFunctorIsDefined)
  {
    _solverOracleFunctorPtr->reducedCostFixingAndEnumeration(BcFormulation(_probConfPtr), enumerationMode,
                                                             threshold.val());
  }
}

void Problem::checkEnumeratedSolutions(const std::vector<Solution *> & solPts, std::vector<bool> & solIsEnumerated)
{
    if (_solverOracleFunctorIsDefined)
    {
        _solverOracleFunctorPtr->checkEnumeratedSolutions(BcFormulation(_probConfPtr), solPts,
                                                          solIsEnumerated);
    }
}

void Problem::getEnumeratedSolutions(const int & maxNumOfSolutions, Solution * solPtr, std::vector<double> & redCosts)
{
    if (_solverOracleFunctorIsDefined)
    {
        BcSolution bcSol(solPtr);
        _solverOracleFunctorPtr->getEnumeratedSolutions(BcFormulation(_probConfPtr), maxNumOfSolutions, bcSol,
                                                        redCosts);
    }
}

bool Problem::setDebugSolution(const std::vector<std::vector<int> > & orderedSolutions, bool vertexBased)
{
    if (_solverOracleFunctorIsDefined)
    {
        return _solverOracleFunctorPtr->setDebugSolution(orderedSolutions, vertexBased);
    }
    return false;
}

void Problem::getDebugSolution(Solution * solPtr)
{
    if (_solverOracleFunctorIsDefined)
    {
        BcSolution bcSol(solPtr);
        _solverOracleFunctorPtr->getDebugSolution(BcFormulation(_probConfPtr), bcSol);
    }
}

bool Problem::isProperSolution(Solution * solPtr)
{
  if (_solverOracleFunctorIsDefined)
    {
      BcSolution bcSol(solPtr);
      return _solverOracleFunctorPtr->isProperSolution(BcFormulation(_probConfPtr), bcSol);
    }
  return true;
}

bool Problem::solSatisfiesCurrentSpRelaxation(Solution * solPtr)
{
  if (_solverOracleFunctorIsDefined)
    {
      BcSolution bcSol(solPtr);
      return _solverOracleFunctorPtr->solSatisfiesCurrentSpRelaxation(BcFormulation(_probConfPtr), bcSol);
    }
  return true;
}

void Problem::callColGenTerminationCallBack(bool afterRedCostFixing, int nodeOrder, int nodeDepth, int cutSepRound,
                                            double dualBound, double elapsedTime, bool masterConverged)
{
  if (_solverOracleFunctorIsDefined)
    _solverOracleFunctorPtr->columnGenerationTerminated(BcFormulation(_probConfPtr), afterRedCostFixing, nodeOrder,
                                                        nodeDepth, cutSepRound, dualBound, elapsedTime,
                                                        masterConverged);
}

/// for the moment, will draw separately for every graph
/// TO DO : draw based on packing sets
bool Problem::drawPrimalSolutionToDotFile(std::vector<MastColumn *> & colsInMasterSolution,
                                          const std::string & filename) const
{
    if (_solverOracleFunctorIsDefined)
    {
        std::vector<std::pair<BcSolution, double> > masterSolution;
        for (std::vector<MastColumn *>::const_iterator colPtrIt = colsInMasterSolution.begin();
             colPtrIt != colsInMasterSolution.end(); ++colPtrIt)
            masterSolution.push_back(std::make_pair(BcSolution((*colPtrIt)->spSol()), (*colPtrIt)->val()));
        return _solverOracleFunctorPtr->drawPrimalSolutionToDotFile(BcFormulation(_probConfPtr), masterSolution, filename);
    }
    return false;
}

bool Problem::improveCurrentSpRelaxation(std::vector<MastColumn *> & colsInMasterSolution,
                                         const bool & masterConverged)
{
  if (_solverOracleFunctorIsDefined)
    {
      std::vector<std::pair<BcSolution, double> > masterSolution;
      for (std::vector<MastColumn *>::const_iterator colPtrIt = colsInMasterSolution.begin();
           colPtrIt != colsInMasterSolution.end(); ++colPtrIt)
        masterSolution.push_back(std::make_pair(BcSolution((*colPtrIt)->spSol()), (*colPtrIt)->val()));
      return _solverOracleFunctorPtr->improveCurrentSpRelaxation(BcFormulation(_probConfPtr), masterSolution,
                                                                 masterConverged);
    }
  return false;
}

bool Problem::lightenCurrentSpRelaxation(const int & masterConverged, const int & callMode)
{
  if (_solverOracleFunctorIsDefined)
    {
      return _solverOracleFunctorPtr->lightenCurrentSpRelaxation(BcFormulation(_probConfPtr),
                                                                 masterConverged, callMode);
    }
  return false;
}


int Problem::getMessageIdToCutGeneration() const
{
  if (_solverOracleFunctorIsDefined)
    {
      return _solverOracleFunctorPtr->getMessageIdToCutGeneration();
    }
  return PricingSolverCutsMessage::noMessage;
}

int Problem::getNumberOfEnumeratedSolutions() const
{
    if (_solverOracleFunctorIsDefined)
    {
        return _solverOracleFunctorPtr->getNumberOfEnumeratedSolutions();
    }
    return -1;
}

bool Problem::customizedSolver(int & maxLevelOfRestriction,
                               const SolutionStatus & requiredSolStat,
                               Double & objVal,
                               Double & primalBound,
                               Double & dualBound,
                               VarPtrSet & inPrimalSol,
                               ConstrPtrSet & inDualSol,
                               BcSolution & primalSolPtr,
                               BcDualSolution & dualSolPtr)
{
    bool status(false);

    if (_solverOracleFunctorIsDefined)
    {
        int phaseOfStageTemp = maxLevelOfRestriction;

        status = (*_solverOracleFunctorPtr)(BcFormulation(_probConfPtr),
                                            objVal.val(),
                                            primalBound.val(),
                                            dualBound.val(),
                                            primalSolPtr,
                                            dualSolPtr,
                                            maxLevelOfRestriction);

        if (status)
        {
            double recomputedObjVal = 0;
            Solution * solPtr = (Solution *) primalSolPtr;

            for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin();
                 it != solPtr->solVarValMap().end(); ++it)
            {
                if (!zeroTest(it->second))
                {
                    recomputedObjVal += it->first->curCost() * it->second;
                    it->first->val(it->second);
                    inPrimalSol.insert(it->first);
                }
            }

            solPtr = solPtr->nextSolPtr();
            while (solPtr != NULL)
            {
                recordSolution(solPtr);
                solPtr = solPtr->nextSolPtr();
            }

            if (bapcodInit().param().CheckOracleOptimality && (phaseOfStageTemp == 0))
            {
                _primalFormulationPtr->updateObjectiveInFormulation();

                Double tempObjVal;
                Double tempPrimBound;
                Double tempDualBound;
                VarPtrSet tempInPrimalSol;
                ConstrPtrSet tempInDualSol;

                ((MIPform*)_primalFormulationPtr)->solve(param().MasterMipSolverBarrierConvergenceTolerance(),
                                                         _rightHandSideZeroTol,
                                                         _reducedCostTolerance,
                                                         'p',
                                                         false,
                                                         SolutionStatus::Optimum,
                                                         tempObjVal,
                                                         tempPrimBound,
                                                         tempDualBound,
                                                         tempInPrimalSol,
                                                         tempInDualSol);

                if ( ( (tempObjVal - objVal)/(objVal+0.000001) > 0.001
                       || (tempObjVal - objVal)/(objVal+0.000001) < - 0.001)
                     && ((tempObjVal - objVal) > 0.001 || (tempObjVal - objVal) < - 0.001))
                {
                    _primalFormulationPtr->printForm();
                    std::cerr << "Problem::customizedSolver OPTIMALITY ERROR" << endl;
                    std::cerr << "Solution of the mip used to CheckOracleOptimality " << endl;

                    for (VarPtrSet::iterator it = tempInPrimalSol.begin(); it != tempInPrimalSol.end(); it++)
                    {
                        std::cerr << "Var " << (*it)->name() << " has val " << (*it)->val() << std::endl;
                    }

                    std::cerr << "Solution of the oracle " << endl;

                    for (VarPtrSet::iterator it = inPrimalSol.begin(); it != inPrimalSol.end(); it++)
                    {
                        std::cerr << "Var " << (*it)->name() << " has val " << (*it)->val() << std::endl;
                    }

                    std::cerr << "objVal of CheckOracleOptimality Mip = " << tempObjVal << endl;
                    std::cerr << "objVal of Oracle = " << objVal << std::endl;
                    exit(1);
                }

                if (printL(1))
                    std::cout << "CheckOracleOptimality MIP ObjVal =" << tempObjVal << endl;
            }

            if (bapcodInit().param().CheckSpOracleFeasibility)
            {
                for (std::list<Variable*>::const_iterator varIt = _probConfPtr->pcVarPtrList().begin();
                     varIt != _probConfPtr->pcVarPtrList().end(); varIt++)
                {
                    if ((*varIt)->isTypeOf(VcId::SubProbVariableMask))
                    {
                        SubProbVariable* spVarPtr = static_cast<SubProbVariable*>(*varIt);
                        if(spVarPtr->costrhs() == 0 && spVarPtr->masterConstrMember2coefMap().empty())
                        {
                            spVarPtr->resetCurCostByValue(0);
                            _primalFormulationPtr->resetObjCoef(spVarPtr);
                        }
                        else
                        {
                            spVarPtr->resetCurCostByValue(10);
                            _primalFormulationPtr->resetObjCoef(spVarPtr);
                        }
                    }
                }

                _primalFormulationPtr->updateObjectiveInFormulation();

                Solution * solPtr = (Solution *) primalSolPtr;
                for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin();
                     it != solPtr->solVarValMap().end(); ++it)
                {
                    if (!zeroTest(it->second))
                    {
                        SubProbVariable* spVarPtr = static_cast<SubProbVariable*>(it->first);
                        spVarPtr->resetCurCostByValue(0);
                        _primalFormulationPtr->resetObjCoef(spVarPtr);
                    }
                }

                _primalFormulationPtr->updateObjectiveInFormulation();

                Double tempObjVal;
                Double tempPrimBound;
                Double tempDualBound;
                VarPtrSet tempInPrimalSol;
                ConstrPtrSet tempInDualSol;

                ((MIPform*)_primalFormulationPtr)->solve(param().MasterMipSolverBarrierConvergenceTolerance(),
                                                         _rightHandSideZeroTol,
                                                         _reducedCostTolerance,
                                                         'p',
                                                         false,
                                                         SolutionStatus::Optimum,
                                                         tempObjVal,
                                                         tempPrimBound,
                                                         tempDualBound,
                                                         tempInPrimalSol,
                                                         tempInDualSol);

                if (tempObjVal > 0.01)
                {
                    _primalFormulationPtr->printForm();
                    std::cerr << "Problem::customizedSolver FEASIBILITY ERROR" << endl;
                    std::cerr << "Solution of the mip used to CheckSpOracleFeasibility " << endl;

                    for (VarPtrSet::iterator it = tempInPrimalSol.begin(); it != tempInPrimalSol.end(); it++)
                    {
                        std::cerr << "Var " << (*it)->name() << " has val " << (*it)->val() << std::endl;
                    }

                    std::cerr << "Solution of the oracle " << std::endl;

                    for (VarPtrSet::iterator it = inPrimalSol.begin(); it != inPrimalSol.end(); it++)
                    {
                        std::cerr << "Var " << (*it)->name() << " has val " << (*it)->val() << std::endl;
                    }

                    std::cerr << "objVal of CheckSpOracleFeasibility Mip = " << tempObjVal << std::endl;
                    std::cerr << "objVal of CheckSpOracleFeasibility Mip should be 0 "  << std::endl;
                    exit(1);
                }

                cout << "CheckSpOracleFeasibility MIP ObjVal =" << tempObjVal << endl;
            }
        }

        inDualSol.clear(); // clear the dual sol at each iteration
        for (ConstrPtr2DoubleMap::const_iterator it = ((DualSolution *) dualSolPtr)->dualSolConstrValMap().begin();
             it != ((DualSolution *) dualSolPtr)->dualSolConstrValMap().end(); ++it)
        {
            if (!zeroTest(it->second))
            {
                it->first->val(it->second);
                inDualSol.insert(it->first);
            }
        }

        return status;
    }
    else
    {

        cout << " customizedSolver NOT DEFINED (old customizedSolver template is now oboslete)"
             << endl;
        exit(1);
    }
    return status;
}

void Problem::getActiveRyanAndFosterBranchingConstraints
     (std::list<BcRyanAndFosterBranchConstr> & ryanAndFosterBranchConstrList) const
{
  ryanAndFosterBranchConstrList.clear();
  for (ConstrIndexManager::const_iterator constrPtrIt = _probConstrManager.begin(Active, 'd');
       constrPtrIt != _probConstrManager.end(Active, 'd'); ++constrPtrIt)
    {
      if ((*constrPtrIt)->isTypeOf(VcId::RyanAndFosterInstSubProbBranchConstrMask))
        ryanAndFosterBranchConstrList.push_back(
          BcRyanAndFosterBranchConstr(dynamic_cast<RyanAndFosterInstSubProbBranchConstr *>(*constrPtrIt)));
    }
}

bool Problem::setupCustomizedSolver(const BcSolverOracleInfo * infoPtr)
{
  bool status(false);
  if (_solverOracleFunctorIsDefined)

    status = _solverOracleFunctorPtr->setupNode(BcFormulation(_probConfPtr), infoPtr);

  return status;
}

BcSolverOracleInfo * Problem::recordSolverOracleInfo()
{
  if (_solverOracleFunctorIsDefined)
    return _solverOracleFunctorPtr->recordSolverOracleInfo(BcFormulation(_probConfPtr));
  return NULL;
}

/*
 * @param flag
 * @param ifPrint
 *
 * @return True if a feasible solution to the problem was found, false otherwise
 */
bool Problem::solveProbLP(const char & flag, const bool & ifPrint)
{
  bool foundSol(false);
  bapcodInit().check(_primalFormulationPtr == NULL,
                     "Problem::solveProbLP(): _solMode == lp or mipSolver => requires  defined formulation");

  foundSol = _primalFormulationPtr->LPform::solve(param().MasterMipSolverBarrierConvergenceTolerance()
                                                  , _rightHandSideZeroTol
                                                  , _reducedCostTolerance
                                                  , flag
                                                  , ifPrint
                                                  , _requiredStatus
                                                  , _objVal
                                                  , _primalBound
                                                  , _dualBound
                                                  , _inPrimalSol
                                                  , _inDualSol
                                                  , _LPpreprocessorOn
                                                  , _LPprobingOn
                                                  , false
                                                  , _LPsolverSelect);

  if (printL(1))
    std::cout  << "Problem::solveProbLP(): " << name() << " _objVal = " << _objVal << std::endl;

  setProbStatus(_primalFormulationPtr->status());

//  if (param().solverName() == "CLP_SOLVER")
//  {
//    /// Ruslan : sometimes CLP solver returns the objective value which is not exact (I do not know why),
//    /// so, we recompute the objective value here using the primal solution
//    std::cout << " init objVal = " << std::setprecision(10) << _objVal << std::endl;
//    _primalBound = 0.0;
//    for (auto * varPtr : _inPrimalSol)
//    {
//        _primalBound += varPtr->val() * varPtr->curCost();
//        std::cout << " " << varPtr->name() << "=" << std::setprecision(10) << varPtr->val();
//    }
//    _objVal = _primalBound;
//    std::cout << " objVal = " << std::setprecision(10) << _objVal << std::endl;
//
//  }

  if (_probStatus.count(SolutionStatus::Optimum) || _probStatus.count(SolutionStatus::OptimumUnscalInfeas))
    _objVal = _dualBound = _primalBound;

  if (printL(5))
    std::cout  << "Problem::solveProbLP(): probStatus() after _primalFormulationPtr->LPform::solve()" << probStatus()
               << std::endl << " _requiredStatus= " << _requiredStatus << std::endl;

  return foundSol;
}

bool Problem::solveProbMIP(const char & flag, const bool & ifPrint)
{
  std::cerr << "Problem::solveProbMIP() should not be called: problem is not a MIP"
            << std::endl;
  exit(1);
}

int Problem::solveProb(int & maxLevelOfRestriction, const char & flag, const bool & ifPrint)
{
  int solverReturnStatus(0);
  resetSolution(flag);

  switch (_solMode.status())
  {
    case SolutionMethod::none: break;
    case SolutionMethod::lpSolver:
    case SolutionMethod::mipSolver:
    {
      if (solveProbLP(flag, ifPrint))
        solverReturnStatus = 1;
      break;
    }
    case SolutionMethod::customSolver:
    case SolutionMethod::custom2mipSolver:
    {
      if (printL(5))
        std::cout  << "Problem::solveProb(): to enter customizedSolver()  " << std::endl;

        if (_primalSolPtr == NULL)
            _primalSolPtr = new Solution(_probConfPtr);

        if (_dualSolPtr == NULL)
            _dualSolPtr = new DualSolution(_probConfPtr);

        BcSolution bcSol(_primalSolPtr);
        BcDualSolution bcDualSol(_dualSolPtr);

        if (customizedSolver(maxLevelOfRestriction,
                             _requiredStatus,
                             _objVal,
                             _primalBound,
                             _dualBound,
                             _inPrimalSol,
                             _inDualSol,
                             bcSol,
                             bcDualSol))
        solverReturnStatus = 1;

      break;
    }
    default:
    {
      bapcodInit().check(true, "Problem solMode undefined");
      break;
    }
  }

  if (printL(1))
  {
    printSolVal();
  }

  setStatusAfterSol();

  return solverReturnStatus;
}

void Problem::setStatusAfterSol()
{
  if (printL(3))
    printSol();

  if (printL(5))
    printDualSol();

  bool statusIsSet = true;
  /// added by Ruslan to correctly set the problem status after solving a MIP
  if ((_objStatus == BcObjStatus::minFloat) || (_objStatus == BcObjStatus::minInt))
    {
      /// minimization problem
      if ((_primalBound == BapcodInfinity) && (_dualBound == BapcodInfinity))
        setProbStatus(SolutionStatus::Infeasible);
      else if ((_primalBound == -BapcodInfinity) && (_dualBound == -BapcodInfinity))
        setProbStatus(SolutionStatus::Unbounded);
      else
        statusIsSet = false;
    }
  else
    {
      /// maximization problem
      if ((_primalBound == BapcodInfinity) && (_dualBound == BapcodInfinity))
        setProbStatus(SolutionStatus::Unbounded);
      else if ((_primalBound == -BapcodInfinity) && (_dualBound == -BapcodInfinity))
        setProbStatus(SolutionStatus::Infeasible);
      else
        statusIsSet = false;
    }

  if (!statusIsSet)
    {
      if (_primalBound == _dualBound)
        setProbStatus(SolutionStatus::Optimum);
      else if (_primalBound == _objVal)
        setProbStatus(SolutionStatus::PrimalFeasSolFound);
      else if (_objVal == _dualBound)
        setProbStatus(SolutionStatus::DualFeasSolFound);
    }

  if (printL(5))
    std::cout << "Problem::setStatusAfterSol(): probStatus()=" << probStatus()
              << ", _requiredStatus= " << _requiredStatus << std::endl;
}

void MipProblem::setStatusAfterSol()
{
  Problem::setStatusAfterSol();
  if (!(probStatus().intersects(_mipRequiredStatus)))
  {
    if (printL(5))
      std::cout << "MipProblem::setStatusAfterSol(): mipProbStatus APPARENTLY DOES NOT SATISFY REQUIRED STATUS, "
                   "test for PrimalFeasSolFound" << std::endl;

  }
  if (printL(5))
    std::cout  << "MipProblem::setStatusAfterSol(): mipProbStatus()" << probStatus() << std::endl
               << "_mipRequiredStatus= " << _mipRequiredStatus << std::endl;
}

std::ostream& Problem::printSolVal(std::ostream& os) const
{
  os << "printSol(Problem name= " << name() << "), objStatus= " << (int) _objStatus << std::endl;
  os << "   objVal = " << _objVal << std::endl;
  os << "   partialSolutionValue = " << _partialSolutionValue << std::endl;
  os << "   totalValue = " << _objVal + _partialSolutionValue << std::endl;

  if (_formulationPtr != NULL)
    os << "   status = " << _formulationPtr->status() << std::endl;

  return (os);
}

std::ostream& Problem::printActiveDynamicConstraints(std::ostream& os) const
{
  std::cout << "Active master branching constraints : " << std::endl;
  for (ConstrIndexManager::const_iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Active, 'd');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    if ((*constrPtrIt)->isTypeOf(VcId::InstMasterBranchingConstrMask))
      {
        const InstMasterBranchingConstr * brConstrPtr = static_cast<const InstMasterBranchingConstr *>(*constrPtrIt);
        brConstrPtr->shortPrint();
        std::cout << ", treatOrderId = " << brConstrPtr->treatOrderId()
                  << ", dualVal = " << brConstrPtr->val() << std::endl;
      }
  std::cout << "Active master cuts : " << std::endl;
  for (ConstrIndexManager::const_iterator constrPtrIt = _probConstrManager.begin(VcIndexStatus::Active, 'd');
       constrPtrIt != _probConstrManager.end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    if ((*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask)
        && !(*constrPtrIt)->isTypeOf(VcId::InstMasterBranchingConstrMask))
      {
        const InstMasterConstr * constrPtr = static_cast<const InstMasterConstr *>(*constrPtrIt);
        std::cout << constrPtr->name() << ", treatOrderId = " << constrPtr->treatOrderId()
                  << ", dualVal = " << constrPtr->val() << std::endl;
      }
  return (os);
}

std::ostream& Problem::printDetailedPrimalSol(std::ostream& os) const
{
  os << "Problem " << _probConfPtr->name() << " solution with value " << _objVal << " : " << std::endl;
  for (VarPtrSet::const_iterator sPtr = _inPrimalSol.begin(); sPtr != _inPrimalSol.end(); ++sPtr)
    {
      os << (*sPtr)->name() << " = " <<  (*sPtr)->val();
      if ((*sPtr)->isTypeOf(VcId::MastColumnMask))
        {
          MastColumn * colPtr = static_cast<MastColumn *>(*sPtr);
          os << ", spId = " << colPtr->cgSpConfPtr()->id().first();
          os << ", treatOrderId = " << colPtr->treatOrderId();
          if (colPtr->spSol()->enumeratedFlag())
            os << ", enumerated";
          if (colPtr->spSol() != NULL)
            {
#ifdef BCP_RCSP_IS_FOUND
              if ((colPtr->spSol()->rcspSolPtr() != nullptr) && !colPtr->spSol()->rcspSolPtr()->arcIds.empty())
#else
              if (!colPtr->spSol()->orderedIds().empty())
#endif
              {
                  colPtr->spSol()->printOrderedSolution(os);
                }
              else
                {
                  os << ", spSol = (";
                  for (VarPtr2DoubleMap::const_iterator mapIt = colPtr->spSol()->solVarValMap().begin();
                       mapIt != colPtr->spSol()->solVarValMap().end(); ++mapIt)
                    {
                      if (mapIt->first->genVarPtr()->defaultName() == "TLCCV")
                          continue;
                      if (mapIt->first->genVarPtr()->defaultName() == "R1CV")
                          continue;
                      if (mapIt != colPtr->spSol()->solVarValMap().begin())
                        os << ", ";
                      os << mapIt->first->name() << " = " << mapIt->second;
                    }
                  os << ")" << std::endl;
                }
            }
        }
      else
        os << std::endl;
    }

  return (os);
}

std::ostream& Problem::printSol(std::ostream& os) const
{
  printSolVal(os);

  std::list<Variable*> masterColumns;

  for (VarPtrSet::const_iterator sPtr = _inPrimalSol.begin(); sPtr != _inPrimalSol.end(); ++sPtr)
  {
    os << "primalSol[" << (*sPtr)->name() << "] = " <<  (*sPtr)->val() << std::endl;
    if ((*sPtr)->isTypeOf(VcId::MastColumnMask))
    {
      masterColumns.push_back((*sPtr));
    }
  }

  for(std::list<Variable*>::iterator it = masterColumns.begin() ; it != masterColumns.end() ; it++)
  {
    os << (*it)->name() << ": " << endl;
    (*it)->spSol()->print(os);
    os << endl;
  }

  return (os);
}

std::ostream& Problem::printDualSol(std::ostream& os, bool compact) const
{
    if (compact)
    {
        os << "Dualsol :";
        os << std::setprecision(12);
        for (ConstrPtrSet::const_iterator sPtr = inDualSol().begin(); sPtr != inDualSol().end(); ++sPtr)
            os << " " << (*sPtr)->name() << "=" << (*sPtr)->valOrSepPointVal();
        os << std::endl;
        os << std::setprecision(6);
        return os;
    }

    os << "printDualSol(Problem name= " << name() << "), objStatus= " << (int) _objStatus << std::endl;

    for (ConstrPtrSet::const_iterator sPtr = inDualSol().begin(); sPtr != inDualSol().end(); ++sPtr)
        os << "dualSol[" << (*sPtr)->name() << "] = " << std::setprecision(10) << (*sPtr)->val()  << std::endl;

    return (os);
}

std::ostream& Problem::printProb(std::ostream& os) const
{
  os << "printProb(Problem name= " << name() << ") , objStatus= " << (int) _objStatus << std::endl;

  os << "  _minCost = " <<  _minCost << std::endl;
  os << "  _maxCost = " <<  _maxCost << std::endl;

  if (printL(7))
  {
    for (VarIndexManager::const_iterator itv = _probVarManager.begin(Active, 's');
         itv != _probVarManager.end(Active, 's'); ++itv)
    {
      (*itv)->print(os);
    }
    for (VarIndexManager::const_iterator itv = _probVarManager.begin(Active, 'd');
         itv != _probVarManager.end(Active, 'd'); ++itv)
    {
      (*itv)->print(os);
    }
    for (VarIndexManager::const_iterator itv = _probVarManager.begin(Active, 'a');
         itv != _probVarManager.end(Active, 'a'); ++itv)
    {
      (*itv)->print(os);
    }

    for (ConstrIndexManager::const_iterator itc = _probConstrManager.begin(Active, 's');
         itc != _probConstrManager.end(Active, 's'); ++itc)
    {
      (*itc)->print(os);
    }
    for (ConstrIndexManager::const_iterator itc = _probConstrManager.begin(Active, 'd');
         itc != _probConstrManager.end(Active, 'd'); ++itc)
    {
      (*itc)->print(os);
    }

    if (_formulationPtr != NULL)
      _formulationPtr->printMatrix(os);
  }

  return (os);
}


bool Problem::primalSolIsFeasible()
{
  /// we use tmpVal value to store
  recordCurRhsInConstrTempVal();

  for (VarPtrSet::const_iterator varPtrIt = _inPrimalSol.begin(); varPtrIt != _inPrimalSol.end(); ++varPtrIt)
    {
      for (ConstVarConstrPtr2Double::const_iterator mapIt = (*varPtrIt)->member2coefMap().begin();
           mapIt != (*varPtrIt)->member2coefMap().end(); ++mapIt)
        mapIt->first->tmpVal(mapIt->first->tmpVal() - mapIt->second * (*varPtrIt)->val());
    }

  return checkIfConstrTempValIsFeasible();
}

bool MipProblem::primalSolIsFeasible()
{
  if (!solIsInt()) return false;

  return Problem::primalSolIsFeasible();
}

bool Problem::primalSolIsFeasible(const VarPtr2DoubleMap & curPrimalSol)
{
    recordCurRhsInConstrTempVal();

    for (VarPtr2DoubleMap::const_iterator varPt = curPrimalSol.begin(); varPt != curPrimalSol.end(); ++varPt)
    {
        for (ConstVarConstrPtr2Double::const_iterator mapIt = varPt->first->member2coefMap().begin();
             mapIt != varPt->first->member2coefMap().end(); ++mapIt)
        {
            mapIt->first->tmpVal(mapIt->first->tmpVal() - mapIt->second * varPt->second);
        }
    }

    return checkIfConstrTempValIsFeasible();
}

bool MipProblem::primalSolIsFeasible(const VarPtr2DoubleMap & curPrimalSol)
{
  if (!solIsInt(curPrimalSol))
      return false;

  return Problem::primalSolIsFeasible(curPrimalSol);
}


void Problem::recordCurRhsInConstrTempVal()
{
  ConstrIndexManager::iterator constrPtrIt;
  for (constrPtrIt = _probConstrManager.begin(Active, 's');
       constrPtrIt != _probConstrManager.end(Active, 's'); ++constrPtrIt)
    {
      (*constrPtrIt)->tmpVal((*constrPtrIt)->curRhs());
    }
  for (constrPtrIt = _probConstrManager.begin(Active, 'd');
       constrPtrIt != _probConstrManager.end(Active, 'd'); ++constrPtrIt)
    {
      (*constrPtrIt)->tmpVal((*constrPtrIt)->curRhs());
    }
  return;
}

bool Problem::checkIfConstrTempValIsFeasible()
{
  ConstrIndexManager::iterator constrPtrIt;

  for (constrPtrIt = _probConstrManager.begin(Active, 's');
       constrPtrIt != _probConstrManager.end(Active, 's'); ++constrPtrIt)
    {
      if (((*constrPtrIt)->sense() == 'E') && ((*constrPtrIt)->tmpVal() != 0))
        return false;
      if (((*constrPtrIt)->sense() == 'G') && ((*constrPtrIt)->tmpVal() > 0))
        return false;
      if (((*constrPtrIt)->sense() == 'L') && ((*constrPtrIt)->tmpVal() < 0))
        return false;
    }
  for (constrPtrIt = _probConstrManager.begin(Active, 'd');
      constrPtrIt != _probConstrManager.end(Active, 'd'); ++constrPtrIt)
    {
      if (((*constrPtrIt)->sense() == 'E') && ((*constrPtrIt)->tmpVal() != 0))
        return false;
      if (((*constrPtrIt)->sense() == 'G') && ((*constrPtrIt)->tmpVal() > 0))
        return false;
      if (((*constrPtrIt)->sense() == 'L') && ((*constrPtrIt)->tmpVal() < 0))
        return false;
    }

  return true;
}

bool Problem::solIsInt()
{
  for (VarPtrSet::const_iterator varPtrIt = _inPrimalSol.begin(); varPtrIt != _inPrimalSol.end(); ++varPtrIt)
    {
      if (((*varPtrIt)->type() == 'B') || ((*varPtrIt)->type() == 'I'))
        {
          if ((*varPtrIt)->val().fractional())
            return false;
        }
    }

  return (true);
}

bool Problem::solIsInt(const VarPtr2DoubleMap & curPrimalSol)
{
  for (VarPtr2DoubleMap::const_iterator varPt = curPrimalSol.begin(); varPt != curPrimalSol.end(); ++varPt)
    {
      if ((varPt->first->type() == 'B') || (varPt->first->type() == 'I'))
        {
          if ((varPt->second).fractional())
            return false;
        }
    }

  return (true);
}


const Double & Problem::objVal() const
{
  return _objVal;
}

const Double & Problem::primalBound() const
{
  return _primalBound;
}

void Problem::setPrimalLpBound(const Double & bd)
{
    _primalBound  = bd;

  return;
}

const Double & Problem::dualBound() const
{
  return _dualBound;
}

void Problem::setDualBound(const Double & bd)
{
  _dualBound  = bd;

  return;
}

const Double & Problem::partialSolutionValue() const
{
  return _partialSolutionValue;
}

void Problem::updatePartialSolution(Variable * varPtr, const Double & value)
{
  _partialSolutionValue += varPtr->costrhs() * value;
  VarPtr2DoubleMap::iterator mapIt = _partialSolution.find(varPtr);
  if (mapIt != _partialSolution.end())
    mapIt->second += value;
  else
    _partialSolution[varPtr] = value;
  if (varPtr->isTypeOf(VcId::MastColumnMask))
    varPtr->incrParticipation(15);

  return;
}

std::ostream& Problem::printPartialSolution(std::ostream& os) const
{
  os << "Problem::printPartialSolution: _partialSolutionValue = " << _partialSolutionValue << std::endl;

  for (VarPtr2DoubleMap::const_iterator mapIt = _partialSolution.begin(); mapIt != _partialSolution.end(); ++mapIt)
  {
    os << "    var "  << mapIt->first->name() << " is used " << mapIt->second <<  std::endl;
  }

  return (os);
}

const VarPtr2DoubleMap & Problem::partialSolution() const
{
  return (_partialSolution);
}

std::ostream& Problem::print(std::ostream& os) const
{
  os << "Problem formPtr = " << _formulationPtr << std::endl;
  if (printL(6))
    printProb(os);

  printForm(os);
  printSol(os);
  printDualSol(os);

  return (os);
}

std::ostream& Problem::printForm(std::ostream& os) const
{

  if  (_formulationPtr != NULL)
    _formulationPtr->printForm(os);

  return (os);
}

const ConstrIndexManager & Problem::probConstrSet() const
{
  return _probConstrManager;
}

ConstrIndexManager & Problem::probConstrSet()
{
  return _probConstrManager;
}

const ConstrPtrSet & Problem::probNonLinearConstrSet() const
{
  return _probNonLinearConstrSet;
}

const VarIndexManager & Problem::probVarSet() const
{
  return (_probVarManager);
}

VarIndexManager & Problem::probVarSet()
{
  return (_probVarManager);
}


void Problem::clearRecordedSol()
{
  while(!_recordedSolPtr.empty())
    {
      delete _recordedSolPtr.back();
      _recordedSolPtr.pop_back();
    }

  _recordedSolPtr.clear();

  return;
}

const SolutionMethod & Problem::solMode() const
{
  return _solMode;
}

void Problem::solverOracleFunctorPtr(BcSolverOracleFunctor * solverOracleFunctPointer)
{
  _solverOracleFunctorPtr = solverOracleFunctPointer;
  _solverOracleFunctorIsDefined = true;
}

const BcSolverOracleFunctor * Problem::solverOracleFunctorPtr()
{
    return _solverOracleFunctorPtr;
}

bool Problem::prepareSolverOracleFunctor()
{
  if (_solverOracleFunctorIsDefined)
    return _solverOracleFunctorPtr->prepareSolver();
  return true;
}

bool Problem::solverOracleFunctorDefined()
{
  return _solverOracleFunctorIsDefined;
}

bool Problem::getEnumeratedStatus()
{
    if (_solverOracleFunctorIsDefined)
        return _solverOracleFunctorPtr->getEnumeratedStatus();
    return false;
}

void Problem::masterHeuristicFunctorPtr(BcMasterHeuristicFunctor * masterHeuristicFunctPointer)
{
  _masterHeuristicFunctorPtr = masterHeuristicFunctPointer;
}

void Problem::fracSolBasedHeuristicFunctorPtr(BcFracSolBasedHeuristicFunctor * heurFunctorPtr)
{
  _fracSolBasedHeuristicFunctorPtr = heurFunctorPtr;
}

void Problem::divingFixingFunctorPtr(BcDivingFixingFunctor * divingFixingFunctPointer)
{
  _divingFixingFunctorPtr = divingFixingFunctPointer;
}

void Problem::enumSolBasedHeuristicFunctorPtr(BcEnumSolBasedHeuristicFunctor * enumSolBasedHeuristicFunctor)
{
    _enumSolBasedHeuristicFunctor = enumSolBasedHeuristicFunctor;
}

bool Problem::divingFixingFunctorDefined()
{
  return (_divingFixingFunctorPtr != NULL);
}

bool Problem::enumSolBasedHeuristicFunctorDefined()
{
    return (_enumSolBasedHeuristicFunctor != NULL);
}


void Problem::runDivingColCutGenTerminatedFunctor(const SolutionVarInfoPtrList & primalSol,
                                                  const VarPtrSet & tabuVarSet,
                                                  std::vector<MastColumn *> & columnsToAdd)
{
  if (_divingFixingFunctorPtr == NULL)
    return;

  std::vector<std::pair<BcSolution, double> > fixedSolution;
  for (VarPtr2DoubleMap::const_iterator mapIt = _partialSolution.begin(); mapIt != _partialSolution.end(); ++mapIt)
    if (mapIt->first->isTypeOf(VcId::MastColumnMask))
      {
        MastColumn * colPtr = static_cast<MastColumn *>(mapIt->first);
        Solution * spSolPtr = colPtr->spSol();
        fixedSolution.push_back(std::make_pair(BcSolution(spSolPtr), mapIt->second));
      }

  std::vector<std::pair<BcSolution, double> > masterSolution;
  for (SolutionVarInfoPtrList::const_iterator solIt = primalSol.begin(); solIt != primalSol.end(); ++solIt)
    if ((*solIt)->varPtr->isTypeOf(VcId::MastColumnMask) && !tabuVarSet.count((*solIt)->varPtr))
      {
        MastColumn * colPtr = static_cast<MastColumn *>((*solIt)->varPtr);
        Solution * spSolPtr = colPtr->spSol();
        masterSolution.push_back(std::make_pair(BcSolution(spSolPtr), (*solIt)->value));
      }

  resetSolution();
  _primalSolPtr = new Solution(_probConfPtr);

  BcSolution bcSol(_primalSolPtr);
  _divingFixingFunctorPtr->colCutGenTerminated(BcFormulation(_probConfPtr), fixedSolution, masterSolution, bcSol);
  Solution * solPtr = (Solution *) _primalSolPtr;
  do
    {
      ColGenSpConf * cgSpConfPtr = dynamic_cast<ColGenSpConf *>(solPtr->probConfPtr());
      if (cgSpConfPtr != NULL)
        {
          MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false);
          cgSpConfPtr->clearColPtrList4Insertion();
          if (!colPtr->suitableForResidualProb())
            {
              if (printL(5))
                std::cout << "Problem::runDivingColCutGenTerminatedFunctor(): generated solution is not feasible"
                          << std::endl;
              return;
            }
          columnsToAdd.push_back(colPtr);
        }
      solPtr = solPtr->nextSolPtr();
    }
  while (solPtr != NULL);
}

bool Problem::runEnumSolBasedHeuristicFunctor(const std::vector<Solution *> & enumSolPts,
                                              const Bound & incPrimalIpBound)
{
    if (!enumSolBasedHeuristicFunctorDefined())
        return false;

    std::vector<BcSolution> enumSolution;
    enumSolution.reserve(enumSolPts.size());
    for (auto * solPtr : enumSolPts)
        enumSolution.push_back(BcSolution(solPtr));
    std::vector<std::pair<int, double> > solution;
    Solution * addSolutionPtr = new Solution(_probConfPtr);
    BcSolution addBcSolution(addSolutionPtr);
    if ((*_enumSolBasedHeuristicFunctor)(BcFormulation(_probConfPtr), incPrimalIpBound.val(), enumSolution,
                                         solution,addBcSolution))
    {
        _objVal = 0;
        for (auto & pair : solution)
        {
            if (pair.first >= 0 && pair.first < enumSolPts.size())
            {
                Solution * solPtr = enumSolPts[pair.first];
                ColGenSpConf * cgSpConfPtr = static_cast<ColGenSpConf *>(solPtr->probConfPtr());
                MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false, 1);
                colPtr->val(pair.second);
                _inPrimalSol.insert(colPtr);
                _objVal += colPtr->costrhs() * pair.second;
            }
        }
        Solution * nextSolPtr = addSolutionPtr->nextSolPtr();
        while (nextSolPtr != nullptr)
        {
            ColGenSpConf * cgSpConfPtr = static_cast<ColGenSpConf *>(nextSolPtr->probConfPtr());
            MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(nextSolPtr, false, 1);
            colPtr->val(nextSolPtr->multiplicity());
            _inPrimalSol.insert(colPtr);
            _objVal += colPtr->costrhs() * nextSolPtr->multiplicity();
            nextSolPtr = nextSolPtr->nextSolPtr();
        }
        return true;
    }
    addSolutionPtr->deleteSolutionsChain();
    delete addSolutionPtr;
    return false;
}

Solution * Problem::runDivingFixingFunctor(const std::list<VariableSolInfo> & partialSol,
                                           const SolutionVarInfoPtrList & primalSol, const VarPtrSet & tabuVarSet)
{
  if (_divingFixingFunctorPtr == NULL)
    return NULL;

  BcDivingFixingFunctorInput in;
  BcDivingFixingFunctorOutput out;

  std::vector<SolutionVarInfo *> solInfoVector;
  for (SolutionVarInfoPtrList::const_iterator solIt = primalSol.begin(); solIt != primalSol.end(); ++solIt)
  {
    if ((*solIt)->varPtr->isTypeOf(VcId::MastColumnMask))
    {
        MastColumn * colPtr = static_cast<MastColumn *>((*solIt)->varPtr);
        Solution * spSolPtr = colPtr->spSol();
        if (tabuVarSet.count((*solIt)->varPtr))
        {
            in.colsInMasterSolAndTabouList.push_back(std::make_pair(BcSolution(spSolPtr), (*solIt)->value));
        }
        else
        {
            in.colsInMasterSolution.push_back(std::make_pair(BcSolution(spSolPtr), (*solIt)->value));
        }

        solInfoVector.push_back(*solIt);
    }
    else if (((*solIt)->varPtr->isTypeOf(VcId::InstMasterVarMask)))
    {
        InstanciatedVar * ivarPtr = static_cast<InstanciatedVar *>((*solIt)->varPtr);
        if (tabuVarSet.count((*solIt)->varPtr))
        {
            in.pureMastVarsInMasterSolAndTabouList.push_back(std::make_pair(BcVar(ivarPtr), (*solIt)->value));
        }
        else
        {
            in.pureMastVarsInMasterSolution.push_back(std::make_pair(BcVar(ivarPtr), (*solIt)->value));
        }
    }
  }

  for (std::list<VariableSolInfo>::const_iterator solIt = partialSol.begin(); solIt != partialSol.end(); ++solIt)
  {
    if (solIt->varPtr->isTypeOf(VcId::MastColumnMask))
    {
        MastColumn * colPtr = static_cast<MastColumn *>(solIt->varPtr);
        Solution * spSolPtr = colPtr->spSol();
        in.colsInFixedSolution.push_back(std::make_pair(BcSolution(spSolPtr), solIt->value));
    }
    else if (solIt->varPtr->isTypeOf(VcId::InstMasterVarMask))
    {
        InstanciatedVar * ivarPtr = static_cast<InstanciatedVar *>(solIt->varPtr);
        in.pureMastVarsInFixedSolution.push_back(std::make_pair(BcVar(ivarPtr), solIt->value));
    }
  }

  in.spPtr = BcFormulation(_probConfPtr);

  (*_divingFixingFunctorPtr)(in, out);

  Solution * solPtr = new Solution(_probConfPtr);
  if (printL(0))
    std::cout << "Functor selected variables : ";
  bool firstWasPrint = false;
  for (std::vector<std::pair<int, int> >::iterator pairIt = out.colsToFix.begin(); pairIt != out.colsToFix.end();
       ++pairIt)
    if ((pairIt->first >= 0) && (pairIt->first < in.colsInMasterSolution.size()))
      {
        SolutionVarInfo * solInfoPtr = solInfoVector[pairIt->first];
        if (printL(0))
        {
            if (firstWasPrint)
                std::cout << "; ";
            else
                firstWasPrint = true;
            std::cout << solInfoPtr->varPtr->name() << " (value = " << solInfoPtr->value
                      << ", rnd.value = " << pairIt->second << ", cost = " << solInfoPtr->varPtr->costrhs() << ")"
                      << ", red.cost = " << solInfoPtr->reducedCost << ")";
        }
        solPtr->includeVar(solInfoPtr->varPtr, pairIt->second, true);
      }
  for (std::vector<std::pair<BcVar, int> >::iterator pairIt = out.pureMastVarsToFix.begin();
       pairIt != out.pureMastVarsToFix.end(); ++pairIt)
  {
      BcVar var = pairIt->first;
      InstanciatedVar * varPtr = (InstanciatedVar *) var;
      if (printL(0))
      {
          if (firstWasPrint)
              std::cout << "; ";
          else
              firstWasPrint = true;
          std::cout << varPtr->name() << "( rnd.value = " << pairIt->second << ", cost = " << varPtr->costrhs() << ")";
      }
      solPtr->includeVar(varPtr, pairIt->second, true);
   }
    if (printL(0))
    {
        if (!firstWasPrint)
            std::cout << "none";
        std::cout << std::endl;
    }

  if (solPtr->solVarValMap().empty())
    {
      delete solPtr;
      solPtr = NULL;
    }

  return solPtr;
}

bool Problem::runFracSolBasedHeuristicFunctor(const std::list<VariableSolInfo> & partialSol,
                                              const SolutionVarInfoPtrList & primalSol)
{
  if (_fracSolBasedHeuristicFunctorPtr == NULL)
    return false;

  if (printL(5))
    std::cout << "Problem::runFracSolBasedHeuristicFunctor(): frac.sol.based heuristic functor is run" << std::endl;

  std::vector<std::pair<BcSolution, double> > fixedSolution;
  BcFracSolBasedHeuristicFunctorInput input;
  input.spPtr = BcFormulation(probConfPtr());
  for (SolutionVarInfoPtrList::const_iterator solIt = primalSol.begin(); solIt != primalSol.end(); ++solIt)
    {
      if ((*solIt)->varPtr->isTypeOf(VcId::MastColumnMask))
        {
          MastColumn * colPtr = static_cast<MastColumn *>((*solIt)->varPtr);
          Solution * spSolPtr = colPtr->spSol();
          input.colsInMasterSolution.push_back(std::make_pair(BcSolution(spSolPtr), (*solIt)->value));
        }
      else if ((*solIt)->varPtr->isTypeOf(VcId::InstanciatedVarMask))
        {
          InstanciatedVar * varPtr = static_cast<InstanciatedVar *>((*solIt)->varPtr);
          input.pureMastVarsInMasterSolution.push_back(std::make_pair(BcVar(varPtr), (*solIt)->value));
        }
    }
    
  for (std::list<VariableSolInfo>::const_iterator solIt = partialSol.begin(); solIt != partialSol.end(); ++solIt)
    if (solIt->varPtr->isTypeOf(VcId::MastColumnMask))
      {
        MastColumn * colPtr = static_cast<MastColumn *>(solIt->varPtr);
        Solution * spSolPtr = colPtr->spSol();
        input.colsInFixedSolution.push_back(std::make_pair(BcSolution(spSolPtr), solIt->value));
      }
    else if (solIt->varPtr->isTypeOf(VcId::InstanciatedVarMask))
      {
        InstanciatedVar * varPtr = static_cast<InstanciatedVar *>(solIt->varPtr);
        input.pureMastVarsInFixedSolution.push_back(std::make_pair(BcVar(varPtr), solIt->value));
      }

  resetSolution();
  _primalSolPtr = new Solution(_probConfPtr);
  BcSolution bcSol(_primalSolPtr);

  bool status = (*_fracSolBasedHeuristicFunctorPtr)(input, bcSol);
  if (!status)
    {
      if (printL(5))
        std::cout << "Problem::runFracSolBasedHeuristicFunctor(): frac.sol.based heur. functor "
                  << "did not return a solution" << std::endl;
      return false;
    }

  _objVal = 0;
  Solution * solPtr = (Solution *) _primalSolPtr;

  /// we first copy master variables from solChainPtr to SolPtr
  for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
       mapIt != solPtr->solVarValMap().end(); ++mapIt)
    {
      if (!mapIt->second.isZero())
        {
          if ((mapIt->second < mapIt->first->curLb()) || (mapIt->second > mapIt->first->curUb()))
            {
              if (printL(5))
                std::cout << "Problem::runFracSolBasedHeuristicFunctor(): generated solution is not feasible"
                          << std::endl;
            }
          mapIt->first->val(mapIt->second);
          _inPrimalSol.insert(mapIt->first);
          _objVal += mapIt->first->costrhs() * mapIt->first->val();
        }
    }
  
  /// now we create columns from subproblem solutions and also add them to
  solPtr = solPtr->nextSolPtr();
  while (solPtr != NULL)
    {
      ColGenSpConf * cgSpConfPtr = dynamic_cast<ColGenSpConf *>(solPtr->probConfPtr());
      if ((cgSpConfPtr != NULL) && (solPtr->multiplicity() > 0))
        {
          MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false);
          cgSpConfPtr->clearColPtrList4Insertion();
          if (!colPtr->suitableForResidualProb())
            {
              if (printL(5))
                std::cout << "Problem::runFracSolBasedHeuristicFunctor(): generated solution is not feasible"
                          << std::endl;
              return false;
            }
          if (_inPrimalSol.count(colPtr))
            {
              colPtr->val(colPtr->val() + solPtr->multiplicity());
            }
          else
            {
              colPtr->val(solPtr->multiplicity());
              _inPrimalSol.insert(colPtr);
            }
          _objVal += colPtr->costrhs() * colPtr->val();
        }
      solPtr = solPtr->nextSolPtr();
    }
  _primalSolPtr->deleteSolutionsChain();

  return true;
}

bool Problem::runMasterHeuristicFunctor()
{
  if (_masterHeuristicFunctorPtr == NULL)
    return false;

  if (printL(5))
    std::cout << "Problem::runMasterHeuristicFunctor(): Master heuristic functor is run" << std::endl;

  std::vector<std::pair<BcSolution, double> > fixedSolution;
  for (VarPtr2DoubleMap::const_iterator mapIt = _partialSolution.begin(); mapIt != _partialSolution.end(); ++mapIt)
    if (mapIt->first->isTypeOf(VcId::MastColumnMask))
      {
        MastColumn * colPtr = static_cast<MastColumn *>(mapIt->first);
        Solution * spSolPtr = colPtr->spSol();
        fixedSolution.push_back(std::make_pair(BcSolution(spSolPtr), mapIt->second));
      }

  resetSolution();
  _primalSolPtr = new Solution(_probConfPtr);
  BcSolution bcSol(_primalSolPtr);

  bool status = (*_masterHeuristicFunctorPtr)(BcFormulation(_probConfPtr), fixedSolution, bcSol);
  if (!status)
    {
      if (printL(5))
        std::cout << "Problem::runMasterHeuristicFunctor(): master heuristic function did not return a solution"
                  << std::endl;
      return false;
    }

  _objVal = 0;
  Solution * solPtr = (Solution *) _primalSolPtr;
  do
    {
      ColGenSpConf * cgSpConfPtr = dynamic_cast<ColGenSpConf *>(solPtr->probConfPtr());
      if ((cgSpConfPtr != NULL) && (solPtr->multiplicity() > 0))
        {
          MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false);
          cgSpConfPtr->clearColPtrList4Insertion();
          if (!colPtr->suitableForResidualProb())
            {
              if (printL(5))
                std::cout << "Problem::runMasterHeuristicFunctor(): generated solution is not feasible" << std::endl;

              return false;
            }

          if (_inPrimalSol.count(colPtr))
            {
              colPtr->val(colPtr->val() + solPtr->multiplicity());
            }
          else
            {
              colPtr->val(solPtr->multiplicity());
              _inPrimalSol.insert(colPtr);
            }
          _objVal += colPtr->costrhs() * colPtr->val();
        }
      else
        {
          for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
               mapIt != solPtr->solVarValMap().end(); ++mapIt)
            {
              if (!mapIt->second.isZero())
                {
                  if ((mapIt->second < mapIt->first->curLb()) || (mapIt->second > mapIt->first->curUb()))
                    {
                      if (printL(5))
                        std::cout << "Problem::runMasterHeuristicFunctor(): generated solution is not feasible" << std::endl;
                      return false;
                    }
                  mapIt->first->val(mapIt->second);
                  _inPrimalSol.insert(mapIt->first);
                  _objVal += mapIt->first->costrhs() * mapIt->first->val();
                }
            }
        }

      solPtr = solPtr->nextSolPtr();
    }
  while (solPtr != NULL);
  
  _primalSolPtr->deleteSolutionsChain();
  
  return true;
}

void Problem::retrieveBasis(LpBasisRecord * curBasisPtr, const bool markInVars, const bool markInConstrs)
{
  if (curBasisPtr == NULL)
  {
    if (printL(5))
      std::cout  << "Problem::retrieveBasis():  UNDEFINED POINTER TO BASIS" << std::endl;
    curBasisPtr = new LpBasisRecord();
  } else
  {
    curBasisPtr->clear();
  }

  _formulationPtr->retrieveBasis(*curBasisPtr, markInVars, markInConstrs);

  if (printL(5))
    std::cout  << "Problem::retrieveBasis():  CALL TO RETRIEVE BASIS " << *curBasisPtr << std::endl;

  return;
}

void Problem::reloadMemorizedBasis(LpBasisRecord * basisPtr)
{
  if (basisPtr != NULL)
    {
      if (printL(5))
        std::cout  << "Problem::reloadMemorizedBasis():  CALL TO RELOAD BASIS " << *basisPtr << std::endl;

      _formulationPtr->loadBasis(*basisPtr);
  }
}

MipProblem::MipProblem(const int & ref,
                       const double & rightHandSideZeroTol,
                       const double & reducedCostTolerance,
                       const BcObjStatus::MinMaxIntFloat & minmaxStatus,
                       const SolutionMethod & solMode,
                       const std::string & name,
                       const SolutionStatus & LPrequiredStatus,
                       const bool & LPpreprocessorOn,
                       const bool & LPprobingOn,
                       const SolutionStatus & MIPrequiredStatus,
                       const bool & MIPpreprocessorOn,
                       const bool & MIPprobingOn,
                       const bool & MIPautomaticCuttingPlanesOn,
                       const char & MIPsolverSelect) :
    Problem(ref,
            rightHandSideZeroTol,
            reducedCostTolerance,
            minmaxStatus,
            solMode,
            name,
            LPrequiredStatus,
            LPpreprocessorOn,
            LPprobingOn,
            MIPsolverSelect),
    _mipRequiredStatus(MIPrequiredStatus),
    _mipProbStatus(SolutionStatus::UnSolved),
    _MIPpreprocessorOn(MIPpreprocessorOn),
    _MIPprobingOn(MIPprobingOn),
    _MIPautomaticCuttingPlanesOn(MIPautomaticCuttingPlanesOn),
    _MIPsolverSelect(MIPsolverSelect)
{
  return;
}

bool MipProblem::solveProbMIP(const char & flag, const bool & ifPrint)
{
  bool foundSol(false);
  bapcodInit().check(_primalFormulationPtr == NULL,
                     "MipProblem::solveProb(): _solMode == lp or mipSolver => requires  defined formulation");

  foundSol = ((MIPform *)_primalFormulationPtr)->solve(param().MasterMipSolverBarrierConvergenceTolerance()
                                          , _rightHandSideZeroTol
                                          , _reducedCostTolerance
                                          , flag
                                          , ifPrint
                                          , _mipRequiredStatus
                                          , _objVal
                                          , _primalBound
                                          , _dualBound
                                          , _inPrimalSol
                                          , _inDualSol
                                          , _MIPpreprocessorOn
                                          , _MIPprobingOn
                                          , _MIPautomaticCuttingPlanesOn
                                          , _MIPsolverSelect);

  if (printL(1))
  {
    printSolVal();
  }

  setProbStatus(_primalFormulationPtr->status());

  if (_primalFormulationPtr->status().count(SolutionStatus::Optimum))
    _dualBound = _primalBound;

  if (printL(5))
    std::cout  << "MipProblem::solveProbMIP(): probStatus() after _primalFormulationPtr->solve()" << probStatus()
               << ", _requiredStatus= " << _requiredStatus << std::endl;
  return foundSol;
}

int MipProblem::solveProb(int & maxLevelOfRestriction, const char & flag, const bool & ifPrint)
{
  int solverReturnStatus(0);
  resetSolution(flag);

  switch (_solMode.status()) {
  case SolutionMethod::none: break;
  case SolutionMethod::lpSolver:
  {
    solverReturnStatus = solveProbLP(flag, ifPrint);
    break;
  }
  case SolutionMethod::mipSolver:
  {
    solverReturnStatus = solveProbMIP(flag, ifPrint);
    break;
  }
  case SolutionMethod::customSolver:
  {
    if (printL(5))
      std::cout  << "Problem::solveProb(): to enter customizedSolver()  "
            << std::endl;

      if (_primalSolPtr == NULL)
        _primalSolPtr = new Solution(_probConfPtr);
      if (_dualSolPtr == NULL)
        _dualSolPtr = new DualSolution(_probConfPtr);
      BcSolution bcSol(_primalSolPtr);
      BcDualSolution bcDualSol(_dualSolPtr);

      solverReturnStatus = customizedSolver(maxLevelOfRestriction,
                                            _requiredStatus,
                                            _objVal,
                                            _primalBound,
                                            _dualBound,
                                            _inPrimalSol,
                                            _inDualSol,
                                            bcSol,
                                            bcDualSol);

    if (printL(1))
      std::cout  << "MipProblem::solveProb(): " << name()
      << " _objVal = " << _objVal
            << std::endl;

    break;
  }
  case SolutionMethod::custom2mipSolver:
  {
    if (printL(5))
      std::cout  << "Problem::solveProb(): to enter customizedSolver()  " << std::endl;

    if (maxLevelOfRestriction >= 1)
      {
          if (_primalSolPtr == NULL)
              _primalSolPtr = new Solution(_probConfPtr);
          if (_dualSolPtr == NULL)
              _dualSolPtr = new DualSolution(_probConfPtr);
          BcSolution bcSol(_primalSolPtr);
          BcDualSolution bcDualSol(_dualSolPtr);

          solverReturnStatus = customizedSolver(maxLevelOfRestriction,
                                                _requiredStatus,
                                                _objVal,
                                                _primalBound,
                                                _dualBound,
                                                _inPrimalSol,
                                                _inDualSol,
                                                bcSol,
                                                bcDualSol);
      }
    else
      {
        solverReturnStatus = solveProbMIP(flag, ifPrint);
      }

    if (printL(1))
      std::cout  << "MipProblem::solveProb(): " << name() << " _objVal = " << _objVal << std::endl;

    break;
  }
  default:
  {
    bapcodInit().check(true, "Problem solMode undefined");
    break;
  }
  }

  setStatusAfterSol();

  return (solverReturnStatus);
}

std::ostream& MipProblem::print(std::ostream& os) const
{
  os << "MipProblem "  << std::endl;
  Problem::print(os);
  return (os);
}

void MipProblem::resetSolution(const char & flag)
{
  _mipProbStatus = SolutionStatus::UnSolved;
  Problem::resetSolution(flag);
  return;
}

BapcodInit & Problem::bapcodInit() const
{
  return _probConfPtr->bapcodInit();
}

void Problem::updateInNonZeroRedCostVarsSet(Variable * varPtr)
{
    if (!varPtr->reducedCost().isZero())
        _nonZeroRedCostVars.insert(varPtr);
    else
        _nonZeroRedCostVars.erase(varPtr);
}

void Problem::removeVarsNotInProblemFromPrimalSolution()
{
    auto varPtrIt = _inPrimalSol.begin();
    while (varPtrIt != _inPrimalSol.end())
    {
        if ((*varPtrIt)->inCurProb())
        {
            varPtrIt++;
        }
        else
        {
            (*varPtrIt)->val(0);
            _inPrimalSol.erase(varPtrIt++);
        }
    }
}

void Problem::removeVarsNotInProblemFromNonZeroRedCostVars()
{
    VarPtrSet::const_iterator varPtrIt = _nonZeroRedCostVars.begin();
    while (varPtrIt != _nonZeroRedCostVars.end())
    {
        if ((*varPtrIt)->inCurProb())
        {
            varPtrIt++;
        }
        else {
            (*varPtrIt)->reducedCost(0);
            _nonZeroRedCostVars.erase(varPtrIt++);
        }
    }
}

void Problem::getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions)
{
    auto * cgSpConf = dynamic_cast<ColGenSpConf *>(_probConfPtr);
    if (cgSpConf != nullptr && _solverOracleFunctorIsDefined && _solverOracleFunctorPtr->getEnumeratedStatus())
    {
        int numEnumSolutions = _solverOracleFunctorPtr->getNumberOfEnumeratedSolutions();
        auto * solPtr = new Solution(cgSpConf);
        std::vector<double> reducedCosts;
        BcSolution bcSol(solPtr);
        _solverOracleFunctorPtr->getEnumeratedSolutions(BcFormulation(cgSpConf), numEnumSolutions, bcSol,
                                                        reducedCosts);
        int solId = 0;
        while (solPtr != nullptr)
        {
            solPtr->resetCost();
            if (!solPtr->solVarValMap().empty())
                enumSolutions.push_back(std::make_tuple(reducedCosts[solId] + cgSpConf->fixedDualCost(), solPtr->cost(),
                                                        BcSolution(solPtr)));
            solPtr = solPtr->nextSolPtr();
            solId += 1;
        }
    }
}
