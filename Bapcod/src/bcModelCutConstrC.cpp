/**
 *
 * This file bcModelCutConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcModelingLanguageC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcProbConfigC.hpp"
#include "bcModelNetworkFlow.hpp"

BcCutConstrArray::BcCutConstrArray():
  BcConstrArray(), _genericCutConstrPtr(NULL)
{

}

BcCutConstrArray::BcCutConstrArray(const BcFormulation & formulation,
			                       const std::string & name,
                                   const char & type,
                                   const double & rootPriorityLevel,
			                       const double & nonRootPriorityLevel,
			                       const bool & toBeUsedInPreprocessing):
  BcConstrArray(), _genericCutConstrPtr(NULL)
{

  if (printL(5))
    std::cout << " BcCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
              << " BcCutConstrArray =  " << name << std::endl;

  _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr(name);

  if (_genericCutConstrPtr == NULL)
    {
      if (printL(5))
        std::cout << " BcCutConstrArray() : need to create cut  "
                  << std::endl;

      _genericCutConstrPtr
        = formulation.probConfPtr()->modelPtr()->createGenericCut(formulation.probConfPtr(), name, type,
                                                                  SelectionStrategy::MostViolated,
                                                                  nonRootPriorityLevel, rootPriorityLevel, 'G', 0.0,
                                                                  toBeUsedInPreprocessing);
    }
  _genericConstrPtr = _genericCutConstrPtr;

}


BcCutConstrArray::~BcCutConstrArray()
{
    return;
}

void BcCutConstrArray::setRootPriorityLevel(const double & rootPriorityLevel_)
{
  if (_genericCutConstrPtr == NULL)
  {
    if (printL(5))
      std::cout << "BaPCod info :  Model _genericCutConstrPtr" << std::endl;
    return;
  }
  _genericCutConstrPtr->rootPriorityLevel(rootPriorityLevel_);
}


void BcCutConstrArray::attach(BcCutSeparationFunctor * separationRoutinePtr)
{
  _genericCutConstrPtr->cutSeparationFunctorPtr(separationRoutinePtr);
}


/// BcCutSeparationFunctor parameters are:
/// @param formPtr is the current formulation on which the cut separation is called
/// @param primalSol is the current solution  on which the cut separation is called; in the case of a master
///                  formulation, primalSol is infact a chain of solution (uste next() to follow the chain);
///                  the first one includes global cost and pure master varaibles, the following one containt
///                  the solution for each subproblme type; for a given subproblem (type), all the solutions associated
///                  with this subproblem (type) have been aggregated inot a signel solution my summing the values taken
///                  in each solution.
/// @param maxViolation as input value equal to the cut violation tolerance of BaPCod; as output value, the user should
///                     provide the maxViolation of all inequality found (althought this information is not used by
///                     bapcod at the moment
/// @param cutList is a container in which to push the constraint generated as cut
/// @return the return value is the number of cut founds (can be less than the size of the list of some of the
///         constraint do not defined violated cut but the user wishes to add them in the pool anyway.

int BcCutSeparationFunctor::operator() (BcFormulation formPtr,
				                        BcSolution & primalSol,
                                        double & maxViolation,
				                        std::list< BcConstr > & cutList)
{
  std::cout << " BcCutSeparationFunctor::operator()  SHOULD NOT BE CALLED" << std::endl;
  return 0;
}

BcCustomNonLinearCut::BcCustomNonLinearCut(CustomNonLinearCut * cutPtr):
  BcConstr(cutPtr), _custNonLinCutPtr(cutPtr)
{
}

const BcCustomNonLinearCutInfo * BcCustomNonLinearCut::cutInfoPtr() const
{
  return _custNonLinCutPtr->_cutInfoPtr;
}

BcCustomNonLinearCutArrayFunctor::BcCustomNonLinearCutArrayFunctor(const BcFormulation & formulation,
                                                                   const std::string & name,
                                                                   const char & type,
                                                                   const SelectionStrategy & priorityRule,
                                                                   const double & priorityLevel):
  BcCutConstrArray(), _numberOfCuts(0)
{
  if (printL(5))
    std::cout << " BcCustomNonLinearCutArrayFunctor() : ProbConfig = "
              << formulation.probConfPtr()->name()
              << " BcCustomNonLinearCutArrayFunctor = " << name << std::endl;

  _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr(name);

  if (_genericCutConstrPtr == NULL)
    {
      if (printL(5))
        std::cout << " BcCustomNonLinearCutArrayFunctor() : need to create cut"
                  << std::endl;

      _genericCutConstrPtr =
        formulation.probConfPtr()->modelPtr()->createGenericCustomNonLinearCutConstr(this, formulation.probConfPtr(),
                                                                                     name, type, priorityRule,
                                                                                     priorityLevel, priorityLevel);
    }
  _genericConstrPtr = _genericCutConstrPtr;
}

BcCustomNonLinearCutArrayFunctor::~BcCustomNonLinearCutArrayFunctor()
{
}

BcCustomNonLinearCut BcCustomNonLinearCutArrayFunctor::createNewCut(BcCustomNonLinearCutInfo * cutInfoPtr)
{
  GenericCustomNonLinearCutConstr * genCustNonLinCutConstrPtr
    = static_cast<GenericCustomNonLinearCutConstr *>(_genericCutConstrPtr);
  std::string name(_genericConstrPtr->defaultName());
  MultiIndex newCutId(_numberOfCuts++);
  newCutId.appendRef2name(name, _genericCutConstrPtr->multiIndexNames());

  CustomNonLinearCut * newCut = new CustomNonLinearCut(MultiIndex(_numberOfCuts++), genCustNonLinCutConstrPtr,
                                                       _genericConstrPtr->probConfPtr(), name, cutInfoPtr);
  return BcCustomNonLinearCut(newCut);
}

int BcCustomNonLinearCutArrayFunctor::cutSeparationRoutine(BcFormulation formPtr,
                                                           BcSolution & projectedSol,
                                                           std::list<std::pair<double, BcSolution> > & columnsInSol,
                                                           const double & violationTolerance,
                                                           std::list<BcCustomNonLinearCut> & cutList)
{
  return 0;
}

double BcCustomNonLinearCutArrayFunctor::getCoefficient(const BcCustomNonLinearCut & cut,
                                                        const BcSolution & spSol)
{
  return 0.0;
}

BcSoftConflictsCut::BcSoftConflictsCut(SoftConflictsCut * cutPtr) :
  BcConstr(cutPtr), _cutPtr(cutPtr)
{
}

double BcSoftConflictsCut::getDualVal() const
{
    return _cutPtr->valOrSepPointVal();
}

int BcSoftConflictsCut::type() const
{
    return _cutPtr->cutType();
}

const std::vector<std::pair<BcVar, BcVar> > BcSoftConflictsCut::conflicts() const
{
    std::vector<std::pair<BcVar, BcVar> > conflictVector;
    conflictVector.reserve(_cutPtr->conflicts().size());
    std::vector<std::pair<SubProbVariable *, SubProbVariable *> >::const_iterator pairIt;
    for (pairIt = _cutPtr->conflicts().begin(); pairIt != _cutPtr->conflicts().end(); ++pairIt)
    {
        BcVar firstVar(pairIt->first);
        BcVar secondVar(pairIt->second);
        conflictVector.push_back(std::make_pair(firstVar, secondVar));
    }

    return conflictVector;
};

BcSoftConflictsCutArrayFunctor::BcSoftConflictsCutArrayFunctor(const BcFormulation & formulation,
                                                               const std::string & name,
                                                               const char & type,
                                                               const SelectionStrategy & priorityRule,
                                                               const double & priorityLevel) :
  BcCutConstrArray(), _numberOfCuts(0)
{
    if (printL(5))
        std::cout << " BcSoftConflictsCutArrayFunctor() : ProbConfig = "
                  << formulation.probConfPtr()->name()
                  << " BcSoftConflictsCutArrayFunctor = " << name << std::endl;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr(name);

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcSoftConflictsCutArrayFunctor() : need to create cut" << std::endl;

        _genericCutConstrPtr = formulation.probConfPtr()->modelPtr()
                                ->createGenericSoftConflictsCutConstr(this, formulation.probConfPtr(),
                                                                      name, type, priorityRule,
                                                                      priorityLevel, priorityLevel);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcSoftConflictsCutArrayFunctor::~BcSoftConflictsCutArrayFunctor()
{

}

BcConstr BcSoftConflictsCutArrayFunctor::createNewCut(const double & rhs,
                                                      std::vector<std::pair<BcVar, BcVar> > & conflicts,
                                                      const int cutType)
{
    GenericSoftConflictsCutConstr * genSoftConflictsCutConstrPtr
            = static_cast<GenericSoftConflictsCutConstr *>(_genericCutConstrPtr);
    std::string name(_genericConstrPtr->defaultName());
    MultiIndex newCutId(_numberOfCuts);
    newCutId.appendRef2name(name, _genericCutConstrPtr->multiIndexNames());

    std::vector<std::pair<SubProbVariable *, SubProbVariable *> > conflictsVector;
    conflictsVector.reserve(conflicts.size());
    for (std::vector<std::pair<BcVar, BcVar> >::iterator pairIt = conflicts.begin(); pairIt != conflicts.end();
         ++pairIt)
    {
        InstanciatedVar * firstVarPtr = (InstanciatedVar *)(pairIt->first);
        InstanciatedVar * secondVarPtr = (InstanciatedVar *)(pairIt->second);
        if (firstVarPtr->isTypeOf(VcId::SubProbVariableMask) && secondVarPtr->isTypeOf(VcId::SubProbVariableMask))
        {
            conflictsVector.push_back(std::make_pair(static_cast<SubProbVariable *>(firstVarPtr),
                                                     static_cast<SubProbVariable *>(secondVarPtr)));
        } else {
            std::cerr << "BaPCod error: a non-subproblem variable is in the conflicts passed to create "
                      << "a soft conflicts cut" << std::endl;
            exit(1);
        }
        if (firstVarPtr->probConfPtr() != secondVarPtr->probConfPtr())
        {
            std::cerr << "BaPCod error: variables in a conflict does not come from the same subproblem when creating "
                      << "a soft conflicts cut" << std::endl;
            exit(1);
        }
    }

    SoftConflictsCut * newCut = new SoftConflictsCut(MultiIndex(_numberOfCuts++), genSoftConflictsCutConstrPtr,
                                                     _genericConstrPtr->probConfPtr(), name, rhs, cutType,
                                                     conflictsVector);
    return BcConstr(newCut);
}

int BcSoftConflictsCutArrayFunctor::cutSeparationRoutine(BcFormulation formPtr,
                                                         std::list<std::pair<double, BcSolution> > & columnsInFixedSol,
                                                         std::list<std::pair<double, BcSolution> > & columnsInSol,
                                                         std::list<BcConstr> & cutList)
{
    return 0;
}

int BcSoftConflictsCutArrayFunctor
    ::cutSeparationBasedOnFixedSol(BcFormulation formPtr,
                                   std::list<std::pair<double, BcSolution> > & columnsInOldFixedSol,
                                   std::list<std::pair<double, BcSolution> > & columnsInNewFixedSol,
                                   std::list<BcConstr> & cutList)
{
    return 0;
}
