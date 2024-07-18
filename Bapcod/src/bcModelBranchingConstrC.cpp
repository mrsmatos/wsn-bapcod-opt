/**
 *
 * This file bcModelBranchingConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"

BcRyanAndFosterBranchConstr::BcRyanAndFosterBranchConstr(RyanAndFosterInstSubProbBranchConstr * iconstrPtr) :
  _iconstrPtr(iconstrPtr), _bcVarI(iconstrPtr->_variPtr), _bcVarJ(iconstrPtr->_varjPtr)
{
}

const BcVar & BcRyanAndFosterBranchConstr::bcVarI() const
{
  return _bcVarI;
}

const BcVar & BcRyanAndFosterBranchConstr::bcVarJ() const
{
  return _bcVarJ;
}

bool BcRyanAndFosterBranchConstr::together() const
{
  return (_iconstrPtr->sense() == 'E');
}

BcAggrSubProbVarBranching::BcAggrSubProbVarBranching(const BcFormulation & formulation,
		                                             const std::string & genVarName,
                                                     const double & highestPriorityFraction,
		                                             const double & priorityLevel,
                                                     const int & numIgnoredIndices,
		                                             const bool & toBeUsedInPreprocessing):
  _genAggrSubProbVarBranchingConstrPtr(NULL)
{
  if (printL(5))
    std::cout << " BcAggrSubProbVarBranching() : ProbConfig =  " << formulation.probConfPtr()->name()
              << " BcAggrSubProbVarBranching =  " << genVarName << std::endl;

  GenericBranchingConstr * genBranchConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(genVarName);
  if (genBranchConstrPtr != NULL)
    _genAggrSubProbVarBranchingConstrPtr = dynamic_cast<GenAggrSubProbVarBranchingConstr *>(genBranchConstrPtr);

  if (_genAggrSubProbVarBranchingConstrPtr == NULL)
    {
      if (printL(5))
        std::cout << " BcAggrSubProbVarBranching() : need to create branching  " << std::endl;

      Model * modelPtr = formulation.probConfPtr()->modelPtr();

      _genAggrSubProbVarBranchingConstrPtr =
        new GenAggrSubProbVarBranchingConstr(modelPtr, formulation.probConfPtr(), genVarName, genVarName,
                                             highestPriorityFraction, numIgnoredIndices,
                                             SelectionStrategy::MostFractional, priorityLevel, priorityLevel,
                                             toBeUsedInPreprocessing);
      _genAggrSubProbVarBranchingConstrPtr->defaultFlag('d');
    }
}

BcAggrSubProbVarBranching::~BcAggrSubProbVarBranching()
{
}

const BcAggrSubProbVarBranching & BcAggrSubProbVarBranching
                                  ::attach(BcVarBranchingPriorityFunctor * priorityRoutinePtr)
{
  _genAggrSubProbVarBranchingConstrPtr->branchingVarPriorityFunctorPtr(priorityRoutinePtr);
  return *this;
}


BcBranchingConstrArray::BcBranchingConstrArray(const BcFormulation & formulation,
					                           const std::string & name,
					                           const SelectionStrategy & priorityRule,
					                           const double & priorityLevel,
					                           const bool & toBeUsedInPreprocessing):
  BcConstrArray(), _genericBranchingConstrPtr(NULL)
{

  if (printL(5))
    std::cout << " BcBranchingConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
              << " BcBranchingConstrArray =  " << name << std::endl;

  _genericBranchingConstrPtr = formulation.probConfPtr()->getGenericBranchingConstr(name);

  if (_genericBranchingConstrPtr == NULL)
    {
      if (printL(5))
        std::cout << " BcBranchingConstrArray() : need to create branching  "
                  << std::endl;

      Model * modelPtr = formulation.probConfPtr()->modelPtr();

      _genericBranchingConstrPtr = modelPtr->createGenericBranching(formulation.probConfPtr(),
                                                                    name,
                                                                    'F',
                                                                    priorityRule,
                                                                    priorityLevel,
                                                                    priorityLevel,
                                                                    'G',
                                                                    0,
                                                                    toBeUsedInPreprocessing);
    }
  _genericConstrPtr = _genericBranchingConstrPtr;
}


BcBranchingConstrArray::~BcBranchingConstrArray()
{
}

GenericBranchingConstr * BcBranchingConstrArray::genericBranchingConstrPtr() const {return _genericBranchingConstrPtr;}


const BcBranchingConstrArray & BcBranchingConstrArray
                               ::attach(BcDisjunctiveBranchingConstrSeparationFunctor * separationRoutinePtr)
{
  _genericBranchingConstrPtr->branchingSeparationFunctorPtr(separationRoutinePtr);
  return *this;
}


std::ostream& operator<<(std::ostream& os, const BcBranchingConstrArray & that)
{
  return that.genericBranchingConstrPtr()->print(os);
}
