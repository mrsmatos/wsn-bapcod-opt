/**
 *
 * This file bcModelObjectiveC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcModelObjectiveC.hpp"

#include "bcModelingLanguageC.hpp"

#include "bcModelC.hpp"

BcObjective::BcObjective():
  _modelPtr(NULL)
{
}

BcObjective::BcObjective(BcModel & model):
  _modelPtr((Model *)model)
{
}

BcObjective::BcObjective(BcFormulation & bcForm):
  _modelPtr(NULL)
{
    if (bcForm.isDefined())
        _modelPtr = bcForm.probConfPtr()->modelPtr();
}

void BcObjective::setMinMaxStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense)
{
    setStatus(newObjectiveSense);
}

void BcObjective::setStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense)
{
    if (newObjectiveSense == BcObjStatus::maxFloat || newObjectiveSense == BcObjStatus::maxInt)
    {
        std::cerr << "BaPCod error : only minimization objective function is possible for the moment" << std::endl;
        exit(1);
    }
    _modelPtr->objectiveSense(newObjectiveSense);
}

const BcObjective & BcObjective::add(const BcVar & var, const double & coef)
{
  if (var._ivarPtr == NULL)
    { 
      if (printL(6))
        std::cout << "BaPCod info :  Model BcVar == NULL" << std::endl;
    }
  else
    {
      var._ivarPtr->costrhs(coef);
    }
  return *this;
}

const BcObjective & BcObjective::operator=(const BcVar & var)
{
  return BcObjective::add(var,1);
}

const BcObjective & BcObjective::operator+(const BcVar & var)
{
  return BcObjective::add(var,1);
}

const BcObjective & BcObjective::operator-(const BcVar & var)
{
  return BcObjective::add(var, -1);
}

const BcObjective & BcObjective::operator+=(const BcVar & var)
{
  return BcObjective::add(var,1);
}

const BcObjective & BcObjective::operator-=(const BcVar & var)
{
  return BcObjective::add(var, -1);
}


const BcObjective & BcObjective::operator=(const BcVarCoef & coef)
{
  return BcObjective::operator+=(coef);
}


const BcObjective & BcObjective::operator+(const BcVarCoef & coef)
{
  return BcObjective::operator+=(coef);
}

const BcObjective & BcObjective::operator-(const BcVarCoef & coef)
{
  return BcObjective::operator-=(coef);
}

const BcObjective & BcObjective::operator+=(const BcVarCoef & coef)
{
  return BcObjective::add(coef.first, coef.second);
}

const BcObjective & BcObjective::operator-=(const BcVarCoef & coef)
{
  return BcObjective::add(coef.first, - coef.second);
}

const BcObjective & BcObjective::operator=(const BcRowExpression & exp)
{
  return BcObjective::operator+=(exp);
}
  
const BcObjective & BcObjective::operator+(const BcRowExpression & exp)
{
  return BcObjective::operator+=(exp);
}

const BcObjective & BcObjective::operator+=(const BcRowExpression & exp)
{
  for(std::list < BcVarCoef >::const_iterator it = exp.begin();
      it != exp.end(); ++it)
    {
      BcObjective::add(it->first, it->second * exp._mult);
    }
  return *this;
}

const BcObjective & BcObjective::operator-(const BcRowExpression & exp)
{
  return BcObjective::operator-=(exp);
}

const BcObjective & BcObjective::operator-=(const BcRowExpression & exp)
{
  for(std::list < BcVarCoef >::const_iterator it = exp.begin();
      it != exp.end(); ++it)
    {
      BcObjective::add(it->first, - it->second * exp._mult);
    }
  return *this;
}

const BcObjective & BcObjective::operator=(BcVarIndex & modVarInd)
{
  return this->operator=( (BcVar) modVarInd); /// casting
}

const BcObjective & BcObjective::operator+(BcVarIndex & modVarInd)
{
  return this->operator+( (BcVar) modVarInd); /// casting
}

const BcObjective & BcObjective::operator-(BcVarIndex & modVarInd)
{
  return this->operator-( (BcVar) modVarInd); /// casting
}

const BcObjective & BcObjective::operator+=(BcVarIndex & modVarInd)
{
  return this->operator+=( (BcVar) modVarInd); /// casting
}

const BcObjective & BcObjective::operator-=(BcVarIndex & modVarInd)
{
  return this->operator-=( (BcVar) modVarInd); /// casting
}

const BcObjective & BcObjective::operator<=(const double & ub)
{
  if ( _modelPtr != NULL) /// global objective
    {
      if ((_modelPtr->objectiveSense() == BcObjStatus::minInt) ||
          (_modelPtr->objectiveSense() == BcObjStatus::minFloat))

        {
          if (_modelPtr->masterConfPtr() != NULL)
            {
              _modelPtr->masterConfPtr()->updatePrimalIncBound(Bound(ub, _modelPtr->objectiveSense()));
            }
          else
            {
              _modelPtr->initialPrimalBound(Bound(ub, _modelPtr->objectiveSense()));
            }
        }
      else
        {
          if (_modelPtr->masterConfPtr() != NULL)
            {
              _modelPtr->masterConfPtr()->updateDualIncBound(Bound(ub, _modelPtr->objectiveSense()));
            }
          if (ub > _modelPtr->initialDualBound())
            _modelPtr->initialDualBound(Bound(ub, _modelPtr->objectiveSense()));
        }
    }
  return *this;
}

void BcObjective::setArtCostValue(const double & cost)
{
    if ( _modelPtr != NULL) /// global objective
    {
        _modelPtr->setArtCostValue(cost);
    }
}

const BcObjective & BcObjective::addCumulatively(const BcRowExpression & exp)
{
    for(std::list < BcVarCoef >::const_iterator it = exp.begin();
        it != exp.end(); ++it)
    {
        if(it->first.isDefined())
        {
            add(it->first, it->first.originalCost() + it->second * exp._mult);
        }
    }
    return *this;
}

const BcObjective & BcObjective::operator==(const double & guessVal)
{
    setArtCostValue(guessVal);
  return *this;
}

const BcObjective & BcObjective::operator>=(const double & lb)
{
  if ( _modelPtr != NULL) /// global objective
    {
      if ((_modelPtr->objectiveSense() == BcObjStatus::maxInt) ||
          (_modelPtr->objectiveSense() == BcObjStatus::maxFloat))
        {
          if (_modelPtr->masterConfPtr() != NULL)
            {
              _modelPtr->masterConfPtr()->updatePrimalIncBound(Bound(lb, _modelPtr->objectiveSense()));
            }
          else
            {
              _modelPtr->initialPrimalBound(Bound(lb, _modelPtr->objectiveSense()));
            }
        }
      else
        {
          if (_modelPtr->masterConfPtr() != NULL)
            {
              _modelPtr->masterConfPtr()->updateDualIncBound(Bound(lb, _modelPtr->objectiveSense()));
            }

          if (lb > _modelPtr->initialDualBound())
            _modelPtr->initialDualBound(Bound(lb, _modelPtr->objectiveSense()));
        }
    }
  return *this;
}


BcObjectiveArray::BcObjectiveArray(BcModel & modelPointer, const std::string & name):
  _modelPtr(modelPointer), _curObjectivePtr(modelPointer)
{
  _modelPtr->_objectiveSense = BcObjStatus::minInt;
}

void BcObjectiveArray::setMinMaxStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense)
{
  if (_modelPtr != NULL)
      _modelPtr->objectiveSense(newObjectiveSense);
  _curObjectivePtr.setMinMaxStatus(newObjectiveSense);
}

BcObjective & BcObjectiveArray::createElement(const MultiIndex & multiIndexId)
{
  return _curObjectivePtr;
}
