/**
 *
 * This file bcModelColGenSpC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

BcColGenSpArray::BcColGenSpArray(BcModel & modelPointer, const std::string & name):
  BcFormulationArray(modelPointer, name),
  _fixedCost(0), _defaultLb(InstMastConvexityConstr::lowerBoundWhenInactive),
  _defaultUb(InstMastConvexityConstr::upperBoundWhenInactive)
{
    _modelPtr->createMasterConf(_name, MultiIndex());
}

void BcColGenSpArray::setFixedCost(const Double & fixedCost)
{
  _fixedCost = fixedCost;
}


BcFormulation & BcColGenSpArray::getElement(const MultiIndex & multiIndexId)
{
  if (printL(6))
    std::cout << " BcColGenSpArray::getElement(const MultiIndex &) is called for id" << multiIndexId << std::endl;

  if (_curForm.isDefined() &&  (_curForm.id() == multiIndexId))
	 return _curForm;

  return _curForm = BcFormulation(_modelPtr->getColGenSubproblem(_name, multiIndexId));
}


BcFormulation & BcColGenSpArray::createElement(const MultiIndex & multiIndexId)
{
    if (printL(6))
        std::cout << " BcColGenSpArray::createElement(const MultiIndex & indexArray)  IS called for id" << multiIndexId
                  << std::endl;

    if (_curForm.isDefined() && (_curForm.id() == multiIndexId))
        return _curForm;

    _curForm = BcFormulation(_modelPtr->getColGenSubproblem(_name, multiIndexId));

    if (_curForm.isDefined())
        return _curForm;

    _curForm = BcFormulation(_modelPtr->createColGenSubproblem(_name, multiIndexId, false,
                                                               _defaultUb, _defaultLb, _fixedCost));

    if (printL(6))
    {
        std::cout << " BcColGenSpArray::operator() spConfPtr->name = " << _curForm.name() << std::endl;
        std::cout << " BcColGenSpArray::operator() spConfPtr->id = " << _curForm.id() << std::endl;
    }

    return _curForm;
}


const BcColGenSpArray & BcColGenSpArray::operator<=(const double & defaultUb)
{
  _defaultUb = defaultUb;
  return *this;
}

const BcColGenSpArray & BcColGenSpArray::operator>=(const double & defaultLb)
{
  _defaultLb = defaultLb;
  return *this;
}


