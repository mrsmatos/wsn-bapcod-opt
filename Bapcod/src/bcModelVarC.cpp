/**
 *
 * This file bcModelVarC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcErrorC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcProbConfigC.hpp"
#include "bcApplicationException.hpp"

using namespace std;

BcVar::BcVar() :
    _ivarPtr(NULL)
{
}

BcVar::BcVar(InstanciatedVar * ivarPtr) :
    _ivarPtr(ivarPtr)
{
}

BcVar::BcVar(BcVarIndex & varInd) :
    _ivarPtr(NULL)
{
  _ivarPtr = varInd.first.getElement(varInd.second);
}

BcVar::~BcVar()
{
}

bool BcVar::isDefined() const
{
  return (_ivarPtr != NULL);
}

const std::string & BcVar::name() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::name ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->name();
}

const std::string & BcVar::genericName() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::genericName ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->genVarPtr()->defaultName();
}

const BcFormulation BcVar::formulation() const
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    return BcFormulation(NULL);
  }
  return BcFormulation(_ivarPtr->probConfPtr());
}

const MultiIndex & BcVar::id() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::id ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->id();
}

double BcVar::curLb() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::curLb ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->curLb();
}

double BcVar::curUb() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::curUb ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->curUb();
}

double BcVar::curCost() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "BcVar::curCost ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->curCost();
}

double BcVar::originalCost() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->costrhs();
}

bool BcVar::inCurForm() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->inCurProb();
}

void BcVar::curVal(double value)
{
    if (_ivarPtr == NULL)
    {
        if (printL(6))
            std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    }
    else
        _ivarPtr->val(value);
}

double BcVar::curVal() const
{
    if (_ivarPtr == NULL)
    {
        std::cout << "ERROR Model BcVar == NULL" << std::endl;
        exit(0);
    }
    return _ivarPtr->val();
}

double BcVar::solVal() const
{
  if (_ivarPtr == NULL)
  {
    std::cout << "ERROR Model BcVar == NULL" << std::endl;
    exit(0);
  }
  return _ivarPtr->solVal();
}

BcVar & BcVar::localLb(double vlb)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
  {
    _ivarPtr->lb(vlb);
  }
  return *this;
}

BcVar & BcVar::localUb(double vub)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
  {
    _ivarPtr->ub(vub);
  }
  return *this;
}

BcVar & BcVar::globalLb(double vlb)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
  {
    _ivarPtr->globalLb(vlb);
  }
  return *this;
}

BcVar & BcVar::globalUb(double vub)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
  {
    _ivarPtr->globalUb(vub);
    if (_ivarPtr->ub() > vub)
    {
      _ivarPtr->ub(vub);
    }
  }
  return *this;
}

void BcVar::type(const char& flag)
{
    if (_ivarPtr == NULL)
    {
        if (printL(1))
            std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    } else
    {
        switch (flag) {
            case 'B':
            {
                _ivarPtr->ub(1);
                _ivarPtr->type('B');
                break;
            }
            case 'I':
            case 'C':
            {
                _ivarPtr->type(flag);
                break;
            }
        }
    }
}

void BcVar::sense(const char & flag)
{
    if (_ivarPtr == NULL)
    {
        if (printL(6))
            std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    } else
    {
        switch (flag) {
            case 'P':
            case 'N':
            case 'F':
            {
                _ivarPtr->type(flag);
                break;
            }
        }
    }
}

void BcVar::setImplicit()
{
    if (_ivarPtr == NULL)
    {
        if (printL(6))
            std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    } else
    _ivarPtr->kind('I');
}

BcVar & BcVar::branchingPriority(double vpriority)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
  {
    if (vpriority <= 0)
    {
      _ivarPtr->isCandForBranching(false);
    }
    _ivarPtr->priority(vpriority);
  }
  return *this;
}

/**
 * todo: explain meaning of U and D
 * @param vdirective 'U' or 'D' branch treated in priority
 */
BcVar & BcVar::branchingDirective(char vdirective)
{
  if (_ivarPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
  } else
    _ivarPtr->directive(vdirective);
  return *this;
}

BcVarCoef BcVar::addCoef(double coef)
{
  return BcVarCoef(*this, coef);
}

BcRowExpression BcVar::addCoef(const BcVar& var)
{
  BcRowExpression exp;
  exp.addCoef(var);
  return exp;
}

BcRowExpression BcVar::addCoef(const BcVarCoef& varCoef)
{
  BcRowExpression exp;
  exp.addCoef(varCoef);
  return exp;
}

BcRowExpression BcVar::addCoef(BcVarIndex& varIndex)
{
  BcRowExpression exp;
  exp.addCoef(varIndex);
  return exp;
}

bool BcVar::operator<(const BcVar & that) const
{
  if (_ivarPtr == NULL)
    return false;
  if (that._ivarPtr == NULL)
    return true;

  if (_ivarPtr->genVarPtr()->defaultName() < that._ivarPtr->genVarPtr()->defaultName())
    return true;
  if (_ivarPtr->genVarPtr()->defaultName() > that._ivarPtr->genVarPtr()->defaultName())
    return false;

  if (_ivarPtr->id() < that._ivarPtr->id())
    return true;
  if (that._ivarPtr->id() < _ivarPtr->id())
    return false;

  return _ivarPtr->operator<(*(that._ivarPtr));
}


BcVarArray::BcVarArray(const BcFormulation & formulation,
                       const std::string & name ,
                       int firstIndexMax, int secondIndexMax, int thirdIndexMax) :
   _genericVarPtr(NULL), _curInstVarPtr(NULL)
{

  _genericVarPtr = formulation.probConfPtr()->getGenericVar(name);

  if (printL(5))
    std::cout << " BcVarArray::BcVarArray() : gvName =  " << name
              << " Exists ? =  " << (_genericVarPtr != NULL)
              << " probConf exists ? " << (formulation.probConfPtr() != NULL) << std::endl;

  if (_genericVarPtr == NULL)
  {
    _genericVarPtr =
            formulation.probConfPtr()->modelPtr()->createGenericVar(formulation.probConfPtr(),
                                                                    BcVarConstrType::local2Formulation,
                                                                    name,
                                                                    MultiIndexNames(),
                                                                    'I',
                                                                    0, /// cost
                                                                    BapcodInfinity, /// ub
                                                                    SelectionStrategy::MostFractional,
                                                                    1.0,
                                                                    0.1,
                                                                    's',
                                                                    'P',
                                                                    firstIndexMax, secondIndexMax, thirdIndexMax);
  }

}

BcVarArray::BcVarArray(GenericVar * genericVarPtr) :
    _genericVarPtr(genericVarPtr), _curInstVarPtr(NULL)
{
}

BcVarArray::~BcVarArray()
{
}

const std::string & BcVarArray::genericName() const
{
  if (_genericVarPtr == NULL)
  {
    std::cerr << "ERROR Model _genericVarPtr == NULL" << std::endl;
    exit(1);
  }
  return _genericVarPtr->defaultName();
}

BcFormulation BcVarArray::formulation() const
{
  if (_genericVarPtr == NULL)
  {
    if (printL(5))
        std::cout << "BaPCod info : Model _genericVarPtr == NULL" << std::endl;
    return BcFormulation(NULL);
  }
  return BcFormulation(_genericVarPtr->probConfPtr());
}

const BcModel& BcVarArray::model() const
{
  if (_genericVarPtr == NULL)
  {
    throw GlobalException("ModelVar::model(): Model _genericVarPtr == NULL");
  }
  return _genericVarPtr->bcModel();
}

BcModel& BcVarArray::model()
{
  if (_genericVarPtr == NULL)
  {
    throw GlobalException("ModelVar::model(): Model _genericVarPtr == NULL");
  }
  return _genericVarPtr->bcModel();
}

BcVar & BcVarArray::createElement(const MultiIndex & indexArray)
{
  if (_genericVarPtr->dimension() == -1)
  {
    _genericVarPtr->dimension(indexArray.endPosition);  
  }
  else 
  {    
    std::stringstream ss;
    ss <<  "BcVarArray::createElement Error : In a BcVarArray that has dimension " 
       << _genericVarPtr->dimension() << ", you can not have an element with " 
       << indexArray.endPosition << " indices";
    std::string s = ss.str();
    _genericVarPtr->bapcodInit().require(
                _genericVarPtr->dimension() == indexArray.endPosition , s.c_str());
  }
  
  ///already in cache memory
  if (_curInstVarPtr._ivarPtr != NULL)
  {
    if (_curInstVarPtr._ivarPtr->id() == indexArray)
      return _curInstVarPtr;
  }

  if (_genericVarPtr == NULL)
  {
    if (printL(5))
      std::cout << "BaPCod info : Model _genericVarPtr == NULL" << std::endl;
    _curInstVarPtr = BcVar(NULL);
    return _curInstVarPtr;
  }

  _curInstVarPtr._ivarPtr = _genericVarPtr->getVarPtr(indexArray);

  if (_curInstVarPtr._ivarPtr != NULL)
  {
    if (printL(5))
      std::cout << "BaPCod info : Model Var with index " << indexArray << " already exists " << std::endl;
    return _curInstVarPtr;
  }

  /// create var
  _curInstVarPtr._ivarPtr = _genericVarPtr->modelPtr()->createVariable(_genericVarPtr->probConfPtr(),
                                                                       _genericVarPtr, indexArray);

   // When the model has been already prepared, we need to add the new created variables
   // in the problem using the following block
   if (_genericVarPtr->probConfPtr()->isPrepared())
   {
       _genericVarPtr->probConfPtr()->probPtr()->addVar(_curInstVarPtr._ivarPtr, 1, 2);
   }

  return _curInstVarPtr;
}

bool BcVarArray::isDefinedAt(const MultiIndex & indexArray)
{

  if (_genericVarPtr == NULL)
    {
      if (printL(5))
        std::cout << "BaPCod info : BcVarArray::isDefinedAt : Model _genericVarPtr == NULL" << std::endl;
      return false;
    }

  if(_genericVarPtr->dimension() != indexArray.endPosition)
    {
      if (printL(5))
	    std::cout << "BaPCod info : BcVarArray::isDefinedAt : "
                  << "In BcVarArray there can not be an element with more indices than the dimension." << std::endl;
	  return false;
    }
    
  _curInstVarPtr._ivarPtr = _genericVarPtr->getVarPtr(indexArray);
  
  if (_curInstVarPtr._ivarPtr == NULL)
    {
      return false;
    }

  return true;
}

BcVar & BcVarArray::getElement(const MultiIndex& indexArray)
{
  if (_genericVarPtr->dimension() != indexArray.endPosition)
  {
          std::cerr << "BcVarArray::getElement(). "
                    << "In BcVarArray there can not be an element with more indices than the dimension. " << std::endl
                    << "      BcVarArray: " << _genericVarPtr->defaultName() << std::endl
                    << "      dimension : " << _genericVarPtr->dimension() << std::endl
                    << "      nbIndices : " << indexArray.endPosition << std::endl;
          exit(1);
  }
  ///already in cache memory
  if (_curInstVarPtr._ivarPtr != NULL)
  {
    if (_curInstVarPtr._ivarPtr->id() == indexArray)
      return _curInstVarPtr;
  }

  if (_genericVarPtr == NULL)
  {
    if (printL(5))
      std::cout << "BaPCod info : Model _genericVarPtr == NULL" << std::endl;
    _curInstVarPtr = BcVar(NULL);
    return _curInstVarPtr;
  }

  _curInstVarPtr._ivarPtr = _genericVarPtr->getVarPtr(indexArray);
  if (printL(5) && _curInstVarPtr._ivarPtr == NULL)
  {
    std::cout << "BaPCod info : Model Var " << _genericVarPtr->defaultName() << " has no index " << indexArray
              << std::endl;
  }

  return _curInstVarPtr;
}

BcVarIndex BcVarArray::operator[](const int & index)
{
  return BcVarIndex(*this, MultiIndex(index));
}

const BcVarArray & BcVarArray::ruleForBranchingVarSelection(const SelectionStrategy & branchingPrRule)
{
  if (_genericVarPtr != NULL)
    _genericVarPtr->priorityRule(branchingPrRule);
  return *this;
}

std::ostream& operator<<(std::ostream& os, const BcVarArray & that)
{
  return that.genericVarPtr()->print(os);
}

const BcVarArray& BcVarArray::defineIndexNames(const MultiIndexNames& multiIndexNames)
{
  if (_genericVarPtr != NULL)
    _genericVarPtr->multiIndexNames(multiIndexNames);

  return *this;
}

const BcVarArray & BcVarArray::priorityForMasterBranching(double priority)
{
  if (_genericVarPtr != NULL)
  {
    if (priority > _genericVarPtr->priorityLevel())
      _genericVarPtr->priorityLevel(priority);

    _genericVarPtr->genericBranchingOnAggregateVarPL(priority);
  }
  return *this;
}

const BcVarArray & BcVarArray::priorityForSubproblemBranching(double priority)
{
  if (_genericVarPtr != NULL)
  {
    if (priority > _genericVarPtr->priorityLevel())
      _genericVarPtr->priorityLevel(priority);

    _genericVarPtr->compBoundSetBranchingPL(priority);
  }
  return *this;
}

const BcVarArray & BcVarArray::priorityForRyanFosterBranching(double priority)
{
  if (_genericVarPtr != NULL)
  {
    if (priority > _genericVarPtr->priorityLevel())
      _genericVarPtr->priorityLevel(priority);

    _genericVarPtr->ryanFosterBranchingPL(priority);
  }
  return *this;
}


const BcVarArray & BcVarArray::branchingDirective(const char & vdefaultDirective)
/// 'U' or 'D' branch treated in priority
{
  if (_genericVarPtr != NULL)
    _genericVarPtr->defaultDirective(vdefaultDirective);
  return *this;
}

const BcVarArray & BcVarArray::globalUb(double vdefaultUb)
{
  if (_genericVarPtr != NULL)
  {
    _genericVarPtr->defaultGlobalUb(vdefaultUb);
  }
  return *this;
}

const BcVarArray & BcVarArray::globalLb(double vdefaultLB)
{
  if (_genericVarPtr != NULL)
  {
    _genericVarPtr->defaultGlobalLb(vdefaultLB);
    _genericVarPtr->defaultSense('F');
  }
  return *this;
}

const BcVarArray & BcVarArray::localUb(double vdefaultUb)
{
  if (_genericVarPtr != NULL)
  {
    _genericVarPtr->defaultUb(vdefaultUb);
  }
  return *this;
}

const BcVarArray & BcVarArray::localLb(double vdefaultLB)
{
  if (_genericVarPtr != NULL)
  {
    _genericVarPtr->defaultLb(vdefaultLB);
    _genericVarPtr->defaultSense('F');
  }
  return *this;
}

const BcVarArray& BcVarArray::defaultVal(double defaultValue)
{
  if (_genericVarPtr != NULL)
  {
    _genericVarPtr->defaultVal(defaultValue);
  }
  return *this;
}

const BcVarArray& BcVarArray::type(const char& defaultType)
{
  if ((_genericVarPtr != NULL) && ((defaultType == 'B') || (defaultType == 'I') || (defaultType == 'C')))
    _genericVarPtr->defaultType(defaultType);
  if (defaultType == 'B')
  {
    _genericVarPtr->defaultUb(1);
    _genericVarPtr->defaultSense('P');
  }
  if (defaultType == 'C')
    ruleForBranchingVarSelection(SelectionStrategy::NotConsideredForSelection);
  return *this;
}

const BcVarArray & BcVarArray::sense(const char & defaultSense)
{
    if ((_genericVarPtr != NULL) && ((defaultSense == 'N') || (defaultSense == 'F') || (defaultSense == 'P')))
        _genericVarPtr->defaultSense(defaultSense);
    return *this;
}

const BcVarArray & BcVarArray::setImplicit()
{
    if (_genericVarPtr != NULL)
        _genericVarPtr->defaultKind('I');
    return *this;
}

BcVarCoef::BcVarCoef(const BcVar & var, double coef) :
    std::pair<BcVar, double>(var, coef)
{
}

BcVarCoef::BcVarCoef(const BcVar & var) :
    std::pair<BcVar, double>(var, Double(1.0))
{
}

BcVarCoef::BcVarCoef(BcVarIndex & modVarInd, double coef) :
    std::pair<BcVar, double>((BcVar) modVarInd, coef)
{
}

BcVarCoef::BcVarCoef(BcVarIndex & modVarInd) :
    std::pair<BcVar, double>((BcVar) modVarInd, double(1.0))
{
}

BcRowExpression BcVarCoef::addCoef(const BcVar& var)
{
  BcRowExpression exp;
  exp.operator+(var);
  return exp;
}

BcRowExpression BcVarCoef::addCoef(const BcVarCoef& varCoef)
{
  BcRowExpression exp;
  exp.operator+(varCoef);
  return exp;
}

BcRowExpression BcVarCoef::addCoef(BcVarIndex& varIndex)
{
  BcRowExpression exp;
  exp.operator+(varIndex);
  return exp;
}

BcVarIndex::BcVarIndex(const BcVarArray & modvar, const MultiIndex & indexArray) :
    std::pair<BcVarArray, MultiIndex>(modvar, indexArray)
{
}

BcVarIndex::BcVarIndex(const BcVarArray & modvar) :
    std::pair<BcVarArray, MultiIndex>(modvar, MultiIndex(0))
{
}

bool BcVarIndex::isDefined()
{
  return first.getElement(second).isDefined();
}

BcVarIndex::operator BcVar()
{
  return first.getElement(second);
}

BcVarIndex & BcVarIndex::operator[](const int & index)
{
  second.operator+=(index);
  return *this;
}

void BcVarIndex::type(const char& flag)
{
  ((BcVar) (*this)).type(flag);
}

const MultiIndex & BcVarIndex::id()
{
  return first.getElement(second).id();
}

double BcVarIndex::curLb()
{
  return first.getElement(second).curLb();
}

double BcVarIndex::curUb()
{
  return first.getElement(second).curUb();
}

double BcVarIndex::curCost()
{
  if (first.isDefinedAt(second))
  {
    return first.curInstVarPtr().curCost();
  }
  else
  {
    std::cout << "Error at BcVarIndex ::curCost() " << std::endl;
    exit(0);
  }
}


double BcVarIndex::solVal()
{
    return first.getElement(second).solVal();
}

double BcVarIndex::curVal()
{
    return first.getElement(second).curVal();
}

bool BcVarIndex::inCurForm()
{
  return first.getElement(second).inCurForm();
}

//std::ostream& operator<<(std::ostream& os, BcVarIndex & that)
//{
//  return os << that.curVal();
//}

const BcVar & BcVarIndex::operator=(double value)
{
    return first.getElement(second).operator=(value);
}

BcVar & BcVarIndex::operator<=(double vub)
{
  return first.getElement(second).operator<=(vub);
}

BcVar & BcVarIndex::operator>=(double vlb)
{
  return first.getElement(second).operator>=(vlb);
}

BcVar & BcVarIndex::branchingPriority(double vpriority)
{
    return first.getElement(second).branchingPriority(vpriority);
}

BcVar & BcVarIndex::branchingDirective(char directive)
{
    return first.getElement(second).branchingDirective(directive);
}

BcVarCoef BcVarIndex::addCoef(double coef)
{
  return first.getElement(second).addCoef(coef);
}

BcVarCoef BcVarIndex::operator*(double coef)
{
  return addCoef(coef);
}

BcRowExpression BcVarIndex::addCoef(const BcVar& var)
{
  BcRowExpression exp;
  exp.addCoef(var);
  return exp;
}

BcRowExpression BcVarIndex::operator+(const BcVar & var)
{
  return addCoef(var);
}

BcRowExpression BcVarIndex::addCoef(const BcVarCoef& varCoef)
{
  BcRowExpression exp;
  exp.addCoef(varCoef);
  return exp;
}

BcRowExpression BcVarIndex::operator+(const BcVarCoef & varCoef)
{
  return addCoef(varCoef);
}

BcRowExpression BcVarIndex::addCoef(BcVarIndex& varIndex)
{
  BcRowExpression exp;
  exp.addCoef(varIndex);
  return exp;
}

BcRowExpression BcVarIndex::operator+(BcVarIndex & varIndex)
{
  return addCoef(varIndex);
}

BcRowExpression::BcRowExpression() :
    _mult(1)
{
}

const BcRowExpression& BcRowExpression::addCoef(const BcVarCoef& coef)
{
  this->push_back(coef);
  return *this;
}

const BcRowExpression& BcRowExpression::removeCoef(BcVarCoef& coef)
{
  coef.second *= -1;
  this->push_back(coef);
  return *this;
}

const BcRowExpression& BcRowExpression::addCoef(const BcVar& coef)
{
  this->push_back(BcVarCoef(coef, 1));
  return *this;
}

const BcRowExpression& BcRowExpression::removeCoef(const BcVar& coef)
{
  this->push_back(BcVarCoef(coef, -1));
  return *this;
}

const BcRowExpression& BcRowExpression::addCoef(BcVarIndex& coef)
{
  return addCoef((BcVar) coef);
}

const BcRowExpression& BcRowExpression::removeCoef(BcVarIndex& coef)
{
  return removeCoef((BcVar) coef);
}

BcRowExpression BcRowExpression::operator*(const double mult) const
{
  BcRowExpression that(*this);
  that._mult = mult;
  return that;
}

BcRowExpression operator*(const double mult, const BcRowExpression & exp)
{
  return (exp * mult);
}

