/**
 *
 * This file bcModelConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcModelConstrC.hpp"

BcConstr::BcConstr(InstanciatedConstr * iconstrPtr) :
    _iconstrPtr(iconstrPtr)
{
}

bool BcConstr::isDefined() const
{
  return (_iconstrPtr != NULL);
}

const std::string & BcConstr::genericName() const
{
  if (_iconstrPtr == NULL)
  {
    std::cerr << "ERROR Model BcConstr == NULL" << std::endl;
    exit(1);
  }
  return _iconstrPtr->genConstrPtr()->defaultName();
}

void BcConstr::nicePrint() const
{
    if (_iconstrPtr == NULL)
    {
        std::cerr << "ERROR Model BcConstr == NULL" << std::endl;
        exit(1);
    }
    return _iconstrPtr->nicePrint(std::cout);
}

const std::string & BcConstr::name() const
{
  if (_iconstrPtr == NULL)
  {
    std::cerr << "ERROR Model BcConstr == NULL" << std::endl;
    exit(1);
  }
  return _iconstrPtr->name();
}

const MultiIndex & BcConstr::id() const
{
  if (_iconstrPtr == NULL)
  {
    std::cerr << "ERROR Model BcConstr == NULL" << std::endl;
    exit(1);
  }
  return _iconstrPtr->id();
}

double BcConstr::curRhs() const
{
  if (_iconstrPtr == NULL)
  {
    std::cerr << "ERROR Model constrPtr == NULL" << std::endl;
    exit(1);
  }
  return _iconstrPtr->curRhs();
}

double BcConstr::curDualVal() const
{
  if (_iconstrPtr == NULL)
  {
    std::cerr << "ERROR Model constrPtr == NULL" << std::endl;
    exit(1);
  }
  return _iconstrPtr->valOrSepPointVal();
}

const BcFormulation BcConstr::formulation() const
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
    return BcFormulation(NULL);
  }
  return BcFormulation(_iconstrPtr->probConfPtr());
}

bool BcConstr::inCurProb() const
{
  if (_iconstrPtr == NULL)
  {
    std::cout << "ERROR Model BcConstr == NULL" << std::endl;
    exit(0);
  }
  return _iconstrPtr->inCurProb();
}

const BcConstr & BcConstr::add(const BcVarCoef & coef)
{

  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    if (coef.first._ivarPtr == NULL)
    {
      if (printL(6))
        std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    }
    else
    {
      _iconstrPtr->genConstrPtr()->modelPtr()->addCoefficient(_iconstrPtr, coef.first._ivarPtr, coef.second);
    }
  }
  return *this;
}
const BcConstr & BcConstr::remove(const BcVarCoef & coef)
{

  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    if (coef.first._ivarPtr == NULL)
    {
      if (printL(6))
        std::cout << "BaPCod info : Model BcVar == NULL" << std::endl;
    }
    else
    {
      _iconstrPtr->genConstrPtr()->modelPtr()->addCoefficient(_iconstrPtr, coef.first._ivarPtr, - coef.second);
    }
  }
  return *this;
}

const BcConstr & BcConstr::add(const BcRowExpression & exp)
{
  for (std::list<BcVarCoef>::const_iterator it = exp.begin(); it != exp.end(); ++it)
  {
    BcConstr::add((exp._mult * it->second) * it->first);
  }
  return *this;
}

const BcConstr & BcConstr::remove(const BcRowExpression & exp)
{
  for (std::list<BcVarCoef>::const_iterator it = exp.begin(); it != exp.end(); ++it)
  {
    BcConstr::remove(*it);
  }
  return *this;
}

bool BcConstr::operator<(const BcConstr & that) const
{
  if (_iconstrPtr == NULL)
    return false;
  if (that._iconstrPtr == NULL)
    return true;

  if (_iconstrPtr->genConstrPtr()->defaultName() < that._iconstrPtr->genConstrPtr()->defaultName())
    return true;
  if (_iconstrPtr->genConstrPtr()->defaultName() > that._iconstrPtr->genConstrPtr()->defaultName())
    return false;

  if (_iconstrPtr->id() < that._iconstrPtr->id())
    return true;
  if (that._iconstrPtr->id() < _iconstrPtr->id())
    return false;

  return _iconstrPtr->operator<(*(that._iconstrPtr));
}

void BcConstr::dualVal(double dualVal)
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    _iconstrPtr->incumbentVal(dualVal);
  }
}

void BcConstr::sense(const char& sense)
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    _iconstrPtr->sense(sense);
  }
}

void BcConstr::rhs(double rhs)
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    _iconstrPtr->costrhs(rhs);
  }
}


void BcConstr::type(const char& type)
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  else
  {
    if ((type == 'E') || (type == 'I'))
    {
      _iconstrPtr->kind(type);
    }
    else
    {
      _iconstrPtr->type(type);
    }
  }
}


const BcConstr & BcConstr::operator-(BcVarIndex & modVarInd)
{
  return this->operator-((BcVar) modVarInd); /// casting
}

const BcConstr & BcConstr::operator+(BcVarIndex & modVarInd)
{
  return this->operator+((BcVar) modVarInd); /// casting
}

const BcConstr & BcConstr::operator+=(BcVarIndex & modVarInd)
{
  return this->operator+=((BcVar) modVarInd); /// casting
}

const BcConstr & BcConstr::operator-=(BcVarIndex & modVarInd)
{
  return this->operator-=((BcVar) modVarInd); /// casting
}

const BcModel& BcConstrArray::model() const
{
  if (_genericConstrPtr == NULL)
  {
    throw GlobalException("ModelConstr::model(): Model _genericConstrPtr == NULL");
  }
  return _genericConstrPtr->bcModel();
}

BcModel& BcConstrArray::model()
{
  if (_genericConstrPtr == NULL)
  {
    throw GlobalException("ModelConstr::model(): Model _genericConstrPtr == NULL");
  }
  return _genericConstrPtr->bcModel();
}

BcConstr& BcConstrArray::createElement(const MultiIndex & indexArray)
{
    if(_genericConstrPtr->dimension() == -1)
    {
        _genericConstrPtr->dimension(indexArray.endPosition);
    }
    else
    {
        std::stringstream ss;
        ss <<  "BcConstrArray::createElement Error : In a BcConstrArray that has dimension "
           << _genericConstrPtr->dimension() << ", you can not have an element with "
           << indexArray.endPosition << " indices";
        std::string s = ss.str();
        _genericConstrPtr->bapcodInit().require(_genericConstrPtr->dimension() == indexArray.endPosition,
                                                s.c_str());
    }

    // already in cache memory
    if (_curInstConstrPtr._iconstrPtr != NULL)
    {
        if (_curInstConstrPtr._iconstrPtr->id() == indexArray)
        {
            return _curInstConstrPtr;
        }
    }

    if (_genericConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << "BaPCod info : Model _genericConstrPtr == NULL" << std::endl;
        _curInstConstrPtr = BcConstr(NULL);
        return _curInstConstrPtr;
    }

    _curInstConstrPtr._iconstrPtr = _genericConstrPtr->getConstrPtr(indexArray);

    if (_curInstConstrPtr._iconstrPtr != NULL)
    {
        if (printL(5))
        {
            std::cout << "BaPCod info : Model Constr with index " << indexArray << " already exists" << std::endl;
            return _curInstConstrPtr;
        }
    }

    /// create constr
    _curInstConstrPtr._iconstrPtr = _genericConstrPtr->modelPtr()->createConstraint(_genericConstrPtr->probConfPtr(),
                                                                                    _genericConstrPtr, indexArray);

    /// create membership
    if (_genericConstrPtr->addConstrFunctorPtr() != NULL)
    {
        (*(_genericConstrPtr->addConstrFunctorPtr()))(indexArray);
    }

    return _curInstConstrPtr;
}

BcConstr::operator InstanciatedConstr *() const
{
  if (_iconstrPtr == NULL)
  {
    if (printL(6))
      std::cout << "BaPCod info : Model BcConstr == NULL" << std::endl;
  }
  return _iconstrPtr;
}

const BcConstrArray & BcConstrArray::defineIndexNames(const MultiIndexNames & newMultiIndexNames)
{
  if (_genericConstrPtr != NULL)
    _genericConstrPtr->multiIndexNames(newMultiIndexNames);
  return *this;
}

BcConstrArray::BcConstrArray() :
    _genericConstrPtr(NULL), _curInstConstrPtr(NULL)
{
}

BcConstrArray::BcConstrArray(const BcFormulation & formulation, const std::string & name) :
    _genericConstrPtr(NULL), _curInstConstrPtr(NULL)
{

  if (printL(6))
    std::cout << " BcConstrArray() : ProbConfig =  " << formulation.name()
              << " BcConstrArray =  " << name << std::endl;
  
  _genericConstrPtr = formulation._probConfPtr->getGenericConstr(name);

  if (_genericConstrPtr == NULL)
  {
    _genericConstrPtr
        = formulation._probConfPtr->modelPtr()->createGenericConstr(formulation._probConfPtr,
			                                                        BcVarConstrType::local2Formulation,
			                                                        name,
                                                    			    MultiIndexNames(),
                                                                    'E',
                                                    			    1,
                                                    			    0,
                                                    			    true,
                                                                    's');
  }
}

const BcFormulation BcConstrArray::formulation() const
{
  if (_genericConstrPtr == NULL)
  {
    std::cout << "ERROR Model _genericConstrPtr == NULL" << std::endl;
    return BcFormulation(NULL);
  }
  return BcFormulation(_genericConstrPtr->probConfPtr());
}


BcConstrArray::~BcConstrArray()
{
}

const std::string & BcConstrArray::genericName() const
{
  if (_genericConstrPtr == NULL)
  {
    std::cout << "ERROR Model _genericConstrPtr == NULL" << std::endl;
    exit(0);
  }
  return _genericConstrPtr->defaultName();
}

void BcConstrArray::toBeUsedInPreprocessing(bool flag)
{
  if (_genericConstrPtr == NULL)
  {
    std::cout << "ERROR Model _genericConstrPtr == NULL" << std::endl;
    exit(0);
  }

  _genericConstrPtr->toBeUsedInPreprocessing(flag);
}

void BcConstrArray::considerAsEqualityInPreprocessing(bool flag)
{
  if (_genericConstrPtr == NULL)
  {
    std::cout << "ERROR Model _genericConstrPtr == NULL" << std::endl;
    exit(0);
  }

  _genericConstrPtr->considerAsEqualityInPreprocessing(flag);
}

BcConstrIndex BcConstrArray::operator[](const int & index)
{
  return BcConstrIndex(*this, MultiIndex(index));
}

bool BcConstrArray::isDefinedAt(const MultiIndex & indexArray) 
{
  
  if (_genericConstrPtr == NULL)
  {
    if (printL(5))
      std::cout << "BaPCod info : Model _genericConstrPtr == NULL" << std::endl;
    return false;
  }

  if(_genericConstrPtr->dimension() != indexArray.endPosition)
  {
     if (printL(5))
     {
       std::cout << "BaPCod info : : In BcConstrArray there can not be an element"
                 << " with more indices than the dimension. " << std::endl;

       std::cout << "      BcConstrArray : " << _genericConstrPtr->defaultName() <<  std::endl;
       std::cout << "          Dimension : " << _genericConstrPtr->dimension() << std::endl;
       std::cout << "  Number of indices : " << indexArray.endPosition << std::endl;   
     }
     
     return false;
  }
  
  _curInstConstrPtr._iconstrPtr = _genericConstrPtr->getConstrPtr(indexArray);
  
  if (_curInstConstrPtr._iconstrPtr == NULL)
  {
    return false;
  }
  
  return true;
  
}

BcConstr& BcConstrArray::getElement(const MultiIndex& indexArray) 
{
  if (_genericConstrPtr->dimension() != indexArray.endPosition)
  {
     std::cerr << "Error : In BcConstrArray there can not be an element"
               << " with more indices than the dimension. " << std::endl;
     
     std::cerr << "      BcConstrArray : " << _genericConstrPtr->defaultName() <<  std::endl;
     std::cerr << "          Dimension : " << _genericConstrPtr->dimension() << std::endl;
     std::cerr << "  Number of indices : " << indexArray.endPosition << std::endl;
     
     exit(1);
  }

  if (_genericConstrPtr == NULL)
  {
    if (printL(5))
        std::cout << "BaPCod info : Model _genericConstrPtr == NULL" << std::endl;
    _curInstConstrPtr._iconstrPtr = NULL;
    return _curInstConstrPtr;
  }

  _curInstConstrPtr._iconstrPtr = _genericConstrPtr->getConstrPtr(indexArray);
  if (printL(5) && (_curInstConstrPtr._iconstrPtr == NULL))
  {
    std::cout << "BaPCod info : Model Constr " << _genericConstrPtr->defaultName()
              << " has no index " << indexArray << std::endl;
  }

  return _curInstConstrPtr;
}


void BcConstrArray::sense(char sense)
{
  _genericConstrPtr->defaultSense(sense);
}

void BcConstrArray::rhs(double rhs)
{
  _genericConstrPtr->defaultCostRhs(rhs);
  if (printL(5))
  {
    std::cout << "BcConstrArray name " << _genericConstrPtr->defaultName()
              << " rhs " << _genericConstrPtr->defaultCostRhs() << std::endl;
  }
}

void BcConstrArray::type(char type)
{
  _genericConstrPtr->defaultType(type);
}

void BcConstrArray::flag(char flag)
{
  _genericConstrPtr->defaultFlag(flag);
}

void BcConstrArray::dualVal(double defaultDualVal)
{
  _genericConstrPtr->defaultVal(defaultDualVal);
}

std::ostream & BcConstrArray::print(std::ostream& os) const
{
    return _genericConstrPtr->print(os);
}

const BcConstrArray & BcConstrArray::attach(BcAddConstrFunctor * addConstrRoutinePtr)
{
    if (_genericConstrPtr == NULL)
    {
        std::cout << "ERROR Model _genericConstrPtr == NULL" << std::endl;
        exit(0);
    }

    _genericConstrPtr->addConstrFunctorPtr(addConstrRoutinePtr);

    return *this;
}

std::ostream& operator<<(std::ostream& os, const BcConstrArray & that)
{
  return that.print(os);
}

std::ostream& operator<<(std::ostream& os, const BcConstr & that)
{
  return ((InstanciatedConstr *)that)->print(os);
}

BcConstrIndex::BcConstrIndex(const BcConstrArray & modConstr,
                             const MultiIndex & indexArray) :
    std::pair<BcConstrArray, MultiIndex>(modConstr, indexArray)
{
}

BcConstrIndex & BcConstrIndex::operator[](const int & index)
{
  second.operator+=(index);
  return *this;
}

bool BcConstrIndex::isDefined()
{
  return first.getElement(second).isDefined();
}

//const BcConstr & BcConstrIndex::operator==(const char & type)
//{
//  /// type = 'C' for core (required for the IP formulation,
//  /// type = 'F' for facultative (only helpfull for the LP formulation),
//  /// type = 'S' for constraints defining a subsystem in column generation for extended formulation approach
//  /// type = 'X' for constraints defining a subproblem convexity constraint in the master
//  /// type = 'I' for constraints defining a implicit constraint used for preprocessing but not place in the formulation
//  /// type = 'E' for constraints defining a explicit constraint inserted in the formulation
//
//  return ((BcConstr) (*this)).operator==(type);
//}

const BcConstr & BcConstrIndex::operator+(const BcVarCoef & coef)
{
  return first.getElement(second).operator+(coef);
}

const BcConstr & BcConstrIndex::operator-(const BcVarCoef & coef)
{
  return first.getElement(second).operator-(coef);
}

const BcConstr & BcConstrIndex::operator+=(const BcVarCoef & coef)
{
  return first.getElement(second).operator+=(coef);
}

const BcConstr & BcConstrIndex::operator-=(const BcVarCoef & coef)
{
  return first.getElement(second).operator-=(coef);
}

const BcConstr & BcConstrIndex::operator+(const BcVar & var)
{
  return first.getElement(second).operator+(var);
}

const BcConstr & BcConstrIndex::operator-(const BcVar & var)
{
  return first.getElement(second).operator-(var);
}

const BcConstr & BcConstrIndex::operator+=(const BcVar & var)
{
  return first.getElement(second).operator+=(var);
}

const BcConstr & BcConstrIndex::operator-=(const BcVar & var)
{
  return first.getElement(second).operator-=(var);
}

const BcConstr & BcConstrIndex::operator+(BcVarIndex & modVarInd)
{
  return first.getElement(second).operator+(modVarInd);
}

const BcConstr & BcConstrIndex::operator-(BcVarIndex & modVarInd)
{
  return first.getElement(second).operator-(modVarInd);
}

const BcConstr & BcConstrIndex::operator+=(BcVarIndex & modVarInd)
{
  return first.getElement(second).operator+=(modVarInd);
}

const BcConstr & BcConstrIndex::operator-=(BcVarIndex & modVarInd)
{
  return first.getElement(second).operator-=(modVarInd);
}

const BcConstr & BcConstrIndex::operator=(const BcRowExpression & exp)
{
  return BcConstrIndex::operator+=(exp);
}

const BcConstr & BcConstrIndex::operator+(const BcRowExpression & exp)
{
  return BcConstrIndex::operator+=(exp);
}

const BcConstr & BcConstrIndex::operator+=(const BcRowExpression & exp)
{
  return first.getElement(second).operator+=(exp);
}

const BcConstr & BcConstrIndex::operator-(const BcRowExpression & exp)
{
  return BcConstrIndex::operator-=(exp);
}

const BcConstr & BcConstrIndex::operator-=(const BcRowExpression & exp)
{
  return first.getElement(second).operator-=(exp);
}

const MultiIndex & BcConstrIndex::id()
{
  return first.getElement(second).id();
}

bool BcConstrIndex::inCurProb()
{
  return first.getElement(second).inCurProb();
}
void BcConstrIndex::sense(const char & sense)
{
  ((BcConstr) (*this)).sense(sense);
}
void BcConstrIndex::rhs(double rhs)
{
  ((BcConstr) (*this)).rhs(rhs);
}

void BcConstrIndex::type(const char& type)
{
  ((BcConstr) (*this)).type(type);
}

void BcConstrIndex::dualVal(double dualVal)
{
  ((BcConstr) (*this)).dualVal(dualVal);
}

BcConstrIndex::operator BcConstr()
{
  return first.getElement(second);
}

void BcAddConstrFunctor::operator()(const MultiIndex & indexArray)
{
    if (printL(0))
        std::cout << " BcAddConstrFunctor::operator()  SHOULD NOT BE CALLED" << std::endl;
    return;
}

//const BcModel& BcAddConstrFunctor::model() const
//{
//    return _model;
//}
//
//BcModel& BcAddConstrFunctor::model()
//{
//    return _model;
//}

