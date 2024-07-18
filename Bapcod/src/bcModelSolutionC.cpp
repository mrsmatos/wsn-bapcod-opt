/**
 *
 * This file bcModelSolutionC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcModelingLanguageC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcBapcodInit.hpp"
#include "bcModelRCSPSolver.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif 
/**
 * Generic code
 */
using namespace std;


BcSolution::BcSolution(Solution * SolutionPointer) :
  _solutionPtr (SolutionPointer)
{}

BcSolution BcSolution::clone()
{
    return BcSolution(this->solutionPtr()->clone());
}

 BcSolution::BcSolution(const BcFormulation & formulation) :
   _solutionPtr(new Solution(formulation.probConfPtr()))
 {}

BcSolution::~BcSolution() 
{
}

BcSolution & BcSolution::includeVarVal(const BcVar & var)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::includeVarVal: undefined solution";
      exit(1);
    }
  _solutionPtr->includeVar(var._ivarPtr, var._ivarPtr->val(), false); 
  return *this;
}

BcSolution & BcSolution::updateVarVal(const BcVar & var)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::updateVarVal: undefined solution";
      exit(1);
    }
  _solutionPtr->includeVar(var._ivarPtr, var._ivarPtr->val(), true); 
  return *this;
}

BcSolution & BcSolution::updateVarVal(const BcVar & var, const double & val)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::updateVarVal: undefined solution";
      exit(1);
    }
  _solutionPtr->includeVar(var._ivarPtr, val, true);
  return *this;
}


BcSolution & BcSolution::appendSol(BcSolution & sol)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::appendSol: undefined solution";
      exit(1);
    }

  if (sol._solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::appendSol: undefined new solution";
      exit(1);
    }

  sol._solutionPtr->previousSolPtr(_solutionPtr);

  return *this;
}

#ifdef BCP_RCSP_IS_FOUND
void BcSolution::setRCSPsolution(bcp_rcsp::Solution * rcspSolPtr)
{
    if (_solutionPtr == NULL)
    {
        std::cerr << "initializeOrderedSolution: undefined solution";
        exit(1);
    }
    Problem * probPtr = _solutionPtr->probConfPtr()->probPtr();
    const BcRCSPFunctor * functorPtr = dynamic_cast<const BcRCSPFunctor *>(probPtr->solverOracleFunctorPtr());
    if (functorPtr != nullptr)
        functorPtr->addToSolution(rcspSolPtr, *this);
}

const bcp_rcsp::Solution * BcSolution::getRCSPsolution()
{
    if (_solutionPtr == NULL)
        return NULL;
    return _solutionPtr->rcspSolPtr();
}
#endif

void BcSolution::initializeOrderedSolution(const std::vector<double> & resConsumption)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "initializeOrderedSolution: undefined solution";
      exit(1);
    }
  _solutionPtr->addToResConsumption(resConsumption);
}

void BcSolution::addToOrderedSolution(const int & id, const std::vector<double> & resConsumption,
                                      const bool & doNotAddSameConsecutiveIds)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "BcSolution::addToOrderedSolution: undefined solution";
      exit(1);
    }
  if (!doNotAddSameConsecutiveIds || _solutionPtr->orderedIds().empty() || (_solutionPtr->orderedIds().back() != id))
    {
      _solutionPtr->addToOrderedIds(id);
      _solutionPtr->addToResConsumption(resConsumption);
    }
}

const std::vector<int> & BcSolution::orderedIds() const
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "BcSolution::orderedIds: undefined solution";
      exit(1);
    }
#ifdef BCP_RCSP_IS_FOUND
  if (_solutionPtr->rcspSolPtr() != nullptr)
    return _solutionPtr->rcspSolPtr()->arcIds;
#endif //BCP_RCSP_IS_FOUND
  return _solutionPtr->orderedIds();
}

#ifdef BCP_RCSP_IS_FOUND
const std::vector<double> & BcSolution::waitingTimes() const
{
  if (_solutionPtr == NULL || _solutionPtr->rcspSolPtr() == nullptr)
    {
      std::cerr << "BcSolution::waitingTimes: undefined solution";
      exit(1);
    }
  return _solutionPtr->rcspSolPtr()->waitingTimes;
}
#endif

const std::vector<std::vector<double> > & BcSolution::resConsumption() const
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "BcSolution::resConsumption: undefined solution";
      exit(1);
    }
#ifdef BCP_RCSP_IS_FOUND
  if (_solutionPtr->rcspSolPtr() != nullptr)
    return _solutionPtr->rcspSolPtr()->resConsumption;
#endif //BCP_RCSP_IS_FOUND
  return _solutionPtr->resConsumption();
}

void BcSolution::setEnumeratedFlag(const bool flag)
{
  if (_solutionPtr == NULL)
  {
    std::cerr << "BcSolution::setEnumeratedFlag: undefined solution";
    exit(1);
  }
  _solutionPtr->enumeratedFlag(flag);
}

std::ostream& BcSolution::print(std::ostream& os) const
{
  if (_solutionPtr == NULL)
    return os << "undefined solution";
  return _solutionPtr->print(os);
}

void BcSolution::printOrderedSolution(std::ostream& os) const
{
  if (_solutionPtr == NULL)
    return;
  _solutionPtr->printOrderedSolution(os);
}

bool BcSolution::defined() const {
    return (_solutionPtr != NULL);
}

double BcSolution::cost() const
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::cost: undefined solution";
      exit(1);
    }

  return _solutionPtr->cost();
}

double BcSolution::resetCost()
{
    if (_solutionPtr == NULL)
    {
        std::cerr << "SolutionPtr::cost: undefined solution";
        exit(1);
    }

    return _solutionPtr->resetCost();
}

void BcSolution::extractVar(const std::string & genericName, std::set< BcVar > & varSet)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVar: undefined solution";
      exit(1);
    }
  _solutionPtr->extractVarWithGenericName(genericName, varSet);
}

std::set< BcVar > BcSolution::extractVar(const std::string & genericName)
{
    if (_solutionPtr == NULL)
    {
        std::cerr << "SolutionPtr::extractVar: undefined solution";
        exit(1);
    }
    std::set<BcVar> varSet;
    if (genericName == "")
        extractVar(varSet);
    else
        extractVar(genericName, varSet);
    return varSet;
}

void BcSolution::extractVar(const std::string & genericName, const BcFormulation & formulation,
                            std::set< BcVar > & varSet)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVar: undefined solution";
      exit(1);
    }
  if (!formulation.isDefined())
    {
      std::cerr << "SolutionPtr::extractVar: undefined formulation";
      exit(1);
    }
  _solutionPtr->extractVarWithGenericName(genericName, formulation.probConfPtr(), varSet);
}

void BcSolution::extractVar(const std::string & genericName, int firstIndex, std::set< BcVar > & varSet)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVar: undefined solution";
      exit(1);
    }
  _solutionPtr->extractVarWithGenericName(genericName, firstIndex, varSet);
}

double BcSolution::getVarVal(BcVar vptr)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVal: undefined solution";
      exit(1);
    }
  return _solutionPtr->solVal((InstanciatedVar *) vptr);
}

void BcSolution::clear()
{
  if (_solutionPtr == NULL)
    {
      std::cerr << "BcSolution::clear: undefined solution";
      exit(1);
    }

  return _solutionPtr->clear();

}

void BcSolution::deleteSolution()
{
    if (_solutionPtr != NULL)
    {
        delete _solutionPtr;
        _solutionPtr = NULL;
    }
}

void BcSolution::deleteSolutionsChain()
{
  Solution * solPtr = _solutionPtr;
  if (solPtr != NULL)
    {
      Solution * thisSolPtr = solPtr;
      solPtr = solPtr->nextSolPtr();
      delete thisSolPtr;
    }
  _solutionPtr = NULL;
}

void BcSolution::getVar(std::set<BcVar> & varSet)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::getVar: undefined solution";
      exit(1);
    }

  _solutionPtr->getVar(varSet);
}

void BcSolution::extractVar(std::set<BcVar> & varSet)
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVar: undefined solution";
      exit(1);
    }

  _solutionPtr->extractVar(varSet);
}

BcFormulation BcSolution::formulation() const 
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVal: undefined solution";
      exit(1);
    }
  return BcFormulation(_solutionPtr->probConfPtr());
}

void BcSolution::formulation(const BcFormulation & formulation)
{
  if (_solutionPtr == NULL)
  {
    std::cerr << "SolutionPtr::extractVal: undefined solution";
    exit(1);
  }
  _solutionPtr->probConfigPtr(formulation.probConfPtr());
}


int BcSolution::getMultiplicity() const
{
  if (_solutionPtr == NULL) 
    {
      std::cerr << "SolutionPtr::extractVal: undefined solution";
      exit(1);
    }
  return _solutionPtr->multiplicity();
}

BcSolution BcSolution::next() const
{
  if (_solutionPtr == NULL) 
    {
      return BcSolution(NULL);
    }
  if (printL(5)) 
    {
      std::cout << "SolutionPtr::next() of solution = " << *_solutionPtr;
      if (_solutionPtr->nextSolPtr() != NULL) 
	std::cout << "SolutionPtr::next() is solution = " << *(_solutionPtr->nextSolPtr());
    }

  return BcSolution(_solutionPtr->nextSolPtr());
}


BcSolution::operator Solution *() const
{
  if (_solutionPtr == NULL) 
    {
      return NULL;
    }
  return _solutionPtr;
}

const Solution * BcSolution::solutionPtr() const
{

  return _solutionPtr;
}

Solution * BcSolution::solutionPtr()
{
  return _solutionPtr;
}


BcDualSolution::BcDualSolution(DualSolution * SolutionPointer) :
        _dualSolutionPtr (SolutionPointer)
{}


BcDualSolution::BcDualSolution(const BcFormulation & formulation) :
        _dualSolutionPtr (new DualSolution(formulation.probConfPtr()))
{}

BcDualSolution::operator DualSolution *() const
{
    return _dualSolutionPtr;
}

bool BcDualSolution::isDefined() const
{
    return (_dualSolutionPtr != NULL);
}

BcDualSolution & BcDualSolution::includeConstr(const BcConstr & constr)
{
    if (_dualSolutionPtr == NULL)
    {
        std::cerr << "DualSolutionPtr::operator +=: undefined solution";
        exit(1);
    }
    _dualSolutionPtr->includeConstr(constr._iconstrPtr, constr._iconstrPtr->val(), false);
    return *this;
}

BcDualSolution & BcDualSolution::updateConstrVal(const BcConstr & constr)
{
    if (_dualSolutionPtr == NULL)
    {
        std::cerr << "DualSolutionPtr::operator +=: undefined solution";
        exit(1);
    }
    _dualSolutionPtr->includeConstr(constr._iconstrPtr, constr._iconstrPtr->val(), true);
    return *this;
}

BcDualSolution & BcDualSolution::appendSol(BcDualSolution & sol)
{
    if (_dualSolutionPtr == NULL)
    {
        std::cerr << "DualSolutionPtr::appendSol undefined solution";
        exit(1);
    }

    if (sol._dualSolutionPtr == NULL)
    {
        std::cerr << "DualSolutionPtr::appendSol undefined new solution";
        exit(1);
    }

    sol._dualSolutionPtr->previousSolPtr(_dualSolutionPtr);

    return *this;
}


std::ostream& BcDualSolution::print(std::ostream& os) const
{
    if (_dualSolutionPtr == NULL) return os << "undefined solution";
    return _dualSolutionPtr->print(os);
}


