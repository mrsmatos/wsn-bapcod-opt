/**
 *
 * This file bcSolutionC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcDoubleC.hpp"
#include "bcModelC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcApplicationException.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcModelVarC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif /* BCP_RCSP_IS_FOUND */

using namespace std;

Solution::Solution(ProbConfig * probConfPtr, Solution * newpreviousSolPtr) :
        _probConfPtr(probConfPtr),
        _ref(0),
        _cost(0),
        _multiplicity(1),
        _previousSolPtr(NULL),
        _nextSolPtr(NULL),
#ifdef BCP_RCSP_IS_FOUND
        _rcspSolPtr(NULL),
#endif //BCP_RCSP_IS_FOUND
        _enumeratedFlag(false)
{
  if (_probConfPtr != NULL)
  {
      _ref = _probConfPtr->pcSolCount();
      _probConfPtr->increasePCSolCount();
  }

  if (newpreviousSolPtr != NULL)
  {
    previousSolPtr(newpreviousSolPtr); /// the function inserts in place and set the next pointer in the previous solution;
  }
}

Solution::Solution(const Solution & that) :
        _probConfPtr(that._probConfPtr),
        _ref(0),
        _cost(that.cost()),
        _multiplicity(that._multiplicity),
        _previousSolPtr(NULL),
        _nextSolPtr(NULL),
#ifdef BCP_RCSP_IS_FOUND
        _rcspSolPtr(NULL),
#endif //BCP_RCSP_IS_FOUND
        _solVarValMap(that._solVarValMap),
        _orderedIds(that._orderedIds),
        _resConsumption(that._resConsumption),
        _enumeratedFlag(that._enumeratedFlag)
{
    if (_probConfPtr != NULL) {
        _ref = _probConfPtr->pcSolCount();
        _probConfPtr->increasePCSolCount();
    }

#ifdef BCP_RCSP_IS_FOUND
    if (that._rcspSolPtr != NULL)
    {
        _rcspSolPtr = new bcp_rcsp::Solution(that._rcspSolPtr);
    }
#endif /* BCP_RCSP_IS_FOUND */

  for (VarPtr2DoubleMap::const_iterator varit = _solVarValMap.begin(); varit != _solVarValMap.end(); ++varit)
  {
    if ((*varit).first->isTypeOf(VcId::MastColumnMask))
    {
      (*varit).first->incrParticipation(4);
    }
  }

  return;
}

Solution::Solution(ProbConfig * probConfPtr, const VarPtr2DoubleMap & solVarValMap):
        _probConfPtr(probConfPtr),
        _ref(0),
        _cost(0),
        _multiplicity(1),
        _previousSolPtr(NULL),
        _nextSolPtr(NULL),
#ifdef BCP_RCSP_IS_FOUND
        _rcspSolPtr(NULL),
#endif //BCP_RCSP_IS_FOUND
        _solVarValMap(solVarValMap),
        _enumeratedFlag(false)
{
    if (_probConfPtr != NULL)
    {
        _ref = _probConfPtr->pcSolCount();
        _probConfPtr->increasePCSolCount();
    }

  for (VarPtr2DoubleMap::const_iterator varIt = _solVarValMap.begin(); varIt != _solVarValMap.end(); ++varIt)
  {
    if ((*varIt).first->isTypeOf(VcId::MastColumnMask))
    {
      (*varIt).first->incrParticipation(4);
    }
  }

  resetCost();
}

Solution::~Solution()
{
  for (VarPtr2DoubleMap::const_iterator varit = _solVarValMap.begin(); varit != _solVarValMap.end(); ++varit)
  {
    if ((*varit).first->isTypeOf(VcId::MastColumnMask))
    {
       (*varit).first->decrParticipation(4);
    }
  }
#ifdef BCP_RCSP_IS_FOUND
  delete _rcspSolPtr;
#endif //BCP_RCSP_IS_FOUND
}

#ifdef BCP_RCSP_IS_FOUND
void Solution::setRcspSolPtr(const bcp_rcsp::Solution * rcspSolPtr)
{
    _rcspSolPtr = rcspSolPtr;
    if (_rcspSolPtr != nullptr)
        _enumeratedFlag = _rcspSolPtr->enumeratedFlag;
}

const bcp_rcsp::Solution * Solution::copyRcspSolPtr() const
{
    if (_rcspSolPtr == nullptr)
        return nullptr;

    return new bcp_rcsp::Solution(* _rcspSolPtr);
}
#endif /* BCP_RCSP_IS_FOUND */

void Solution::enumeratedFlag(bool flag)
{
    _enumeratedFlag = flag;

#ifdef BCP_RCSP_IS_FOUND
    if ((_rcspSolPtr != nullptr) && (_rcspSolPtr->enumeratedFlag != _enumeratedFlag))
    {
        bcp_rcsp::Solution * rcspEnumSolPtr = new bcp_rcsp::Solution(* _rcspSolPtr);
        rcspEnumSolPtr->enumeratedFlag = _enumeratedFlag;
        delete _rcspSolPtr;
        _rcspSolPtr = rcspEnumSolPtr;
    }
#endif /* BCP_RCSP_IS_FOUND */
}


void Solution::previousSolPtr(Solution * prevSolPtr)
{
  if (prevSolPtr != NULL)
  {
    /// attach next to prevSolPtr at the end of this solution chain
    if (prevSolPtr->_nextSolPtr != NULL)
    {
      Solution * tmpEnd = this;
      while (tmpEnd->_nextSolPtr != NULL)
      {
        tmpEnd = tmpEnd->_nextSolPtr;
      }
      tmpEnd->_nextSolPtr = prevSolPtr->_nextSolPtr;
      prevSolPtr->_nextSolPtr->_previousSolPtr = tmpEnd;
    }
    /// attach this solution chain as next to prevSolPtr
    prevSolPtr->_nextSolPtr = this;
    _previousSolPtr = prevSolPtr;
  }
  return;
}

void Solution::includeVarSet(const VarPtrSet & varPtrSet)
{
  for (VarPtrSet::const_iterator varit = varPtrSet.begin(); varit != varPtrSet.end(); ++varit)
  {
    VarPtr2DoubleMap::iterator mapit = _solVarValMap.find(*varit);

    if (mapit != _solVarValMap.end())
    {
      mapit->second += (*varit)->val();
    }
    else
    {
      if((*varit)->isTypeOf(VcId::MastColumnMask))
        {
          (*varit)->incrParticipation(6);
        }
      _solVarValMap[*varit] = (*varit)->val();
    }
  }

  return;
}

const Double & Solution::solVal(Variable * vptr)
{
  return _solVarValMap[vptr];
}

void Solution::getVar(std::set< BcVar > & varSet)
{
    varSet.clear();

    for (VarPtr2DoubleMap::const_iterator it = solVarValMap().begin(); it != solVarValMap().end(); ++it)
    {
        it->first->solVal(it->second);

        if (printL(6))
            std::cout << "Solution::getVarSet()" << it->first->name() << " = " << it->first->solVal() << std::endl;

        InstanciatedVar * ivPtr = dynamic_cast<InstanciatedVar *>(it->first);
        if (ivPtr != NULL)
        {
            varSet.insert(BcVar(ivPtr));
        }
    }
}

void Solution::extractVar(std::set< BcVar > & varSet)
{
  varSet.clear();

  Solution * solPtr = this;

  while (solPtr != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin(); it != solPtr->solVarValMap().end(); ++it)
    {
      it->first->solVal(it->second);

      if (printL(6))
        std::cout << "Solution::getVarSet()" << it->first->name() << " = " << it->first->solVal() << std::endl;

      InstanciatedVar * ivPtr = dynamic_cast<InstanciatedVar *>(it->first);
      if (ivPtr != NULL)
      {
        varSet.insert(BcVar(ivPtr));
      }
    }
    solPtr = solPtr->_nextSolPtr;
  }
}

void Solution::extractVarWithGenericName(const std::string & name, std::set< BcVar > & varSet)
{
  varSet.clear();

  Solution * solPtr = this;

  while (solPtr != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin(); it != solPtr->solVarValMap().end(); ++it)
    {

      if (printL(6))
        std::cout << "Solution::extractVarWithGenericName()" << it->first->name() << " = " << it->first->solVal()
                  << std::endl;

      if (it->first->isTypeOf(VcId::InstanciatedVarMask))
      {
        if (it->first->genVarPtr()->defaultName() == name)
        {
          it->first->solVal(it->second);
          varSet.insert(BcVar(static_cast<InstanciatedVar*>(it->first)));
        }
      }
    }

    solPtr = solPtr->_nextSolPtr;
  }
}

void Solution::extractVarWithGenericName(const std::string & name, ProbConfig * probConfigPtr,
                                         std::set< BcVar > & varSet)
{
  varSet.clear();

  Solution * solPtr = this;
  while (solPtr != NULL)
  {
    if (solPtr->_probConfPtr == probConfigPtr)
    {
      for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin(); it != solPtr->solVarValMap().end();
           ++it)
      {
        if (printL(6))
          std::cout << "Solution::extractVarWithGenericName()" << it->first->name() << " = " << it->first->solVal()
                    << std::endl;

        if (it->first->isTypeOf(VcId::InstanciatedVarMask))
        {
          if (it->first->genVarPtr()->defaultName() == name)
          {
            it->first->solVal(it->second);
            varSet.insert(BcVar(static_cast<InstanciatedVar*>(it->first)));
          }
        }

      }
    }
    solPtr = solPtr->_nextSolPtr;
  }
}

void Solution::extractVarWithGenericName(const std::string & name, int firstIndex, std::set< BcVar > & varSet)
{
  varSet.clear();
  std::set< BcVar > _cacheVarPtrSet; /// chache for extractVarWithGenericName()

  Solution * solPtr = this;
  while (solPtr != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin(); it != solPtr->solVarValMap().end(); ++it)
    {
      if (printL(6))
        std::cout << "Solution::extractVarWithGenericName()" << it->first->name() << " = " << it->first->solVal()
                  << std::endl;

      InstanciatedVar * ivPtr = dynamic_cast<InstanciatedVar *>(it->first);
      if (ivPtr != NULL)
      {
        if ((ivPtr->genVarPtr()->defaultName() == name) && (ivPtr->id().first() == firstIndex))
        {
          ivPtr->solVal(it->second);
          varSet.insert(BcVar(ivPtr));
        }
      }
    }
    solPtr = solPtr->_nextSolPtr;
  }
}

void Solution::includeVars(const VarPtr2DoubleMap & varValMap, const bool & cumulativeVal)
{
  for (VarPtr2DoubleMap::const_iterator mapIt = varValMap.begin();
      mapIt != varValMap.end(); ++mapIt)
    includeVar(mapIt->first, mapIt->second, cumulativeVal);
}

void Solution::includeVar(Variable * varPtr, const Double & val, const bool & cumulativeVal)
{
  if (printL(6))
    std::cout <<"Var " <<varPtr->name() <<" with val " << val << std::endl;

    VarPtr2DoubleMap::iterator mapit = _solVarValMap.find(varPtr);

    if (mapit != _solVarValMap.end())
    {
      if (cumulativeVal)
      {
        mapit->second += val;
      }
      
      else
      {
        _solVarValMap[varPtr] = val;
      }
    }
    else
    {
      if(varPtr->isTypeOf(VcId::MastColumnMask))
        {
          varPtr->incrParticipation(8);
        }
      
      _solVarValMap[varPtr] = val;
    }

  return;
}

const VarPtr2DoubleMap & Solution::solVarValMap() const
{
  return (_solVarValMap);
}

const Double & Solution::resetCost()
{
  _cost = 0;

  for (VarPtr2DoubleMap::const_iterator it = _solVarValMap.begin(); it != _solVarValMap.end(); ++it)
  {
    if (printL(5))
      std::cout << " Solution::resetCost  sol[" << it->first->name() << "] = " << it->second << " cost = "
                << it->first->costrhs() << std::endl;

    _cost += it->second * it->first->costrhs();
  }

  return (_cost);
}

void Solution::deleteSolutionsChain()
{
  Solution * nextSolPtr = _nextSolPtr;
  while (nextSolPtr != NULL)
    {
      Solution * tempSolPtr = nextSolPtr->_nextSolPtr;
      delete nextSolPtr;
      nextSolPtr = tempSolPtr;
    }
  _nextSolPtr = NULL;
}

/// Virtual copy constructor
Solution * Solution::clone() const
{
  return (new Solution(*this));
}

bool Solution::operator<(const Solution & that) const
{
  if (_probConfPtr < that._probConfPtr)
    return (true);

  if (_probConfPtr > that._probConfPtr)
    return (false);

  if (cost() < that.cost())
    return (true);

  if (cost() > that.cost())
    return (false);

  return ((ref() < that.ref()));
}

void Solution::printDetailedSolution(std::ostream& os) const
{
    for (VarPtr2DoubleMap::const_iterator mapIt = solVarValMap().begin(); mapIt != solVarValMap().end(); ++mapIt)
    {
        os << mapIt->first->name() << " = " << mapIt->second;
        if (mapIt->first->isTypeOf(VcId::MastColumnMask))
        {
            MastColumn *colPtr = static_cast<MastColumn *>(mapIt->first);
            os << ", spId = " << colPtr->cgSpConfPtr()->id().first();
            os << ", treatOrderId = " << colPtr->treatOrderId();
            if (colPtr->spSol()->enumeratedFlag())
                os << ", enumerated";
            if (colPtr->spSol() != NULL) {
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
                         mapIt != colPtr->spSol()->solVarValMap().end(); ++mapIt) {
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
        } else {
            os << std::endl;
        }
    }
}

void Solution::shortPrint(std::ostream& os) const
{
  for (VarPtr2DoubleMap::const_iterator it = solVarValMap().begin(); it != solVarValMap().end(); ++it)
    {
      if (it != solVarValMap().begin())
        os << ", ";
      os << it->first->name() << "=" << it->second ;
    }
}

void Solution::printOrderedSolution(std::ostream& os) const
{
  const NetworkFlow * netFlowPtr = probConfPtr()->networkFlowPtr();

  if (
#ifdef BCP_RCSP_IS_FOUND
      ((_rcspSolPtr == nullptr) || (probConfPtr()->rcspGraphPtr() == nullptr))
      && 
#endif //BCP_RCSP_IS_FOUND
      (_orderedIds.empty() || (netFlowPtr == NULL)))
        return;

  std::vector<std::vector<double> > resConsumption;
  std::vector<int> vertexIds;
  std::vector<int> arcIds;
#ifdef BCP_RCSP_IS_FOUND
  os << "   Ordered solution " << (_rcspSolPtr->enumeratedFlag ? " (enum.) " : " (not. enum.)") << " : ";
  arcIds = _rcspSolPtr->arcIds;
  if ((_rcspSolPtr != nullptr) && (probConfPtr()->rcspGraphPtr() != nullptr))
  {
      if (!probConfPtr()->rcspGraphPtr()->obtainVertexIds(_rcspSolPtr, vertexIds))
      {
          os << "could not retrieve" << std::endl;
          return;
      }
      std::vector<int> stdResourceIndices;
      for (int resIndex = 0; resIndex < probConfPtr()->rcspGraphPtr()->resources.size(); ++resIndex)
      {
           if (probConfPtr()->rcspGraphPtr()->resources[resIndex].type != bcp_rcsp::ResourceData::Selection)
               stdResourceIndices.push_back(resIndex);
      }
      for (auto & vertResCons : _rcspSolPtr->resConsumption)
      {
          resConsumption.emplace_back();
          resConsumption.back().reserve(stdResourceIndices.size());
          for (auto resIndex: stdResourceIndices)
              resConsumption.back().push_back(vertResCons[resIndex]);
      }
  }
  else
#endif //BCP_RCSP_IS_FOUND
  {
      arcIds = _orderedIds;
      for (std::vector<int>::const_iterator arcIdIt = _orderedIds.begin(); arcIdIt != _orderedIds.end(); ++arcIdIt)
      {
          if (arcIdIt == _orderedIds.begin())
            vertexIds.push_back(netFlowPtr->netArcPtr(*arcIdIt)->tailVertexPtr()->id());
          vertexIds.push_back(netFlowPtr->netArcPtr(*arcIdIt)->headVertexPtr()->id());
      }
      resConsumption = _resConsumption;
  }

  if (vertexIds.empty())
  {
      os << "empty" << std::endl;
      return;
  }

    std::vector<std::vector<double> >::const_iterator resConsIt = resConsumption.begin();
    auto arcIdIt = arcIds.begin();
    int prevVertId = vertexIds.front();
    os << prevVertId;
    std::vector<double>::const_iterator prevResIt, resIt;
    if (!resConsIt->empty())
    {
        resIt = resConsIt->begin();
        os << "(" << *resIt;
        for (++resIt; resIt != resConsIt->end(); ++resIt)
            os << "," << *resIt;
        os << ")";
    }
    std::vector<int>::const_iterator vertIdIt = vertexIds.begin();
    std::vector<std::vector<double> >::const_iterator prevResConsIt = resConsIt;
    for (++vertIdIt, ++resConsIt; vertIdIt != vertexIds.end(); ++vertIdIt, ++resConsIt)
    {
        int thisVertId = *vertIdIt;
        bool resConsChanged = false;
        if (!resConsIt->empty())
        {
            for (resIt = resConsIt->begin(), prevResIt = prevResConsIt->begin(); resIt != resConsIt->end();
                 resIt++, prevResIt++)
                if (*resIt != *prevResIt)
                    resConsChanged = true;
        }
        else
        {
            resConsChanged = true;
        }
        //if (resConsChanged) : now show every vertex, as the models with idle arcs are rarely used
        {
            if (resConsIt->empty())
            {
                os << " -> (" << *arcIdIt << ") ";
                ++arcIdIt;
            }
            os << " -> " << thisVertId;
            if (!resConsIt->empty())
            {
                resIt = resConsIt->begin();
                os << "(" << *resIt;
                for (++resIt; resIt != resConsIt->end(); ++resIt)
                    os << "," << *resIt;
                os << ")";
            }
        }
        prevResConsIt = resConsIt;
        prevVertId = thisVertId;
    }
    os << std::endl;
}

std::ostream& Solution::print(std::ostream& os) const
{
  os << "Solution: ref = " << ref();
  os << "   solCost = " << _cost << std::endl;
  os << "   multiplicity = " << _multiplicity << std::endl;
  os << "   nbVar = " << _solVarValMap.size() << std::endl;

  printOrderedSolution(os);

  for (VarPtr2DoubleMap::const_iterator it = solVarValMap().begin();
      it != solVarValMap().end(); ++it)
  {
    os << "Solution includes var[" << it->first->name() << "("<< it->first->costrhs() << ")] = " << it->second
       << std::endl;
  }
  
  os << "-------------------------" << std::endl;

  if (_nextSolPtr != NULL)
    return _nextSolPtr->print(os);
  return (os);
}

std::ostream& Solution::printVar(std::ostream& os) const
{
  for (VarPtr2DoubleMap::const_iterator it = solVarValMap().begin(); it != solVarValMap().end(); ++it)
  {
    os << "Solution includes var[" << it->first->name() << "] = " << it->second << std::endl;

    it->first->print(os);
  }

  return (os);
}

void Solution::clear()
{
  for (VarPtr2DoubleMap::const_iterator mapIt = _solVarValMap.begin(); mapIt != _solVarValMap.end(); ++mapIt)
    {
      if (mapIt->first->isTypeOf(VcId::MastColumnMask))
        {
          mapIt->first->decrParticipation(4);
        }
    }
  _orderedIds.clear();
  _resConsumption.clear();
  _enumeratedFlag = false;
  _solVarValMap.clear();
  _cost = 0;
}

DualSolution::DualSolution(const DualSolution & that) :
    _probConfPtr(that._probConfPtr), _ref(0), _rhs(that._rhs), _dualSolConstrValMap(that._dualSolConstrValMap),
    _previousSolPtr(NULL), _nextSolPtr(NULL)
{
    if (_probConfPtr != NULL)
    {
        _ref = _probConfPtr->pcSolCount();
        _probConfPtr->increasePCSolCount();
    }
}

DualSolution::DualSolution(ProbConfig * probConfPtr) :
    _probConfPtr(probConfPtr), _ref(0), _rhs(0), _previousSolPtr(NULL), _nextSolPtr(NULL)
{
    if (_probConfPtr != NULL)
    {
        _ref = _probConfPtr->pcSolCount();
        _probConfPtr->increasePCSolCount();
    }
}

void DualSolution::includeConstr(const ConstrPtrSet & constrPtrSet)
{
    _rhs = 0;

    for (ConstrPtrSet::const_iterator it = constrPtrSet.begin(); it != constrPtrSet.end(); ++it)
    {
        includeConstr(*it, (*it)->val(), true);
        _rhs += (*it)->val() * (*it)->curRhs();
    }

    return;
}

void DualSolution::includeConstr(Constraint * constrPtr, const Double & val, const bool & cumulativeVal)
{
    if (cumulativeVal)
    {
        ConstrPtr2DoubleMap::iterator mapit = _dualSolConstrValMap.find(constrPtr);

        if (mapit != _dualSolConstrValMap.end())
        {
            mapit->second += val;
        }
        else
        {
            _dualSolConstrValMap[constrPtr] = val;
        }
    }
    else
    {
        _dualSolConstrValMap[constrPtr] = val;
    }

    return;
}

/// Virtual copy constructor
DualSolution * DualSolution::clone() const
{
    return (new DualSolution(*this));
}

const Double & DualSolution::computeTrueRhs()
{
    _rhs = 0;

    for (ConstrPtr2DoubleMap::const_iterator it = _dualSolConstrValMap.begin();
         it != _dualSolConstrValMap.end(); ++it)
    {
        _rhs -= it->second * it->first->costrhs() * it->first->sign();

        if (printL(6))
            std::cout << "      constr[" << it->first->name() << "] = " << it->second
                      << " and rhs = " << it->first->rhs() << " DualSol._rhs = " << _rhs << std::endl;
    }

    return (_rhs);
}

std::ostream& DualSolution::print(std::ostream& os) const
{
    os << "DualSolution, ref = " << ref() << " rhs = " << _rhs  << std::endl;
    for (ConstrPtr2DoubleMap::const_iterator it = _dualSolConstrValMap.begin();
         it != _dualSolConstrValMap.end(); ++it)
    {
        os << "      constr[" << it->first->name() << "] has val = " << it->second
           << " and rhs = " << it->first->rhs() << std::endl;
    }

    if (_nextSolPtr != NULL)
        return _nextSolPtr->print(os);
    return (os);
}

std::ostream& DualSolution::printConstr(std::ostream& os) const
{
    for (ConstrPtr2DoubleMap::const_iterator it = _dualSolConstrValMap.begin();
         it != _dualSolConstrValMap.end(); ++it)
    {
        os << "      constr[" << it->first->name() << "] = " << it->second
           << std::endl;

        it->first->print(os);
    }

    return (os);
}

void DualSolution::previousSolPtr(DualSolution * prevSolPtr)
{
    if (prevSolPtr != NULL)
    {
        /// attach next to prevSolPtr at the end of this solution chain
        if (prevSolPtr->_nextSolPtr != NULL)
        {
            DualSolution * tmpEnd = this;
            while (tmpEnd->_nextSolPtr != NULL)
            {
                tmpEnd = tmpEnd->_nextSolPtr;
            }
            tmpEnd->_nextSolPtr = prevSolPtr->_nextSolPtr;
            prevSolPtr->_nextSolPtr->_previousSolPtr = tmpEnd;
        }
        prevSolPtr->_nextSolPtr = this;
        _previousSolPtr = prevSolPtr;
    }
    return;
}
