/**
 *
 * This file bcCoefAboundC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcCoefAboundC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"

#include "bcProbConfigC.hpp"
#include "bcVarConstrC.hpp"
#include "bcGlobalException.hpp"
#include "bcColGenSpConfC.hpp"

/// Methods of class ProbCoef
ProbCoef::ProbCoef() :
    rowRef(-1), colRef(-1), coef(0)
{
    return;
}


ProbCoef::ProbCoef(const int & rowNb, const int & colNb, const Double & coefVal) :
    rowRef(rowNb), colRef(colNb), coef(coefVal)
{
  return;
}


ProbCoef::ProbCoef(const ProbCoef & pc) :
    rowRef(pc.rowRef), colRef(pc.colRef), coef(pc.coef)
{
}

bool ProbCoef::operator<(const ProbCoef & b) const
{
  if (colRef < b.colRef)
    return (true);

  if (colRef > b.colRef)
    return (false);

  if (rowRef < b.rowRef)
    return (true);

  return (false);
}

bool ProbCoef::operator==(const ProbCoef & b) const
{
  if (colRef != b.colRef)
    return (false);

  if (rowRef != b.rowRef)
    return (false);

  return (true);
}
std::ostream& ProbCoef::print(std::ostream& os) const
{
  os << "ProbCoef: rowRef= " << rowRef << ", colRef= " << colRef << ", coef= " << coef << std::endl;

  return (os);
}

/// Methods of class ProbSetCoef

ProbSetCoef::ProbSetCoef() :
    setRef(-1), setType(' '), colRef(-1), coef(0), setPriorityIndex(0)
{

  return;
}

ProbSetCoef::ProbSetCoef(const int & setNb, const char & type, const int & colNb, const Double & coefVal,
                         const Double & pri) :
    setRef(setNb), setType(type), colRef(colNb), coef(coefVal), setPriorityIndex(pri)
{
  return;
}

bool ProbSetCoef::operator<(const ProbSetCoef & b) const
{
  if (setRef < b.setRef)
    return (true);

  if (setRef > b.setRef)
    return (false);

  if (colRef < b.colRef)
    return (true);

  return (false);
}
bool ProbSetCoef::operator==(const ProbSetCoef & b) const
{
  if (setRef != b.setRef)
    return (false);

  if (colRef != b.colRef)
    return (false);

  return (true);
}
std::ostream& ProbSetCoef::print(std::ostream& os) const
{
  os << "ProbSetCoef: setRef= " << setRef << ", setType= " << setType << ", colRef= " << colRef
     << ", coef= " << coef << ", priIndex= " << setPriorityIndex << std::endl;

  return (os);
}

/// Methods of class ProbBound

ProbBound::ProbBound() :
    ref(-1), sense(' ')
{
  return;
}

ProbBound::ProbBound(const int & Lref, const char & Lsense, const Double & Lbound) :
    ref(Lref), sense(Lsense), bound(Lbound)
{
  return;
}

bool ProbBound::operator<(const ProbBound & b) const
{
  if (ref < b.ref)
    return (true);

  if (ref > b.ref)
    return (false);

  if (sense > b.sense)
    return (true);

  return (false);
}
bool ProbBound::operator==(const ProbBound & b) const
{
  if (ref != b.ref)
    return (false);

  if (sense != b.sense)
    return (false);

  return (true);
}
std::ostream& ProbBound::print(std::ostream& os) const
{
  return (os << "ref= " << ref << ", sense= " << sense << ", bound= " << bound << std::endl);
}

/// Methods of class ProbType

ProbType::ProbType() :
    ref(-1), type(' ')
{
  return;
}

ProbType::ProbType(const int & Lref, const char & Ltype) :
    ref(Lref), type(Ltype)
{
  return;
}

bool ProbType::operator<(const ProbType & b) const
{
  if (ref < b.ref)
    return (true);

  return (false);
}

bool ProbType::operator==(const ProbType & b) const
{
  if (ref == b.ref)
    return (true);

  return (false);
}
std::ostream& ProbType::print(std::ostream& os) const
{
  return os << "ref= " << ref << ", type= " << type << std::endl;
}

/// Methods of class ProbIntC

ProbIntC::ProbIntC() :
    ref(-1), type(' ')
{
  return;
}

ProbIntC::ProbIntC(const int & Lref, const char & Ltype, const Double & Lval) :
    ref(Lref), type(Ltype), val(Lval)
{
  return;
}

bool ProbIntC::operator<(const ProbIntC & b) const
{
  return (ref < b.ref);
}

bool ProbIntC::operator==(const ProbIntC & b) const
{
  return (ref == b.ref);
}

std::ostream& ProbIntC::print(std::ostream& os) const
{
  return (os << "ref= " << ref << ", type= " << type << ", val= " << val << std::endl);
}

/// Methods of class ComponentBound

ComponentBound::ComponentBound(InstanciatedVar * varPtr, const Double & val,
			                   const char & sign, const Double & card, const Double & compCard) :
  _varPtr(varPtr), _val(val), _sign(sign), _cardinality(card),
  _complementCard(compCard)
{
  return;
}

ComponentBound::ComponentBound(Variable * varPtr, const Double & val,
			                   const char & sign, const Double & card, const Double & compCard) :
  _val(val), _sign(sign), _cardinality(card), _complementCard(compCard)
{
  _varPtr = static_cast<InstanciatedVar *>(varPtr);

  if (_varPtr == NULL)
      std::cout << "ComponentBound::ComponentBound() require an InstanciatedVar" << std::endl;
}

void ComponentBound::complement()
{
  if (_sign == 'G')
    {
      _sign = 'L';
      _val = _val - 1;
    }
  else
    {
      _sign = 'G';
      _val = _val + 1;
    }

  Double tmp;
  tmp = _cardinality;
  _cardinality = _complementCard;
  _complementCard = tmp;

  return;
}

bool ComponentBound::operator==(const ComponentBound & that) const
{
  return (!((VcRefST(_varPtr, that._varPtr)) || (VcRefST(that._varPtr, _varPtr)) || (_sign != that._sign)
          || (_val < that._val) || (that._val < _val)));
}

bool ComponentBound::satisfiedBy(Variable * varPtr, const Double & fracVal) const
{
  if (VcRefST(_varPtr, varPtr))
    return (true);

  if (VcRefST(varPtr, _varPtr))
    return (true);

  if (_sign == 'G')
    return (fracVal >= _val);
  return (fracVal <= _val);
}

bool ComponentBound::satisfiedBy(const Double & val) const
{
  if (_sign == 'G')
    {
      if (val < _val)
        {
          return (false);
        }
      else
        return (true);
    }
  if (val > _val)
    return (false);

  return (true);
}

std::ostream& ComponentBound::print(std::ostream& os) const
{
  if (_varPtr == NULL)
    {
      os << "empty ComponentBound" << std::endl;
    }
  else
    {
      os << "var:" << _varPtr->name() << " " << (_sign == 'G' ? " >= " : " <= ")
          << " bound:" << _val << " c = " << _cardinality << " cc = " << _complementCard << std::endl;
    }

  return (os);
}

/// Methods of class ComponentSequence

ComponentSequence::ComponentSequence(ColGenSpConf * cgSpConfPtr) :
  std::vector<ComponentBound>(), 
  _cgSpConfPtr(cgSpConfPtr),
  _allBranchesGenerated(false), 
  _fracWeightRU(0), 
  _fracWeightRD(0),
  _activeSense('G'), /// Default is 'G' for empty class sequence or other component sequence
  _directPred(NULL),
  _additionalNbOfCpBd(0)
{
  clear();

  return;
}

ComponentSequence::ComponentSequence(const ComponentSequence & cs) :
  std::vector<ComponentBound>(cs),
  _cgSpConfPtr(cs._cgSpConfPtr),
  _allBranchesGenerated(cs._allBranchesGenerated),
  _fracWeight(cs._fracWeight),
  _fracWeightRU(cs._fracWeightRU), 
  _fracWeightRD(cs._fracWeightRD),
  _activeSense(cs._activeSense),
  _directPred(cs._directPred),
  _additionalNbOfCpBd(cs._additionalNbOfCpBd),
  _priorityFactor(cs._priorityFactor)
{
  return;
}

void ComponentSequence::reset() 
{
  _allBranchesGenerated = false;
  _fracWeight = 0;
  _fracWeightRU = 0;
  _fracWeightRD = 0;
  _activeSense = 'G'; /// Default is 'G' for all non-empty class sequence
  _directPred = NULL;
  _additionalNbOfCpBd = 0;
  _priorityFactor = 0;
  return;
}

ComponentSequence * ComponentSequence::directPred() const
{
  return _directPred;
}

ComponentSequence::InclusionStatus compareCbS(const ComponentSequence & listA, const ComponentSequence & listB)
{
  if (listA._cgSpConfPtr != listB._cgSpConfPtr)
    return (ComponentSequence::different);

  /// moved here from operator==(const CompSetInstMastBranchConstr & a, const CompSetInstMastBranchConstr & b)
  if (listA.empty() && listB.empty())
    return (ComponentSequence::identical);

  ComponentSequence::const_iterator cbAPt = listA.begin();
  ComponentSequence::const_iterator cbBPt = listB.begin();
  for (; (cbAPt != listA.end()) && (cbBPt != listB.end()); cbAPt++, cbBPt++)
    {
      if (*cbAPt != *cbBPt)
        return (ComponentSequence::different);
    }
  /// listA is a subset of listB
  if (cbAPt == listA.end())
    {
      /// listA is a superset of listB
      if (cbBPt == listB.end())
        /// Both list are identical
        return (ComponentSequence::identical);

      return (ComponentSequence::superclass);
    }
  /// listA is a superset of listB
  if (cbBPt == listB.end())
    return (ComponentSequence::subclass);

  return (ComponentSequence::different);
}

void ComponentSequence::complement()
{
  if (!empty())
    {
      rbegin()->complement();
      roundFracWeight();
      return;
    }
  /// Else Branching on Q
  if (_activeSense == 'L')
    _activeSense = 'G';
  else
    _activeSense = 'L';

  return;
}

const Double & ComponentSequence::classCardinality() const
{
  if (activeSense() == 'G')
    return (_fracWeightRU);

  return (_fracWeightRD);
}

char ComponentSequence::activeSense() const
{
  if (!empty())
    return ('G');

  return (_activeSense);
}

void ComponentSequence::fracWeight(const Double & fw)
{
  _fracWeight = fw;

  return;
}

void ComponentSequence::roundFracWeight()
{
  if (empty())
    {
      /// When branching on cardinality
      _fracWeightRU = Dceil(fracWeight());
      _fracWeightRD = _fracWeightRU - 1;
    }
  else
    {
      /// Round up to fw + 1 if fw is integer
      _fracWeightRU = Dfloor(fracWeight() + 1);
      _fracWeightRD = _fracWeightRU - 1;
    }

  return;
}

const Double & ComponentSequence::fracWeight() const
{
  if (empty())
    return (_fracWeight);

  return (rbegin()->cardinality());
}

void ComponentSequence::allCompLbVarPts(std::vector<InstanciatedVar *> & varPts) const
{
  varPts.clear();
  varPts.push_back(NULL);
  /// Go up to the last lower bound fixing : sign == 'G'
  for (std::vector<ComponentBound>::const_reverse_iterator it = this->rbegin();
       it != this->rend(); ++it)
    {
      if (it->sign() != 'L')
        varPts.push_back(it->varPtr());
    }
}

InstanciatedVar * ComponentSequence::lastCompLbVarPtr() const
{
  if (this->empty())
    return NULL;

  /// Go up to the last lower bound fixing : sign == 'G'
  std::vector<ComponentBound>::const_reverse_iterator it = this->rbegin();
  while (it->sign() == 'L')
    {
      ++it;
      if (it == this->rend())
        return NULL;
    }

  return it->varPtr();
}

bool ComponentSequence::satisfiedBy(Solution * solPtr) const
{
  if (solPtr == NULL)
    {
      for (ComponentSequence::const_iterator it = begin(); it != end(); ++it)
        if (it->sign() == 'G')
          return (false);

      return (true);
    }
  else
    {
      for (ComponentSequence::const_iterator it = begin(); it != end(); ++it)
        {
            VarPtr2DoubleMap::const_iterator valIt = solPtr->solVarValMap().find(it->varPtr());
            if (valIt != solPtr->solVarValMap().end())
            {
                if (!(it->satisfiedBy(valIt->second)))
                {
                    return false;
                }
            }
            else if (it->sign() == 'G')
                return (false);
        }
      return (true);
    }

  return (false);
}

std::ostream& ComponentSequence::print(std::ostream& os) const
{
  os << "     ComponentSequence " << std::endl
     << " cgSpConf = " << (_cgSpConfPtr != NULL ? _cgSpConfPtr->name() : "undefined") << std::endl;

  os << " ";
  for (ComponentSequence::const_iterator it = begin(); it != end(); ++it)
    os << *it;

  os << " activeSense = " << activeSense() << std::endl;
  os << " totalNbCol = " << _fracWeight << std::endl;
  os << " fracWeight = " << _fracWeight << std::endl;
  os << " fracWeightRU = " << _fracWeightRU << std::endl;
  os << " fracWeightRD = " << _fracWeightRD << std::endl;
  os << " classCardinality() = " << classCardinality() << std::endl;
  os << " directPred = " << _directPred << std::endl;

  return os;
}
