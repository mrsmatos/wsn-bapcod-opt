/**
 *
 * This file bcSpVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcMasterConfC.hpp"

using namespace std;

/**
 * 
 * @param priority higher priority index means chosen first for branching
 */
SubProbVariable::SubProbVariable(MasterConf * masterConfPtr,
                                 const IndexCell & id, 
				                 GenericVar * genVarPtr,
                                 ProbConfig * probConfigPtr,
                                 const std::string & name,
                                 const Double & costrhs, 
				                 const char & sense,
                                 const char & type, 
				                 const char & kind,
                                 const Double & upperBound,
                                 const Double & lowerBound, 
				                 const char & flag,
                                 const char & directive,
                                 const Double & priority, 
				                 const Double & val,
                                 const Double & globalUb,
                                 const Double & globalLb,
                                 const bool & presetMembership) :
InstanciatedVar(id, 
		        genVarPtr,
		        probConfigPtr,
		        std::string(name + "_OspV"),
		         costrhs,
		         sense,
		         type,
		         kind,
		         upperBound,
		         lowerBound,
		         flag,
		         directive,
		         priority,
		         val,
		         globalUb,
		         globalLb,
		         presetMembership),
		        _masterConfPtr(masterConfPtr),
		        _cgSpConfPtr(dynamic_cast<ColGenSpConf *>(probConfigPtr)),
		        _memorisedCurGlobUb(globalUb), _memorisedCurGlobLb(globalLb)
{
  bapcodInit().check(_cgSpConfPtr == NULL,
                     "SubProbVariable(): probConfigPtr should be of type ColGenSpConf *");

  if (printL(6))
    std::cout << "new SubProbVariable : lowerBound = " << lowerBound << this;
}

SubProbVariable::SubProbVariable(InstanciatedVar * iv, MasterConf * masterConfPtr) :
    InstanciatedVar(iv->id(),
                    iv->genVarPtr(),
                    iv->probConfPtr(),
		            std::string(iv->name() + "_CspV"),
                    iv->costrhs(), 
		            iv->sense(),
		            iv->type(),
		            iv->kind(),
                    iv->ub(), 
		            iv->lb(),
		            iv->flag(),
		            iv->directive(),
                    iv->priority(), 
		            iv->val(),
		            iv->globalUb(),
		            iv->globalLb(),
                    iv->presetMembership()),
    _masterConfPtr(masterConfPtr),
    _cgSpConfPtr(dynamic_cast<ColGenSpConf *>(iv->probConfPtr())),
    _memorisedCurGlobUb(iv->globalUb()),
    _memorisedCurGlobLb(iv->globalLb())
{
  bapcodInit().check(_cgSpConfPtr == NULL,
                     "SubProbVariable(): probConfigPtr should be of type ColGenSpConf *");

  buildMembershipHasBeenPerformed(iv->buildMembershipHasBeenPerformed());

  for (ConstVarConstrPtr2Double::const_iterator vcit = iv->member2coefMap().begin();
       vcit != iv->member2coefMap().end(); ++vcit)
    {
      if (vcit->first->isTypeOf(VcId::MasterConstrMask))
      {
          MasterConstr * mcPtr = dynamic_cast<MasterConstr *>(vcit->first);
          mcPtr->includeSubProbVarAsMember(this, vcit->second);
          includeMasterConstrAsMember(vcit->first, vcit->second);
      }
      else
      {
          includeMember(vcit->first, vcit->second, false);
      }
  }
}

SubProbVariable::~SubProbVariable()
{
  clearMembership();
  return;
}

MasterConf * SubProbVariable::masterConfPtr() const
{
  return (_masterConfPtr);
}

ColGenSpConf * SubProbVariable::cgSpConfPtr() const
{
  return (_cgSpConfPtr);
}

bool SubProbVariable::infeasible() const
{
  if (printL(6))
  {
    std::cout << " SubProbVariable::infeasible(): name() = " << name()
              << " curGlobLb() = " << curGlobLb() << " curGlobUb() = " << curGlobUb() << std::endl;
  }
  if (curGlobLb() > curGlobUb())
    return (true);

  return (Variable::infeasible());
}

bool SubProbVariable::membCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SubProbVariable::membCount this =  " << name() << ", that = "
              << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::MasterConstrMask))
    {
      if (printL(6))
	std::cout << "SubProbVariable::membCount MasterConstr" << std::endl;

      if (membershipUpToDate() && vcPtr->membershipUpToDate())
	{
	  if (printL(6))
	    std::cout
	      << "SubProbVariable::membCount membershipUpToDate Mast Constr "
	      << vcPtr->name() << std::endl;

	  return (masterConstrMember2coefMap().count(vcPtr));
	}
      /// else need to recompute count

      /// Already computed before
      if (masterConstrMember2coefMap().count(vcPtr))
	{
	  if (printL(6))
	    std::cout << "SubProbVariable::membCount count Mast Constr "
		      << vcPtr->name() << std::endl;

	  return true;
	}

      /// Already excluded before
      if (nonMemberSet().count(vcPtr))
	{
	  if (printL(6))
	    std::cout << "SubProbVariable:::membCount nonMember Mast Constr "
		      << vcPtr->name() << std::endl;

	  return false;
	}

      /// Do not know if already computed before
      LpCoef coef(computeCoef(vcPtr));
      if (coef.first)
	{
	  if (printL(6))
	    std::cout << "SubProbVariable:::membCount compute Mast Constr "
		      << vcPtr->name() << " coef = " << coef << std::endl;

	  includeMasterConstrAsMember(vcPtr, coef.second);
	  if ((vcPtr)->isTypeOf(VcId::MasterConstrMask))
	    {
	      dynamic_cast<MasterConstr*>(vcPtr)->includeSubProbVarAsMember(this,
									   coef.second);
	    }

	  return true;
	}
      else
	{
	  if (printL(6))
	    std::cout << "SubProbVariable:::membCount recordNonMember Mast Constr "
		      << vcPtr->name() << std::endl;

	  recordNonMember(vcPtr);
	  return false;
	}
    }

  return InstanciatedVar::membCount(vcPtr);
}

const Double & SubProbVariable::membCoef(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SubProbVariable:::membCoef  this =  " << name() << ", that = "
              << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::MasterConstrMask))
  {
    if (printL(6))
      std::cout << "SubProbVariable:::membCoef MasterConstr" << std::endl;

    if (membershipUpToDate() && vcPtr->membershipUpToDate())
    {
      if (printL(6))
        std::cout
            << "SubProbVariable:::membCoef membershipUpToDate Mast Constr "
            << vcPtr->name() << std::endl;

      return upToDateMasterConstrMembCoef(vcPtr);
    }

    /// Already computed before
    ConstVarConstrPtr2Double::const_iterator constrMemberIt = masterConstrMember2coefMap().find(vcPtr);
    if (constrMemberIt != masterConstrMember2coefMap().end())
    {
      if (printL(6))
        std::cout << "SubProbVariable:::membCoef count Mast Constr " << vcPtr->name() << std::endl;

      return constrMemberIt->second;
    }

    /// Already excluded before
    if (nonMemberSet().count(vcPtr))
    {
      if (printL(6))
        std::cout << "SubProbVariable:::membCoef nonMember Mast Constr " << vcPtr->name() << std::endl;

      return Double::staticZero;
    }

    /// Do not know if already computed before
    LpCoef coef(computeCoef(vcPtr));
    if (coef.first)
    {
      if (printL(6))
        std::cout << "SubProbVariable:::membCoef compute Mast Constr "
                  << vcPtr->name() << " coef = " << coef << std::endl;

	  if ((vcPtr)->isTypeOf(VcId::InstMasterConstrMask))
	    {
	      dynamic_cast<MasterConstr*>(vcPtr)->includeSubProbVarAsMember(this, coef.second);
	    }
        return includeMasterConstrAsMember(vcPtr, coef.second);
    }
    else
    {
      if (printL(6))
        std::cout << "SubProbVariable:::membCoef recordNonMember Mast Constr " << vcPtr->name() << std::endl;

      recordNonMember(vcPtr);
      return Double::staticZero;
    }
  }
  
  return InstanciatedVar::membCoef(vcPtr);
}

const Double & SubProbVariable::includeMasterConstrAsMember(VarConstr * mcPtr,
                                                            const Double & coef)
{
  if (printL(6))
    std::cout << "SubProbVariable::includeMasterConstrAsMember this " << name()
              << " that " << mcPtr->name() << " coef = " << coef << std::endl;
        
  ConstVarConstrPtr2Double::iterator it = _masterConstrMember2coefMap.find(
      mcPtr);
  if (it != _masterConstrMember2coefMap.end())
  {
    it->second += coef;
  }
  else
  {
    _masterConstrMember2coefMap[mcPtr] = coef;
  }
  
  if (printL(5))
    {
      std::cout << "SubProbVariable::includeMasterConstrAsMember var " << name()
                << " included in constr " << mcPtr->name() << " with coef = "
                << _masterConstrMember2coefMap[mcPtr] << std::endl;
    }

  return _masterConstrMember2coefMap[mcPtr];
}

const Double & SubProbVariable::includeMasterColAsMember(MastColumn * colPtr,
                                                         const Double & coef)
{
  MapMastColumnPtr2Double::iterator it = _masterColumnMember2coefMap.find(colPtr);
  if (it != _masterColumnMember2coefMap.end())
  {
    it->second += coef;
    return it->second;
  }
  else
  {
    _masterColumnMember2coefMap[colPtr] = coef;
  }
  return coef;
}

const Double & SubProbVariable::upToDateMasterConstrMembCoef(ConstVarConstrConstPtr vcPtr) const
{
  ConstVarConstrPtr2Double::const_iterator constrMemberIt = _masterConstrMember2coefMap.find(vcPtr);

  if (constrMemberIt != _masterConstrMember2coefMap.end())
  {
    return constrMemberIt->second;
  }
  else
  {
    return Double::staticZero;
  }
}

const Double & SubProbVariable::includeMember(VarConstr * vcPtr,
                                              const Double & coef,
                                              const bool & cumulativeCoef)
{
  if (printL(6))
    std::cout << "SubProbVariable::includeMember this =  " << name()
              << ", that = " << vcPtr->name() << ", coef = " << coef
              << std::endl;

  if (vcPtr->isTypeOf(VcId::MasterConstrMask))
    {
      MasterConstr * mcPtr = dynamic_cast<MasterConstr *>(vcPtr);
      includeMasterConstrAsMember(vcPtr, coef);
      return mcPtr->includeSubProbVarAsMember(this, coef);
    }

  return Variable::includeMember(vcPtr, coef, cumulativeCoef);
}

void SubProbVariable::eraseMasterConstrAsMember(VarConstr * mcPtr)
{
  if (printL(6))
    std::cout << "SubProbVariable::eraseMasterConstrAsMember() " << mcPtr->name() << std::endl;

  _masterConstrMember2coefMap.erase(mcPtr);
}

void SubProbVariable::setMembership()
{
  InstanciatedVar::setMembership();
  return;
}

/// Also check that the membership in master constraint
void SubProbVariable::enumerativeSetMembership()
{
  InstanciatedVar::enumerativeSetMembership();

  /** 
   * Only if not present membership as SP var 
   * in Master constraints are managed in buildMembership
   */
  if (!presetMembership())
  {
    if ((_masterConfPtr != NULL) && (_masterConfPtr->probPtr() != NULL))
      {
        Variable::setMembership(_masterConfPtr->probPtr()->probConstrSet());
      }
  }

  return;
}

void SubProbVariable::clearMembership()
{
  for (ConstVarConstrPtr2Double::iterator cPt = _masterConstrMember2coefMap.begin();
       cPt != _masterConstrMember2coefMap.end(); cPt++)
    {
        if (cPt->first->isTypeOf(VcId::MasterConstrMask))
        {
            MasterConstr * mcPtr = dynamic_cast<MasterConstr * >(cPt->first);
            mcPtr->eraseSubProbVarAsMember(this);
        }
    }

  _masterConstrMember2coefMap.clear();
  InstanciatedVar::clearMembership();

  return;
}

void SubProbVariable::recallMemorisedBounds()
{
  Variable::recallMemorisedBounds();

  return;
}

/// Measures the variable min contribution to the satisfaction of a specific  constraint
const Double SubProbVariable::lhsMinContrib(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::MasterConstrMask))
  {
    Double coef(vcPtr->membCoef(this));
    if (printL(6))
      std::cout << "lhsMinContrib(" << name() << ") coef = " << coef
                << " curGlobLb() = " << curGlobLb() << " curGlobUb() = "
                << curGlobUb() << std::endl;

    if (coef > 0)
      return (minGlobCurLb() * coef);
    else if (coef < 0)
      return (maxGlobCurUb() * coef);
    else
      bapcodInit().check(1, "SubProbVariable::lhsMinContrib: var should not be in membership map ",
                         ProgStatus::run);
  }

  return (VarConstr::lhsMinContrib(vcPtr));
}

/// Measures the variable max contribution to the satisfaction of a specific  constraint
const Double SubProbVariable::lhsMaxContrib(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::MasterConstrMask))
  {
    Double coef(vcPtr->membCoef(this));
    if (coef > 0)
      return (maxGlobCurUb() * coef);
    else if (coef < 0)
      return (minGlobCurLb() * coef);
    else
      bapcodInit().check(
          1,
          "SubProbVariable::lhsMaxContrib: var should not be in membership map ",
          ProgStatus::run);
  }

  return (VarConstr::lhsMaxContrib(vcPtr));
}

void SubProbVariable::resetBoundsAndCostToDefaults()
{
  _memorisedCurGlobLb = _globalLb;
  _memorisedCurGlobUb = _globalUb;
  _curLb = _memorisedCurLb = _lowerBound;
  _curUb = _memorisedCurUb = _upperBound;
  _memorisedCurCost = _costrhs;
}

const Double SubProbVariable::maxGlobCurUb() const
{
  if (_cgSpConfPtr->upperBoundMastConstrPtr() != NULL)
  {
    return ((_cgSpConfPtr->upperBoundMastConstrPtr()->curRhs() * curUb()));
  }
  else
    return (BapcodInfinity);
}

const Double SubProbVariable::minGlobCurLb() const
{
  if (_cgSpConfPtr->lowerBoundMastConstrPtr() != NULL)
  {
    if (printL(6))
      std::cout << "SubProbVariable::minGlobCurLb(): var " << name()
                << " cur Sp lb = "
                << _cgSpConfPtr->lowerBoundMastConstrPtr()->curRhs()
                << " curLb = " << curLb() << std::endl;
    return ((_cgSpConfPtr->lowerBoundMastConstrPtr()->curRhs() * curLb()));
  }

  return (0);
}

const Double & SubProbVariable::curCost() const
{
  if (printL(6))
    std::cout << " SubProbVariable::curCost() " << name() << " _costrhs = "
              << costrhs() << "  _memorisedCurCost = " << _memorisedCurCost
              << std::endl;

  return (_memorisedCurCost);
}

const Double & SubProbVariable::costrhs() const
{
  if (printL(6))
    std::cout << " SubProbVariable::costrhs() " << name() << " _costrhs = "
              << _costrhs << "  _memorisedCurCost = " << _memorisedCurCost << std::endl;

  return _costrhs;
}

void SubProbVariable::updateCurCostWithConstraint(const Constraint * constrPtr)
{
    if (param().SafeDualBoundScaleFactor() > 0)
        _memorisedCurCost += ceil(constrPtr->valOrSepPointVal()._val * param().SafeDualBoundScaleFactor());
    else
        _memorisedCurCost += constrPtr->valOrSepPointVal();
}

void SubProbVariable::calculateScaledCurCostForSafeDualBound(const bool & inPurePhaseOne)
{
   long long int scaledCurCost = 0;
   long long int scaleFactor = param().SafeDualBoundScaleFactor();

    if (!inPurePhaseOne)
        scaledCurCost = floor(costrhs()._val * scaleFactor);

    if (_masterConfPtr->probPtr()->inDualSol().size() <= _masterConstrMember2coefMap.size())
    {
        for (ConstrPtrSet::const_iterator cPt = _masterConfPtr->probPtr()->inDualSol().begin();
             cPt != _masterConfPtr->probPtr()->inDualSol().end(); cPt++)
            if (((*cPt)->type() != 'S') && (*cPt)->membCount(this))
            {
                scaledCurCost += ceil((*cPt)->valOrSepPointVal()._val * (*cPt)->membCoef(this)._val * scaleFactor);
            }
    }
    else
    {
        for (ConstVarConstrPtr2Double::iterator cPt = _masterConstrMember2coefMap.begin();
             cPt != _masterConstrMember2coefMap.end(); cPt++)
        {
            if ((cPt->first->type() != 'S') && cPt->first->inCurProb())
            {
                if (cPt->first->isTypeOf(VcId::ConstraintMask))
                {
                    Constraint * constrPtr = static_cast<Constraint *>(cPt->first);
                    if (_masterConfPtr->probPtr()->inDualSol().count(constrPtr))
                    {
                         scaledCurCost += ceil(constrPtr->valOrSepPointVal()._val * cPt->second._val * scaleFactor);
                    }
                }
            }
        }
    }
    _memorisedCurCost = (double)scaledCurCost;
}

void SubProbVariable::resetCost(const bool & inPurePhaseOne)
{
    int localprintlevel(6);

    /// Compute SP var reduced cost
    bapcodInit().require(_masterConfPtr != NULL,
                         "SubProbVariable::resetCost(): masterConfPtr should be well defined");

    if (param().SafeDualBoundScaleFactor() > 0)
    {
        calculateScaledCurCostForSafeDualBound(inPurePhaseOne);
        return;
    }

    if (printL(localprintlevel))
        std::cout << "SubProbVariable::resetCost() spCost[" << name()
                  << "] = _memorisedCurCost = " << _memorisedCurCost << std::endl;

    if (inPurePhaseOne)
    {
        if (printL(localprintlevel))
            std::cout << " SubProbVariable::costrhs() " << name() << " staticZero " << std::endl;
        _memorisedCurCost = Double::staticZero;
    }
    else
        _memorisedCurCost = costrhs();


    if (_masterConfPtr->probPtr()->inDualSol().size() <= _masterConstrMember2coefMap.size())
    {
        for (ConstrPtrSet::const_iterator cPt = _masterConfPtr->probPtr()->inDualSol().begin();
             cPt != _masterConfPtr->probPtr()->inDualSol().end(); cPt++)
            if (((*cPt)->type() != 'S') && (*cPt)->membCount(this))
            {
                _memorisedCurCost += (*cPt)->valOrSepPointVal() * (*cPt)->membCoef(this);

                if (printL(localprintlevel))
                    std::cout << "SubProbVariable[" << name() << "] in const[" << (*cPt)->name() << "] of val["
                              << (*cPt)->valOrSepPointVal() << "] has coef[" << (*cPt)->membCoef(this)
                              << "]  spCost= " << _memorisedCurCost << std::endl;
            }
    }
    else
    {
        for (ConstVarConstrPtr2Double::iterator cPt = _masterConstrMember2coefMap.begin();
             cPt != _masterConstrMember2coefMap.end(); cPt++)
        {
            if ((cPt->first->type() != 'S') && cPt->first->inCurProb())
            {
                if (printL(localprintlevel))
                    std::cout << "SubProbVariable[" << name() << "] in constr[" << cPt->first->name() << std::endl;
                if (cPt->first->isTypeOf(VcId::ConstraintMask))
                {
                    Constraint * constrPtr = static_cast<Constraint *>(cPt->first);

                    if (_masterConfPtr->probPtr()->inDualSol().count(constrPtr))
                    {
                        _memorisedCurCost += constrPtr->valOrSepPointVal() * cPt->second;

                        if (printL(localprintlevel))
                            std::cout << "SubProbVariable[" << name() << "] in const[" << constrPtr->name()
                                      << "] of val[" << constrPtr->valOrSepPointVal() << "] has coef[" << cPt->second
                                      << "]  spCost= " << _memorisedCurCost << std::endl;
                    }
                }
            }
        }
    }

    if (printL(localprintlevel))
        std::cout << "SubProbVariable::resetCost()  spCost[" << name()
                  << "] = _memorisedCurCost = " << _memorisedCurCost << std::endl;
    return;
}

std::ostream& SubProbVariable::print(std::ostream& os) const
{
  os << "SubProbVariable" << std::endl;
  InstanciatedVar::print(os);
  os << "    masterConstrMember2coefMap" << std::endl;

  for (ConstVarConstrPtr2Double::const_iterator itm = _masterConstrMember2coefMap.begin();
      itm != _masterConstrMember2coefMap.end(); ++itm)
  {
    os << "   coef[" << itm->first->name() << "] = " << itm->second << std::endl;
  }
  return (os);
}

bool SubProbVariable::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SubProbVariableMask, vcIdentifier);
}

/// Methods of class InstSubProbBranchingConstr 

InstSubProbBranchingConstr::InstSubProbBranchingConstr(const IndexCell& id,
                                                       GenericConstr * genConstrPtr,
                                                       ProbConfig* probConfigPtr,
                                                       const std::string& name,
                                                       const Double& costrhs,
                                                       const char& sense,
                                                       const char& type,
                                                       const char& kind,
                                                       const char& flag,
                                                       const Double& val,
                                                       const Double& upperBound,
                                                       const Double& lowerBound,
                                                       const char & directive,
                                                       const Double & priority) :
    InstanciatedConstr(id, genConstrPtr, probConfigPtr, name, costrhs, sense,
                       type, kind, flag, val, upperBound, lowerBound, directive,
                       priority),
    BranchingConstrBaseType(probConfigPtr)
{
  presetMembership(false);

  return;
}

/// this function is needed to determine whether a column violates a supbroblem branching constraint or not
bool InstSubProbBranchingConstr::violated(Variable * varPtr, const Double & useLevel)
{
  if (varPtr->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * colPtr = static_cast<MastColumn *>(varPtr);
      Solution * solPtr = colPtr->spSol();

      Double coeff(0.0);
      for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
           mapIt != solPtr->solVarValMap().end(); ++mapIt)
        coeff += membCoef(mapIt->first) * mapIt->second;
      return computeViolation(coeff * useLevel).positive();
    }
  return Constraint::violated(varPtr, useLevel);
}

std::ostream& InstSubProbBranchingConstr::print(std::ostream& os) const
{
  os << "InstSubProbBranchingConstr" << std::endl;
  InstanciatedConstr::print(os);

  return (os);
}

bool InstSubProbBranchingConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstSubProbBranchingConstrMask, vcIdentifier);
}
