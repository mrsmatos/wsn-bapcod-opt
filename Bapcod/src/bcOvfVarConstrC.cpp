/**
 *
 * This file bcOvfVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcNodeC.hpp"
#include "bcModelParameterC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcMasterConfC.hpp"

/**
 * Generic code
 */
using namespace std;

/// Methods of class OvfVar

OvfVar::OvfVar(ProbConfig * originatingPconfPtr, Variable * originatingVarPtr, int cnt):
  Variable(originatingPconfPtr->modelPtr(),
	    string("O") + cnt + originatingVarPtr->name(),
	    originatingVarPtr->costrhs(),
	    originatingVarPtr->sense(),
	    originatingVarPtr->type(),
	    originatingVarPtr->kind(),
	    originatingVarPtr->ub(),
	    originatingVarPtr->lb(),
	    originatingVarPtr->flag(),
	    originatingVarPtr->directive(),
	    /// Lower priority index means chosen first for branching
	    originatingVarPtr->priority(),
	    0,
	    originatingVarPtr->ub(),
	    originatingVarPtr->lb(),
	   false,
	    -1),
  _originatingPconfPtr(originatingPconfPtr),
  _originatingVarPtr(originatingVarPtr),
  _cnt(cnt)
{
  presetMembership(false);
  if (printL(6))
    std::cout << "OvfVar::OvfVar() new var name = " << name() << std::endl;

  _originatingVarPtr->includeOvfPtr(this);

}

OvfVar::OvfVar(ProbConfig * originatingPconfPtr, int cnt):
        Variable(originatingPconfPtr->modelPtr(),
                 "",
                 (Double) 0,
                 'P',
                 'C',
                 'E',
                 (Double) BapcodInfinity,
                 (Double) 0,
                 's',
                 'U',
                /// Lower priority index means chosen first for branching
                 (Double) 0,
                 0,
                 (Double) BapcodInfinity,
                 (Double) 0,
                 false,
                 -1),
        _originatingPconfPtr(originatingPconfPtr),
        _originatingVarPtr(NULL),
        _cnt(cnt)
{
}

bool OvfVar::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "OvfVar::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  /// Invert call
  if (vcPtr->isTypeOf(VcId::SpecificOvfConstrMask))
    return static_cast<SpecificOvfConstr *>(vcPtr)->computeCount(this);

  if (vcPtr->isTypeOf(VcId::OvfConstrMask))
  {
    OvfConstr * ovcPtr = static_cast<OvfConstr * >(vcPtr);

    if (ovcPtr->originatingPconfPtr()->configType() == ProbConfig::master)
    {
      return(ovcPtr->originatingConstrPtr()->membCount(_originatingVarPtr));
    }
    else if (ovcPtr->originatingPconfPtr()->configType() == ProbConfig::colGenSp)
    {
      if ((ovcPtr->originatingPconfPtr() == _originatingPconfPtr)
          && (ovcPtr->cnt() == _cnt))
        return(ovcPtr->originatingConstrPtr()->membCount(_originatingVarPtr));
      else
        return(false);
    }
  }
  else
  {
    bapcodInit().check(0, "OvfVar::computeCount: parameter should be and OvfConstr");
  }

  return(false);
}

const LpCoef OvfVar::computeCoef(ConstVarConstrConstPtr  vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpecificOvfConstrMask))
    return vcPtr->computeCoef(this);

  if (vcPtr->isTypeOf(VcId::OvfConstrMask))
  {
    OvfConstr * ovcPtr = static_cast<OvfConstr * >(vcPtr);

    if (ovcPtr->originatingPconfPtr()->configType() == ProbConfig::master)
    {
      ovcPtr->originatingConstrPtr()->membCoef(_originatingVarPtr);
    }
    else if (ovcPtr->originatingPconfPtr()->configType() == ProbConfig::colGenSp)
    {
      if ((ovcPtr->originatingPconfPtr() == _originatingPconfPtr)
          && (ovcPtr->cnt() == _cnt))
        return(ovcPtr->originatingConstrPtr()->membCoef(_originatingVarPtr));
      else
        return(Double::staticZero);
    }
  }
  else
  {
    bapcodInit().check(0, "OvfVar::computeCoef: parameter should be and OvfConstr");
  }

  return(LpCoef::ZeroCoef);
}

bool OvfVar::membCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "OvfVar::membCount this " << name() << " that "  << vcPtr->name() << std::endl;

  return Variable::membCount(vcPtr);

}

const Double & OvfVar::membCoef(ConstVarConstrConstPtr  vcPtr)
{
  if (printL(6))
    std::cout << "OvfVar::membCoef this " << name() << " that "  << vcPtr->name() << std::endl;

  return Variable::membCoef(vcPtr);

}

void OvfVar::setMembership()
{
    if (printL(6))
        std::cout <<  name() << std::endl;

    if (_originatingVarPtr != NULL)
    {
        if (_originatingVarPtr->presetMembership())
        {
            presetMembership(true);

            if (!_originatingVarPtr->buildMembershipHasBeenPerformed())
            {

                if (_originatingVarPtr->isTypeOf(VcId::InstanciatedVarMask))
                {
                    _originatingVarPtr->genVarConstrPtr()->
                            buildMembership(static_cast<InstanciatedVar *>(_originatingVarPtr));
                    _originatingVarPtr->buildMembershipHasBeenPerformed(true);
                }
            }


            buildMembershipHasBeenPerformed(true);
        }
    }

    Variable::setMembership();

    return;
}

std::ostream&  OvfVar::print(std::ostream& os) const
{
  if (_originatingVarPtr != NULL)
    os <<   "OvfVar whose originating var is "  << _originatingVarPtr->name() << std::endl;

  return(os);
}

bool OvfVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::OvfVarMask, vcIdentifier);
}

/// Methods of class SpSetupOvfVar

SpSetupOvfVar::SpSetupOvfVar(ProbConfig * originatingPconfPtr,  int cnt):
  OvfVar(originatingPconfPtr,  cnt)
{
  name(name() + "suv" + originatingPconfPtr->ref());
  const ColGenSpConf * const cgSpConfPtr = dynamic_cast<const ColGenSpConf * const >(originatingPconfPtr);

  bapcodInit().check(cgSpConfPtr == NULL,
                     "SpSetupOvfVar::SpSetupOvfVar(): originatingPconfPt should be of type ColGenSpConf");
  costrhs(cgSpConfPtr->fixedCost());
  sense('P');
  type('B');
  ub(1);
  lb(0);
  flag('s');
  directive('U');
  priority(1);

  resetCost(false);
  recallMemorisedBounds();
  if (printL(6))
    std::cout << "SpSetupOvfVar::SpSetupOvfVar() new var name = " << name() << std::endl;

  print(cout);
}

bool SpSetupOvfVar::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SpSetupOvfVar::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::SpVarLbOvfConstrMask))
    {
      SpVarLbOvfConstr * ovlcPtr = static_cast<SpVarLbOvfConstr * >(vcPtr);
      return(((_originatingPconfPtr == ovlcPtr->originatingPconfPtr())
	      && (ovlcPtr->originatingVarPtr() == _originatingVarPtr)
	      && (_cnt == ovlcPtr->cnt())));
    }

  if (vcPtr->isTypeOf(VcId::SpVarUbOvfConstrMask))
    {
      SpVarUbOvfConstr * ovucPtr = static_cast<SpVarUbOvfConstr * >(vcPtr);
      return(((_originatingPconfPtr == ovucPtr->originatingPconfPtr())
	      && (ovucPtr->originatingVarPtr() == _originatingVarPtr)
	      && (_cnt == ovucPtr->cnt())));
    }

  if (vcPtr->isTypeOf(VcId::SpVarUbOvfConstrMask))
    {
      SpLbOvfConstr * ovlbcPtr = static_cast<SpLbOvfConstr * >(vcPtr);
      return((_originatingPconfPtr == ovlbcPtr->originatingPconfPtr()));
    }

  return(false);
}

const LpCoef SpSetupOvfVar::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpVarLbOvfConstrMask))
    {
      SpVarLbOvfConstr * ovlcPtr = static_cast<SpVarLbOvfConstr * >(vcPtr);
      if ((_originatingPconfPtr == ovlcPtr->originatingPconfPtr())
	  && (ovlcPtr->originatingVarPtr() == _originatingVarPtr)
	  && (_cnt == ovlcPtr->cnt()))
	return(LpCoef::UnitCoef);
      else
	return(LpCoef::ZeroCoef);
    }

  if (vcPtr->isTypeOf(VcId::SpVarUbOvfConstrMask))
    {
      SpVarUbOvfConstr * ovucPtr = static_cast<SpVarUbOvfConstr * >(vcPtr);
      if ((_originatingPconfPtr == ovucPtr->originatingPconfPtr())
	  && (ovucPtr->originatingVarPtr() == _originatingVarPtr)
	  && (_cnt == ovucPtr->cnt()))
	return(LpCoef::MinusOneCoef);
      else
	return(LpCoef::ZeroCoef);
    }
  if (vcPtr->isTypeOf(VcId::SpVarUbOvfConstrMask))
    {
      SpLbOvfConstr * ovlbcPtr = static_cast<SpLbOvfConstr * >(vcPtr);
      if (_originatingPconfPtr == ovlbcPtr->originatingPconfPtr())
	return(LpCoef::UnitCoef);
      else
	return(LpCoef::ZeroCoef);
    }
  return(LpCoef::ZeroCoef);
}

void SpSetupOvfVar::setMembership()
{
  if (printL(6))
    std::cout <<  name() << std::endl;

  Variable::setMembership();

  return;
}


std::ostream&  SpSetupOvfVar::print(std::ostream& os) const
{
  os <<   "SpSetupOvfVar "  << std::endl;
  OvfVar::print(os);

  return(os);
}

bool SpSetupOvfVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SpSetupOvfVarMask, vcIdentifier);
}


/// Methods of class OvfConstr

OvfConstr::OvfConstr(ProbConfig * originatingPconfPtr, Constraint * originatingConstrPtr, int cnt):
  Constraint(originatingPconfPtr->modelPtr(),
	         (originatingConstrPtr != NULL ? std::string( std::string("O") + cnt + originatingConstrPtr->name())
	                                       : std::string( std::string("O") + cnt)),
             (originatingConstrPtr != NULL?  originatingConstrPtr->costrhs() : (Double) 0 ),
             (originatingConstrPtr != NULL?  originatingConstrPtr->sense() : 'E' ),
             (originatingConstrPtr != NULL?  originatingConstrPtr->type() : ' ' ),
             (originatingConstrPtr != NULL?  originatingConstrPtr->kind() : 'E' ),
             (originatingConstrPtr != NULL?  originatingConstrPtr->flag() : 's' ),
             -1,
             0.0,
             BapcodInfinity,
             - BapcodInfinity,
             (originatingConstrPtr != NULL?  originatingConstrPtr->directive() : 'U'),
             (originatingConstrPtr != NULL?  originatingConstrPtr->priority() : Double(1)),
             (originatingConstrPtr != NULL?  originatingConstrPtr->presetMembership() : true),
             (originatingConstrPtr != NULL?  originatingConstrPtr->toBeUsedInPreprocessing() : true)),
             _originatingPconfPtr(originatingPconfPtr),
             _originatingConstrPtr(originatingConstrPtr),
             _cnt(cnt)
{
  presetMembership(false);
  if (printL(6))
    std::cout << "OvfConstr::OvfConstr() new constr name = " << name() << " rhs = "
	          << costrhs() << " curRhs() = " << curRhs() << std::endl;
  if (_originatingConstrPtr != NULL)
    {
      _originatingConstrPtr->includeOvfPtr(this);
    }
}

bool OvfConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "OvfConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  /// Invert call
  if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    return vcPtr->computeCount(this);

  if (vcPtr->isTypeOf(VcId::OvfVarMask))
  {
    OvfVar * ovvPtr = static_cast<OvfVar * >(vcPtr);

    if (_originatingPconfPtr->configType() == ProbConfig::master)
    {
      if (printL(6))
      {
        std::cout << "OvfConstr::membCount masterConstr " << _originatingConstrPtr->name()
                  << " original Var "  << ovvPtr->originatingVarPtr()->name() << std::endl;
      }

      return _originatingConstrPtr->membCount(ovvPtr->originatingVarPtr());
    }
    else if ((ovvPtr->originatingPconfPtr() == _originatingPconfPtr) && (ovvPtr->cnt() == _cnt))
    {
      if (printL(6))
      {
        std::cout << "OvfConstr::membCount spmakeConstr " << _originatingConstrPtr->name()
                  << " original Var "  << ovvPtr->originatingVarPtr()->name() << std::endl;
      }

      return(_originatingConstrPtr->membCount(ovvPtr->originatingVarPtr()));
    }
    else
    {
      return(false);
    }
  }
  else
    bapcodInit().check(0, "OvfConstr::computeCount: parameter should be and OvfVar");

  return(false);
}


const LpCoef OvfConstr::computeCoef(ConstVarConstrConstPtr  vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    /// Invert call
    return vcPtr->computeCoef(this);

  if (vcPtr->isTypeOf(VcId::OvfVarMask))
  {
    OvfVar * ovvPtr = static_cast<OvfVar *>(vcPtr);

    if (_originatingPconfPtr->configType() == ProbConfig::master)
    {
      return(_originatingConstrPtr->membCoef(ovvPtr->originatingVarPtr()));
    }
    /// Sp conf
    else if ((ovvPtr->originatingPconfPtr() == _originatingPconfPtr) &&
        (ovvPtr->cnt() == _cnt))
    {
      return(_originatingConstrPtr->membCoef(ovvPtr->originatingVarPtr()));
    }
    else
    {
      return(LpCoef::ZeroCoef);
    }
  }
  else
  {
    bapcodInit().check(0, "OvfConstr::computeCoef: parameter should be and OvfVar");
  }

  return(LpCoef::ZeroCoef);
}

bool OvfConstr::membCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "OvfConstr::membCount this " << name() << " that "  << vcPtr->name() << std::endl;

  return Constraint::membCount(vcPtr);
}

const Double & OvfConstr::membCoef(ConstVarConstrConstPtr  vcPtr)
{
  if (printL(6))
    std::cout << "OvfConstr::membCoef this " << name() << " that "  << vcPtr->name() << std::endl;

  return Constraint::membCoef(vcPtr);
}

void OvfConstr::setMembership()
{
  if (_originatingConstrPtr != NULL)
  {
    if (_originatingConstrPtr->presetMembership())
    {
      presetMembership(true);

      if (!_originatingConstrPtr->buildMembershipHasBeenPerformed())
      {
        _originatingConstrPtr->genVarConstrPtr()->buildMembership(
                static_cast<InstanciatedConstr * >(_originatingConstrPtr));
        _originatingConstrPtr->buildMembershipHasBeenPerformed(true);
      }

      int printLevel = 6;

      for (ConstVarConstrPtr2Double::iterator itm = _originatingConstrPtr->member2coefMap().begin();
          itm != _originatingConstrPtr->member2coefMap().end(); ++itm)
      {
        if (printL(printLevel))
          std::cout << "OvfVar::setMembership(): constr " << name()
                    << " _originatingConstrPtr " << _originatingConstrPtr->name() << " in var " << itm->first->name()
                    <<  std::endl;

        for (std::list< OvfVar * >::const_iterator ovfVarIt = itm->first->ovfVarList().begin();
            ovfVarIt != itm->first->ovfVarList().end();
            ++ovfVarIt)
        {
          if ( (*ovfVarIt)->cnt() == _cnt) //added by Ruslan
            includeMember(*ovfVarIt, itm->second, false);
        }
      }


      if (_originatingConstrPtr->isTypeOf(VcId::InstMasterConstrMask))
      {
        InstMasterConstr * mcPtr = static_cast< InstMasterConstr * >(_originatingConstrPtr);
        if (printL(printLevel))
          std::cout << "OvfVar::setMembership(): constr " << name() << " mcPtr->subProbVarMember2coefMap().size() "
                    <<  mcPtr->subProbVarMember2coefMap().size() <<  std::endl;
        for(MapSubProbVariablePtr2Double::iterator itm = mcPtr->subProbVarMember2coefMap().begin();
            itm != mcPtr->subProbVarMember2coefMap().end(); ++itm)
        {
          if (printL(printLevel))
            std::cout << "OvfVar::setMembership(): constr " << name()
                      << " _originatingConstrPtr " << _originatingConstrPtr->name()
                      << " in sp var " << itm->first->name() <<  std::endl;

          for (std::list< OvfVar * >::const_iterator ovfVarIt = itm->first->ovfVarList().begin();
               ovfVarIt != itm->first->ovfVarList().end(); ++ovfVarIt)
            includeMember(*ovfVarIt, itm->second, false);
        }
      }

      buildMembershipHasBeenPerformed(true);
    }
  }
  Constraint::setMembership();
}

const Double OvfConstr::curRhs() const
{
  return Constraint::curRhs();
}


std::ostream & OvfConstr::print(std::ostream& os) const
{
  if (_originatingConstrPtr != NULL)
    os <<   "OvfConstr whose originating constr is "  << _originatingConstrPtr->name() << std::endl;

  return(os);
}

bool OvfConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::OvfConstrMask, vcIdentifier);
}

/// Methods of class SpVarLbOvfConstr
SpVarLbOvfConstr::SpVarLbOvfConstr(ProbConfig * originatingPconfPtr,
				                   Variable * originatingVarPtr,
				                   int cnt):
  SpecificOvfConstr(originatingPconfPtr, cnt),
  _originatingVarPtr(originatingVarPtr)
{
  name(name() + "vlb"+ originatingPconfPtr->ref());
  costrhs(0);
  sense('L');
  type('E');
  flag( 's');
}

bool SpVarLbOvfConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SpVarLbOvfConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    {
      SpSetupOvfVar * ovsvPtr = static_cast<SpSetupOvfVar * >(vcPtr);
      return ((ovsvPtr->cnt() == _cnt) && (ovsvPtr->originatingPconfPtr() == _originatingPconfPtr));
    }

  if (vcPtr->isTypeOf(VcId::OvfVarMask))
    {
      OvfVar * ovvPtr = static_cast<OvfVar * >(vcPtr);
      return ((ovvPtr->cnt() == _cnt) && (ovvPtr->originatingVarPtr() == _originatingVarPtr)
	          && (ovvPtr->originatingPconfPtr() == _originatingPconfPtr));
    }

  return(false);
}

const LpCoef SpVarLbOvfConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    {
      SpSetupOvfVar * ovsvPtr = static_cast<SpSetupOvfVar * >(vcPtr);
      if ((ovsvPtr->cnt() == _cnt) && (ovsvPtr->originatingPconfPtr() == _originatingPconfPtr))
	    return(_originatingVarPtr->curLb());

      return(LpCoef::ZeroCoef);
    }

  if (vcPtr->isTypeOf(VcId::OvfVarMask))
    {
      OvfVar * ovvPtr = static_cast<OvfVar * >(vcPtr);
      if ((ovvPtr->cnt() == _cnt)
	  && (ovvPtr->originatingVarPtr() == _originatingVarPtr)
	  && (ovvPtr->originatingPconfPtr() == _originatingPconfPtr) )
        return(LpCoef::MinusOneCoef);

      return(LpCoef::ZeroCoef);
    }

  return(LpCoef::ZeroCoef);
}

void SpVarLbOvfConstr::setMembership()
{
  Constraint::setMembership();

  return;
}


std::ostream&  SpVarLbOvfConstr::print(std::ostream& os) const
{
  os <<   "SpVarLbOvfConstr "  << std::endl;
  OvfConstr::print(os);

  return(os);
}

bool SpVarLbOvfConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SpVarLbOvfConstrMask, vcIdentifier);
}


/// Methods of class SpVarUbOvfConstr

SpVarUbOvfConstr::SpVarUbOvfConstr(ProbConfig * originatingPconfPtr,
				                   Variable * originatingVarPtr,
				                   int cnt):
  SpecificOvfConstr(originatingPconfPtr, cnt),
  _originatingVarPtr(originatingVarPtr)
{
  name(name() + "vub" + originatingPconfPtr->ref());
  costrhs(0);
  sense('G');
  type('E');
  flag( 's');
}

bool  SpVarUbOvfConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SpVarUbOvfConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

 if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    {
      SpSetupOvfVar * ovsvPtr = static_cast<SpSetupOvfVar * >(vcPtr);
      return ((ovsvPtr->cnt() == _cnt) && (ovsvPtr->originatingPconfPtr() == _originatingPconfPtr));
    }

  if (vcPtr->isTypeOf(VcId::OvfVarMask))
    {
      OvfVar * ovvPtr = static_cast<OvfVar * >(vcPtr);
      return ((ovvPtr->cnt() == _cnt) && (ovvPtr->originatingVarPtr() == _originatingVarPtr)
	          && (ovvPtr->originatingPconfPtr() == _originatingPconfPtr) );
    }

  return(false);
}

const LpCoef SpVarUbOvfConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpSetupOvfVarMask))
    {
      SpSetupOvfVar * ovsvPtr = static_cast<SpSetupOvfVar * >(vcPtr);
      if  ((ovsvPtr->cnt() == _cnt) && (ovsvPtr->originatingPconfPtr() == _originatingPconfPtr))
	    return(_originatingVarPtr->curUb());

      return(LpCoef::ZeroCoef);
    }
  if (vcPtr->isTypeOf(VcId::OvfVarMask))
    {
      OvfVar * ovvPtr = static_cast<OvfVar * >(vcPtr);
      if ((ovvPtr->cnt() == _cnt) && (ovvPtr->originatingVarPtr() == _originatingVarPtr)
	      && (ovvPtr->originatingPconfPtr() == _originatingPconfPtr))
        return(LpCoef::MinusOneCoef);

      return(LpCoef::ZeroCoef);
    }

  return(LpCoef::ZeroCoef);
}

void SpVarUbOvfConstr::setMembership()
{
  Constraint::setMembership();

  return;
}

std::ostream&  SpVarUbOvfConstr::print(std::ostream& os) const
{
  os <<   "SpVarUbOvfConstr "  << std::endl;
  OvfConstr::print(os);

  return(os);
}

bool SpVarUbOvfConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SpVarUbOvfConstrMask, vcIdentifier);
}


/// Methods of class SpLbOvfConstr

SpLbOvfConstr::SpLbOvfConstr(ProbConfig * originatingPconfPtr):
    SpecificOvfConstr(originatingPconfPtr, 0)
{
  name(name() + "su");
  const ColGenSpConf * const cgSpConfPtr =
    dynamic_cast<const ColGenSpConf * const >(originatingPconfPtr);

  bapcodInit().check(cgSpConfPtr == NULL, "SpLbOvfConstr::SpLbOvfConstr(): "
                                                     "originatingPconfPt should be of type ColGenSpConf *");

  costrhs(*(cgSpConfPtr->lowerBoundPtr()));
  sense('G');
  type('E');
  flag( 's');
}

bool SpLbOvfConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "SpLbOvfConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  SpSetupOvfVar * ovsvPtr = dynamic_cast<SpSetupOvfVar * >(vcPtr);
  if (ovsvPtr)
    return((ovsvPtr->originatingPconfPtr() == _originatingPconfPtr));

  return(false);
}

const LpCoef SpLbOvfConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::SpLbOvfConstrMask))
    {
      SpSetupOvfVar * ovsvPtr = static_cast<SpSetupOvfVar * >(vcPtr);
      if (ovsvPtr->originatingPconfPtr() == _originatingPconfPtr)
	    return(LpCoef::UnitCoef);

      return(LpCoef::ZeroCoef);
    }

  return(LpCoef::ZeroCoef);
}

void SpLbOvfConstr::setMembership()
{
  Constraint::setMembership();

  return;
}

std::ostream&  SpLbOvfConstr::print(std::ostream& os) const
{
  os <<   "SpLbOvfConstr "  << std::endl;
  OvfConstr::print(os);

  return(os);
}

bool SpLbOvfConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SpLbOvfConstrMask, vcIdentifier);
}


bool SpecificOvfConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SpecificOvfConstrMask, vcIdentifier);
}
