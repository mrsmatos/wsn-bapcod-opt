/**
 *
 * This file bcMastVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcGenVarConstrC.hpp"
#include "bcIndexC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcModelC.hpp"
#include "bcNodeC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcStabilizationInfo.hpp"
#include "bcModelCutConstrC.hpp"

using namespace std;
using namespace VcIndexStatus;

/**
 * Sense need to define unbounded variable in the master
 * (for otherwise there is a dual variable associated with this bound)
 */
InstMasterVar::InstMasterVar(const IndexCell & id,
			                 GenericVar * genVarPtr,
			                 MasterConf * masterConfPtr,
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
					masterConfPtr,
					std::string(name) + "_OmastV",
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
					presetMembership)
{
}

InstMasterVar::InstMasterVar(InstanciatedVar * iv) :
  InstanciatedVar(iv->id(),
		  iv->genVarPtr(),
		  iv->probConfPtr(),
		  iv->name(),
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
		  iv->presetMembership())
{
}

InstMasterVar::~InstMasterVar()
{
}

std::ostream & InstMasterVar::print(std::ostream& os) const
{
  os << "InstMasterVar: " << std::endl;
  InstanciatedVar::print(os);

  return (os);
}

void InstMasterVar::enumerativeSetMembership()
{
}

bool InstMasterVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstMasterVarMask, vcIdentifier);
}

bool InstMasterVar::membCount(ConstVarConstrConstPtr vcPtr)
{
    if (printL(6))
        std::cout << "InstMasterVar::membCount() this =  " << name() << ", that = " << vcPtr->name() << std::endl;

    return InstanciatedVar::membCount(vcPtr);
}

const Double & InstMasterVar::membCoef(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "InstMasterVar::membCoef() this =  " << name() << ", that = " << vcPtr->name() << std::endl;

  return InstanciatedVar::membCoef(vcPtr);
}

bool InstMasterVar::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(7))
    std::cout << "InstMasterVar::computeCount this " << name() << " that " << vcPtr->name() << std::endl;

  return (InstanciatedVar::computeCount(vcPtr));
}

const LpCoef InstMasterVar::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  return InstanciatedVar::computeCoef(vcPtr);
}

/*
 * Methods of class MasterConstr
 */
MasterConstr::MasterConstr(MasterConf * masterConfPtr) :
    _masterConfPtr(masterConfPtr)
{
  return;
}

MasterConstr::~MasterConstr()
{
}

bool MasterConstr::countColumn(ConstVarConstrConstPtr vcPtr)
{

  if (vcPtr->isTypeOf(VcId::MastColumnMask) && this->isTypeOf(VcId::InstMasterConstrMask))
    {
      MastColumn * mastColumnPtr = static_cast<MastColumn *> (vcPtr);
      return mastColumnPtr->membCount(static_cast<InstMasterConstr *>(this));
    }
  return false;
}

Double MasterConstr::coefColumn(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::MastColumnMask) && this->isTypeOf(VcId::InstMasterConstrMask))
    {
      MastColumn * mastColumnPtr = static_cast<MastColumn *> (vcPtr);
      return (mastColumnPtr->membCoef(static_cast<InstMasterConstr *>(this)));
    }

  return (Double::staticZero);
}

const Double & MasterConstr::includeSubProbVarAsMember(SubProbVariable * spvPtr, const Double & coef)
{
  if (printL(6))
  {
    std::cout << "MasterConstr::includeSubProbVarAsMember  spVar " << spvPtr->name() << " coef = " << coef
              << " in constr " << (dynamic_cast<Constraint*> (this))->name() << std::endl;
  }

  MapSubProbVariablePtr2Double::iterator it = _subProbVarMember2coefMap.find(spvPtr);

  if (it != _subProbVarMember2coefMap.end())
  {
    it->second += coef;
    return it->second;
  }
  else
  {
    _subProbVarMember2coefMap[spvPtr] = coef;
  }
  return coef;
}

const Double & MasterConstr::includePureMastVarAsMember(InstMasterVar * pmvPtr, const Double & coef)
{
  if (printL(6))
    std::cout << "MasterConstr::includeSubProbVarAsMember  spVar " << pmvPtr->name() << " coef = " << coef << std::endl;

  ConstVarConstrPtr2Double::iterator it = _pureMastVarMember2coefMap.find(pmvPtr);

  if (it != _pureMastVarMember2coefMap.end())
  {
    it->second += coef;
    return it->second;
  }
  else
  {
    _pureMastVarMember2coefMap[pmvPtr] = coef;
  }
  return coef;
}

void MasterConstr::eraseSubProbVarAsMember(SubProbVariable * spvPtr)
{
  if (printL(6))
    std::cout << "MasterConstr::eraseSubProbVarAsMember() " << spvPtr->name() << std::endl;

  _subProbVarMember2coefMap.erase(spvPtr);
}

void MasterConstr::erasePureMastVarAsMember(VarConstr * spvPtr)
{
  if (printL(6))
    std::cout << "MasterConstr::erasepureMastVarAsMember() " << spvPtr->name()
    << std::endl;

  _pureMastVarMember2coefMap.erase(spvPtr);
}

void MasterConstr::clearSubProbVarMember()
{
  if (printL(6))
    std::cout << "MasterConstr::clearSubProbVarMember() " << std::endl;

  Constraint * cPtr = dynamic_cast<Constraint *> (this);
  if (cPtr != NULL)
  {
    for (MapSubProbVariablePtr2Double::iterator spvPt = _subProbVarMember2coefMap.begin();
         spvPt != _subProbVarMember2coefMap.end(); ++spvPt)
    {
      spvPt->first->eraseMasterConstrAsMember(cPtr);
    }
  }
  _subProbVarMember2coefMap.clear();

  return;
}

void MasterConstr::clearPureMastVarMember()
{
  _pureMastVarMember2coefMap.clear();

  return;
}

const Double & MasterConstr::upToDateSubProbVarMembCoef(SubProbVariable * vcPtr) const
{
  MapSubProbVariablePtr2Double::const_iterator varConstrIt = _subProbVarMember2coefMap.find(vcPtr);

  if (varConstrIt != _subProbVarMember2coefMap.end())
  {
    return varConstrIt->second;
  }
  return Double::staticZero;
}

bool MasterConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::MasterConstrMask, vcIdentifier);
}

/*
 * Methods of class InstMasterConstr
 */
InstMasterConstr::InstMasterConstr(InstanciatedConstr * ic) :
  MasterConstr(dynamic_cast<MasterConf *> (ic->probConfPtr())),
  InstanciatedConstr(ic->id(),
		     ic->genConstrPtr(),
		     ic->probConfPtr(),
		     ic->name(),
		     ic->costrhs(),
		     ic->sense(),
		     ic->type(),
		     ic->kind(),
		     ic->flag(),
		     ic->val(),
		     ic->ub(),
		     ic->lb(),
		     ic->directive(),
		     ic->priority(),
		     ic->presetMembership(),
		     ic->toBeUsedInPreprocessing())
{
  name(ic->name() + "_CmastC");

  if (printL(6))
    std::cout << "InstMasterConstr::InstMasterConstr(upcasting instanciatedConstr) " << name() << std::endl;

  buildMembershipHasBeenPerformed(ic->buildMembershipHasBeenPerformed());

  bool cumulativeCoef(false);

  for (ConstVarConstrPtr2Double::iterator vcit = ic->member2coefMap().begin();
       vcit != ic->member2coefMap().end(); ++vcit)
  {
    includeMember(vcit->first, vcit->second, cumulativeCoef);
  }
}

/**
 *
 *
 * @param priority For branching on constraints
 */
InstMasterConstr::InstMasterConstr(const IndexCell & id,
                                   GenericConstr * genConstrPtr,
                                   ProbConfig * probConfigPtr,
                                   const std::string & name,
                                   const Double & costrhs,
                                   const char & sense,
                                   const char & type,
                                   const char & kind,
                                   const char & flag,
                                   const Double & val,
                                   const Double & upperBound,
                                   const Double & lowerBound,
                                   const char & directive,
                                   const Double & priority,
                                   const bool & presetMembership,
                                   const bool & toBeUsedInPreprocessing,
                                   const bool & considerAsEqualityInPreprocessing) :
  MasterConstr(dynamic_cast<MasterConf *> (probConfigPtr)),
  InstanciatedConstr(id, genConstrPtr, probConfigPtr,
                     std::string(name) + "_OmastC",
                     costrhs, sense, type, kind, flag, val,
                     upperBound, lowerBound, directive,
                     priority, presetMembership, toBeUsedInPreprocessing,
                     considerAsEqualityInPreprocessing),
  _treatOrderId(-1)
{
  if (printL(6))
    std::cout << "InstMasterConstr::InstMasterConstr() " << name << " toBeUsedInPreprocessing "
	          << toBeUsedInPreprocessing << std::endl;

  if (_flag == 'd') /// dynamic constraints
    {
      /// TO DO : the sense of the dynamic constraints can be wrong in the constructor (it can be changed by user after)
      ///  so, adding of artificial variables should be moved (probably to MasterConf::castAndAddConstraint())
      if ((probConfigPtr->param().mastInitMode().status() == MasterInitMode::localAndGlobAc)
	      || (probConfigPtr->param().mastInitMode().status() == MasterInitMode::localArtCol)
	      || (probConfigPtr->param().mastInitMode().status() == MasterInitMode::incSolColAndLac))
	    addLocalArtVar(probConfigPtr->modelPtr()->objectiveSense());

      if ((probConfigPtr->param().colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
          || (probConfigPtr->param().colGenDualPriceSmoothingAlphaFactor() > 0))
        createStabInfo(probConfigPtr->modelPtr()->objectiveSense());
    }
}

InstMasterConstr::InstMasterConstr(InstMasterConstr * imcPtr,
                                   GenericConstr * genBrConstrPtr,
                                   const std::string & name,
                                   const Double & rhs,
                                   const char & sense,
                                   const char & type,
                                   const char & kind,
                                   const char & flag) :
    MasterConstr(*imcPtr),
    InstanciatedConstr(imcPtr, genBrConstrPtr, name, rhs, sense, type, kind, flag),
    _treatOrderId(-1)
{
}

InstMasterConstr::InstMasterConstr(const InstMasterConstr &ic) :
  MasterConstr(ic), InstanciatedConstr(ic), _treatOrderId(-1)
{
  if (printL(6))
    std::cout << "InstMasterConstr::InstanciatedConstr(that) " << name()
              << " presetMembership = " << VarConstr::presetMembership() << std::endl;
}

InstMasterConstr::~InstMasterConstr()
{
  clearMembership();

  if (_flag == 'd')
    {
      /// if associated local artificial variables were not added to the problem
      /// we delete them, otherwise there will be memory leak
      /// (if they were added to the problem, they will be deleted when erased from the
      ///  problem's var. manager)
      if ((_posLocalArtVarPtr != NULL) && (_posLocalArtVarPtr->problemPtr() == NULL))
        delete _posLocalArtVarPtr;
      _posLocalArtVarPtr = NULL;
      if ((_negLocalArtVarPtr != NULL) && (_negLocalArtVarPtr->problemPtr() == NULL))
        delete _negLocalArtVarPtr;
      _negLocalArtVarPtr = NULL;
    }
  return;
}

bool InstMasterConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(7))
    std::cout << "InstMasterConstr::computeCount this " << name() << " that " << vcPtr->name() << std::endl;

  return InstanciatedConstr::computeCount(vcPtr);
}

const LpCoef InstMasterConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  return InstanciatedConstr::computeCoef(vcPtr);
}

bool InstMasterConstr::membCount(ConstVarConstrConstPtr vcPtr)
{
    if (printL(6))
        std::cout << "InstMasterConstr::membCount() this =  " << name() << ", that = " << vcPtr->name() << std::endl;


    if (vcPtr->isTypeOf(VcId::SubProbVariableMask))
    {
        SubProbVariable * spvPtr = static_cast<SubProbVariable *> (vcPtr);
        if (membershipUpToDate() && spvPtr->membershipUpToDate())
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCount() membershipUpToDate SP var "
                          << spvPtr->name() << " count ? " << subProbVarMember2coefMap().count(spvPtr) << std::endl;

            return (subProbVarMember2coefMap().count(spvPtr));
        }
        /// Already computed before
        if (subProbVarMember2coefMap().count(spvPtr))
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCount() count SP var " << spvPtr->name() << std::endl;

            return true;
        }

        /// Already excluded before
        if (nonMemberSet().count(spvPtr))
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCount() nonMember SP var " << spvPtr->name() << std::endl;

            return false;
        }

        /// Do not know if already computed before
        LpCoef coef(computeCoef(spvPtr));

        if (coef.first)
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCount() compute SP var " << spvPtr->name() << " coef = " << coef
                          << std::endl;
            includeSubProbVarAsMember(spvPtr, coef.second);
            spvPtr->includeMasterConstrAsMember(this, coef.second);
            return true;
        }
        else
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCount() recordNonMember SP var " << spvPtr->name() << std::endl;

            recordNonMember(spvPtr);
            return false;

        }
    }

    return InstanciatedConstr::membCount(vcPtr);
}

const Double & InstMasterConstr::membCoef(ConstVarConstrConstPtr vcPtr)
{

    if (printL(6))
        std::cout << "InstMasterConstr::membCoef() this =  " << name() << ", that = " << vcPtr->name() << std::endl;

    if (vcPtr->isTypeOf(VcId::SubProbVariableMask))
    {
        SubProbVariable * spvPtr = static_cast<SubProbVariable *> (vcPtr);
        if (membershipUpToDate() && spvPtr->membershipUpToDate())
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCoef() membershipUpToDate SP var " << spvPtr->name() << std::endl;

            return upToDateSubProbVarMembCoef(spvPtr);
        }

        /// Already computed before
        MapSubProbVariablePtr2Double::const_iterator probVarIt = subProbVarMember2coefMap().find(spvPtr);

        if (probVarIt != subProbVarMember2coefMap().end())
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCoef() count SP var " << spvPtr->name() << std::endl;

            return probVarIt->second;
        }

        /// Already excluded before
        if (nonMemberSet().count(spvPtr))
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCoef() nonMember SP var " << spvPtr->name() << std::endl;

            return Double::staticZero;
        }

        /// Do not know if already computed before
        LpCoef coef(computeCoef(spvPtr));
        if (coef.first)
        {
            if (printL(7))
                std::cout << "InstMasterConstr::membCoef() compute SP var " << spvPtr->name()
                          << " coef = " << coef << std::endl;

            spvPtr->includeMasterConstrAsMember(this, coef.second);
            return includeSubProbVarAsMember(spvPtr, coef.second);
        }
        else
        {
            recordNonMember(spvPtr);
            return Double::staticZero;
        }
    }

    return InstanciatedConstr::membCoef(vcPtr);
}

const Double & InstMasterConstr::includeMember(VarConstr * vcPtr, const Double & coef, const bool & cumulativeCoef)
{
    if (printL(6))
        std::cout << "InstMasterConstr::includeMember this =  " << name()
                  << ", that = " << vcPtr->name() << ", coef = " << coef << std::endl;

    if (vcPtr->isTypeOf(VcId::SubProbVariableMask))
    {
        SubProbVariable * spvPtr = static_cast<SubProbVariable *> (vcPtr);
        spvPtr->includeMasterConstrAsMember(this, coef);
        return includeSubProbVarAsMember(spvPtr, coef);
    }
    else if (vcPtr->isTypeOf(VcId::InstMasterVarMask))
    {
        InstMasterVar * instMasterVarPtr = static_cast<InstMasterVar *> (vcPtr);
        includePureMastVarAsMember(instMasterVarPtr, coef);
    }

    return InstanciatedConstr::includeMember(vcPtr, coef, cumulativeCoef);
}

void InstMasterConstr::createStabInfo(const BcObjStatus::MinMaxIntFloat & minmax)
{
  /// we do not stabilize implicit, type S, convexity constraints
  if ((kind() == 'I') || (type() == 'S') || isTypeOf(VcId::InstMastConvexityConstrMask))
    return;

  _stabInfoPtr = new VarConstrStabInfo(this);

  LocalArtificialVar * artVarPtr(NULL);

  if (param().StabilFuncOuterAngle() > 0)
  {
    artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::PosOuterId, minmax,
                                       "ltap", _modelPtr->artVarCost());
    if (printL(6))
      std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                << " instMasterConstr name  " << name() << ", sense = " << sense()
                << ", objStatus =   " << minmax << std::endl;
    _stabInfoPtr->posOuterArtVarPtr(artVarPtr);

    artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::NegOuterId, minmax,
                                       "ltan", _modelPtr->artVarCost());
    if (printL(6))
      std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                << " instMasterConstr name  " << name() << ", sense = " << sense()
                << ", objStatus =   " << minmax << std::endl;
    _stabInfoPtr->negOuterArtVarPtr(artVarPtr);
  }

  if (param().StabilFuncInnerAngle() > 0)
  {
    artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::PosInnerId, minmax, "lgap", _modelPtr->artVarCost());
    if (printL(6))
      std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
      << " instMasterConstr name  " << name() << ", sense = " << sense()
      << ", objStatus =   " << minmax << std::endl;
    _stabInfoPtr->posInnerArtVarPtr(artVarPtr);

    artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::NegInnerId, minmax, "lgan", _modelPtr->artVarCost());
    if (printL(6))
      std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
      << " instMasterConstr name  " << name() << ", sense = " << sense()
      << ", objStatus =   " << minmax << std::endl;
    _stabInfoPtr->negInnerArtVarPtr(artVarPtr);
  }
}

bool InstMasterConstr::addLocalArtVar(const BcObjStatus::MinMaxIntFloat & minmax)
{
  if (printL(6))
    std::cout << " InstMasterConstr::addLocalArtVar TRYING to add a localArtVar in instMasterConstr name  "
              << name() << " subProbVarMember2coefMap().empty() "
              << subProbVarMember2coefMap().empty() << std::endl;

  /// Ignore implicit constraints
  if (kind() == 'I')
    return false;

  /// Ignore type S constraints
  if (type() == 'S')
    return false;

  LocalArtificialVar * artVarPtr(NULL);

  if ((posLocalArtVarPtr() == NULL) && (sense() != 'L'))
    {
      artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::PosLocalId, minmax, "lap",
                                         _modelPtr->artVarCost());
      if (printL(5))
        std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                  << " instMasterConstr name  " << name() << ", sense = " << sense() << ", objStatus =   "
                  << minmax << std::endl;
      posLocalArtVarPtr(artVarPtr);
      /// art. vars. associated to dynamic constraints are not stored in master conf.,
      /// as they may be deleted in the middle of Branch-and-Price
      if (_flag == 's')
        _masterConfPtr->addNonStabilizationLocalArtVar(artVarPtr);
    }
  if ((negLocalArtVarPtr() == NULL) && (sense() != 'G'))
    {
      artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::NegLocalId, minmax, "lan",
                                         _modelPtr->artVarCost());
      if (printL(5))
        std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                  << " instMasterConstr name  " << name() << ", sense = " << sense() << ", objStatus =   "
                  << minmax << std::endl;
      negLocalArtVarPtr(artVarPtr);
      /// art. vars. associated to dynamic constraints are not stored in master conf.,
      /// as they may be deleted in the middle of Branch-and-Price
      if (_flag == 's')
        _masterConfPtr->addNonStabilizationLocalArtVar(artVarPtr);
    }
  return (true);
}


void InstMasterConstr::setMembership()
{
    if (!buildMembershipHasBeenPerformed())
      {
        genVarConstrPtr()->buildMembership(this);
        buildMembershipHasBeenPerformed(true);
      }

    if (isTypeOf(VcId::NonLinearInstMastConstrMask))
      {
        /// TO DO;
      }
    else
      {
        bool cumulativeCoef(true);

        for (MapSubProbVariablePtr2Double::const_iterator spvpPt =
             subProbVarMember2coefMap().begin();
             spvpPt != subProbVarMember2coefMap().end(); ++spvpPt)
          {
            for (MapMastColumnPtr2Double::const_iterator mcPt =
                spvpPt->first->masterColumnMember2coefMap().begin();
                mcPt != spvpPt->first->masterColumnMember2coefMap().end();
                ++mcPt)
              {
                if (!param().UseColumnsPool() && (_flag == 'd')
                    && (mcPt->first->vcIndexStatus() != VcIndexStatus::Active)
                    && (mcPt->first->vcIndexStatus() != VcIndexStatus::Inactive))
                  continue;

                if (printL(6))
                  {
                    Double coeff(0.0);
                    coeff = includeMember(mcPt->first, mcPt->second * spvpPt->second, cumulativeCoef);
                    std::cout << " InstMasterConstr::setMembership["
                    << name() << " column = " << mcPt->first->name()
                    << ", coeff = " << coeff << std::endl;
                  }
                else
                    includeMember(mcPt->first, mcPt->second * spvpPt->second, cumulativeCoef);
              }

          }
      }

    Constraint::setMembership();

    return;
}

void InstMasterConstr::enumerativeSetMembership()
{
  if (!presetMembership())
  {
    if (_masterConfPtr != NULL)
    {
      for (std::vector<ColGenSpConf *>::const_iterator spcPt =
              _masterConfPtr->colGenSubProbConfPts().begin();
              spcPt != _masterConfPtr->colGenSubProbConfPts().end(); spcPt++)
        if ((*spcPt)->probPtr() != NULL)
        {
          Constraint::setMembership((*spcPt)->probPtr()->probVarSet());
        }
    }
  }

  return;
}

void InstMasterConstr::clearMembership()
{
  clearSubProbVarMember();
  VarConstr::clearMembership();

  return;
}

std::ostream & InstMasterConstr::print(std::ostream& os) const
{
  os << "InstMasterConstr" << std::endl;
  InstanciatedConstr::print(os);

  return (os);
}

void InstMasterConstr::nicePrint(std::ostream& os) const
{
	os << "Master Constraint " << name() << " :";
	for (ConstVarConstrPtr2Double::const_iterator mapIt = _pureMastVarMember2coefMap.begin();
		 mapIt != _pureMastVarMember2coefMap.end(); ++mapIt)
	{
		if (mapIt->second >= 0.0)
			os << "+";
		os << mapIt->second << "*" << mapIt->first->name();
	}
	for (MapSubProbVariablePtr2Double::const_iterator mapIt = _subProbVarMember2coefMap.begin();
		 mapIt != _subProbVarMember2coefMap.end(); ++mapIt)
	{
		if (mapIt->second >= 0.0)
			os << "+";
		os << mapIt->second << "*" << mapIt->first->name();
	}

	if (_sense == 'G')
		os << " >= ";
	else if (_sense == 'L')
		os << " <= ";
	else
		os << " == ";
	os << _costrhs << std::endl;
}

bool InstMasterConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstMasterConstrMask, vcIdentifier);
}

/*
 * Methods of class NonLinearInstMastConstr
 */

NonLinearInstMastConstr::NonLinearInstMastConstr(InstanciatedConstr * ic) :
    InstMasterConstr(ic)
{
  return;
}

NonLinearInstMastConstr::NonLinearInstMastConstr(const IndexCell & id,
                                                 GenericConstr * genConstrPtr, ProbConfig * probConfigPtr,
                                                 const std::string & name, const Double & costrhs, const char & sense,
                                                 const char & type, const char & kind, const char & flag,
                                                 const Double & val, const Double & upperBound,
                                                 const Double & lowerBound, const char & directive,
                                                 const Double & priority) :
    InstMasterConstr(id, genConstrPtr, probConfigPtr, name, costrhs, sense, type, kind, flag, val,
                     upperBound, lowerBound, directive, priority)
{
}

NonLinearInstMastConstr::NonLinearInstMastConstr(const NonLinearInstMastConstr &ic) :
InstMasterConstr(ic)
{
}

NonLinearInstMastConstr::~NonLinearInstMastConstr()
{
}

bool NonLinearInstMastConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "NonLinearInstMastConstr::computeCount this " << name() << " that " << vcPtr->name() << std::endl;

  if (MasterConstr::countColumn(vcPtr))
    return (true);

  /// Artificial var?
  return (InstMasterConstr::computeCount(vcPtr));
}

const LpCoef NonLinearInstMastConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (MasterConstr::countColumn(vcPtr))
    return MasterConstr::coefColumn(vcPtr);

  /// Artificial var?
  return InstMasterConstr::computeCoef(vcPtr);
}


std::ostream & NonLinearInstMastConstr::print(std::ostream& os) const
{
  os << "NonLinearInstMastConstr" << std::endl;
  InstMasterConstr::print(os);

  return (os);
}

bool NonLinearInstMastConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::NonLinearInstMastConstrMask, vcIdentifier);
}

/*
 * Methods of class InstMastConvexityConstr
 */
InstMastConvexityConstr::InstMastConvexityConstr(const IndexCell & id,
        GenericConstr * genConstrPtr,
        MasterConf * masterConfPtr,
        ColGenSpConf * cgSpConfPtr,
        const std::string & name,
        const Double & costrhs,
        const char & sense,
        const char & type,
        const char & kind,
        const char & flag,
        const int & index,
        const Double & val,
        const Double & upperBound,
        const Double & lowerBound) :
InstMasterConstr(id,
genConstrPtr,
masterConfPtr,
name,
costrhs,
sense,
type, /// should be 'X'
kind,
flag,
val,
upperBound,
lowerBound),
_cgSpConfPtr(cgSpConfPtr),
_memorisedCsBcConstr(NULL),
_locallyValidRhsDefined(false),
_locallyValidRhs(costrhs)
{
  //  printout("InstMastConvexityConstr::InstMastConvexityConstr");

  if (type != 'X')
    std::cout << "InstMastConvexityConstr ERROR wrong type" << name << std::endl;

  if (printL(7))
    std::cout << "InstMastConvexityConstr " << name << std::endl;
}

InstMastConvexityConstr::~InstMastConvexityConstr()
{
}

ColGenSpConf * InstMastConvexityConstr::cgSpConfPtr() const
{
  return (_cgSpConfPtr);
}

bool InstMastConvexityConstr::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "InstMastConvexityConstr::computeCount this " << name() << " that " << vcPtr->name() << std::endl;

  return (InstMasterConstr::computeCount(vcPtr));
}

const LpCoef InstMastConvexityConstr::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  return InstMasterConstr::computeCoef(vcPtr);
}

void InstMastConvexityConstr::enumerativeSetMembership()
{
  InstanciatedConstr::enumerativeSetMembership();

  return;
}

void InstMastConvexityConstr::newLocalRhs(const Double & localRhs)
{
  _locallyValidRhs = localRhs;
}

const Double & InstMastConvexityConstr::newLocalRhs() const
{
  return _locallyValidRhs;
}

void InstMastConvexityConstr::defineLocalRhs(const Double & localRhs)
{
  bool tighterBound(false);

  if (printL(5))
  {
    std::cout << "InstMastConvexityConstr::defineLocalRhs() " << name()
            << " sense() " << sense()
            << " localRhs " << localRhs
            << " _locallyValidRhs " << _locallyValidRhs
            << " curRhs() " << curRhs()
            << std::endl;
  }
  switch (sense())
  {
    case 'G':
    {
      tighterBound = (localRhs > _locallyValidRhs);
      break;
    }
    case 'L':
    {
      tighterBound = (localRhs < _locallyValidRhs);
      break;
    }
  }
  if (tighterBound)

  {
    _locallyValidRhsDefined = true;
    _locallyValidRhs = localRhs;

    if (printL(5))
    {
      std::cout << "InstMastConvexityConstr:: AFTER defineLocalRhs()  " << name()
              << " sense() " << sense()
              << " localRhs " << localRhs
              << " _locallyValidRhs " << _locallyValidRhs
              << " curRhs() " << curRhs()
              << std::endl;
    }
  }
  return;
}

void InstMastConvexityConstr::resetLocalRhs()
{
  _locallyValidRhs = curRhs();
  _locallyValidRhsDefined = false;
}

const Double & InstMastConvexityConstr::costrhs() const
{
  return (InstMasterConstr::costrhs());
}

void InstMastConvexityConstr::addMember(VarConstr * vcPtr)
{
  if (vcPtr->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * mastColumnPtr = static_cast<MastColumn *> (vcPtr);
      mastColumnPtr->addMember(this);
      return;
    }

  InstanciatedConstr::addMember(vcPtr);

  return;
}

void InstMastConvexityConstr::setMembership()
{
  VarConstr::setMembership();
  return;
}

void InstMastConvexityConstr::clearMembership()
{
  VarConstr::clearMembership();

  return;
}

std::ostream & InstMastConvexityConstr::print(std::ostream& os) const
{
  os << "InstMastConvexityConstr" << std::endl;
  InstMasterConstr::print(os);

  return (os);
}

bool InstMastConvexityConstr::addLocalArtVar(const BcObjStatus::MinMaxIntFloat & minmax)
{
    if (!param().LocArtVarInConvexityConstr())
        return false;

    LocalArtificialVar * artVarPtr(NULL);

    if ((posLocalArtVarPtr() == NULL) && (sense() != 'L'))
    {
        artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::PosLocalId, minmax, "lap",
                                           _modelPtr->artVarCost());
        if (printL(6))
            std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                      << " instMasterConstr name  " << name() << ", sense = " << sense() << ", objStatus =   "
                      << minmax << std::endl;
        posLocalArtVarPtr(artVarPtr);
        _masterConfPtr->addNonStabilizationLocalArtVar(artVarPtr);
    }
    if ((negLocalArtVarPtr() == NULL) && (sense() != 'G'))
    {
        artVarPtr = new LocalArtificialVar(this, LocalArtificialVar::NegLocalId, minmax, "lan",
                                           _modelPtr->artVarCost());
        if (printL(6))
            std::cout << " InstMasterConstr::addLocalArtVar add localArtVar " << artVarPtr->name()
                      << " instMasterConstr name  " << name() << ", sense = " << sense() << ", objStatus =   "
                      << minmax << std::endl;
        negLocalArtVarPtr(artVarPtr);
        _masterConfPtr->addNonStabilizationLocalArtVar(artVarPtr);
    }
    return true;
}

bool InstMastConvexityConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstMastConvexityConstrMask, vcIdentifier);
}

/*
 * Methods of class InstMasterBranchingConstr
 */

InstMasterBranchingConstr::InstMasterBranchingConstr(const IndexCell & id,
        GenericConstr * genBrConstrPtr, ProbConfig * probConfPtr,
        const std::string & name, const Double & rhs, const char & sense,
        const char & type, const char & kind, const char & flag, const Double& val,
        const Double& upperBound, const Double& lowerBound, const char & directive,
        const Double & priority) :
    InstMasterConstr(id, genBrConstrPtr, probConfPtr, name, rhs, sense, type, kind, flag, val,
                     upperBound, lowerBound, directive, priority),
    BranchingConstrBaseType(probConfPtr)
{
  presetMembership(false);
}

InstMasterBranchingConstr::InstMasterBranchingConstr(InstanciatedConstr * icPtr) :
    InstMasterConstr(icPtr),
    BranchingConstrBaseType(icPtr->probConfPtr())
{
  presetMembership(false);
}

InstMasterBranchingConstr::InstMasterBranchingConstr(InstMasterConstr * imcPtr,
        GenericConstr * genBrConstrPtr, const std::string & name,
        const Double & rhs, const char & sense, const char & type,
        const char & kind, const char & flag) :
    InstMasterConstr(imcPtr, genBrConstrPtr, name, rhs, sense, type, kind, flag),
    BranchingConstrBaseType(imcPtr->probConfPtr())
{
  presetMembership(false);
}

InstMasterBranchingConstr::InstMasterBranchingConstr(const InstMasterBranchingConstr & that) :
    InstMasterConstr(that), BranchingConstrBaseType(that.probConfPtr())
{
  presetMembership(false);
}

InstMasterBranchingConstr::~InstMasterBranchingConstr()
{
  return;
}


std::ostream & InstMasterBranchingConstr::print(std::ostream& os) const
{
  os << "InstMasterBranchingConstr" << std::endl;
  InstMasterConstr::print(os);

  return (os);
}

bool InstMasterBranchingConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstMasterBranchingConstrMask, vcIdentifier);
}

/*
 * Methods of class BasicConstrInstMastBranchingConstr
 */

BasicConstrInstMastBranchingConstr
::BasicConstrInstMastBranchingConstr(const IndexCell & id, GenericConstr * genBrConstrPtr,
                                     ProbConfig * probConfPtr, InstanciatedConstr * instConstrPtr,
                                     const std::string & description, const std::string & name,
                                     const Double & rhs, const char & sense, const char & type,
                                     const char & kind, const char & flag) :
    InstMasterBranchingConstr(id, genBrConstrPtr, probConfPtr, name, rhs, sense, type, kind, flag),
    _instConstrPtr(instConstrPtr), _description(description)
{
  presetMembership(false);
  if (printL(6))
    std::cout << "BasicConstrInstMastBranchingConstr() " << name
              << " presetMembership =  " << presetMembership() << std::endl;
  return;
}

BasicConstrInstMastBranchingConstr
::BasicConstrInstMastBranchingConstr(InstMasterConstr * imcPtr, GenericConstr * genBrConstrPtr,
                                     const std::string & description, const std::string & name,
                                     const Double & rhs, const char & sense,
                                     const char & type, const char & kind, const char & flag) :
    InstMasterBranchingConstr(imcPtr, genBrConstrPtr, name, rhs, sense, type, kind, flag),
    _instConstrPtr(imcPtr)
{
  presetMembership(false);
  if (printL(6))
    std::cout << "BasicConstrInstMastBranchingConstr() " << name
              << " presetMembership =  " << presetMembership() << std::endl;

  return;
}

BasicConstrInstMastBranchingConstr::~BasicConstrInstMastBranchingConstr()
{
  return;
}

std::vector<std::string> BasicConstrInstMastBranchingConstr::forDotPrint() const
{
  std::stringstream ss;
  shortPrint(ss);
  std::vector<std::string> dotStrings(1, std::string(ss.str()));
  return dotStrings;
}

void BasicConstrInstMastBranchingConstr::shortPrint(std::ostream& os) const
{
  os << _description;
  switch(_sense)
    {
    case 'G':
      os << " >= ";
      break;
    case 'L':
      os << " <= ";
      break;
    case 'E':
      os << " == ";
      break;
    default:
      os << " ?= ";
      break;
    }
  os << _costrhs << " ";
}

void BasicConstrInstMastBranchingConstr::setMembership()
{
  if (printL(6))
    std::cout << "BasicConstrInstMastBranchingConstr::setMembership() genVarConstrPtr() ="
              << genVarConstrPtr()->defaultName() << "  constr=" << _instConstrPtr->name() << std::endl;

  if(!buildMembershipHasBeenPerformed())
    {
      presetMembership(true);
      buildMembershipHasBeenPerformed(true);
    }

  /// here we basically copy the membership information of _instConstrPtr
  for (ConstVarConstrPtr2Double::const_iterator mapIt = _instConstrPtr->member2coefMap().begin();
       mapIt != _instConstrPtr->member2coefMap().end(); ++mapIt)
    {
      includeMember(mapIt->first, mapIt->second, false);
    }

  InstMasterConstr * instMastConstrPtr = static_cast<InstMasterConstr *>(_instConstrPtr);

  /// compute columns membership (it was not set for _instConstrPtr)
  for (MapSubProbVariablePtr2Double::const_iterator mapIt = instMastConstrPtr->subProbVarMember2coefMap().begin();
       mapIt != instMastConstrPtr->subProbVarMember2coefMap().end(); ++mapIt)
    {
      SubProbVariable * spVarPtr = mapIt->first;
      includeMember(spVarPtr, mapIt->second, false);
      for (MapMastColumnPtr2Double::const_iterator colMapIt = spVarPtr->masterColumnMember2coefMap().begin();
           colMapIt != spVarPtr->masterColumnMember2coefMap().end(); ++colMapIt)
        {
          /// commented by Ruslan : at the moment a new branching constraint is added to the problem,
          /// columns are not yet reset, so they do not have the status they will have at the current node
//          if (!param().UseColumnsPool()
//              && (colMapIt->first->vcIndexStatus() != VcIndexStatus::Active)
//              && (colMapIt->first->vcIndexStatus() != VcIndexStatus::Inactive))
//            continue;
          includeMember(colMapIt->first, mapIt->second * colMapIt->second, true);
        }
    }

  Constraint::setMembership();

  return;
}

std::ostream & BasicConstrInstMastBranchingConstr::print(std::ostream& os) const
{
  os << "BasicConstrInstMastBranchingConstr" << std::endl;
  InstMasterBranchingConstr::print(os);
  return (os);
}

bool BasicConstrInstMastBranchingConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::BasicConstrInstMastBranchingConstrMask, vcIdentifier);
}

void MasterColSolution::push_back(Variable * varPtr, const ValueRecord & rec)
{
  std::list< std::pair < MastColumn *, ValueRecord > >::push_back(make_pair(static_cast<MastColumn *> (varPtr), rec));
  return;
}

/***************************************************************************
 **************   TODO: Methods for CustomNonLinearCut   *******************
 ***************************************************************************/

CustomNonLinearCut::CustomNonLinearCut(const IndexCell& id,
                                       GenericCustomNonLinearCutConstr * genConstrPtr,
                                       ProbConfig* probConfigPtr,
                                       const std::string & name,
                                       const BcCustomNonLinearCutInfo * cutInfoPtr):
  InstMasterConstr(id, genConstrPtr, probConfigPtr, name, genConstrPtr->defaultCostRhs(),
                   genConstrPtr->defaultSense(), genConstrPtr->defaultType(),
                   genConstrPtr->defaultKind(), genConstrPtr->defaultFlag()),
  Base4NonLinearConstraint(), _cutInfoPtr(cutInfoPtr),
  _genUserNonLinearCutConstrPtr(genConstrPtr)
{
}

CustomNonLinearCut::~CustomNonLinearCut()
{
  delete _cutInfoPtr;
}

void CustomNonLinearCut::setMembership()
{
    if(!buildMembershipHasBeenPerformed())
    {
        genVarConstrPtr()->buildMembership(this);
        buildMembershipHasBeenPerformed(true);
    }

    bool cumulativeCoef(false);
    
    VarIndexManager::const_iterator it;
    for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Active, 'd');
         it != problemPtr()->probVarSet().end(VcIndexStatus::Active, 'd'); ++it)
      if ((*it)->isTypeOf(VcId::MastColumnMask))
        {
          LpCoef lpCoeff = _genUserNonLinearCutConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }
    for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
         it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
      if ((*it)->isTypeOf(VcId::MastColumnMask))
        {
          LpCoef lpCoeff = _genUserNonLinearCutConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }
  /// if column pool is not used, we do not generate membership of the constraint
  /// in and unsuitable column, as generated column will be active
  /// only at the node it was generated and in the subtree rooted at this node
  if (param().UseColumnsPool())
    for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
         it != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); ++it)
      if ((*it)->isTypeOf(VcId::MastColumnMask))
        {
          LpCoef lpCoeff = _genUserNonLinearCutConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }
    
    Constraint::setMembership();

    return;
}

bool CustomNonLinearCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::CustomNonLinearCutConstrMask, vcIdentifier);
}

/***************************************************************************
 **************   TODO: Methods for SoftConlictCut   *******************
 ***************************************************************************/

SoftConflictsCut::SoftConflictsCut(const IndexCell& id,
								   GenericSoftConflictsCutConstr * genConstrPtr,
				                   ProbConfig * probConfigPtr,
				                   const std::string & name,
								   const Double & rhs,
								   const int cutType,
				                   const std::vector<std::pair<SubProbVariable *, SubProbVariable *> > & conflicts):
	InstMasterConstr(id, genConstrPtr, probConfigPtr, name, rhs, 'L', genConstrPtr->defaultType(),
		             genConstrPtr->defaultKind(), genConstrPtr->defaultFlag()),
	Base4NonLinearConstraint(), _cutType(cutType), _genSoftConflictConstrPtr(genConstrPtr), _conflicts(conflicts)
{
}

SoftConflictsCut::~SoftConflictsCut()
{
}

void SoftConflictsCut::setMembership()
{
	if(!buildMembershipHasBeenPerformed())
	{
		genVarConstrPtr()->buildMembership(this);
		buildMembershipHasBeenPerformed(true);
	}

	bool cumulativeCoef(false);

	VarIndexManager::const_iterator it;
	for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Active, 'd');
		 it != problemPtr()->probVarSet().end(VcIndexStatus::Active, 'd'); ++it)
		if ((*it)->isTypeOf(VcId::MastColumnMask))
		{
			LpCoef lpCoeff = _genSoftConflictConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
			if (lpCoeff.first)
				includeMember(*it, lpCoeff.second, cumulativeCoef);
		}
	for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
		 it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
		if ((*it)->isTypeOf(VcId::MastColumnMask))
		{
			LpCoef lpCoeff = _genSoftConflictConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
			if (lpCoeff.first)
				includeMember(*it, lpCoeff.second, cumulativeCoef);
		}
	/// if column pool is not used, we do not generate membership of the constraint
	/// in and unsuitable column, as generated column will be active
	/// only at the node it was generated and in the subtree rooted at this node
	if (param().UseColumnsPool())
		for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
			 it != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); ++it)
			if ((*it)->isTypeOf(VcId::MastColumnMask))
			{
				LpCoef lpCoeff = _genSoftConflictConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
				if (lpCoeff.first)
					includeMember(*it, lpCoeff.second, cumulativeCoef);
			}

	Constraint::setMembership();

	return;

}

bool SoftConflictsCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::SoftConflictsCutConstrMask, vcIdentifier);
}

const std::vector<std::pair<SubProbVariable *, SubProbVariable *> > & SoftConflictsCut::conflicts() const
{
  return _conflicts;
};
