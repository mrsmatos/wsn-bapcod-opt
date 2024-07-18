/**
 *
 * This file bcInstanciatedVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcModelC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcVarConstrC.hpp"
#include "bcStabilizationInfo.hpp"

/**
 * Generic code
 */
using namespace std;

/**
 * Methods of class InstanciatedVarConstr
 *
 */
InstanciatedVarConstr::InstanciatedVarConstr(const IndexCell & id, GenericVarConstr* genericVarConstrPtr,
                                             ProbConfig* probConfigPtr):
  _id(id),
  _genericVarConstrPtr(genericVarConstrPtr), 
  _probConfigPtr(probConfigPtr)
{
  if (printL(6))
    std::cout <<  "InstanciatedVarConstr::InstanciatedVarConstr()  id = " <<  this->id() << std::endl;
}


InstanciatedVarConstr::InstanciatedVarConstr(InstanciatedVarConstr * imcPtr, GenericVarConstr * genericVarConstrPtr):
  _id(imcPtr->id()),
  _genericVarConstrPtr(genericVarConstrPtr), 
  _probConfigPtr(imcPtr->_probConfigPtr)
{
  if (printL(6)) 
    std::cout <<  "InstanciatedVarConstr::InstanciatedVarConstr(InstanciatedVarConstr *)  id = " <<  id() << std::endl;
}

InstanciatedVarConstr::InstanciatedVarConstr(const InstanciatedVarConstr & that):
  _id(that.id()),
  _genericVarConstrPtr(that._genericVarConstrPtr), 
  _probConfigPtr(that._probConfigPtr)
{
  if (printL(6)) 
    std::cout <<  "InstanciatedVarConstr::InstanciatedVarConstr(copy) id = " <<  id() << std::endl;
}

const IndexCell & InstanciatedVarConstr::id() const 
{
  return(_id);
}

bool InstanciatedVarConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstanciatedVarConstrMask, vcIdentifier);
}

/**
 * Methods of class InstanciatedVar
 *
 */

/** 
 * 
 * 
 * @param id 
 * @param genericVarConstrPtr 
 * @param probConfigPt 
 * @param name 
 * @param costrhs 
 * @param sense 
 * @param type 
 * @param kind 
 * @param upperBound 
 * @param lowerBound 
 * @param flag 
 * @param directive 
 * @param priority // higher priority means chosen first for branching
 * @param val 
 * @param globalUb 
 * @param globalLb 
 * @param presetMembership 
 */
InstanciatedVar::InstanciatedVar(const IndexCell& id, 
                                 GenericVar * genericVarConstrPtr, 
                                 ProbConfig * probConfigPtr, 
                                 const std::string& name, 
                                 const Double& costrhs, 
                                 const char& sense, 
                                 const char& type, 
                                 const char& kind, 
                                 const Double& upperBound, 
                                 const Double& lowerBound, 
                                 const char& flag, 
                                 const char& directive, 
                                 const Double& priority, 
                                 const Double& val,
                                 const Double& globalUb, 
                                 const Double& globalLb,
				                 const bool & presetMembership):
  Variable(genericVarConstrPtr->modelPtr(), name, costrhs, sense, type, kind, upperBound, lowerBound, flag, directive,
	       priority, val, globalUb, globalLb, presetMembership, -1),
  InstanciatedVarConstr(id, genericVarConstrPtr, probConfigPtr),
  _genVarPtr(genericVarConstrPtr)
{
  if (printL(6))
    std::cout << "Instanciatedvar() " << name << " sense = " << sense << " lowerBound = " << lowerBound
	          << " upperBound = " << upperBound << " globalLb = " << globalLb << " globalUb = " << globalUb
	          << " presetMembership = " << VarConstr::presetMembership() << std::endl;

  if (probConfPtr() != NULL) 
      probConfPtr()->insertInstVar(this);

  if (_genVarPtr->priorityRule() == SelectionStrategy::NotConsideredForSelection)
    isCandForBranching(false);
  
  genericVarConstrPtr->recordInstanciation(this);
}

/** 
 * @param ivar
 */
InstanciatedVar::InstanciatedVar(const InstanciatedVar & ivar):
  Variable(ivar), InstanciatedVarConstr(ivar.id(), ivar.genVarConstrPtr(), ivar.probConfPtr()),
  _genVarPtr(ivar.genVarPtr())
{
  if (printL(6))
    std::cout <<  " : recordInstanciation/2 " << name() << " of type " <<  _genericVarConstrPtr << std::endl;
  
  _genericVarConstrPtr->recordInstanciation(this);
}

InstanciatedVar::~InstanciatedVar()
{
  GenericVar * gvPtr = dynamic_cast<GenericVar *>(_genericVarConstrPtr);
  bapcodInit().require(gvPtr != NULL, "InstanciatedVar::~InstanciatedVar(): genvar undefined");
  gvPtr->deleteInstanciation(this);
}

const Double InstanciatedVar::fracPartRelativeToTarget() const
{
  return Dfrac(tmpVal(), _genVarPtr->target4MostFractionalBranchingPriority());
}

const Double & InstanciatedVar::costrhs() const
{
  if (printL(6))
    std::cout << "this->name = " << name() << std::endl;
  
  return genVarConstrPtr()->genericCost(this);
}

void InstanciatedVar::costrhs(const Double & newCostRhs)
{
  if (printL(6))
    std::cout << "this->name = " << name() << std::endl;
  Variable::costrhs(newCostRhs);
  return;
}

bool InstanciatedVar::computeCount(ConstVarConstrConstPtr vcPtr) 
{
  if (printL(7))
    std::cout << "InstanciatedVar::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  if (presetMembership() && vcPtr->presetMembership() && vcPtr->membershipUpToDate()) 
    return (member2coefMap().count(vcPtr));

  if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask)) 
    {
      InstanciatedConstr * icPtr = static_cast<InstanciatedConstr * >(vcPtr);
      return (genVarConstrPtr()->genericCount(icPtr, this));
    }
    
  return false;
}


const LpCoef InstanciatedVar::computeCoef(ConstVarConstrConstPtr vcPtr) 
{
  if (presetMembership() && vcPtr->presetMembership() && vcPtr->membershipUpToDate())
    return upToDateMembCoef(vcPtr);

  if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask)) 
    {
      InstanciatedConstr * icPtr = static_cast<InstanciatedConstr * >(vcPtr);
      return (genVarConstrPtr()->genericCoef(icPtr, this));
    }
  
  return LpCoef::ZeroCoef;
}

void InstanciatedVar::setMembership()
{
  if (printL(6))
    std::cout <<  name() << std::endl;

  if (!buildMembershipHasBeenPerformed())
    {
      genVarConstrPtr()->buildMembership(this);
      buildMembershipHasBeenPerformed(true);
    }

  Variable::setMembership();

  return;
}

void InstanciatedVar::enumerativeSetMembership()
{
  Variable::enumerativeSetMembership();

  return;
}

bool InstanciatedVar::consecutive2varWhenBrOnCBS(InstanciatedVar * varPtr)
{
  return _genericVarConstrPtr->consecutiveVarWhenBrOnCBS(varPtr,this);
}

GenericVarConstr * InstanciatedVar::genVarConstrPtr() const 
{
  return(_genericVarConstrPtr);
}

void InstanciatedVar::createStabInfo(const BcObjStatus::MinMaxIntFloat & minmax)
{
  /// we do not stabilize implicit, type S varaints
  if ((kind() == 'I') || (type() == 'S'))
    return;

  _stabInfoPtr = new VarConstrStabInfo(this);
}

std::ostream& InstanciatedVar::print(std::ostream& os) const
{
  os << "InstanciatedVar" << std::endl;
  os << "   id = " << _id << std::endl;
  if (genVarConstrPtr() != NULL) 
    os << "   genericVarConstr = " << (long)genVarConstrPtr() << " " << genVarConstrPtr();
  
  if (probConfPtr() != NULL) 
    os << "   probConfig name = " << probConfPtr()->name() <<  std::endl;
  
  Variable::print(os);
  
  return(os);
}

const std::string & InstanciatedVar::genericName()
{
  return genVarPtr()->defaultName();
}
 
bool InstanciatedVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstanciatedVarMask, vcIdentifier);
}

void InstanciatedVar::setArcMembership(const int arcId, const double value)
{
  /// for the moment, there a variables may be associated to at most one arc,
  /// TO DO : modify this function with additional cumulative argument which says whether
  /// the value should be cumulated with already existing one or no
  _arcIdToCoeff[arcId] = value;
}


/**
 * Methods of class InstanciatedConstr
 */

InstanciatedConstr::InstanciatedConstr(const IndexCell& id,
                                       GenericConstr* genericVarConstrPtr,
                                       ProbConfig* probConfigPt,
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
                                       const Double & priority,
                                       const bool & presetMembership,
                                       const bool & toBeUsedInPreprocessing,
                                       const bool & considerAsEqualityInPreprocessing) :
  Constraint(genericVarConstrPtr->modelPtr(), name, costrhs, sense,
             type, kind, flag, -1, val, upperBound, lowerBound, directive, priority,
             presetMembership, toBeUsedInPreprocessing, considerAsEqualityInPreprocessing),
  InstanciatedVarConstr(id, genericVarConstrPtr, probConfigPt), _genConstrPtr(genericVarConstrPtr)
{
  if (printL(6))
    std::cout << "InstanciatedConstr::InstanciatedConstr(...) " << name
              << " presetMembership = " << VarConstr::presetMembership() << std::endl;

  _genericVarConstrPtr->recordInstanciation(this);
  if ((probConfPtr() != NULL) && (flag == 's'))
      probConfPtr()->insertInstConstr(this);
}

InstanciatedConstr::InstanciatedConstr(InstanciatedConstr * imcPtr,  
				       GenericConstr * genConstrPtr,
				       const std::string & name, 
				       const Double & rhs, 
				       const char & sense,
				       const char & type, 
				       const char & kind, 
				       const char & flag):
  Constraint(*imcPtr),
  InstanciatedVarConstr(imcPtr, genConstrPtr),
  _genConstrPtr(genConstrPtr)
{
  VarConstr::name(name);
  VarConstr::costrhs(rhs);
  VarConstr::sense(sense);
  VarConstr::type(type);
  VarConstr::kind(kind);
  VarConstr::flag(flag);
  _genericVarConstrPtr->recordInstanciation(this);
}


InstanciatedConstr::InstanciatedConstr(const InstanciatedConstr & that):
  Constraint(that),
  InstanciatedVarConstr(that),
  _genConstrPtr(that.genConstrPtr())
{
  if (printL(6))
    std::cout << "InstanciatedConstr::InstanciatedConstr(that) " << name() 
	          << " presetMembership = " << VarConstr::presetMembership() << std::endl;
}

InstanciatedConstr::~InstanciatedConstr()
{
  GenericConstr * gcPtr = dynamic_cast<GenericConstr *>(_genericVarConstrPtr);
  bapcodInit().require(gcPtr != NULL, "InstanciatedVar::~InstanciatedVar(): genvar undefined");
  gcPtr->deleteInstanciation(this);
}

const std::string & InstanciatedConstr::genericName()
{
  return genConstrPtr()->defaultName();
}

GenericVarConstr * InstanciatedConstr::genVarConstrPtr() const
{
  return(_genericVarConstrPtr);
}

void InstanciatedConstr::setMembership()
{
  if (!buildMembershipHasBeenPerformed())
    {
      genVarConstrPtr()->buildMembership(this);
      buildMembershipHasBeenPerformed(true);
    }

  Constraint::setMembership();

  return;
}

void InstanciatedConstr::enumerativeSetMembership()
{
  Constraint::enumerativeSetMembership();

  return;
}

const Double & InstanciatedConstr::costrhs() const
{
  if (printL(8))
    std::cout << "this->name = " << name() << std::endl;
  
  return(genVarConstrPtr()->genericRhs(this));
}


bool InstanciatedConstr::computeCount(ConstVarConstrConstPtr vcPtr) 
{
  if (printL(7))
    std::cout << "InstanciatedConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  bool countStatus(false);

  if (vcPtr->isTypeOf(VcId::AggregateVariableMask))
    { 
      AggregateVariable * agvPtr = dynamic_cast<AggregateVariable * >(vcPtr);
      return vcPtr->computeCount(this);
    }

  if (membershipUpToDate() && vcPtr->membershipUpToDate()) 
    {
      if (member2coefMap().count(vcPtr)) 
	    countStatus = true;
      
      if (printL(7)) 
	    std::cout << "InstanciatedConstr::computeCount both are already set, countStatus = " << countStatus
                  << std::endl;
    }
  else
    {
      if (vcPtr->isTypeOf(VcId::InstanciatedVarMask))
	    {
	  InstanciatedVar * ivPtr = static_cast<InstanciatedVar * >(vcPtr);
	  bapcodInit().check(genVarConstrPtr() == NULL,
			             "InstanciatedConstr::count(): genericVarConstrPtr should be defined");
	  
	  if (genVarConstrPtr()->genericCount(this, ivPtr)) 
	    countStatus = true;
	}
      else
	{
	  if (Constraint::computeCount(vcPtr)) 
	    countStatus = true;
	}
    }

  return  countStatus;
}

const LpCoef InstanciatedConstr::computeCoef(ConstVarConstrConstPtr vcPtr) 
{
  if (vcPtr->isTypeOf(VcId::AggregateVariableMask))
    {
      return vcPtr->computeCoef(this);
    }

  if (membershipUpToDate() && vcPtr->membershipUpToDate()) 
    {
      if (printL(7)) 
	    std::cout << "InstanciatedConstr::computeCoef both are already set "  << std::endl;
      
      return upToDateMembCoef(vcPtr);
    }
  if (vcPtr->isTypeOf(VcId::InstanciatedVarMask))
    {
      InstanciatedVar * ivPtr = static_cast<InstanciatedVar * >(vcPtr);
      bapcodInit().check(genVarConstrPtr() == NULL,
			             "InstanciatedConstr::count(): _genericVarConstrPtr should be defined");
      
      return genVarConstrPtr()->genericCoef(this, ivPtr);
    }
  return Constraint::computeCoef(vcPtr);
}

std::ostream& InstanciatedConstr::print(std::ostream& os) const
{
  os << "InstanciatedConstr" << std::endl;
  os << "   id = " << _id << std::endl;
  if (genVarConstrPtr() != NULL) 
    os << "   genericVarConstr = " <<  genVarConstrPtr();
  
  if (probConfPtr() != NULL) 
    os << "   probConfig name = " << probConfPtr()->name() <<  std::endl;
  
  Constraint::print(os);
  
  return(os);
}

void InstanciatedConstr::nicePrint(std::ostream& os) const
{
  os << "Constraint " << name() << " :";
  for (ConstVarConstrPtr2Double::const_iterator mapIt = _member2coefMap.begin(); mapIt != _member2coefMap.end();
       ++mapIt)
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

bool InstanciatedConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::InstanciatedConstrMask, vcIdentifier);
}

/**
 * Methods of class NonLinearInstConstr
 *
 */

NonLinearInstConstr::NonLinearInstConstr(const IndexCell& id, 
					 GenericConstr * genConstrPtr, 
					 ProbConfig * probConfigPt, 
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
					 const Double & priority):
  InstanciatedConstr(id, 
		     genConstrPtr,  
		     probConfigPt, 
		     name, 
		     costrhs, 
		     sense, 
		     type, 
		     kind, 
		     flag, 
		     val, 
		     upperBound, 
		     lowerBound, 
		     directive, 
		     priority,
		     false, //presetMembership
		     false //toBeUsedInPreprocessing
		     )
{
  if (printL(5))
    std::cout << "NonLinearInstConstr::NonLinearInstConstr(...) construct" << name << std::endl;
  
  if (probConfPtr() != NULL)
    {
        DynamicGenericConstr * dgcPtr = dynamic_cast<DynamicGenericConstr * > (genVarConstrPtr());
        if (printL(5))
            std::cout << "NonLinearInstConstr::NonLinearInstConstr(...) DynamicGenericConstr ? "
                      << (dgcPtr != NULL) << std::endl;

        if (dgcPtr != NULL)
        {
            /// Ruslan : seems to be not used
            /// anyway, we do not insert into constr prototypes automatically,
            /// as a constraint may by a prototype or may be a "real" one
//          MasterConf * mastConfPtr = dynamic_cast<MasterConf * > (probConfigPt);
//        if (mastConfPtr != NULL)
//	        dgcPtr->insertConstrPrototypes(new NonLinearInstMastConstr(this));
//	      else
//	        dgcPtr->insertConstrPrototypes(this);
        }
      else
        {
	  if (printL(5)) 
	    std::cout << "NonLinearInstConstr::NonLinearInstConstr(...) insert NonLinearInstConstr in probConf" 
		      << std::endl;
	  
          probConfPtr()->insertInstConstr(this);
        }
    }
}

bool NonLinearInstConstr::computeCount(ConstVarConstrConstPtr vcPtr) 
{
  if (printL(6))
    std::cout << "NonLinearInstConstr::computeCount this " << name() << " that "  << vcPtr->name() << std::endl;

  bapcodInit().check(genVarConstrPtr() == NULL,
	                 "NonLinearInstConstr::count(): genericVarConstrPtr should be defined");

  InstanciatedVar * ivPtr = dynamic_cast<InstanciatedVar * >(vcPtr);
  if (ivPtr != NULL)
    return(genVarConstrPtr()->genericCount(this, ivPtr));
    
  return(Constraint::computeCount(vcPtr));
}

const LpCoef NonLinearInstConstr::computeCoef(ConstVarConstrConstPtr vcPtr) 
{
  return LpCoef::ZeroCoef;
}

std::ostream& NonLinearInstConstr::print(std::ostream& os) const
{
  os << "NonLinearInstConstr" << std::endl;
  InstanciatedConstr::print(os);
  
  return(os);
}

bool NonLinearInstConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::NonLinearInstConstrMask, vcIdentifier);
}

/**
 * Methods of class BranchingConstrBaseType
 *
 */

void BranchingConstrBaseType::append2name(const std::string & ad)
{
  Constraint * constrPtr = dynamic_cast<Constraint * >(this);
  if (constrPtr)
    constrPtr->name(constrPtr->name() + ad);
    
  return;
}

std::ostream & BranchingConstrBaseType::print(std::ostream& os) const
{
  os << "BranchingConstrBaseType" << std::endl;

   for (std::set<ProbConfig *>::const_iterator it = _probConfigPtrSet.begin(); it != _probConfigPtrSet.end();
       ++it) 
     os << "    probConfig name = " << (*it)->name() << std::endl;

  return(os);
}

bool BranchingConstrBaseType::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::BranchingConstrBaseTypeMask, vcIdentifier);
}
