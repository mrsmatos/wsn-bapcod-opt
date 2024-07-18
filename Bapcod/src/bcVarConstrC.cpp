/**
 *
 * This file bcVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcFormC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcOvfVarConstrC.hpp"

#include "bcPrintC.hpp"
#include "bcProblemC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcStabilizationInfo.hpp"
#include "bcGlobalException.hpp"

using namespace std;
using namespace VcIndexStatus;

bool VarConstrSort::operator()(const VarConstr * a, const VarConstr * b) const
{
  return (a->VCref() < b->VCref());
}

/**
 * Methods of class VarConstr
 */


VarConstr::VarConstr(const VarConstr & vc, const int & ref) :
  _VCref(ref),
  _name(""),
  _directive(vc.directive()),
  _priority(vc.priority()),
  _modelPtr(vc.modelPtr()),
  _costrhs(vc.costrhs()),
  _sense(vc.sense()),
  _sign(vc.sign()),
  _type(vc.type()),
  _kind(vc.kind()),
  _flag(vc.flag()),
  _inCurProb(false),
  _inCurForm(false),
  _index(-1),
  _multiIndex(),
  _val(vc.val()),
 _tmpVal(vc.tmpVal()),
 _solVal(vc.solVal()),
 _incumbentVal(vc.incumbentVal()),
 _upperBound(vc.ub()),
 _lowerBound(vc.lb()),
 _globalUb(vc.globalUb()),
 _globalLb(vc.globalLb()),
 _memorisedCurUb(vc.ub()),
 _memorisedCurLb(vc.lb()),
 _curUb(vc.curUb()),
 _curLb(vc.curLb()),
 _mult(1),
 _presetMembership(vc.presetMembership()),
_membershipUpToDate(vc.membershipUpToDate()),
 _buildMembershipHasBeenPerformed(false),
_member2coefMap(vc.member2coefMap()),
 _nonMemberSet(vc.nonMemberSet()),
 _problemPtr(vc.problemPtr()),
 _participation(0),
 _violation(vc.violation()),
 _reducedCost(vc.reducedCost()),
 _ovfVarList(),
 _ovfConstrList(),
_stabInfoPtr(NULL)

{
  _modelPtr->garbageCollector().insert(this);

  if (vc._stabInfoPtr != NULL)
    _stabInfoPtr = new VarConstrStabInfo(*(vc._stabInfoPtr), this);

  if (printL(6))
    std::cout << "VarConstr::VarConstr(COPY) "  << name()
	          << " lowerBound = " << _lowerBound
	          << " upperBound = " << _upperBound
	          << " globalLb = " << _globalLb
	          << " globalUb = " << _globalUb
	          << " presetMembership = " << VarConstr::presetMembership()
	          << std::endl;

  if (vc.name() != "")
    {
      name(vc.name() + "_copy");
    }
}

VarConstr::VarConstr(Model * modelPtr,
		             const long & ref,
		             const std::string & cname,
		             const Double & costrhs,
		             const char & sense,
		             const char & type,
		             const char & kind,
		             const char & flag,
		             const int & index,
		             const Double & val,
		             const Double & upperBound,
		             const Double & lowerBound,
		             const Double & globalUb,
		             const Double & globalLb,
		             const char & directive,
		             const Double & priority,
		             const bool & presetMembership) :
		 _VCref(ref),
		 _name(cname),
		 _directive(directive),
		 _priority(priority),
		 _modelPtr(modelPtr),
		 _costrhs(costrhs),
		 _sense(sense),
		 _sign((sense == 'L') ? (((modelPtr->objectiveSense() == BcObjStatus:: minInt) ||
					  (modelPtr->objectiveSense() == BcObjStatus:: minFloat))? -1: 1)
		       : (((modelPtr->objectiveSense() == BcObjStatus:: minInt) ||
			   (modelPtr->objectiveSense() == BcObjStatus:: minFloat))? 1: -1)),
		 _type(type),
		 _kind(kind),
		 _flag(flag),
		 _inCurProb(false),
		 _inCurForm(false),
		 _index(index),
		 _multiIndex(),
		 _val(val),
		 _tmpVal(0),
		 _challengerRoundedValue(0),
		 _solVal(0),
		 _incumbentVal(val),
		 _upperBound(upperBound),
		 _lowerBound(lowerBound),
		 _globalUb(globalUb),
		 _globalLb(globalLb),
		 _memorisedCurUb(upperBound),
		 _memorisedCurLb(lowerBound),
		 _curUb(upperBound),
		 _curLb(lowerBound),
		 _mult(1),
		 _presetMembership(presetMembership),
		 _membershipUpToDate(false),
		 _buildMembershipHasBeenPerformed(false),
		 _member2coefMap(),
		 _nonMemberSet(),
		 _problemPtr(NULL),
         _infoIsUpdated(false),
         _inPreprocessedList(false),
         _participation(0),
		 _violation(0),
		 _reducedCost(0),
		 _ovfVarList(),
		 _ovfConstrList(),
		 _stabInfoPtr(NULL)

{
  _modelPtr->garbageCollector().insert(this);

  if (printL(6))
    std::cout << "VarConstr::VarConstr() name = " << name()
	      << " cname = " << cname
	      << " lowerBound = " << lowerBound
	      << " upperBound = " << upperBound
	      << " globalLb = " << globalLb
	      << " globalUb = " << globalUb
	      << " presetMembership = " << VarConstr::presetMembership()
	      << std::endl;

  return;
}

VarConstr::~VarConstr()
{
  clearMembership();
  _modelPtr->garbageCollector().erase(this);
  if (_stabInfoPtr != NULL)
    delete _stabInfoPtr;

  return;
}

BapcodInit & VarConstr::bapcodInit() const
{
  return _modelPtr->bapcodInit();
}

const ControlParameters & VarConstr::param() const
{
  return _modelPtr->param();
}
ControlParameters & VarConstr::param()
{
  return _modelPtr->param();
}


const Double & VarConstr::costrhs() const
{
  return _costrhs;
}

const Double VarConstr::cost() const
{
  return _costrhs;
}

const Double & VarConstr::rhs() const
{
  return _costrhs;
}

void VarConstr::resetCostFromDefaultCost(const Double & factor)
{
}

  /*
   * Variables have participation incremented in the following methods
   *
   * 0
   * 1  SolutionVarInfo::SolutionVarInfo
   * 2  ProblemSetupInfo::ProblemSetupInfo
   * 3  ProblemFullSetDownAlgorithm::recordProblemInfo
   * 4  Solution::Solution
   * 6  Solution::includeVarSet
   * 8  Solution::includeVar
   * 10 DiveInfo::tabuVariables
   * 11 DiveAlgorithm::runStrongDive
   * 12 DiveAlgorithm::runStrongDive
   * 15 Problem::updatePartialSolution
   * ...
   *
   */

  /*
   * Variables have participation decremented in the following methods
   *
   * 1
   * 2  SolutionVarInfo::~SolutionVarInfo()
   * 3  ProblemSetupInfo::~ProblemSetupInfo
   * 4  Solution::~Solution
   * 5  DiveInfo::~DiveInfo
   * 6  DiveAlgorithm::runStrongDive
   * 9  Problem::resetPartialSolution
   *
   */

  /*
   * Constraints have participation incremented in the following methods
   *
   * 0  Node::Node
   * 1  ProblemFullSetDownAlgorithm::recordProblemInfo
   * 2  BcExtendedArcCut::reserveCut
   * 3  StabilizationInfo::StabilizationInfo
   * 4  ColGenStabilization::addConstrAndAssociatedArtVarsToStabCandList
   * 5  ColGenStabilization::recordStabilizationInfo
   * 6  BranchingConstrGenerator::BranchingConstrGenerator
   * 7  BranchingConstrGenerator::BranchingConstrGenerator
   * ...
   *
   */

  /*
   * Constraints have participation decremented in the following methods
   *
   * 0  Node::clearLocalNodeBrConstrList
   * 1  ProblemSetupInfo::~ProblemSetupInfo
   * 2  BcExtendedArcCut::releaseCut
   * 3  StabilizationInfo::~StabilizationInfo
   * 4  ColGenStabilization::deactivate
   * 5  BranchingConstrGenerator::~BranchingConstrGenerator
   *
   */


void VarConstr::incrParticipation(int where)
{
    _participation++;
}

void VarConstr::decrParticipation(int where)
{
    _participation--;
}

const Double & VarConstr::violation() const
{
  return _violation;
}

void VarConstr::violation(const Double & viol)
{
  _violation = viol;
  _violation.Czero();
}

void VarConstr::costrhs(const Double & newCostrhs)
{
  _costrhs = newCostrhs;
  return;
}

void VarConstr::addToProb(Problem * probPtr)
{
  problemPtr(probPtr);
  setMembership();
}

void VarConstr::deleteFromProb()
{
  clearMembership();
  problemPtr(NULL);
}

const std::string & VarConstr::name() const
{
    return _name;
}

void VarConstr::name(const std::string & cname)
{
    _name = cname;
    return;
}

void VarConstr::name(char * cname)
{
    _name = "";
    _name.append(cname);
    return;
}

const char & VarConstr::sense() const
{
  return (_sense);
}

const int & VarConstr::sign() const
{
  return (_sign);
}

const char & VarConstr::sense(const char & sense)
{
  return (_sense = sense);
}

const char & VarConstr::type() const
{
  return (_type);
}

const char & VarConstr::kind() const
{
  return (_kind);
}

const char & VarConstr::type(const char & s)
{
  return (_type = s);
}

const char & VarConstr::kind(const char & s)
{
  return (_kind = s);
}

const char & VarConstr::flag() const
{
  return (_flag);
}

void VarConstr::flag(const char & cflag)
{
  _flag = cflag;
  return;
}

const int & VarConstr::index() const
{
  return (_index);
}

const int & VarConstr::index(const Double & newIndex)
{
  _index = newIndex;
  return (_index);
}

const Double VarConstr::curRhs() const
{
  return costrhs();
}

const bool VarConstr::inCurProb()
{
  return (_inCurProb);
}

const bool & VarConstr::inCurForm() const
{
  return (_inCurForm);
}

void VarConstr::activate()
{
  _inCurProb = true;
  return;
}

void VarConstr::setInForm()
{
  _inCurForm = true;
  return;
}

void VarConstr::desactivate()
{
  _inCurProb = false;
  _val = 0;

  return;
}

void VarConstr::unsetInForm()
{
  if (printL(7))
    std::cout << " VarConstr::unsetInForm() " << name() << std::endl;

  _inCurForm = false;

  return;
}

const Double & VarConstr::val() const
{
  return (_val);
}

void VarConstr::val(const Double & newVal)
{
  _val = newVal;

  return;
}

const Double & VarConstr::tmpVal() const
{
  return (_tmpVal);
}

void VarConstr::tmpVal(const Double & newVal)
{
  _tmpVal = newVal;

  return;
}

const Double & VarConstr::solVal() const
{
  return (_solVal);
}

void VarConstr::solVal(const Double & newVal)
{
  _solVal = newVal;
  return;
}

const Double VarConstr::fracPart() const
{
  return (Dfrac(tmpVal()));
}


const Double VarConstr::fracPartRelativeToTarget() const
{
  return (Dfrac(tmpVal()));
}

const Double VarConstr::lFracPart() const
{
  return (Lfrac(tmpVal()));
}

const Double VarConstr::uFracPart() const
{
  return (Ufrac(tmpVal()));
}

const Double & VarConstr::incumbentVal() const
{
  return (_incumbentVal);
}

void VarConstr::incumbentVal(const Double & newVal)
{
  _incumbentVal = newVal;

  return;
}

const Double & VarConstr::ub() const
{
  return (_upperBound);
}

void VarConstr::ub(const Double & ub)
{
  _curUb = _memorisedCurUb = _upperBound = ub;

  return;
}
const Double & VarConstr::lb() const
{
  return _lowerBound;
}

void VarConstr::lb(const Double & lb)
{
  _curLb = _memorisedCurLb = _lowerBound = lb;
  return;
}

const Double & VarConstr::globalUb() const
{
  return _globalUb;
}

void VarConstr::globalUb(const Double & ub)
{
  _globalUb = ub;
  return;
}

const Double & VarConstr::globalLb() const
{
  return _globalLb;
}

void VarConstr::globalLb(const Double & lb)
{
  _globalLb = lb;
  return;
}

const Double & VarConstr::mult() const
{
  return _mult;
}

void VarConstr::mult(const Double & m)
{
  _mult = m;
  return;
}

ConstVarConstrPtr2Double & VarConstr::member2coefMap()
{
  return _member2coefMap;
}

const ConstVarConstrPtr2Double & VarConstr::member2coefMap()  const
{
  return _member2coefMap;
}

void VarConstr::addMember(VarConstr * vcPtr)
{
  /// Needed, because not done externally ; it does the includeMember(vcPtr) if necessary
  membCount(vcPtr);
  return;
}

const Double & VarConstr::includeMember(VarConstr * vcPtr, const Double & coef, const bool & cumulativeCoef)
{
  if (printL(6))
    std::cout << "VarConstr::includeMember this =  " << name() << ", that = " << vcPtr->name()
	          << ", coef = " << coef << "  cumulativeCoef " <<  cumulativeCoef << std::endl;

  vcPtr->includeAsMember(this, coef, cumulativeCoef);

  return includeAsMember(vcPtr, coef, cumulativeCoef);
}

const Double & VarConstr::includeAsMember(VarConstr * vcPtr, const Double & coef, const bool & cumulativeCoef)
{
  if (printL(7))
    std::cout << "VarConstr::includeAsMember this =  " << name() << ", that = " << vcPtr->name()
              << ", coef = " << coef << std::endl;

    if (cumulativeCoef)
    {
        ConstVarConstrPtr2Double::iterator position = _member2coefMap.find(vcPtr);

        /// Comparison not already done
        if (position != _member2coefMap.end())
        {
            position->second += coef;
            return position->second;
        }
        else
        {
            _member2coefMap[vcPtr] = coef;
        }
    }
    else
    {
        _member2coefMap[vcPtr] = coef;
    }
    return coef;
}

void VarConstr::eraseAsMember(VarConstr * vcPtr)
{
  _member2coefMap.erase(vcPtr);
}

/**
 * Do behind the schene membership for MasterConstr and SubProblemVariables
 */
void VarConstr::setMembership()
{
  _membershipUpToDate = true;

  return;
}

void VarConstr::clearMembership()
{
  if (printL(6))
    std::cout << "VarConstr::clearMembership() this =  " << name() << std::endl;

  for (ConstVarConstrPtr2Double::iterator itm = _member2coefMap.begin(); itm != _member2coefMap.end(); ++itm)
    {
      if (printL(6))
        std::cout << "VarConstr::clearMembership() that =  " << itm->first->name() << std::endl;

      itm->first->eraseAsMember(this);
    }
  _member2coefMap.clear();
  _nonMemberSet.clear();
  _membershipUpToDate = false;

  return;
}


bool VarConstr::operator<(const VarConstr & that) const
{
  return VCref() < that.VCref();
}

bool VarConstr::operator==(const VarConstr & that) const
{
  if (this->operator<(that))
      return false;
  if (that.operator<(*this))
      return false;
  return true;
}

bool VarConstr::operator!=(const VarConstr & that) const
{
  return !this->operator==(that);
}

/**
 * Measures the variable min contribution
 * to the satisfaction of a specific constraint
 *
 * @param vcPtr
 *
 * @return
 */
const Double VarConstr::lhsMinContrib(ConstVarConstrConstPtr vcPtr)
{
  Double coef(vcPtr->membCoef(this));

  if (coef > 0)
    return (curLb() * coef);
  else
    return (curUb() * coef);
}

/**
 * Measures the variable max contribution
 * to the satisfaction of a specific  constraint
 */
const Double VarConstr::lhsMaxContrib(ConstVarConstrConstPtr vcPtr)
{
  Double coef(vcPtr->membCoef(this));

  if (coef > 0)
    return (curUb() * coef);
  else
    return (curLb() * coef);
}

std::ostream & VarConstr::print(std::ostream& os) const
{
  os << "   ref = " << _VCref << std::endl;
  os << "   name = " << name() << std::endl;
  os << "   costrhs = " << costrhs() << std::endl;
  os << "   redCost = " << reducedCost() << std::endl;
  os << "   sense = " << _sense << std::endl;
  os << "   sign = " << _sign << std::endl;
  os << "   type = " << _type << std::endl;
  os << "   kind = " << _kind << std::endl;
  os << "   directive = " << _directive << std::endl;
  os << "   priority = " << _priority << std::endl;
  os << "   flag = " << _flag << std::endl;
  os << "   active = " << _inCurProb << std::endl;
  os << "   inCurForm = " << _inCurForm << std::endl;
  os << "   index = " << _index << std::endl;
  os << "   val = " << _val << std::endl;
  os << "   incumbentVal = " << _incumbentVal << std::endl;
  os << "   globalLowerBound = " << _globalLb << std::endl;
  os << "   lowerBound = " << _lowerBound << std::endl;
  os << "   upperBound = " << _upperBound << std::endl;
  os << "   globalUpperBound = " << _globalUb << std::endl;
  os << "   mult = " << _mult << std::endl;

  if (_member2coefMap.empty())
	  os << "     whose membership is not defined " << std::endl;
  else
  {
	  os << "     whose membership is: " << std::endl;

	  std::vector<std::pair<VCPTR, Double> > temp;
	  for (ConstVarConstrPtr2Double::const_iterator itm = _member2coefMap.begin(); itm != _member2coefMap.end(); ++itm)
	  {
		  temp.push_back(std::pair<VCPTR, Double>(itm->first,itm->second));
	  }

	  for (std::vector<std::pair<VCPTR, Double> >::const_iterator itm = temp.begin(); itm != temp.end(); ++itm)
		  os << "   coef[" << itm->first->name() << "] = " << itm->second << std::endl;
  }
  return (os);
}

const char & VarConstr::directive() const
{
  return (_directive);
}

void VarConstr::directive(const char & dir)
{
  _directive = dir;

  return;
}

const Double & VarConstr::priority() const
{
  return (_priority);
}

void VarConstr::priority(const Double & prior)
{
  _priority = prior;

  return;
}

void VarConstr::recordNonMember(ConstVarConstrConstPtr vcPtr)
{
}

bool VarConstr::membCount(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr == NULL)
    return false;

  if (printL(7))
    std::cout << "VarConstr::membCount this =  " << name() << ", that = " << vcPtr->name() << std::endl;

  if (membershipUpToDate() && vcPtr->membershipUpToDate())
    {
      if (printL(7))
        std::cout << "membershipUpToDate vc " << vcPtr->name() << std::endl;

      return (member2coefMap().count(vcPtr));
    }
  /// else need to recompute count

  /// Already computed before
  if (member2coefMap().count(vcPtr))
    {
      if (printL(7))
        std::cout << "count vc " << vcPtr->name() << std::endl;

      return true;
    }

  /// Already excluded before
  if (nonMemberSet().count(vcPtr))
    {
      if (printL(7))
        std::cout << "nonMember vc " << vcPtr->name() << std::endl;

      return false;
    }

  LpCoef coef(computeCoef(vcPtr));
  if (coef.first)
    {
      if (printL(7))
        std::cout << "compute vc " << vcPtr->name() << std::endl;

      includeMember(vcPtr, coef.second, false);
      return true;
    }
  else
    {
      if (printL(7))
        std::cout << "recordNonMember vc " << vcPtr->name() << std::endl;

      recordNonMember(vcPtr);
    }

  return false;
}

const Double & VarConstr::membCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr == NULL)
    return Double::staticZero;

  if (printL(7))
    std::cout << "VarConstr::membCoef this =  " << name() << ", that = "
        << vcPtr->name() << std::endl;

  if (membershipUpToDate() && vcPtr->membershipUpToDate())
    {
      if (printL(7))
        std::cout << "membershipUpToDate vc " << vcPtr->name() << std::endl;

      return upToDateMembCoef(vcPtr);
    }

  /// Already computed before
  ConstVarConstrPtr2Double::const_iterator constrMemberIt = member2coefMap().find(vcPtr);
  if (constrMemberIt != member2coefMap().end())
    {
      if (printL(7))
        std::cout << "count vc " << vcPtr->name() << std::endl;

      return constrMemberIt->second;
    }

  /// Already excluded before
  if (nonMemberSet().count(vcPtr))
    {
      if (printL(7))
        std::cout << "nonMember vc " << vcPtr->name() << std::endl;

      return Double::staticZero;
    }

 /// Do not know if already computed before
  LpCoef coef(computeCoef(vcPtr));
  if (coef.first)
    {
      if (printL(7))
        std::cout << "compute vc " << vcPtr->name() << std::endl;

      return includeMember(vcPtr, coef.second, false);
      //return true;
    }
  else
    {
      if (printL(7))
        std::cout << "recordNonMember vc " << vcPtr->name() << std::endl;

      recordNonMember(vcPtr);
      return Double::staticZero;
    }
  return Double::staticZero;
}

const Double & VarConstr::upToDateMembCoef(ConstVarConstrConstPtr vcPtr) const
{
    ConstVarConstrPtr2Double::const_iterator cur = _member2coefMap.find(vcPtr);

    if (cur != _member2coefMap.end())
        return cur->second;
    return Double::staticZero;
}

const Double & VarConstr::reducedCost() const
{
  return (_reducedCost);
}

void VarConstr::reducedCost(const Double & newRedCost)
{
  _reducedCost = newRedCost;
}

const Double & VarConstr::computeReducedCost()
{
  int printlevel = 6;

  long long int scaleFactor = param().SafeDualBoundScaleFactor();

  if (scaleFactor > 0)
      _reducedCost = floor(cost()._val * scaleFactor);
  else
      _reducedCost = cost();

  bapcodInit().require(_problemPtr != NULL, "VarConstr::computeReducedCost(): problemPtr == NULL");

  if (param().useSPVarMembershipCache())
  {
      if (_cachedMemberShip.empty())
      {
          if (printL(printlevel))
          {
              std::cout << "Formulation is assumed static. Computing cached membership of " << "Var[" << name() << "] "
                        << std::endl;
          }
          for (ConstVarConstrPtr2Double::iterator itm = member2coefMap().begin();
               itm != member2coefMap().end(); ++itm)
          {
              _cachedMemberShip.push_back({itm->first, itm->first->membCoef(this)});
              if (printL(printlevel))
              {
                  std::cout << "Var[" << name() << "] in const[" << itm->first->name()
                            << "] of val[" << itm->first->val()
                            << "] has coef[" << itm->first->membCoef(this) << "]" << std::endl;
              }
          }

      }
      for (const auto &memb : _cachedMemberShip)
      {
          if (memb.first->inCurForm())
          {
              if (scaleFactor > 0)
                  _reducedCost += ceil(memb.first->val()._val * memb.second * scaleFactor);
              else
                  _reducedCost += memb.first->val() * memb.second;
              if (printL(printlevel))
              {
                  std::cout << "Var[" << name() << "] in const["
                            << memb.first->name() << "] of val[" << memb.first->val()
                            << "] has coef[" << memb.second << "]  rc= "
                            << _reducedCost << std::endl;
              }
          }
      }
  }
  else
  {
      if (member2coefMap().size() <= _problemPtr->inDualSol().size())
      {
          for (ConstVarConstrPtr2Double::iterator itm = member2coefMap().begin();
               itm != member2coefMap().end(); ++itm)

              if (itm->first->inCurForm())
              {
                  if (scaleFactor > 0)
                      _reducedCost += ceil(itm->first->val()._val * itm->first->membCoef(this)._val * scaleFactor);
                  else
                      _reducedCost += itm->first->val() * itm->first->membCoef(this);
                  if (printL(printlevel))
                  {
                      std::cout << "Var[" << name() << "] in const["
                                << itm->first->name() << "] of val[" << itm->first->val()
                                << "] has coef[" << itm->first->membCoef(this) << "]  rc= "
                                << _reducedCost << " second=" << itm->second << std::endl;
                  }
              }
      }
      else
      {
          for (ConstrPtrSet::const_iterator cPt = _problemPtr->inDualSol().begin();
               cPt != _problemPtr->inDualSol().end(); cPt++)
              if (member2coefMap().count(*cPt))
              {
                  if (scaleFactor > 0)
                      _reducedCost += ceil((*cPt)->val()._val * (*cPt)->membCoef(this)._val * scaleFactor);
                  else
                      _reducedCost += (*cPt)->val() * (*cPt)->membCoef(this);
                  if (printL(printlevel))
                  {
                      std::cout << "Var[" << name() << "] in const[";
                      // bug: method VarConstr::name() cannot be accessed twice in the same line since all calls share the
                      // same temporary static string object. The same string will appear for all calls
                      std::cout << (*cPt)->name() << "] of val[" << (*cPt)->val() << "] has coef["
                                << (*cPt)->membCoef(this) << "]  rc= " << _reducedCost << std::endl;
                  }
              }
      }
  }

  if (scaleFactor > 0)
      _reducedCost._val /= scaleFactor;

  if (printL(printlevel))
    std::cout << "VarConstr::computeReducedCost()   rc[" << name() << "] = "
              << _reducedCost << " while cost is " << cost() << std::endl;

  /// added by Ruslan: reduced cost should be always synchronized with _NonZeroRedCostVars of the problem
  if (isTypeOf(VcId::VariableMask))
    problemPtr()->updateInNonZeroRedCostVarsSet(static_cast<Variable *>(this));

  return _reducedCost;
}

const Double VarConstr::greedyCost()
{
  Double contribution = contrib();
  if (printL(6))
    std::cout << "VarConstr::greedyCost[" << name() << "]: cost() = " << cost()
              << " / contrib() " << contribution << "; ratio = "
              << (cost() / contribution) << std::endl;

  if (contribution.positive())
    return ((cost() / contribution));
  else
    return (BapcodInfinity);
}

bool VarConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::VarConstrMask, vcIdentifier);
}

const Double & VarConstr::valOrSepPointVal() const
{
  if ((_stabInfoPtr != NULL) && (_stabInfoPtr->stabilizationParticipationFlag() == 2))
    return _stabInfoPtr->sepPointVal();
  else
    return _val;
}

void VarConstr::deleteStabInfoPtr()
  {
    if(_stabInfoPtr != NULL)
      delete _stabInfoPtr;
    _stabInfoPtr = NULL;
  }


/**
 * Methods of class Variable
 */

Variable::Variable(Model * modelPtr,
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
		           const bool & presetMembership,
		           const int & index) :
  VarConstr(modelPtr,
	        modelPtr->modelVarCnt(),
            name,
	        costrhs,
	        sense,
	        type,
	        kind,
	        flag,
            index,
	        val,
	        upperBound,
	        lowerBound,
	        globalUb,
	        globalLb,
	        directive,
            priority,
	        presetMembership),
  _memorisedCurCost(costrhs),
  _isCandForBranching((type == 'B') || (type == 'I'))
{
  _modelPtr->increaseModelVarCnt();
  if (_type == 'B')
    {
      if (_upperBound > 1)
	{
	  _upperBound = _memorisedCurUb = _curUb = 1;
	}
      if (_lowerBound < 0)
	{
	  _lowerBound = _globalLb = _memorisedCurLb = _curLb = 0;
	}
    }

  if (printL(6))
    std::cout << "Variable::Variable() " << name << " in [" << _lowerBound
        << ", " << _upperBound << "] " << std::endl;

}

Variable::Variable(const Variable & that) :
    VarConstr(that, that._modelPtr->modelVarCnt()),
    _memorisedCurCost(that._memorisedCurCost),
    _isCandForBranching(that._isCandForBranching)
{
    that._modelPtr->increaseModelVarCnt();
}

Variable::~Variable()
{
  return;
}

void Variable::addToPreprocessedList()
{
  if (!_inPreprocessedList)
    {
      _problemPtr->preprocessedVarsList().push_back(this);
      _inPreprocessedList = true;
    }
}

bool Variable::boundsAreDifferentFromDefaultInForm()
{
  return (_inCurForm && ((_lowerBound != _memorisedCurLb) || (_upperBound != _memorisedCurUb)));
}

void Variable::resetBoundsAndCostToDefaults()
{
  _curLb = _memorisedCurLb = _lowerBound;
  _curUb = _memorisedCurUb = _upperBound;
  _memorisedCurCost = _costrhs;
}

bool Variable::activateVariable(bool toSetInForm)
{
    if (_problemPtr == NULL)
        return false;

    _problemPtr->probVarSet().insert(this, Active);
    activate();
    if (toSetInForm)
        _problemPtr->setVar2Form(this);
    return true;
}

bool Variable::desactivateVariable(const VcIndexStatus::VcStatus & status, bool toUnsetInForm)
{
  if (_problemPtr == NULL)
    return false;

  _problemPtr->probVarSet().insert(this, status);
  desactivate();
  if (toUnsetInForm)
    _problemPtr->unsetVar2Form(this);
  return true;
}

bool Variable::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(7))
    std::cout << " Variable::computeCount this " << name() << " that "
	      << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::AggregateVariableMask))
    {
      AggregateVariable * agvPtr = dynamic_cast<AggregateVariable *>(this);
      return agvPtr->agvComputeCount(vcPtr);
    }

  return vcPtr->computeCount(this);
}

const LpCoef Variable::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (printL(7))
    std::cout << " Variable::computeCoef this " << name() << " that "
        << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::AggregateVariableMask))
    {
      AggregateVariable * agvPtr = dynamic_cast<AggregateVariable *>(this);
      return agvPtr->agvComputeCoef(vcPtr);
    }

  return vcPtr->computeCoef(this);
}

void Variable::addMember(VarConstr * vcPtr)
{
  VarConstr::addMember(vcPtr);

  return;
}

bool Variable::PrioritySort::operator()(const VarPointer & a, const VarPointer & b) const
{
  if (a->priority() > b->priority())
    return (true);

  if (a->priority() < b->priority())
    return (false);

  return (a->VCref() < b->VCref());
}

const Double Variable::cost() const
{
  return curCost();
}


const Double & Variable::curCost() const
{
    if (printL(6))
        std::cout << " Variable::curCost() " << name() << " _costrhs = " << costrhs() << "  _memorisedCurCost = "
                  << _memorisedCurCost << std::endl;

    return _memorisedCurCost;
}


const Double & Variable::costrhs() const
{
    if (printL(6))
        std::cout << " Variable::costrhs() " << name() << " _costrhs = " << _costrhs << "  _memorisedCurCost = "
                  << _memorisedCurCost << std::endl;

    return _costrhs;
}

void Variable::costrhs(const Double & newCostRhs)
{
  VarConstr::costrhs(newCostRhs);

  return;
}

void Variable::setMembership()
{
  if (!presetMembership())
    {
      AggregateVariable * agvPtr = dynamic_cast<AggregateVariable *>(this);
      if (agvPtr != NULL)
        agvPtr->agvSetMembership(this);
      /// it is enough to force enumerativeSetMembership() in constraint or if preset membership is not on
      enumerativeSetMembership();
    }

  VarConstr::setMembership();

  return;
}

void Variable::enumerativeSetMembership()
{
  if (printL(6))
    std::cout << name() << "  presetMembership() = " << presetMembership()
        << std::endl;

  if (presetMembership())
    return;

  if (problemPtr() != NULL)
    {
      setMembership(problemPtr()->probConstrSet());
    }

  return;
}

void Variable::setMembership(const ConstrIndexManager & constrSet)
{
  for (ConstrIndexManager::const_iterator it = constrSet.begin(Active, 's');
       it != constrSet.end(Active, 's'); ++it)
    {
        if (printL(6))
            std::cout << " Variable::setMembership TRY to add constr " << it->name() << std::endl;

      /// Does not work if this of that is a new added VarConstr if (!(*it)->presetMembership())
      addMember(*it);
    }

  for (ConstrIndexManager::const_iterator it = constrSet.begin(Inactive, 's');
       it != constrSet.end(Inactive, 's'); ++it)
    {
      addMember(*it);
    }

  for (ConstrIndexManager::const_iterator it = constrSet.begin(Unsuitable, 's');
       it != constrSet.end(Unsuitable, 's'); ++it)
    {
      addMember(*it);
    }

  for (ConstrIndexManager::const_iterator it = constrSet.begin(Active, 'd');
      it != constrSet.end(Active, 'd'); ++it)
    {
      if (printL(6))
        std::cout << " Variable::setMembership TRY to add constr "
        << it->name() << std::endl;

      /// Does not work if this of that is a new added VarConstr if (!(*it)->presetMembership())
      addMember(*it);
    }

  for (ConstrIndexManager::const_iterator it = constrSet.begin(Inactive, 'd');
      it != constrSet.end(Inactive, 'd'); ++it)
    {
      addMember(*it);
    }

  for (ConstrIndexManager::const_iterator it = constrSet.begin(Unsuitable, 'd');
      it != constrSet.end(Unsuitable, 'd'); ++it)
    {
      addMember(*it);
    }
}


void Variable::clearMembership()
{
  VarConstr::clearMembership();
  return;
}

bool Variable::suitableToFixValue(const Double & value)
{
  return (!value.isZero());
}

bool Variable::consecutive2varWhenBrOnCBS(InstanciatedVar * varPtr)
{
  return true;
}

const Double Variable::contrib(const Double & use)
{
  Double ctb(0);
  Double coef(0);

  for (ConstVarConstrPtr2Double::iterator itm = member2coefMap().begin(); itm != member2coefMap().end(); ++itm)
    {
      if (itm->first->inCurProb())
        {
          Constraint * constrPtr = dynamic_cast<Constraint *>(itm->first);
          coef = constrPtr->membCoef(this);
          switch (constrPtr->sense())
            {
                case 'L':
                {
                    if (coef * use > constrPtr->curRhs())
                        return (0);

                    break;
                }
                case 'E':
                {
                    if (coef * use > constrPtr->curRhs())
                        return (0);

                    ctb += coef * use;
                    break;
                }
                default: // 'G'
                {
                    if (coef * use > constrPtr->curRhs())
                        ctb += Dmax(constrPtr->curRhs(), 0);
                    else
                        ctb += coef * use;
                    break;
                }
            }

          if (printL(6))
            cout << "Var[" << name() << "] of use = " << use << " has coef "
                << constrPtr->membCoef(this) << " in const["
                << constrPtr->name() << "] of curRhs " << constrPtr->curRhs()
                << " ctb = " << ctb << std::endl;
        }
    }

  return (ctb);
}

void Variable::includeOvfPtr(OvfVar * ovfPtr)
{
  _ovfVarList.push_back(ovfPtr);
  return;
}

std::ostream & Variable::print(std::ostream& os) const
{
  os << "Variable base" << std::endl;
  VarConstr::print(os);
  os << "   isCandForBranching = " << _isCandForBranching  << std::endl;
  const AggregateVariable * const agvPtr = dynamic_cast<const AggregateVariable * const >(this);

  if (agvPtr != NULL)
    agvPtr->printAgv(os);

  return (os);
}

void Variable::resetCost(const bool & inPurePhaseOne)
{
  if (inPurePhaseOne)
    _memorisedCurCost = Double::staticZero;
  else
    _memorisedCurCost = costrhs();

  if (printL(6))
    std::cout << "  _memorisedCurCost = " << _memorisedCurCost << std::endl;

  return;
}

void Variable::resetCurCostByValue(const Double & newcost)
{
  if (printL(6))
    std::cout << " Variable::resetCurCostByValue()  var = " << name() << " curCost = " << newcost << std::endl;

  _memorisedCurCost = newcost;
}

void Variable::recallMemorisedBounds()
{
  _curLb = _memorisedCurLb;
  _curUb = _memorisedCurUb;
  if (printL(5))
    std::cout << "Variable::recallMemorisedBounds() " << name() << " in ["
              << _curLb << ", " << _curUb << "] " << std::endl;
  return;
}

bool Variable::infeasible() const
{
  /**
   * Variables:
   * 'P' = positive
   * 'N' = negative
   * 'F' = free
   */

  if (printL(6))
      std::cout << " Variable::infeasible(): name() = " << name() << " curLb() = " << curLb()
                << " curUb() = " << curUb() << std::endl;

  return curUb() < curLb();
}


/// TO DO : check whether reducedCost and incumbentVal are initialized correctly
SolutionVarInfo::SolutionVarInfo(Variable * varPtrV) :
    varPtr(varPtrV),
    value(varPtrV->val()),
    reducedCost(varPtrV->reducedCost()),
    incumbentValue(varPtrV->incumbentVal()),
    cost(varPtrV->curCost()),
    priorityLevel(1.0),
    canRoundDown(true),
    canRoundUp(true)
{
  if (value.fractional())
    {
      canRoundDown = varPtr->suitableToFixValue(Dfloor(varPtr->val()));
      canRoundUp = varPtr->suitableToFixValue(Dceil(varPtr->val()));
    }
  else if (value.isZero())
    {
      canRoundUp = canRoundDown = false;
    }

  if (varPtr->isTypeOf(VcId::MastColumnMask))
    {
      varPtr->incrParticipation(1);
      MastColumn * colPtr = static_cast<MastColumn *>(varPtrV);
      priorityLevel = colPtr->cgSpConfPtr()->priorityLevel();
    }
  else if (varPtr->isTypeOf(VcId::InstMasterVarMask))
  {
      InstanciatedVar * ivPtr = static_cast<InstanciatedVar *>(varPtr);
      priorityLevel = ivPtr->genVarPtr()->priorityLevel();
  }
}

SolutionVarInfo::SolutionVarInfo(const SolutionVarInfo & that):
    varPtr(that.varPtr),
    value(that.value),
    reducedCost(that.reducedCost),
    incumbentValue(that.incumbentValue),
    cost(that.cost),
    priorityLevel(that.priorityLevel),
    canRoundDown(that.canRoundDown),
    canRoundUp(that.canRoundUp)
{
    if (varPtr->isTypeOf(VcId::MastColumnMask))
        varPtr->incrParticipation(1);
}


SolutionVarInfo::~SolutionVarInfo()
{
    if (varPtr->isTypeOf(VcId::MastColumnMask))
        varPtr->decrParticipation(2);
}

bool Variable::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::VariableMask, vcIdentifier);
}


/**
 *  Methods of class AggregateVariable
 */

AggregateVariable::AggregateVariable(Solution * spSol) :
  _spSol((spSol == NULL ? NULL :spSol->clone())), _varPtr(NULL)
{
}

AggregateVariable::AggregateVariable(const AggregateVariable & that) :
    _spSol((that._spSol == NULL ? NULL : that._spSol->clone())), _varPtr(NULL)
{
}

void AggregateVariable::setAggregateVariable(Variable * varPtr)
{

  if (_varPtr != NULL) return; // setAggregateVariable already done

  _varPtr = varPtr;
  if (_spSol != NULL)
    {
      for (VarPtr2DoubleMap::const_iterator sPtr =
          _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end();
          sPtr++)
        {
          sPtr->first->includeAggregateVarAsMember(varPtr);
        }
    }

  return;
}
void AggregateVariable::unsetAggregateVariable()
{
  if (_varPtr == NULL) return;

  if (_spSol == NULL) return;

   for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
        sPtr != _spSol->solVarValMap().end(); sPtr++)
      sPtr->first->eraseAggregateVarAsMember(_varPtr);

  return;
}

AggregateVariable::~AggregateVariable()
{
  unsetAggregateVariable();
  /// Solution belong to subproblem (it s not the case when it is cloned inside a masterColumn)
     if (_spSol != NULL)
     {
       delete _spSol; _spSol = NULL;
     }
}


bool AggregateVariable::operator<(const VarConstr & that) const
{

  if (!that.isTypeOf(VcId::AggregateVariableMask))
    {
      return _varPtr->VarConstr::operator<(that);
    }

  if (this->_spSol == NULL)
    return (true);

  const AggregateVariable & bTemp = dynamic_cast<const AggregateVariable&>(that);

  if (bTemp._spSol == NULL)
    return (false);

  if (printL(6))
    std::cout << "AggregateVariable::LexicographicallyST(): both have SpSol" << std::endl;

  VarPtr2DoubleMap::const_iterator it1(this->_spSol->solVarValMap().begin());
  VarPtr2DoubleMap::const_iterator it2(bTemp._spSol->solVarValMap().begin());

  /// Map may not have same size
  for (;
       (it1 != this->_spSol->solVarValMap().end())
	 && (it2 != bTemp._spSol->solVarValMap().end()); ++it1, ++it2)
    {
      if (printL(6))
	std::cout << "AggregateVariable::LexicographicallyST(): comp v[" << it1->first->ref()
		  << "] = " << it1->second << " with v[" << it2->first->ref()
		  << "] = " << it2->second << std::endl;

      if (VcRefST(it1->first, it2->first))
	return (true);

      if (VcRefST(it2->first, it1->first))
	return (false);

      /// Same item
      if (it2->second < it1->second)
	return (true);

      if (it1->second < it2->second)
	return (false);
    }

  if (printL(6))
    std::cout << "AggregateVariable::LexicographicallyST(): same sp sol so far" << std::endl;

  /// Not same size
  if (this->_spSol->solVarValMap().size() < bTemp._spSol->solVarValMap().size())
    return (true);

  if (this->_spSol->solVarValMap().size() > bTemp._spSol->solVarValMap().size())
    return (false);

  return (false);
}


/// Assumes spSol var have been set before
void AggregateVariable::agvSetMembership(Variable * varPtr)
{
  if (_spSol == NULL)
    return;

  for (VarPtr2DoubleMap::const_iterator varIt = _spSol->solVarValMap().begin();
       varIt != _spSol->solVarValMap().end(); ++varIt)
    {
      if (printL(6))
        std::cout << "AggregateVariable::agvSetMembership(): solVarValMap includes var "
	              << varIt->first->name() << std::endl;

      for (ConstVarConstrPtr2Double::iterator constrIt = varIt->first->member2coefMap().begin();
	       constrIt != varIt->first->member2coefMap().end(); ++constrIt)
        {
          if (printL(6))
            std::cout << "AggregateVariable::agvSetMembership(): var "
		              << varIt->first->name() << " in constr " << constrIt->first->name() << std::endl;
	      varPtr->includeMember(constrIt->first, constrIt->second * varIt->second, true);
        }
    }

  return;
}

bool AggregateVariable::agvComputeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(7))
    std::cout << "AggregateVariable::agvComputeCount() this " << _varPtr->name()
              << " that " << vcPtr->name() << std::endl;

  if (_spSol == NULL)
    return (false);

  bool status(false);
  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end();
       sPtr++)
    {
      if (printL(7))
        std::cout << "AggregateVariable::computeCount(): SP var "
            << sPtr->first->name() << " constr " << vcPtr->name() << std::endl;

      if (vcPtr->membCount(sPtr->first))
        {
          if (printL(7))
            std::cout << "YES" << std::endl;

          status = true;
        }
      else if (printL(7))
        std::cout << "NO" << std::endl;
    }

  if (printL(7))
    std::cout << "AggregateVariable::computeCount() :   " << status << std::endl;

  return status;
}

const LpCoef AggregateVariable::agvComputeCoef(ConstVarConstrConstPtr vcPtr)
{
  Double coef(0);

  if (printL(7))
    std::cout << "AggregateVariable::computeCoef() constrPtr->name() in " << vcPtr->name() << std::endl;

  if (_spSol == NULL)
    return LpCoef::ZeroCoef;

  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
      sPtr != _spSol->solVarValMap().end(); sPtr++)
    {
      coef += vcPtr->membCoef(sPtr->first) * sPtr->second;
      if (printL(7))
        std::cout << "AggregateVariable::computeCoef(): var "
            << sPtr->first->name() << " has val " << sPtr->second
            << " in constr " << vcPtr->name() << " aggregatCoef = " << coef << std::endl;
    }

  if (printL(7))
    std::cout << "AggregateVariable::computeCoef() coef =  " << coef << std::endl;

  return LpCoef(coef);
}

std::ostream & AggregateVariable::printAgv(std::ostream& os) const
{
  os << " AggregateVariable Solution has spSol " << (_spSol != NULL) << std::endl;

  if (printL(6))
    {
      if (_spSol != NULL)
        for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
             sPtr != _spSol->solVarValMap().end(); sPtr++)
          std::cout << "AggregateVariable include var " << sPtr->first->name()
                    << " at coef = " << sPtr->second << std::endl;
    }

  return os;
}

void AggregateVariable::includeVar(Variable * aggVarPtr, Variable * varPtr,
                                   const Double & val, const bool & cumulativeVal)
{
    _spSol->includeVar(varPtr, val, cumulativeVal);
    varPtr->includeAggregateVarAsMember(aggVarPtr);
}

bool AggregateVariable::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::AggregateVariableMask, vcIdentifier);
}

/**
 * Methods of class Constraint
 */

Constraint::Constraint(Model * modelPtr,
		       const std::string & name,
		       const Double & costrhs,
		       const char & sense,
		       const char & type,
		       const char & kind,
		       const char & flag,
		       const int & index,
		       const Double & val,
		       const Double & upperBound,
		       const Double & lowerBound,
		       const char & directive,
		       const Double & priority,
		       const bool & presetMembership,
		       const bool & ltoBeUsedInPreprocessing,
		       const bool & considerAsEqualityInPreprocessing) :
  VarConstr(modelPtr,
	        modelPtr->modelConstrCnt(),
            name,
	        costrhs,
	        sense,
	        type,
	        kind,
	        flag,
	        index,
	        val,
	        upperBound,
	        lowerBound,
	        upperBound,
	        lowerBound,
	        directive,
	        priority,
	        presetMembership),
  _memorisedCurRhs(costrhs),
  _lhsComputed((ltoBeUsedInPreprocessing?false:true)),
  _toBeUsedInPreprocessing(ltoBeUsedInPreprocessing),
  _considerAsEqualityInPreprocessing(considerAsEqualityInPreprocessing),
  _posLocalArtVarPtr(NULL),
  _negLocalArtVarPtr(NULL),
  _minSlack(-BapcodInfinity),
  _maxSlack(BapcodInfinity),
  _curMinSlack(-BapcodInfinity),
  _curMaxSlack(BapcodInfinity),
  _inPropagateList(false)
{
    modelPtr->increaseModelConstrCnt();
    if (printL(6))
        std::cout << "Constraint::Constraint() " << name << " toBeUsedInPreprocessing = " << toBeUsedInPreprocessing()
	              << " presetMembership = " << presetMembership << " presetMembership() = "
	              << VarConstr::presetMembership() << std::endl;
}

Constraint::Constraint(const Constraint & that) :
  VarConstr(that, that._modelPtr->modelConstrCnt()),
  _minSlack(that._minSlack),
  _maxSlack(that._maxSlack),
  _curMinSlack(that._curMinSlack),
  _curMaxSlack(that._curMaxSlack),
  _memorisedCurRhs(that._costrhs),
  _lhsComputed(that._lhsComputed),
  _toBeUsedInPreprocessing(that._toBeUsedInPreprocessing),
  _considerAsEqualityInPreprocessing(that._considerAsEqualityInPreprocessing),
  _posLocalArtVarPtr(that._posLocalArtVarPtr),
  _negLocalArtVarPtr(that._negLocalArtVarPtr)
{
    that._modelPtr->increaseModelConstrCnt();
    if (_posLocalArtVarPtr != NULL)
        _posLocalArtVarPtr->setConstraintPtr(this);
    if (_negLocalArtVarPtr != NULL)
        _negLocalArtVarPtr->setConstraintPtr(this);
    if (that._stabInfoPtr != NULL)
        _stabInfoPtr = new VarConstrStabInfo(*(that._stabInfoPtr), this);
    else
        _stabInfoPtr = NULL;
}

Constraint::~Constraint()
{
  /**
   * Do not decrement _staticCounter as it would perturbate
   * lexicographix ordering of remaining constraints
   */
  if(printL(7))
    cout << "constraint at 0x" << hex << (long)this << dec << " has been deleted" << endl;

  deleteStabInfoPtr();
  return;
}

void Constraint::includeOvfPtr(OvfConstr * ovfPtr)
{
  _ovfConstrList.push_back(ovfPtr);
  return;
}

void Constraint::resetRhs()
{
  _memorisedCurRhs = costrhs();
  if (printL(6))
    std::cout << "Constraint::resetRhs: chg rhs of  " << name() << " index " << index() << std::endl;
  return;
}

void Constraint::curRhs(const Double & rhs)
{
  _memorisedCurRhs = rhs;
}

const Double Constraint::curRhs() const
{
  return _memorisedCurRhs;
}

void Constraint::addToPreprocessedList()
{
  if (!_inPreprocessedList)
    {
      _problemPtr->preprocessedConstrsList().push_back(this);
      _inPreprocessedList = true;
    }
}

void Constraint::resetSlacksAndRhsToDefaults()
{
  _curMaxSlack = _maxSlack;
  _curMinSlack = _minSlack;
  _memorisedCurRhs = _costrhs;
}

bool Constraint::activateConstraint(bool toSetInForm)
{
    if (_problemPtr == NULL)
        return false;

    _problemPtr->probConstrSet().insert(this, Active);
    activate();
    if (toSetInForm)
        _problemPtr->setConstr2Form(this);

    if (_posLocalArtVarPtr != NULL)
        _posLocalArtVarPtr->activateVariable(toSetInForm);
    if (_negLocalArtVarPtr != NULL)
        _negLocalArtVarPtr->activateVariable(toSetInForm);
    if (_stabInfoPtr != NULL)
    {
        Variable * varPtr = _stabInfoPtr->negInnerArtVarPtr();
        if (varPtr != NULL)
            varPtr->activateVariable(toSetInForm);
        varPtr = _stabInfoPtr->negOuterArtVarPtr();
        if (varPtr != NULL)
            varPtr->activateVariable(toSetInForm);
        varPtr = _stabInfoPtr->posInnerArtVarPtr();
        if (varPtr != NULL)
            varPtr->activateVariable(toSetInForm);
        varPtr = _stabInfoPtr->posOuterArtVarPtr();
        if (varPtr != NULL)
            varPtr->activateVariable(toSetInForm);
    }
    return true;
}

bool Constraint::desactivateConstraint(const VcIndexStatus::VcStatus & status, bool toUnsetInForm)
{
  if (_problemPtr == NULL)
    return false;

  _problemPtr->probConstrSet().insert(this, status);
  desactivate();
  if (toUnsetInForm)
    _problemPtr->unsetConstr2Form(this);

  if (_posLocalArtVarPtr != NULL)
    _posLocalArtVarPtr->desactivateVariable(status, toUnsetInForm);
  if (_negLocalArtVarPtr != NULL)
    _negLocalArtVarPtr->desactivateVariable(status, toUnsetInForm);
  if (_stabInfoPtr != NULL)
    {
      Variable * varPtr = _stabInfoPtr->negInnerArtVarPtr();
      if (varPtr != NULL)
        varPtr->desactivateVariable(status, toUnsetInForm);
      varPtr = _stabInfoPtr->negOuterArtVarPtr();
      if (varPtr != NULL)
        varPtr->desactivateVariable(status, toUnsetInForm);
      varPtr = _stabInfoPtr->posInnerArtVarPtr();
      if (varPtr != NULL)
        varPtr->desactivateVariable(status, toUnsetInForm);
      varPtr = _stabInfoPtr->posOuterArtVarPtr();
      if (varPtr != NULL)
        varPtr->desactivateVariable(status, toUnsetInForm);
    }
  return true;
}

void Constraint::setMembership()
{
  if (printL(6))
    std::cout << " Constraint::setMembership   " << name() << " presetMembership() = " << presetMembership()
              << std::endl;

  if (!presetMembership())
  {
    enumerativeSetMembership();
  }
  else
  {
    bool addGlobalArtVars = false;
    if (isTypeOf(VcId::InstMasterConstrMask))
    {
      InstMasterConstr * iMastConstrPtr = static_cast<InstMasterConstr*> (this);
      if (!iMastConstrPtr->subProbVarMember2coefMap().empty())
      {
        addGlobalArtVars = true;
      }
    }
    if (addGlobalArtVars)
    {
      if (printL(6))
        std::cout << " Constraint::setMembership   " << name()
        << " adding artificial variables as members" << std::endl;

      if (problemPtr()->posGlobalArtVarPtr() != NULL)
        problemPtr()->posGlobalArtVarPtr()->addMember(this);
      if (problemPtr()->negGlobalArtVarPtr() != NULL)
        problemPtr()->negGlobalArtVarPtr()->addMember(this);
      if (negLocalArtVarPtr() != NULL)
        negLocalArtVarPtr()->addMember(this);
      if (posLocalArtVarPtr() != NULL)
        posLocalArtVarPtr()->addMember(this);
    }

    if (_stabInfoPtr != NULL)
      _stabInfoPtr->setMembership();
  }

  VarConstr::setMembership();

  return;
}

void Constraint::enumerativeSetMembership()
{
  if (presetMembership())
      return;

  if (problemPtr() != NULL)
      setMembership(problemPtr()->probVarSet());

  return;
}

void Constraint::setMembership(const VarPtrSet & varSet)
{
  if (printL(6))
    std::cout << " Constraint::setMembership " << name() << "  varSet size =   " << varSet.size() << std::endl;

  for (VarPtrSet::const_iterator it = varSet.begin(); it != varSet.end(); ++it)
    {
      if (printL(6))
	    std::cout << " Constraint::setMembership try adding  " << (*it)->name() << std::endl;

	  addMember(*it);
    }
}

void Constraint::setMembership(const VarIndexManager & varSet, const VcIndexStatus::VcStatus& status, char flag)
{
  for (VarIndexManager::const_iterator it = varSet.begin(status, flag); it != varSet.end(status, flag); ++it)
    {
      addMember(*it);
    }
}

void Constraint::setMembership(const VarIndexManager & varSet)
{
  setMembership(varSet, Active, 's');
  setMembership(varSet, Inactive, 's');
  setMembership(varSet, Unsuitable, 's');

  setMembership(varSet, Active, 'd');
  setMembership(varSet, Inactive, 'd');
  setMembership(varSet, Unsuitable, 'd');

  setMembership(varSet, Active, 'a');
  setMembership(varSet, Inactive, 'a');
  setMembership(varSet, Unsuitable, 'a');
}


/// Is there  only to  handle the case of artificial variables
bool Constraint::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << " Constraint::computeCount this " << name() << " that "
        << vcPtr->name() << std::endl;

  /// start changed by Ruslan

  if (vcPtr->isTypeOf(VcId::GlobalArtificialVarMask))
    return true;
  else if (vcPtr->isTypeOf(VcId::LocalArtificialVarMask))
    {
     if ((_negLocalArtVarPtr == vcPtr) || (_posLocalArtVarPtr == vcPtr))
       return true;
     else if ((_stabInfoPtr != NULL) && _stabInfoPtr->computeCount(vcPtr) )
      return true;
    }

  /// end changed by Ruslan

  return (false);
}

/**
 * Is there  only to  handle the case of artificial variables
 */
const LpCoef Constraint::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::GlobalArtificialVarMask))
  {

   if (static_cast<GlobalArtificialVar*>(vcPtr)->senseType() == 'G')
    {
      if (sense() != 'L') /// in 'G' or 'E' constraints
      {
        return LpCoef::UnitCoef;
      }
      else
      {
        return LpCoef::MinusOneCoef;
      }
    }
    else
    {
      if (sense() != 'L') /// in 'G' or 'E' constraints
      {
        return LpCoef::MinusOneCoef;
      }
      else
      {
        return LpCoef::UnitCoef;
      }
    }
  }

  if (vcPtr->isTypeOf(VcId::LocalArtificialVarMask))//classId() == VarConstr::LocalArtificialVarId)
    {
      if (_negLocalArtVarPtr == vcPtr)
        return LpCoef::MinusOneCoef;
      else if (_posLocalArtVarPtr == vcPtr)
        return LpCoef::UnitCoef;
      else if (_stabInfoPtr != NULL)
        {
    	  return _stabInfoPtr->computeCoef(vcPtr);
        }
    }

  return LpCoef::ZeroCoef;
}

const Double & Constraint::computeViolation(const Double & lhs)
{
    switch (sense())
    {
        case 'G':
        {
            violation(Dmax(curRhs() - lhs, 0));
            break;
        }
        case 'L':
        {
            violation(Dmax(lhs - curRhs(), 0));
            break;
        }
        default:
            violation(Dmax(lhs - curRhs(), curRhs() - lhs));
            break;
    }

    return (violation());
}

bool Constraint::violated(Variable * varPtr, const Double & useLevel)
{
  return computeViolation(varPtr->membCoef(this) * useLevel).positive();
}

Double Constraint::computeLhs(const std::list<std::pair<Variable *, Double> > & curSol)
{
  Double curLhs(0);
  for (std::list<std::pair<Variable *, Double> >::const_iterator it = curSol.begin(); it != curSol.end(); ++it)
    {
      curLhs += it->second * membCoef(it->first);
      if (printL(6))
        std::cout << "Constraint::computeLhs(): curSol includes " << it->first->name() << " at val = " << it->second
                  << " curLhs = " << curLhs << std::endl;
    }

  return curLhs;
}

Double Constraint::computeLhs(const VarPtr2DoubleMap & solVarValMap)
{
    Double curLhs(0);
    for (VarPtr2DoubleMap::const_iterator it = solVarValMap.begin(); it != solVarValMap.end(); ++it)
    {
        curLhs += it->second * membCoef(it->first);
        if (printL(6))
            std::cout << "Constraint::computeLhs(): curSol includes " << it->first->name() << " at val = "
                      << it->second << " curLhs = " << curLhs << std::endl;
    }

    return curLhs;
}


Double Constraint::computeLhs(const VarPtrSet & curSol)
{
    Double curLhs(0);
    for (VarPtrSet::const_iterator it = curSol.begin(); it != curSol.end(); ++it)
    {
        curLhs += (*it)->val() * membCoef(*it);
        if (printL(6))
            std::cout << "Constraint::computeLhs(): curSol includes "
                      << (*it)->name() << " at val = " << (*it)->val() << " curLhs = " << curLhs << std::endl;
    }

    return curLhs;
}

Double Constraint::computeLhs(const std::list<Variable *> & curSol)
{
    Double curLhs(0);
    for (std::list<Variable *>::const_iterator it = curSol.begin(); it != curSol.end(); ++it)
    {
        curLhs += (*it)->val() * membCoef(*it);
        if (printL(6))
            std::cout << "Constraint::computeLhs(): curSol includes "
                      << (*it)->name() << " at val = " << (*it)->val() << " _curLhs = "
                      << curLhs << std::endl;
    }

    return curLhs;
}

std::ostream & Constraint::print(std::ostream& os) const
{
  os << "Constraint base" << std::endl;
  VarConstr::print(os);
  if (_stabInfoPtr != NULL)
    _stabInfoPtr->print(os);

  return (os);
}

bool Constraint::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::ConstraintMask, vcIdentifier);
}

LocalArtificialVar * Constraint::posLocalArtVarPtr() const
{
  return (_posLocalArtVarPtr);
}
LocalArtificialVar * Constraint::negLocalArtVarPtr() const
{
  return (_negLocalArtVarPtr);
}

void Constraint::posLocalArtVarPtr(LocalArtificialVar * ptr)
{
  _posLocalArtVarPtr = ptr;
  return;
}

void Constraint::negLocalArtVarPtr(LocalArtificialVar * ptr)
{
  _negLocalArtVarPtr = ptr;
  return;
}

bool Constraint::operator<(const VarConstr& b) const
{
  if (!b.isTypeOf(VcId::ConstraintMask))
    {
      return VarConstr::operator<(b);
    }

  const Constraint & bTemp = static_cast<const Constraint&>(b);

  if (printL(8))
    std::cout << " " << this->name() << " " << std::hex << (long)(this) << std::dec
              << " L< " << bTemp.name() << " " << std::hex << (long)(&bTemp) << std::dec
              <<  std::endl;

  if (printL(7))
    std::cout << "comp rhs " << this->costrhs()
        << " with " << bTemp.costrhs() << std::endl;

  if (this->costrhs() < bTemp.costrhs())
    return (true);

  if (this->costrhs() > bTemp.costrhs())
    return (false);

  ConstVarConstrPtr2Double::const_iterator it1 = this->member2coefMap().begin();
  ConstVarConstrPtr2Double::const_iterator it2 = bTemp.member2coefMap().begin();
  for (;(it1 != this->member2coefMap().end()) && (it2 != bTemp.member2coefMap().end()); ++it1, ++it2)
    {
      if (VcRefST(it1->first, it2->first))
        return (true);

      if (VcRefST(it2->first, it1->first))
        return (false);
      /// Same item
      if (it1->second > it2->second)
        return (true);

      if (it1->second < it2->second)
        return (false);
    }

  /// Not same size
  if (this->member2coefMap().size() > bTemp.member2coefMap().size())
    return (true);

  if (this->member2coefMap().size() < bTemp.member2coefMap().size())
    return (false);

  return (false);
}

/**
 * Methods of class ArtificialVar
 */

ArtificialVar::ArtificialVar(const Double & defaultCost) :
    _defaultCost(defaultCost)
{
  return;
}

ArtificialVar::~ArtificialVar()
{
  return;
}

const Double & ArtificialVar::defaultCost()
{
  return (_defaultCost);
}

void ArtificialVar::setCost(const Double & cost)
{
  _defaultCost = cost;
  return;
}

bool ArtificialVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::ArtificialVarMask, vcIdentifier);
}

/**
 * Methods of class GlobalArtificialVar
 */

GlobalArtificialVar::GlobalArtificialVar(Model * modelPtr,
                                         const Double & cost,
					                     const char & senseType,
					                     const BcObjStatus::MinMaxIntFloat & minmax,
					                     const std::string & lname) :
  Variable(modelPtr,
	       lname,
    	   (( minmax == BcObjStatus::minInt) || ( minmax == BcObjStatus::minFloat)) ? cost : -cost,
    	   'P',
	       'C',
	       'E',
	       BapcodInfinity,
    	   0,
	       'a',
	       'U',
    	   1.0,
	       0,
	       BapcodInfinity,
    	   0,
	       true,
	       -1),
  ArtificialVar((( minmax == BcObjStatus::minInt) || ( minmax == BcObjStatus::minFloat)) ? cost : -cost),
  _senseType(senseType)
{
  std::string rname;
  if (_senseType == 'G')
    rname = "posGlobArtVar";
  else rname = "negGlobArtVar";

  name(rname);

  if (printL(6))
    std::cout << "GlobalArtificialVar::GlobalArtificialVar() " << name()
      	      << " _senseType = " << _senseType
	          << " in [" << _lowerBound << ", " << _upperBound << "] " << std::endl;
}

GlobalArtificialVar::~GlobalArtificialVar()
{
  return;
}

void GlobalArtificialVar::resetCostFromDefaultCost(const Double & factor)
{
  costrhs(defaultCost() * (1+factor));
  _memorisedCurCost = _costrhs; /// added by Ruslan
}

const Double & GlobalArtificialVar::costrhs() const
{
 if (printL(6))
   std::cout << " GlobalArtificialVar::costrhs() " << name() << " _costrhs = " << _costrhs << "  _memorisedCurCost = "
             << _memorisedCurCost << std::endl;

    return _costrhs;
}

const Double & GlobalArtificialVar::curCost() const
{
  if (printL(6))
    std::cout << " GlobalArtificialVar::curCost() " << name() << " _costrhs = " << costrhs() << "  _memorisedCurCost = "
	          << _memorisedCurCost << std::endl;

  return (_memorisedCurCost);
}


bool GlobalArtificialVar::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "GlobalArtificialVar::computeCount() this " << name()
              << " that " << vcPtr->name() << std::endl;

  bapcodInit().require( vcPtr->isTypeOf(VcId::ConstraintMask),
                        "GlobalArtificialVar::count() should not be called with parameter other than constraint");

  return vcPtr->computeCount(this);
}

const LpCoef GlobalArtificialVar::computeCoef(ConstVarConstrConstPtr vcPtr)
{

  bapcodInit().require(vcPtr->isTypeOf(VcId::ConstraintMask),
                       "GlobalArtificialVar::count() should not be called with parameter other than constraint");
  return vcPtr->computeCoef(this);
}

void GlobalArtificialVar::setMembership(const ConstrIndexManager & constrSet, const VcIndexStatus::VcStatus& status,
                                        char flag)
{
  for (ConstrIndexManager::const_iterator it = constrSet.begin(status, flag);
        it != constrSet.end(status, flag); ++it)
      {
        if (printL(6))
          std::cout << " GlobalArtificialVar::setMembership TRY to add constr " << (*it)->name() << std::endl;

        if (senseType() == 'G')
          {
            if ((*it)->sense() != 'L') /// in 'G' or 'E' constraints
              includeMember((*it), 1, false);
            else
              includeMember((*it), -1, false);
          }
        else ///(senseType() == 'L')
          {
            if ((*it)->sense() != 'L')  /// in 'G' or 'E' constraints
              includeMember((*it), -1, false);
            else
              includeMember((*it), 1, false);
          }
      }
}

void GlobalArtificialVar::setMembership(const ConstrIndexManager & constrSet)
{
  if (printL(6))
    std::cout << " GlobalArtificialVar::setMembership constrSet size = " << constrSet.size() << std::endl;

  setMembership(constrSet, Active, 's');
  setMembership(constrSet, Inactive, 's');
  setMembership(constrSet, Unsuitable, 's');

  setMembership(constrSet, Active, 'd');
  setMembership(constrSet, Inactive, 'd');
  setMembership(constrSet, Unsuitable, 'd');

  return;
}

bool GlobalArtificialVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::GlobalArtificialVarMask, vcIdentifier);
}


/**
 * Methods of class LocalArtificialVar
 */

LocalArtificialVar::LocalArtificialVar(Constraint * constrPtr,
                                       const LocalArtClassId locClassId,
                                       const BcObjStatus::MinMaxIntFloat & minmax,
                                       const std::string & name,
                                       const Double & defCost,
                                       const Double & ub) :
  Variable(constrPtr->modelPtr(),
	   std::string(name) + constrPtr->name(),
	   ((( minmax == BcObjStatus::minInt) || ( minmax == BcObjStatus::minFloat))? defCost : - defCost),
	   'P',
	   'I',
	   'E',
	   ub,
	   0,
	   'a', /// changed from 's' to 'a' by Ruslan
	   'U',
	   1.0,
	   0,
	   BapcodInfinity,
	   0,
	   true,
	   -1),
  ArtificialVar((((minmax == BcObjStatus::minInt) || ( minmax == BcObjStatus::minFloat)) ? defCost
                                                                                                    : - defCost)),
  _constrPtr(constrPtr),
  _tpVal(1),
  _locClassId(locClassId)
{

  if (printL(6))
    std::cout << "LocalArtificialVar::LocalArtificialVar() " << name
              << " in [" << _lowerBound << ", " << _upperBound << "] " << " objStatus = " << minmax
	          << " defaultCost() = " << defaultCost() << std::endl;
}

LocalArtificialVar::~LocalArtificialVar()
{
  return;
}


const Double & LocalArtificialVar::defaultCost()
{

  _defaultCost = ArtificialVar::defaultCost();

  return _defaultCost;
}

const Double & LocalArtificialVar::costrhs() const
{
 if (printL(6))
   std::cout << " LocalArtificialVar::costrhs() " << name() << " _costrhs = " << _costrhs << "  _memorisedCurCost = "
             << _memorisedCurCost << std::endl;

    return _costrhs;
}

const Double & LocalArtificialVar::curCost() const
{
  if (printL(6))
    std::cout << " LocalArtificialVar::curCost() " << name() << " _costrhs = " << costrhs() << "  _memorisedCurCost = "
	          << _memorisedCurCost << std::endl;

  return (_memorisedCurCost);
}


void LocalArtificialVar::resetCostFromDefaultCost(const Double & factor) ///changed by Ruslan
{
  costrhs(defaultCost() * (1+factor));
 _memorisedCurCost = _costrhs; /// added by Ruslan
}

bool LocalArtificialVar::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "LocalArtificialVar::computeCount() this " << name()
	      << " that " << vcPtr->name() << std::endl;

  bapcodInit().require(vcPtr->isTypeOf(VcId::ConstraintMask),
                       "LocalArtificialVar::coef() should not be called wipt parameter other than constraint");

  return (vcPtr->computeCount(this));
}

const LpCoef LocalArtificialVar::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  bapcodInit().require(vcPtr->isTypeOf(VcId::ConstraintMask),
                       "LocalArtificialVar::count() should not be called wipt parameter other than constraint");

  return (vcPtr->computeCoef(this));
}

void LocalArtificialVar::setMembership()
{
  if (printL(5))
    std::cout << "setMembership of local artificial variable " << name() << std::endl;

  if (_constrPtr != NULL)
    {
      Double coeff(_constrPtr->costrhs());
      if (coeff == 0)
        coeff = 1;
      switch (_locClassId)
      {
      case PosOuterId:
      case PosInnerId:
        includeMember(_constrPtr, coeff, false);
        break;
      case NegOuterId:
      case NegInnerId:
        includeMember(_constrPtr, -coeff, false);
        break;
      case PosLocalId:
        includeMember(_constrPtr, 1, false);
        break;
      case NegLocalId:
        includeMember(_constrPtr, -1, false);
        break;
      }
    }

  VarConstr::setMembership();

  return;
}

const LocalArtificialVar::LocalArtClassId & LocalArtificialVar::localClassId() const
{
  return _locClassId;
}

void LocalArtificialVar::setConstraintPtr(Constraint * constrPtr)
{
  _constrPtr = constrPtr;
  return;
}

bool LocalArtificialVar::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::LocalArtificialVarMask, vcIdentifier);
}

void MasterVarSolution::push_back(Variable * varPtr, const ValueRecord  & rec)
{
  std::list< std::pair < Variable *, ValueRecord > >::push_back(make_pair(varPtr, rec));
  return;
}
