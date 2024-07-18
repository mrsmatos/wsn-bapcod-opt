/**
 *
 * This file bcGenVarBranchConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcGenBranchingConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"


using namespace std;

/// Methods of class GenVarBranchConstrGenerator

GenVarBranchConstrGenerator::GenVarBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
			                                             Variable * varPtr,
                                                         const char & priorityDir,
			                                             const Double & candidateLhs) :
  BranchingConstrGenerator(gbcPtr, priorityDir, candidateLhs), _varPtr(varPtr)
{
  _description = _varPtr->name();
}

GenVarBranchConstrGenerator::GenVarBranchConstrGenerator(const GenVarBranchConstrGenerator & that) :
  BranchingConstrGenerator(that), _varPtr(that._varPtr)
{
}

GenVarBranchConstrGenerator::~GenVarBranchConstrGenerator()
{
}

void GenVarBranchConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curSol)
{
  _candidateLhs = 0;

  for(SolutionVarInfoPtrList::const_iterator infoIt = curSol.begin(); infoIt != curSol.end();
          infoIt++)
    {
      if (!(*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
        {
          /// OFV type of variable, i.e. pure master variable
          if ((*infoIt)->varPtr == _varPtr)
            _candidateLhs += (*infoIt)->value;

          continue;
        }

      if ((*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
	    {
	      MastColumn * MastColumnPtr = static_cast<MastColumn *> ((*infoIt)->varPtr);

	      if (!MastColumnPtr->spVarCount(_varPtr))
	        continue;

	      _candidateLhs += MastColumnPtr->spVarVal(_varPtr) * (*infoIt)->value;
	    }
    }

  return;
}

bool GenVarBranchConstrGenerator::nextNodeBrConstr(Node * parentNodePtr,
                                                   std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList,
                                                   const ConstrPtrSet & existingMasterBranchingConstr)
{
  /// Node branching constraint defined last is treated first (LIFO)
  nextBranchingConstrPtrList.clear();
  bool success(true);

  int ancestorNodeRef(-1);
  if(parentNodePtr != NULL)
    ancestorNodeRef = parentNodePtr->ref();

  switch(_direction)
    {
    case 'U':
      {
        if(_childNbCounter == 0)
          instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G',
                              nextBranchingConstrPtrList);
        else if(_childNbCounter == 1)
          instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L',
                              nextBranchingConstrPtrList);
        else
          success = false;

        return (success);
      }
    default:
      {
        if(_childNbCounter == 0)
          instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L',
                              nextBranchingConstrPtrList);
        else if(_childNbCounter == 1)
          instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G',
                              nextBranchingConstrPtrList);
        else
          success = false;

        return (success);
      }
    }

  /// Not path to this instruction
  return (success);
}

void GenVarBranchConstrGenerator::instanciateBrConstr(const int & parentNodeNb, const int & childNb, const Double & rhs,
                                                      const char & sense,
                                                      std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
  bapcodInit().require(_varPtr->isTypeOf(VcId::InstanciatedVarMask),
	  "GenVarBranchConstrGenerator::instanciateBrConstr(): varPtr should be of type InstanciatedVar");

  std::string name("BCV");
  name = name + _varPtr->name();

  InstanciatedVar * ivPtr = static_cast<InstanciatedVar *> (_varPtr);

  if (printL(5))
      std::cout << "GenVarBranchConstrGenerator::instanciateBrConstr() " << name << std::endl;


  BranchingConstrBaseType * ptr = new GenVarInstMastBranchConstr(ivPtr->id(), _genericBrConstrPtr,
								                                 _genericBrConstrPtr->modelPtr()->masterConfPtr(),
                                                                 ivPtr, name + "p" + parentNodeNb + "c" + childNb, rhs,
                                                                 sense, ' ', 'E', 'd');
  if(printL(5))
    ptr->print(cout);

  nextBranchingConstrPtrList.push_back(ptr);

  return;
}

std::ostream & GenVarBranchConstrGenerator::print(std::ostream& os) const
{
  BranchingConstrGenerator::print(os);
  os << "GenVarBranchConstrGenerator" << std::endl;
  os << "   ivar = " << _varPtr->name() << std::endl;
  os << "   candidateLhs = " << _candidateLhs << std::endl;

  return (os);
}

void GenVarBranchConstrGenerator::nicePrint(std::ostream& os) const
{
  os << "var " << _varPtr->name() << " (lhs=" << _candidateLhs << ")";
}

bool BrVarPriorityCalcAndComp_FirstFound::operator()(Variable * a, Variable * b)
{
  /// Follows priority
  if (a->priority() > b->priority())
    /// Higher priority var chosen for branching
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_HighestPriority::operator()(Variable * a, Variable * b)
{
  /// Follows priority
  if (a->priority() > b->priority())
    /// Higher priority var chosen for branching
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_FracWeightedPriority::operator()(Variable * a, Variable * b)
{
  /// Follows priority weighted by distance from int value with direction precribed by userclosest value
  if (a->fracPartRelativeToTarget() * a->priority() > b->fracPartRelativeToTarget() * b->priority())
    // Higher frac weighted priority var chosen for branching
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_MostFractional::operator()(Variable * a, Variable * b)
{
  /// distance from int value with direction precribed by userclosest value
  if(a->fracPartRelativeToTarget() > b->fracPartRelativeToTarget())
    return true;
  if((a->fracPartRelativeToTarget() >= b->fracPartRelativeToTarget()) && (a->priority() > b->priority()))
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_LeastFractional::operator()(Variable * a, Variable * b)
{
  /// Closest to integer value with direction toward closest value
  Double viol_a, viol_b;
  char dir;
  calculateVarViolAndDir(a, viol_a, dir);
  calculateVarViolAndDir(b, viol_b, dir);
  if(viol_a < viol_b)
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_Closest2RoundUp::operator()(Variable * a, Variable * b)
{
  /// Closest to rounded up
  if(a->lFracPart() > b->lFracPart())
    return true;
  if((a->lFracPart() >= b->lFracPart()) && (a->priority() > b->priority()))
    return true;
  return false;
}

bool BrVarPriorityCalcAndComp_Closest2RoundDown::operator()(Variable * a, Variable * b)
{
  /// Closest to rounded down
  if(a->uFracPart() > b->uFracPart())
    return true;
  if((a->uFracPart() >= b->uFracPart()) && (a->priority() > b->priority()))
    return true;
  return false;
}

void BrVarPriorityCalcAndComp_FirstFound::calculateVarViolAndDir(Variable * var, Double & viol, char & dir)
{
  /// Follows priority
  viol = var->fracPart();
  dir = var->directive();
}

/// added by Ruslan : please verify that correct
void BrVarPriorityCalcAndComp_HighestPriority::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// Follows priority
  viol = var->fracPart();
  dir = var->directive();
}

void BrVarPriorityCalcAndComp_FracWeightedPriority::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// Follows priority weighted by distance from int value with direction precribed by userclosest value
  viol = var->fracPartRelativeToTarget() * var->priority();
  dir = var->directive();
}

void BrVarPriorityCalcAndComp_MostFractional::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// distance from int value with direction precribed by userclosest value
  viol = var->fracPartRelativeToTarget();
  dir = var->directive();
}

void BrVarPriorityCalcAndComp_LeastFractional::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// Closest to integer value with direction toward closest value
  Double temp = var->lFracPart();
  if(var->tmpVal() < temp + 0.5)
    {
      viol = var->tmpVal() - temp;
      dir = 'L';
    }
  else // (var->tmpVal() >= temp + 0.5)
    {
      viol = Dceil(var->tmpVal()) - var->tmpVal();
      dir = 'U';
    }
}

void BrVarPriorityCalcAndComp_Closest2RoundUp::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// Closest to rounded up
  viol = var->lFracPart();
  dir = 'U';
}

void BrVarPriorityCalcAndComp_Closest2RoundDown::calculateVarViolAndDir(Variable* var, Double & viol, char & dir)
{
  /// Closest to rounded down
  viol = var->uFracPart();
  dir = 'L';
}

/// Methods of class GenVarGenBranchConstr

GenVarGenBranchConstr::GenVarGenBranchConstr(Model * modelPtr, ProbConfig * probConfPtr, GenericVar * genVarPtr,
                                             const SelectionStrategy & priorityRule, const Double & priorityLevel) :
  GenericBranchingConstr(modelPtr, probConfPtr, /*genVarPtr, */"GVGbr", priorityRule, priorityLevel, true),
  _genVarPtr(genVarPtr)
{
}

GenVarGenBranchConstr::~GenVarGenBranchConstr()
{
}

const Double & GenVarGenBranchConstr::genericCostRhs(const InstanciatedVarConstr * const ivarconstrPtr) const
{
  std::cout << "GenVarGenBranchConstr::genericCostRhs(): error should not be called" << std::endl;
  return Double::staticZero;
}

void GenVarGenBranchConstr::branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                              const int & maxNumOfCandidates,
                                                              BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  /// project onto  OVF solution
  VarPtr2DoubleMap curAggregateMastSol;
  for (MasterColSolution::const_iterator listIt = curListOfMasterCol.begin(); listIt != curListOfMasterCol.end();
       ++listIt)
  {
    if (printL(5))
      std::cout << "consider Master Var " << listIt->first->name() << " with val "
                << listIt->second._value << std::endl;

    listIt->first->fillAggregateSol(curAggregateMastSol, listIt->second._value);
  }

  std::list<Variable *> spVarList;

  /// Record only frac spVar Value with Dfrac > param.BapCodCutViolationTolerance
  for (VarPtr2DoubleMap::iterator camIt = curAggregateMastSol.begin(); camIt != curAggregateMastSol.end(); camIt++)
    {
      if (printL(5))
        std::cout << "consider sp var " << (camIt->first)->name() << " with aggregate use " << camIt->second
                  << std::endl;

      if(!camIt->first->isCandForBranching())
        continue;
 
      if (!camIt->first->isTypeOf(VcId::InstanciatedVarMask))
        continue;

      InstanciatedVar * ivarPtr = static_cast<InstanciatedVar *> (camIt->first);
      bapcodInit().check(ivarPtr == NULL, "GenVarGenBranchConstr::branchingSeparationRoutine() "
                                          "column sp var should be of type InstanciatedVar");

      if (ivarPtr->genVarConstrPtr() != genVarPtr())
        /// Not concerned by this GenVarGenBranchConstr
        continue;

      ivarPtr->tmpVal(camIt->second);
      if (printL(5))
        std::cout << "check sp var = " << ivarPtr->name() << " use " << ivarPtr->tmpVal() << std::endl;

      if (ivarPtr->fracPart() > param().BapCodIntegralityTolerance())
        spVarList.push_back(ivarPtr);
    }

  /// create a set of candidate branching variables
  /// pass thourgh the list of variables to select the best branching candidate
  for (std::list<Variable *>::const_iterator it = spVarList.begin(); it != spVarList.end(); ++it)
    {
      /// calculate the new violation and direction
      Double varViolation(0);
      char varDirection(' ');
      brCmpPtr()->calculateVarViolAndDir(*it, varViolation, varDirection);

      generatedBrConstrGeneratorSet.insert(new GenVarBranchConstrGenerator(this, *it, varDirection, varViolation));

      if (printL(5))
        std::cout << "branchingSeparationFindCandidates : new candidate = " << (*it)->name()
                  << " use " << (*it)->tmpVal() << " frac " << varViolation
                  << " generatedBrConstrGeneratorSet.size() = " << generatedBrConstrGeneratorSet.size() << std::endl;

    if ((int) generatedBrConstrGeneratorSet.size() > maxNumOfCandidates)
      {
        if (printL(5))
          std::cout << "branchingSeparationFindCandidates : remove last  candidate of list of size = "
                    << generatedBrConstrGeneratorSet.size() << std::endl;
      
        delete *(--(generatedBrConstrGeneratorSet.end()));
        generatedBrConstrGeneratorSet.erase(--(generatedBrConstrGeneratorSet.end()));
      }
  }

  return;
}


void GenVarGenBranchConstr::branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                              const MasterColSolution & curListOfMasterCol,
                                                              const int & maxNumOfCandidates,
                                                              BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  /// @todo implment as in branchingSeparationRoutine(const std::list< MastColumn * > & curListOfMasterCol,
  if (priorityLevel() <= 0)
    return;

  std::list<InstanciatedVar *> fracCandVarList;
  for (MasterVarSolution::const_iterator it = curListOfMastAndSubprobVar.begin();
       it != curListOfMastAndSubprobVar.end(); ++it)
  {
    if (printL(5))
      std::cout << "consider Master Var " << it->first->name() << " with val " << it->second._value << std::endl;

    if (!it->first->isCandForBranching())
      continue;

    if (!it->first->isTypeOf(VcId::InstanciatedVarMask))
      continue;

    InstanciatedVar *ivarPtr = static_cast<InstanciatedVar *> (it->first);
    bapcodInit().check(ivarPtr == NULL, "GenVarGenBranchConstr::branchingSeparationRoutine() "
                                        "column sp var should be of type InstanciatedVar");

    if (ivarPtr->genVarConstrPtr() != genVarPtr())
      /// Not concerned by this GenVarGenBranchConstr
      continue;
    ivarPtr->tmpVal(it->second._value);
    fracCandVarList.push_back(ivarPtr);
  }

  /// reimplemented by Ruslan with support of maxNumOfCandidates
  /// sorting of candidates is performed in generatedBrConstrGeneratorSet using BranchingSeparationPriorityComp
  for (std::list<InstanciatedVar *>::const_iterator it = fracCandVarList.begin(); it != fracCandVarList.end(); ++it) {
    if (printL(5))
      std::cout << " Current selected Branching Variable " << (*it)->name()
                << " bestBrVarUse =  " << (*it)->tmpVal()
                << " bestBrVarViolation =  " << (*it)->fracPart()
                << " bestBrVarDir =  " << (*it)->directive() << std::endl;

    generatedBrConstrGeneratorSet.insert(new GenVarBranchConstrGenerator(this, *it, (*it)->directive(),
                                                                         (*it)->fracPart()));

    if ((int) generatedBrConstrGeneratorSet.size() > maxNumOfCandidates) {
      if (printL(6))
        std::cout << "GenVarGenBranchConstr::branchingSeparationFindCandidates(): remove last cand."
                " of list of size = " << generatedBrConstrGeneratorSet.size() << std::endl;
      delete *(--(generatedBrConstrGeneratorSet.end()));
      generatedBrConstrGeneratorSet.erase(--(generatedBrConstrGeneratorSet.end()));
    }
  }

  return;
}

bool GenVarGenBranchConstr::genericCount(const InstanciatedConstr * const iconstrPtr,
                                         const InstanciatedVar * const ivarPtr) const
{
  if (printL(5))
    std::cout << "GenVarGenBranchConstr::genericCount() constr=" << iconstrPtr->name() << "  var=" << ivarPtr->name()
              << std::endl;

  if (iconstrPtr->genVarConstrPtr() != this)
    {
      if(printL(7))
        std::cout << "NO: (iconstrPtr->genVarConstrPtr() != this)" << std::endl;
      return (false);
    }

  if (ivarPtr->genVarConstrPtr() != _genVarPtr)
    {
      if(printL(7))
        std::cout << "NO: (ivarPtr->genVarConstrPtr() != _genVarPtr)"
              << std::endl;

      return (false);
    }

  if ((iconstrPtr->id() == ivarPtr->id()))
    {
      if(printL(7))
        std::cout << "YES" << std::endl;

      return (true);
    }

  return (false);
}

const LpCoef GenVarGenBranchConstr::genericCoef(const InstanciatedConstr * const iconstrPtr,
                                                const InstanciatedVar * const ivarPtr) const
{
  if(printL(5))
    std::cout << "GenVarGenBranchConstr::genericCoef() constr=" << iconstrPtr->name()
              << "  var=" << ivarPtr->name() << std::endl;

  if(iconstrPtr->genVarConstrPtr() != this)
    return LpCoef::ZeroCoef;

  if(ivarPtr->genVarConstrPtr() != _genVarPtr)
    return LpCoef::ZeroCoef;

  if(iconstrPtr->id() == ivarPtr->id())
    return LpCoef::UnitCoef;

  return LpCoef::ZeroCoef;
}

void GenVarGenBranchConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  if(printL(6))
      std::cout << "GenVarGenBranchConstr::buildMembership " << iconstrPtr->name() << std::endl;

  iconstrPtr->presetMembership(true);

  return;
}

std::ostream& GenVarGenBranchConstr::print(std::ostream& os) const
{
  return os << "GenVarGenBranchConstr of genVar " << _genVarPtr->defaultName()
            << " with priorityLevel " << priorityLevel() << std::endl;
}

/*
 * Methods of class GenVarInstMastBranchConstr
 */

GenVarInstMastBranchConstr::GenVarInstMastBranchConstr(const IndexCell & id, GenericConstr * genBrConstrPtr,
                                                       ProbConfig * probConfPtr, InstanciatedVar * instVarPtr,
                                                       const std::string & name, const Double & rhs, const char & sense,
                                                       const char & type, const char & kind, const char & flag,
                                                       const Double & val, const Double & upperBound,
                                                       const Double & lowerBound, const char & directive,
                                                       const Double & priority) :
  InstMasterBranchingConstr(id, genBrConstrPtr, probConfPtr, name, rhs, sense, type, kind, flag, val,
                            upperBound, lowerBound, directive, priority),
  _instVarPtr(instVarPtr)
{
}

GenVarInstMastBranchConstr::~GenVarInstMastBranchConstr()
{
}

void GenVarInstMastBranchConstr::setMembership()
{
  if (printL(6))
    std::cout << "GenVarInstMastBranchConstr::setMembership() genVarConstrPtr() ="
              << genVarConstrPtr()->defaultName() << "  var=" << _instVarPtr->name() << std::endl;

  if (!buildMembershipHasBeenPerformed())
    {
      genVarConstrPtr()->buildMembership(this);
      buildMembershipHasBeenPerformed(true);
    }

  bool cumulativeCoef(false);

  includeMember(_instVarPtr, 1, cumulativeCoef);

  if (_instVarPtr->isTypeOf(VcId::SubProbVariableMask))
    {
      SubProbVariable * spVarptr = static_cast<SubProbVariable *>(_instVarPtr);
 
      if (printL(5))
        std::cout << "GenVarInstMastBranchConstr::setMembership()  spVarptr->masterColumnMember2coefMap().size() ="
	              << spVarptr->masterColumnMember2coefMap().size() << std::endl;

      for (MapMastColumnPtr2Double::const_iterator mcPt = spVarptr->masterColumnMember2coefMap().begin();
	       mcPt != spVarptr->masterColumnMember2coefMap().end(); ++mcPt)
        {
          includeMember(mcPt->first, mcPt->second, cumulativeCoef);
        }
    }

  Constraint::setMembership();

  return;
}

std::ostream & GenVarInstMastBranchConstr::print(std::ostream& os) const
{
  os << "GenVarInstMastBranchConstr" << std::endl;
  InstMasterBranchingConstr::print(os);

  return (os);
}

std::vector<std::string> GenVarInstMastBranchConstr::forDotPrint() const
{
  std::stringstream sstream;
  shortPrint(sstream);
  std::string s = sstream.str();

  s = s.substr(0, s.size() - 1); //to remove the tailing one space.

  vector<string> ret;
  ret.push_back(s);
  return ret;
}

void GenVarInstMastBranchConstr::shortPrint(std::ostream& os) const
{
  os << _instVarPtr->name();
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

bool GenVarInstMastBranchConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::GenVarInstMastBranchConstrMask, vcIdentifier);
}
