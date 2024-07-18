/**
 *
 * This file bcRyanFosterBranchConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"

/**
 * Generic code
 */
using namespace std;



//////////////////////////////////////////////////////////////////////
/********************************************************************/
/****** Methods of class RyanAndFosterBranchConstrGenerator *********/
/********************************************************************/
//////////////////////////////////////////////////////////////////////


RyanAndFosterBranchConstrGenerator::RyanAndFosterBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
                                                                       InstanciatedVar * variPtr,
                                                                       InstanciatedVar * varjPtr,
				                                                       const Double & candidateLhs,
				                                                       const char & priorityDir):
  BranchingConstrGenerator(gbcPtr, priorityDir, candidateLhs), _variPtr(variPtr), _varjPtr(varjPtr)
{
  _description = _variPtr->name() + std::string(" and ") + _varjPtr->name() + std::string(" in ")
                 + _variPtr->probConfPtr()->name();
}

RyanAndFosterBranchConstrGenerator
::RyanAndFosterBranchConstrGenerator(const RyanAndFosterBranchConstrGenerator & other):
  BranchingConstrGenerator(other), _variPtr(other._variPtr), _varjPtr(other._varjPtr)
{
}

RyanAndFosterBranchConstrGenerator::~RyanAndFosterBranchConstrGenerator()
{
}

void RyanAndFosterBranchConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curSol)

{
  _candidateLhs = 0;
  for (SolutionVarInfoPtrList::const_iterator infoIt = curSol.begin(); infoIt != curSol.end(); infoIt++)
  {
    if (!(*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
      continue;

    MastColumn * mastColumnPtr = static_cast<MastColumn *>((*infoIt)->varPtr);

    if (!mastColumnPtr->spVarCount(_variPtr))
      continue;

    if (!mastColumnPtr->spVarCount(_varjPtr))
      continue;

    _candidateLhs +=  (*infoIt)->value;
  }

  if (printL(5)) 
    std::cout << "RyanAndFosterBranchConstrGenerato on var pair(" << _variPtr->name() << ",  " << _varjPtr->name()
              << "); _candidateLhs = " << _candidateLhs << std::endl;

  return;
}

void RyanAndFosterBranchConstrGenerator
     ::instanciateBrConstr(const int & parentNodeNb,
                           const int & childNb,
                           const Double & rhs,
                           const char & sense,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
    std::string name("BCrf");
    name = name + _variPtr->ref() + "-" + _varjPtr->ref();

    InstSubProbBranchingConstr * ISPBCptr(NULL);

    switch (sense)
    {
        case 'G':
        {
            /// Add Subproblem Branching Constraint x^g_i =  x^g_j
            ISPBCptr = new RyanAndFosterInstSubProbBranchConstr(_variPtr, _varjPtr, _genericBrConstrPtr,
                                                                _variPtr->probConfPtr(),
                                                                name + "S" + "p" + parentNodeNb + "c" + childNb, 0,
                                                                'E');
            if (printL(5))
                ISPBCptr->print(std::cout);

            nextBranchingConstrPtrList.push_back(ISPBCptr);
            break;
        }
        default:
        {
            /// Add Subproblem Branching Constraint x^g_i +  x^g_j <= 1
            ISPBCptr = new RyanAndFosterInstSubProbBranchConstr(_variPtr, _varjPtr, _genericBrConstrPtr,
                                                                _variPtr->probConfPtr(),
                                                                name + "S" + "p" + parentNodeNb + "c" + childNb, 1,
                                                                'L');
            if (printL(5))
                ISPBCptr->print(cout);

            nextBranchingConstrPtrList.push_back(ISPBCptr);
            break;
        }
    }

    return;
}

bool RyanAndFosterBranchConstrGenerator
     ::nextNodeBrConstr(Node * parentNodePtr,
                        std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList,
                        const ConstrPtrSet & existingMasterBranchingConstr)
{
  /// Node branching constraint defined last is treated first (LIFO)
  nextBranchingConstrPtrList.clear();
  bool success(true);

  int ancestorNodeRef(-1);
  if (parentNodePtr != NULL)
      ancestorNodeRef = parentNodePtr->ref();

  switch (_direction)
    {
    case 'U':
      {
        if (_childNbCounter == 0) 
	        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G',
                                nextBranchingConstrPtrList);
        else if (_childNbCounter == 1) 
	        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L',
                                nextBranchingConstrPtrList);
        else
	        success = false;
        break;
      }
    default:
      {
        if (_childNbCounter == 0) 
	        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L',
                                nextBranchingConstrPtrList);
        else if (_childNbCounter == 1) 
	        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G',
                                nextBranchingConstrPtrList);
        else
            success = false;
        break;
      }
    }
  return success;
}

std::ostream & RyanAndFosterBranchConstrGenerator::print(std::ostream& os) const
{
  BranchingConstrGenerator::print(os);
  os << "RyanAndFosterBranchConstrGenerator" << std::endl;
  os << "   ivar = " <<  _variPtr->name() << std::endl;
  os << "   jvar = " <<  _varjPtr->name() << std::endl;
  os << "   candidateLhs = " <<  _candidateLhs << std::endl;
  
  return(os);
}

void RyanAndFosterBranchConstrGenerator::nicePrint(std::ostream& os) const
{
  os << "Ryan&Foster pair " << _variPtr->name() << " and " << _varjPtr->name()
     << " in " << _variPtr->probConfPtr()->name() << " (lhs=" << _candidateLhs << ")";
}

//////////////////////////////////////////////////////////////////////
/********************************************************************/
/********** Methods of class RyanAndFosterGenBranchConstr ***********/
/********************************************************************/
//////////////////////////////////////////////////////////////////////


RyanAndFosterGenBranchConstr::RyanAndFosterGenBranchConstr(Model * modelPtr,
							                               GenericVar * genVarPtr,
							                               const SelectionStrategy & priorityRule,
							                               const Double & priorityLevel) :
   GenVarGenBranchConstr(modelPtr, genVarPtr->probConfPtr(), genVarPtr, priorityRule, priorityLevel)
{
}

RyanAndFosterGenBranchConstr::~RyanAndFosterGenBranchConstr()
{
}

/// branchingSeparationRoutine()
/// input is the list of fractional columns
/// output is a BranchingConstrGenerator added to the priority sorted set generatedBrConstrGeneratorSet

void RyanAndFosterGenBranchConstr
     ::branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                         const int & maxNumOfCandidates,
                                         BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  std::list < Variable * > spVarList;

  Double candidateLhs(0), fractPart(0);

  const IndexCell2InstancVarPtrMap & varMap = genVarPtr()->indexCell2InstancVarPtrMap();

  for (IndexCell2InstancVarPtrMap::const_iterator var1it = varMap.begin(); var1it != varMap.end(); ++var1it)
    {
      bapcodInit().check((var1it->second) == NULL,
                          "RyanAndFosterGenBranchConstr::branchingSeparationRoutine() "
                          "column sp var should be of type InstanciatedVar");

      if (printL(5))
        std::cout << "check sp var 1 = " << (var1it->second)->name() << std::endl;

      if (!(var1it->second)->isCandForBranching())
        continue;

      if (_genVarPtr != (var1it->second)->genVarConstrPtr() )
        continue;

      IndexCell2InstancVarPtrMap::const_iterator var2it = var1it;
      for (++var2it; var2it != varMap.end(); ++var2it)
        {
          bapcodInit().check ((var2it->second) == NULL,
                              "RyanAndFosterGenBranchConstr::branchingSeparationRoutine() "
                              "column sp var should be of type InstanciatedVar");

          if (printL(5))
            std::cout << "check sp var 2 = " << (var2it->second)->name() << std::endl;

          if (!(var2it->second)->isCandForBranching()) continue;

          if (_genVarPtr != (var2it->second)->genVarConstrPtr() )
            continue;

          /// Compute fractionality of \sum_{g : x_i^g = x_j^g = 1} \lambda_g
          candidateLhs = 0;
          for (MasterColSolution::const_iterator mcIt = curListOfMasterCol.begin(); mcIt != curListOfMasterCol.end();
               mcIt++)
            {
              if (mcIt->first->spVarCount(var1it->second) && mcIt->first->spVarCount(var2it->second))
                candidateLhs += mcIt->second._lfracValue;
            }
          fractPart = Dfrac(candidateLhs);
          if (fractPart.isZero())
            continue;

          generatedBrConstrGeneratorSet.insert(new RyanAndFosterBranchConstrGenerator(this, var1it->second,
                                                                                      var2it->second, candidateLhs));

          if ((int) generatedBrConstrGeneratorSet.size() > maxNumOfCandidates)
            {
              if (printL(6))
                {
                  std::cout << " RyanAndFosterGenBranchConstr::branchingSeparationFindCandidates():"
                            << " remove last candidate of list of size = " << generatedBrConstrGeneratorSet.size()
                            << std::endl;
                }
              delete *(--(generatedBrConstrGeneratorSet.end()));
              generatedBrConstrGeneratorSet.erase(--(generatedBrConstrGeneratorSet.end()));
            }
        }
    }
}

bool RyanAndFosterGenBranchConstr::genericCount(const InstanciatedConstr * const iconstrPtr, 
						                        const InstanciatedVar * const ivarPtr) const
{
  return false;
}

const LpCoef RyanAndFosterGenBranchConstr::genericCoef(const InstanciatedConstr * const iconstrPtr, 
						                               const InstanciatedVar * const ivarPtr) const
{
  return LpCoef::ZeroCoef;
}


std::ostream& RyanAndFosterGenBranchConstr::print(std::ostream& os) const
{
  return(os << "RyanAndFosterGenBranchConstr" << std::endl);
}


//////////////////////////////////////////////////////////////////////////////
/****************************************************************************/
/********** Methods of class RyanAndFosterInstSubProbBranchConstr ***********/
/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////


RyanAndFosterInstSubProbBranchConstr::RyanAndFosterInstSubProbBranchConstr(InstanciatedVar * variPtr,
                                                                           InstanciatedVar * varjPtr,
									                                       GenericConstr * genConstrPtr,
                                                                           ProbConfig * probConfigPt,
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
									                                       const Double & priority):
  InstSubProbBranchingConstr(MultiIndex(variPtr->ref(),varjPtr->ref()), genConstrPtr, probConfigPt, name, costrhs,
			                 sense, type, kind, flag, val, upperBound, lowerBound, directive, priority),
  _variPtr(variPtr), _varjPtr(varjPtr)
{
}

void RyanAndFosterInstSubProbBranchConstr::setMembership()
{
  if (printL(6))
    std::cout << "RyanAndFosterInstSubProbBranchConstr::setMembership() genVarConstrPtr() ="
              << genVarConstrPtr()->defaultName() << std::endl;

  if (!buildMembershipHasBeenPerformed())
    {
      genVarConstrPtr()->buildMembership(this);
      buildMembershipHasBeenPerformed(true);
    }

  bool cumulativeCoef(false);

  includeMember(_variPtr, 1, cumulativeCoef);

  cumulativeCoef = true;

  Double varJCoeff(1);
  if (sense() == 'E')
    varJCoeff = -1;
  includeMember(_varjPtr, varJCoeff, cumulativeCoef);

  Constraint::setMembership();

  return;
}

std::ostream& RyanAndFosterInstSubProbBranchConstr::print(std::ostream& os) const
{
  os << "RyanAndFosterInstSubProbBranchConstr" << std::endl;
  os << "   ivar = " <<  _variPtr->name() << std::endl;
  os << "   jvar = " <<  _varjPtr->name() << std::endl;
  
  InstSubProbBranchingConstr::print(os);
  
  return(os);
}

void RyanAndFosterInstSubProbBranchConstr::shortPrint(std::ostream& os) const
{
  if (_sense == 'E')
    os << _variPtr->name() << " = " << _varjPtr->name() << " ";
  else
    os << _variPtr->name() << " + " << _varjPtr->name() << " <= 1 ";
  os << "in " << genConstrPtr()->probConfPtr()->name() << " ";
}

std::vector<std::string> RyanAndFosterInstSubProbBranchConstr::forDotPrint() const
{
  std::stringstream ss;
  if (_sense == 'E')
    ss << _variPtr->name() << " = " << _varjPtr->name();
  else
    ss << _variPtr->name() << " <> " << _varjPtr->name();

  std::vector<std::string> stringVector(1, std::string(ss.str()));
  return stringVector;
}

bool RyanAndFosterInstSubProbBranchConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::RyanAndFosterInstSubProbBranchConstrMask, vcIdentifier);
}
