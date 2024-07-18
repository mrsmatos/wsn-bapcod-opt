/**
 *
 * This file bcMastColumnC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcModelC.hpp"
#include "bcNodeC.hpp"

#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcStabilizationInfo.hpp"

#include <limits.h>

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif /* BCP_RCSP_IS_FOUND */

using namespace std;
using namespace VcIndexStatus;

/*
 * Methods of class MastColumn
 */

MastColumn::MastColumn(MasterConf * masterConfPtr,
                       ColGenSpConf * cgSpConfPtr,
                       Solution * spSol,
                       const std::string & name) :
  AggregateVariable(spSol),
  Variable(masterConfPtr->modelPtr(), name + masterConfPtr->modelPtr()->modelMastColCnt(),
           0, 'P', 'I', 'E', BapcodInfinity, 0, 'd',
           'U', -1, 0, BapcodInfinity, 0, false, -1),
  _mcref(masterConfPtr->modelPtr()->modelMastColCnt()),
  _treatOrderId(0), _nodeTreatOrder(-1), _cgSpConfPtr(cgSpConfPtr)
{
    masterConfPtr->modelPtr()->increaseModelMastColCnt();

    mult(cgSpConfPtr->mult());

    if (_spSol != NULL)
    {
        costrhs(zero(cgSpConfPtr->fixedCost() + _spSol->cost()));

        if (cgSpConfPtr->param().SplitColIntoDissagregateSpVar())
        {
            VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
            lb(sPtr->first->globalLb());
            ub(sPtr->first->globalUb());
        }
    }

#if(_USING_HASHTABLE)
    computeHashKey();
#endif
}

MastColumn::MastColumn(const MastColumn & that) :
  AggregateVariable(that), Variable(that), _mcref(that.mcref()),
  _treatOrderId(0), _nodeTreatOrder(that._nodeTreatOrder), _cgSpConfPtr(that.cgSpConfPtr())
{
  setAggregateVariable(this);

  sense('P');
}

MastColumn::~MastColumn()
{
  if (_spSol != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
         sPtr != _spSol->solVarValMap().end(); ++sPtr)
    {
      SubProbVariable * spVarptr = dynamic_cast<SubProbVariable *> (sPtr->first);
      if (spVarptr == NULL) continue;
      spVarptr->masterColumnMember2coefMap().erase(this);
    }
  }
}

const Double & MastColumn::curCost() const
{
  if (printL(6))
    std::cout << " MastColumn::curCost() " << name()  << std::endl;

  return Variable::curCost();
}

const Double & MastColumn::costrhs() const
{
  if (printL(6))
    std::cout << " MastColumn::costrhs() " << name() << " _costrhs = " << _costrhs
              << "  _memorisedCurCost = " << _memorisedCurCost << std::endl;

  return Variable::costrhs();
}

void MastColumn::recordInMembershipOfSubProbVar()
{
  if (_spSol != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end();
         ++sPtr)
    {
      SubProbVariable * spVarptr = dynamic_cast<SubProbVariable *> (sPtr->first);
      if (spVarptr == NULL)
        continue;
      spVarptr->includeMasterColAsMember(this, sPtr->second);
    }
  }
  return;
}

const Double MastColumn::spVarVal(Variable * varPtr) const
{
  if (_spSol == NULL)
    return 0;

  VarPtr2DoubleMap::const_iterator varValIt = _spSol->solVarValMap().find(varPtr);
  if (varValIt != _spSol->solVarValMap().end())
  {
    return varValIt->second;
  }
  return 0;
}

bool MastColumn::spVarCount(Variable * varPtr) const
{
  if (_spSol == NULL)
    return (false);

  return (_spSol->solVarValMap().count(varPtr) != 0);
}

const Double MastColumn::contrib(const Double & use)
{
  Double ctb(0);
  Double coef(0);

  for (ConstVarConstrPtr2Double::iterator itm = member2coefMap().begin(); itm != member2coefMap().end(); ++itm)
  {
    if (itm->first->inCurProb())
    {
      if ((itm->first->type() == 'S') || (itm->first->type() == 'X')) /// otherwise this is a subproblem constraint or convexity constraint
      {
        continue;
      }

      Constraint * constrPtr = dynamic_cast<Constraint *> (itm->first);
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
          /// 'G'
        default:
        {
          if (coef * use > constrPtr->curRhs())
            ctb += Dmax(constrPtr->curRhs(), 0);
          else
            ctb += coef * use;
        }
      }
      if (printL(6))
        std::cout << "Var[" << name() << "] of use = " << use << " has coef "
                  << constrPtr->membCoef(this) << " in const[" << constrPtr->name() << "] of curRhs "
                  << constrPtr->curRhs() << " ctb = " << ctb << std::endl;
    }

  }

  return (ctb);
}

ColGenSpConf * MastColumn::cgSpConfPtr() const
{
  return (_cgSpConfPtr);
}

bool MastColumn::suitableToFixValue(const Double & value)
{
  bool returnValue = true;

  if (value.isZero() || (type() == 'C') || _spSol->solVarValMap().empty())
    {
      /// we cannot fix a column to zero or fix an empty column
      returnValue = false;
    }
  else
    {
      /// Here we simply verify whether fixing this column to the given value violates
      /// global bounds of its subproblem variables
      /// This test does not always detect infeasibility when several subproblem variables
      /// of this column participate in a same master constraint
      for (VarPtr2DoubleMap::const_iterator mapIt = _spSol->solVarValMap().begin();
	       mapIt != _spSol->solVarValMap().end(); mapIt++)
      {
        Double val = value * mapIt->second;
        //if ((val < mapIt->first->globalCurLb()) || (val > mapIt->first->globalCurUb()))
        /// checking of lower bound is removed because in the case of positive coeffs
        /// in the master constraints, it can declare unsuitable columns which are
        /// actually suitable
        /// TO DO: reimplement this check taking into account the coefficients of the
        /// column in the master constraint
        if (val > mapIt->first->globalCurUb())
          {
            returnValue = false;
            break;
          }
      }
      /// Here the generic verification is started, but left for the future
      /// For the moment the user can use isProperSolution() callback to check if it is suitable to fix a column
//      std::map<const Constraint * const, Double> remainingRhs;
//      std::map<const Constraint * const, Double>::iterator rhsMapIt;
//      for (VarPtr2DoubleMap::const_iterator mapIt = _spSol->solVarValMap().begin();
//           mapIt != _spSol->solVarValMap().end(); mapIt++)
//        {
//          SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(mapIt->first);
//          for (ConstVarConstrPtr2Double::const_iterator mastMapIt = spVarPtr->masterConstrMember2coefMap().begin();
//               mastMapIt != spVarPtr->masterConstrMember2coefMap().end(); ++mastMapIt)
//            {
//              const Constraint * const mastConstrPtr = static_cast<const Constraint * const>(mastMapIt->first);
//              rhsMapIt = remainingRhs.find(mastConstrPtr);
//              if (rhsMapIt->first->vcIndexStatus() != VcIndexStatus::Active)
//                continue;
//              if (rhsMapIt != remainingRhs.end())
//                rhsMapIt->second -= value * mapIt->second;
//              else
//                remainingRhs[mastConstrPtr] = mastConstrPtr->curRhs() - value * mapIt->second;
//            }
//        }
//      for (rhsMapIt = remainingRhs.begin(); rhsMapIt != remainingRhs.end(); ++rhsMapIt)
//        {
             /// here we need to calculate minSlack and maxSlack of the constraint the same way as in
             /// preprocessing
//        }

    }

  returnValue = returnValue && cgSpConfPtr()->probPtr()->isProperSolution(_spSol);

  if (printL(5)) /// TO DO : put it back to 6
    std::cout << "MastColumn::suitableToFixValue() : column " << name()
              << ", value " << value << ", result is " << returnValue << std::endl;

  return returnValue;
}

int MastColumn::maxValueInCurrentMasterProblem()
{
  int maxValue(*_cgSpConfPtr->upperBoundPtr());

  if (_spSol == NULL)
    return maxValue;

  for (VarPtr2DoubleMap::const_iterator mapIt = _spSol->solVarValMap().begin();
       mapIt != _spSol->solVarValMap().end(); mapIt++)
  {
    if (mapIt->second > 0)
    {
      int candidateValue = floor((double) mapIt->first->globalCurUb() / (double) mapIt->second);
      if (maxValue > candidateValue)
        maxValue = candidateValue;
    }
    if (mapIt->second < 0)
    {
      int candidateValue = floor((double) mapIt->first->globalCurLb() / (double) mapIt->second);
      if (maxValue > candidateValue)
        maxValue = candidateValue;
    }
  }

  if (maxValue < 0)
  {
    std::cerr << "BaPCod WARNING: maxValue is negative in MastColumn::maxValueInCurrentMasterProblem()" << std::endl;
    maxValue = 0;
  }

  return maxValue;
}

#define suitableForResidualProbPrintLevel 5

bool MastColumn::suitableForResidualProb(const Double & useLevel)
{
  /**
   * Check subproblem variable lower and upper bound constraints
   * because sp var bounds are not up to date if (_spSol != NULL)
   */
  if (_spSol != NULL)
  {
    VarPtr2DoubleMap::const_iterator sPtr;
    for (sPtr = _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end(); sPtr++)
    {
      if (printL(suitableForResidualProbPrintLevel))
        std::cout << "MastColumn(" << name() << ") sp var  " << sPtr->first->name() << " in [curLb="
                  << sPtr->first->curLb() << ",curUb = " << sPtr->first->curUb() << "] and in [curGlobLb="
                  << sPtr->first->curGlobLb() << ",curGlobUb = " << sPtr->first->curGlobUb()
                  << "] while  val = " << sPtr->second << " and useLevel = " << useLevel << std::endl;

      if (sPtr->first->curUb() < sPtr->second)
      {
        if (printL(suitableForResidualProbPrintLevel))
          std::cout << "MastColumn[" << name() << "] is unsuitable" << std::endl;
        return (false);
      }

      if (sPtr->first->curLb() > sPtr->second)
      {
        if (printL(suitableForResidualProbPrintLevel))
          std::cout << "MastColumn[" << name() << "] is unsuitable" << std::endl;
        return (false);
      }

      /// commented by Ruslan: global and local bounds for subproblem variables may be
      /// not "linked" (we can have local ub = 1 and global ub = 0) because sometimes
      /// we do not want to change bounds in the subproblem
      /// in this case, we should not check global bounds for a newly created column
      //          if (sPtr->first->curGlobUb() < sPtr->second * useLevel)
      //            {
      //              if (printL(suitableForResidualProbPrintLevel))
      //                std::cout << "MastColumn[" << name() << "] is unsuitable"
      //                    << std::endl;
      //              return (false);
      //            }

      if (printL(suitableForResidualProbPrintLevel))
        std::cout << "MastColumn[" << name() << "] sp var " << sPtr->first->name() << " verifies component bounds"
                  << std::endl;
    }
    /// we check the callback which can declare the column non-proper
    if (!cgSpConfPtr()->probPtr()->isProperSolution(_spSol))
    {
      return false;
    }
  }

  /// Check subproblem branching constraints
  ConstrIndexManager::iterator constrPtrIt;
  for (constrPtrIt = _cgSpConfPtr->probPtr()->probConstrSet().begin(VcIndexStatus::Active, 'd');
       constrPtrIt != _cgSpConfPtr->probPtr()->probConstrSet().end(VcIndexStatus::Active, 'd'); ++constrPtrIt)
    if ((*constrPtrIt)->isTypeOf(VcId::BranchingConstrBaseTypeMask))
    {
      if (printL(suitableForResidualProbPrintLevel))
      {
        Double cCoef(membCoef(*constrPtrIt));
        std::cout << "MastColumn::suitableForResidualProb(): col " << name()
                  << " spConf Branching  Constraint " << *constrPtrIt;

        std::cout << "coef " << cCoef << " useLevel " << useLevel << " violation "
                  << (*constrPtrIt)->computeViolation(cCoef * useLevel) << " violated "
                  << (*constrPtrIt)->violated(this, useLevel) << std::endl;
      }

      if (!(*constrPtrIt)->toBeUsedInPreprocessing())
        continue;

      /// Need to use base class violation subroutine and not branching constraint one
      if ((*constrPtrIt)->violated(this, useLevel))
      {
        if (printL(suitableForResidualProbPrintLevel))
          std::cout << "MastColumn[" << name() << "] is unsuitable" << std::endl;
        return (false);
      }

    }

  if (printL(suitableForResidualProbPrintLevel))
    std::cout << "MastColumn[" << name() << "] is suitable for ColGenSpConf extra Constraints" << std::endl;

  if (printL(6))
    std::cout << "MastColumn[" << name() << "] is suitable" << std::endl;

  return (true);
}

void MastColumn::fillAggregateSol(VarPtr2DoubleMap & curAggregateMastSol,
                                  const Double & val) const
{
  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
       sPtr != _spSol->solVarValMap().end(); sPtr++)
      curAggregateMastSol[sPtr->first] += sPtr->second * val;

  return;
}

void MastColumn::fillMapOfIntSpVar(VarPtr2DoubleMap & curSolMap, const Double & val) const
{
  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end();
       sPtr++)
    if (sPtr->first->isCandForBranching())
      curSolMap[sPtr->first] += sPtr->second * val;

  return;
}

void MastColumn
     ::fillSpIndicatorMap(std::map<Variable *,
                                   std::multiset< std::pair< MastColumn *, ValueRecord >,
                                                  SortMastColPerDecreasingSpVal>, Variable::PrioritySort> & curMap,
                          const ValueRecord & valRec, const VarPtrSet & candVarForCompBoundBranching)
{
  if (printL(6))
    for (VarPtrSet::const_iterator cspPt = candVarForCompBoundBranching.begin();
         cspPt != candVarForCompBoundBranching.end(); cspPt++)
      std::cout << " cand sp var " << (*cspPt)->name() << std::endl;

  /// Fills colPt list sorted in lexicographic order with regards to each spVar
  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin(); sPtr != _spSol->solVarValMap().end();
       sPtr++)
  {
    if (printL(6))
      std::cout << " sp var " << sPtr->first->name() << ", candVarForCompBoundBranching.count() = "
                << candVarForCompBoundBranching.count(sPtr->first) << std::endl;

    if (candVarForCompBoundBranching.count(sPtr->first))
    {
      if (printL(6))
        std::cout << " update sp var " << sPtr->first->name() << std::endl;


      if (curMap.find(sPtr->first) == curMap.end())
      {
        SortMastColPerDecreasingSpVal comparator(sPtr->first);
        std::multiset< std::pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal> ms(comparator);
        curMap[sPtr->first] = ms;
      }

      curMap[sPtr->first].insert(make_pair< MastColumn * const, const ValueRecord&>(this, valRec));
    }
  }
  if (printL(6))
    std::cout << " fillSpIndicatorMap col " << name() << " curMap size " << curMap.size() << std::endl;

  return;
}

void MastColumn::fillSpIndicatorColMultiset(std::multiset< std::pair< MastColumn *, ValueRecord >,
                                                           SortMastColPerDecreasingSpVal> & curMultiSet,
                                            Variable * spVarPtr, const ValueRecord & valRec)
{
  if (printL(6))
    std::cout << " fillSpIndicatorColMultiset col " << name() << std::endl;

  /// Fills colPt list sorted in lexicographic order with regards to each spVar
  if (_spSol != NULL)
  {
    if (printL(6))
      std::cout << " sp var " << spVarPtr->name() << std::endl;

    curMultiSet.insert(make_pair< MastColumn * const, const ValueRecord& >(this, valRec));
  }

  return;
}

bool MastColumn::defaultComputeCount(ConstVarConstrConstPtr vcPtr)
{
  return AggregateVariable::agvComputeCount(vcPtr);
}

const LpCoef MastColumn::defaultComputeCoef(ConstVarConstrConstPtr vcPtr)
{
  return AggregateVariable::agvComputeCoef(vcPtr);
}

/// Assumes spSol var have been set before

void MastColumn::agvSetMembership(Variable * varPtr)
{
  if (_spSol == NULL)
    return;

  bool cumulativeCoef(true);

  for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
       sPtr != _spSol->solVarValMap().end(); ++sPtr)
    {
      if (printL(6))
        std::cout << "MastColumn::agvSetMembership(): solVarValMap includes var "
                  << sPtr->first->name() << std::endl;

      if (sPtr->first->isTypeOf(VcId::SubProbVariableMask))
        {
          SubProbVariable * spvPtr =  static_cast<SubProbVariable *> (sPtr->first);

          for (ConstVarConstrPtr2Double::const_iterator itm = spvPtr->masterConstrMember2coefMap().begin();
               itm != spvPtr->masterConstrMember2coefMap().end(); ++itm)
            {
              if (printL(6))
                std::cout << "MastColumn::agvSetMembership(): var "
                          << sPtr->first->name() << " with val " << sPtr->second
                          << " in constr " << itm->first->name() << " with coef "
                          << itm->second << std::endl;

              /// if column pool is not used, we do not generate membership of the column
              /// in unsuitable dynamic constraints, as generated column will be active
              /// only at the node it was generated and in the subtree rooted at this node
              if (!param().UseColumnsPool() && (itm->first->flag() == 'd')
                  && (itm->first->vcIndexStatus() != VcIndexStatus::Active)
                  && (itm->first->vcIndexStatus() != VcIndexStatus::Inactive))
                continue;
              /// we need to verify membershipUpToDate because variables and constraints may be added simultaneously
              /// (in the column-and-row generation approach) and in this case without verification
              ///  includeMember would be called two times
              if (itm->first->membershipUpToDate())
                {
                  if (printL(6))
                    std::cout << "MastColumn::agvSetMembership(): var "
                              << sPtr->first->name() << " with val " << sPtr->second
                              << " in constr " << itm->first->name() << " with coef "
                              << itm->second << std::endl;

                  varPtr->includeMember(itm->first, itm->second * sPtr->second, cumulativeCoef);
                }

            }
        }
      else
        {
          bapcodInit().require(false,
                               "MastColumn::agvSetMembership() composing var should be of type SubProbVariable");
        }
    }

  return;
}

void MastColumn::setMembership()
{
  MastColumn::agvSetMembership(this);
  bool cumulativeCoef(false);

  if (!_cgSpConfPtr->param().SplitColIntoDissagregateSpVar())
  {

      ///  add Membership in convexity constraint
      if (_cgSpConfPtr->lowerBoundMastConstrPtr() != NULL)
      {
          LpCoef coef(computeCoef(_cgSpConfPtr->lowerBoundMastConstrPtr()));
          if (coef.first)
              includeMember(_cgSpConfPtr->lowerBoundMastConstrPtr(), coef.second, cumulativeCoef);
      }
      if (_cgSpConfPtr->upperBoundMastConstrPtr() != NULL)
      {
          LpCoef coef(computeCoef(_cgSpConfPtr->upperBoundMastConstrPtr()));
          if (coef.first)
              includeMember(_cgSpConfPtr->upperBoundMastConstrPtr(), coef.second, cumulativeCoef);
      }

    ///  add Membership in non linear  Master constraint
    for (ConstrPtrSet::const_iterator constrPt = problemPtr()->probNonLinearConstrSet().begin();
	     constrPt != problemPtr()->probNonLinearConstrSet().end(); ++constrPt)
        {
          if (printL(5))
            std::cout << " MastColumn::setMembership() call addMember for constr "
                      << (*constrPt)->name() << std::endl;

          if (!param().UseColumnsPool() && ((*constrPt)->flag() == 'd')
              && ((*constrPt)->vcIndexStatus() != VcIndexStatus::Active)
              && ((*constrPt)->vcIndexStatus() != VcIndexStatus::Inactive))
            continue;

          addMember(*constrPt);
        }
    }

  VarConstr::setMembership();

  return;
}

bool MastColumn::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "MastColumn::computeCount this " << name() << " that "
              << vcPtr->name() << std::endl;

  if (vcPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    {
        InstMastConvexityConstr * imccPtr = static_cast<InstMastConvexityConstr *> (vcPtr);
        return (imccPtr->cgSpConfPtr() == _cgSpConfPtr);
    }

  if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask))
    {
      InstanciatedConstr * icPtr = static_cast<InstanciatedConstr *> (vcPtr);
      if (printL(6))
        std::cout << "MastColumn::computeCount : InstanciatedConstr "
                  << std::endl;

      const Base4NonLinearGenericConstr * const genNLMCPtr =
        dynamic_cast<const Base4NonLinearGenericConstr * const> (icPtr->genVarConstrPtr());

      if (genNLMCPtr != NULL)
        {
          if (printL(6))
            std::cout
              << "MastColumn::computeCount : Base4NonLinearGenericConstr "
              << icPtr->name() << " - " << name() << std::endl;

          return (genNLMCPtr->genericMastColumnCount(icPtr, this));
        }
    }

  if (vcPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
    {
      Base4NonLinearConstraint * nliConstrPtr =
        dynamic_cast<Base4NonLinearConstraint *> (vcPtr);

      if (printL(6))
        std::cout << "MastColumn::computeCount : Base4NonLinearConstraint "
                  << std::endl;

      return vcPtr->computeCount(this);
    }

  return (defaultComputeCount(vcPtr));
}

const LpCoef MastColumn::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    {
        InstMastConvexityConstr * imccPtr = static_cast<InstMastConvexityConstr *> (vcPtr);
        if (imccPtr->cgSpConfPtr() == _cgSpConfPtr)
            return LpCoef::UnitCoef;
        else
            return LpCoef::ZeroCoef;
    }
  if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask))
    {
      InstanciatedConstr * icPtr = static_cast<InstanciatedConstr *> (vcPtr);

      const Base4NonLinearGenericConstr * const genNLMCPtr =
        dynamic_cast<const Base4NonLinearGenericConstr * const> (icPtr->genVarConstrPtr());

      if (genNLMCPtr != NULL)
        return (genNLMCPtr->genericMastColumnCoef(icPtr, this));
    }
  if (vcPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
    {
      Base4NonLinearConstraint * nliConstrPtr =
        dynamic_cast<Base4NonLinearConstraint *> (vcPtr);

      /// Instead of genericCount, handled directly by count (dyamic_cast<MastColumn *>)
      return (vcPtr->computeCoef(this));
    }

  return (defaultComputeCoef(vcPtr));
}

std::ostream& MastColumn::printColVector(std::ostream& os) const
{
  if (_spSol != NULL)
    for (VarPtr2DoubleMap::const_iterator sPtr = _spSol->solVarValMap().begin();
            sPtr != _spSol->solVarValMap().end(); sPtr++)
    {
      os << "   MC includes spVar[" << sPtr->first->name() << "] = "
              << sPtr->second << std::endl;
    }

  return (os);
}

std::ostream& MastColumn::print(std::ostream& os) const
{
  os << "MastColumn ref = " << _mcref << std::endl;
  if (_cgSpConfPtr != NULL)
    os << "   ColGenSpConf ref = " << _cgSpConfPtr->ref() << std::endl;
  Variable::print(os);
  printColVector(os);
  if (_spSol != NULL)
    _spSol->printOrderedSolution(os);

  return (os);
}

bool MastColumn::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::MastColumnMask, vcIdentifier);
}

bool MastColumn::operator<(const VarConstr& b) const
{
  if (!b.isTypeOf(VcId::MastColumnMask))
  {
    return VarConstr::operator<(b);
  }

  const MastColumn & bTemp = static_cast<const MastColumn&> (b);

  if (printL(6))
    std::cout << "MastColumn::LexicographicallyST()" << this->name()
    << " L< " << bTemp.name()
    << std::endl;

  if (this->cgSpConfPtr()->ref() < bTemp.cgSpConfPtr()->ref())
    return (true);

  if (this->cgSpConfPtr()->ref() > bTemp.cgSpConfPtr()->ref())
    return (false);

  if (printL(6))
    std::cout << "MastColumn::LexicographicallyST(): same ColGenSpConf" << std::endl;

  if (this->spSol() == NULL)
    return (true);

  if (bTemp.spSol() == NULL)
    return (false);

  if (printL(6))
    std::cout << "MastColumn::LexicographicallyST(): both have SpSol" << std::endl;


  VarPtr2DoubleMap::const_iterator it1(this->spSol()->solVarValMap().begin());
  VarPtr2DoubleMap::const_iterator it2(bTemp.spSol()->solVarValMap().begin());

  /// Map may not have same size
  for (; (it1 != this->spSol()->solVarValMap().end())
          && (it2 != bTemp.spSol()->solVarValMap().end()); ++it1, ++it2)
  {
    if (printL(6))
      std::cout << "MastColumn::LexicographicallyST(): comp v[" << it1->first->ref()
                << "] = " << it1->second << " with v[" << it2->first->ref() << "] = " << it2->second << std::endl;

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
    std::cout << "MastColumn::LexicographicallyST(): same sp sol so far" << std::endl;

  /// Not same size
  if (this->spSol()->solVarValMap().size() < bTemp.spSol()->solVarValMap().size())
    return (true);

  if (this->spSol()->solVarValMap().size() > bTemp.spSol()->solVarValMap().size())
    return (false);

    /// added by Ruslan : compare now the ordered vertices (routes)
    /// and their resource consumption
    /// (same "var-val collection" can have different routes
    ///  and even the same routes and different resource consumption
    ///  and thus different coefficients in the limited memory cuts)

    /// the same column may have two copies ("enumerated" and "non-enumerated")
    if (this->spSol()->enumeratedFlag() != bTemp.spSol()->enumeratedFlag() )
        return this->spSol()->enumeratedFlag() < bTemp.spSol()->enumeratedFlag();

#ifdef BCP_RCSP_IS_FOUND
    const bcp_rcsp::Solution * thisRcspSolPtr = this->spSol()->rcspSolPtr();
    const bcp_rcsp::Solution * otherRcspSolPtr = bTemp.spSol()->rcspSolPtr();

    if ((thisRcspSolPtr != NULL) || (otherRcspSolPtr != NULL))
    {
        if (thisRcspSolPtr == NULL)
            return (true);

        if (otherRcspSolPtr == NULL)
            return (false);


        auto thisArcIdIt = thisRcspSolPtr->arcIds.begin();
        auto otherArcIdIt = otherRcspSolPtr->arcIds.begin();
        while ((thisArcIdIt != thisRcspSolPtr->arcIds.end()) && (otherArcIdIt != otherRcspSolPtr->arcIds.end()))
        {
            if (*thisArcIdIt < *otherArcIdIt)
                return true;
            else if (*thisArcIdIt > *otherArcIdIt)
                return false;
            ++thisArcIdIt;
            ++otherArcIdIt;
        }
        if (thisRcspSolPtr->arcIds.size() < otherRcspSolPtr->arcIds.size())
            return true;
        else if (thisRcspSolPtr->arcIds.size() > otherRcspSolPtr->arcIds.size())
            return false;
        auto thisResConsIt = thisRcspSolPtr->resConsumption.begin();
        auto otherResConsIt = otherRcspSolPtr->resConsumption.begin();
        while ((thisResConsIt != thisRcspSolPtr->resConsumption.end())
               && (otherResConsIt != otherRcspSolPtr->resConsumption.end()))
        {
            if (*thisResConsIt < *otherResConsIt)
                return true;
            else if (*thisResConsIt > *otherResConsIt)
                return false;
            ++thisResConsIt, ++otherResConsIt;
        }
        auto thisWaitTimesIt = thisRcspSolPtr->waitingTimes.begin();
        auto otherWaitTimesIt = otherRcspSolPtr->waitingTimes.begin();
        while ((thisWaitTimesIt != thisRcspSolPtr->waitingTimes.end())
               && (otherWaitTimesIt != otherRcspSolPtr->waitingTimes.end()))
        {
            if (*thisWaitTimesIt < *otherWaitTimesIt)
                return true;
            else if (*thisWaitTimesIt > *otherWaitTimesIt)
                return false;
            ++thisWaitTimesIt, ++otherWaitTimesIt;
        }
    }
#endif /* BCP_RCSP_IS_FOUND */

    const std::vector<int> & thisOrder = this->spSol()->orderedIds();
    const std::vector<int> & otherOrder = bTemp.spSol()->orderedIds();
    if (!thisOrder.empty() || !otherOrder.empty())
    {
        std::vector<int>::const_iterator toIt(thisOrder.begin());
        std::vector<int>::const_iterator ooIt(otherOrder.begin());
        while ((toIt != thisOrder.end()) && (ooIt != otherOrder.end()))
        {
            if (*toIt < *ooIt)
                return true;
            else if (*toIt > *ooIt)
                return false;
            ++toIt, ++ooIt;
        }
        if (thisOrder.size() < otherOrder.size())
            return true;
        else if (thisOrder.size() > otherOrder.size())
            return false;
        const std::vector<std::vector<double> > & thisResConsumption = this->spSol()->resConsumption();
        const std::vector<std::vector<double> > & otherResConsumption = bTemp.spSol()->resConsumption();
        std::vector<std::vector<double> >::const_iterator trcIt(thisResConsumption.begin());
        std::vector<std::vector<double> >::const_iterator orcIt(otherResConsumption.begin());
        while ((trcIt != thisResConsumption.end()) && (orcIt != otherResConsumption.end()))
        {
            if (*trcIt < *orcIt)
                return true;
            else if (*trcIt > *orcIt)
                return false;
            ++trcIt, ++orcIt;
        }
        if (thisResConsumption.size() < otherResConsumption.size())
            return true;
        else if (thisResConsumption.size() > otherResConsumption.size())
            return false;
    }

  return (false);
}


bool SortMastColPerDecreasingSpVal::operator()(const std::pair < MastColumn *, ValueRecord > & a,
                                               const std::pair < MastColumn *, ValueRecord > & b) const
{
  return (a.first->spVarVal(_spVarPtr) > b.first->spVarVal(_spVarPtr));
}

/*
 * Methods of class MissingColumn
 */
MissingColumn::MissingColumn(MasterConf * masterConfPtr, ColGenSpConf * cgSpConfPt, const Double & artVarCost) :
    MastColumn(masterConfPtr, cgSpConfPt), ArtificialVar()
{
    bapcodInit().require((cgSpConfPtr() != NULL) && (cgSpConfPtr()->probPtr() != NULL),
                         "MissingColumn() cgSpConfPtr and probPtr should be defined");

    cgSpConfPtr()->misColPtr(this);

    Double MCcost(0);
    if (artVarCost != 0.0)
        MCcost = artVarCost;
    else {
        /// Compute worst case cost
        Double vc = 0;

        for (VarIndexManager::const_iterator sPtr = cgSpConfPtr()->probPtr()->probVarSet().begin(Active, 's');
             sPtr != cgSpConfPtr()->probPtr()->probVarSet().end(Active, 's'); sPtr++) {
            vc = (*sPtr)->curCost();
            if (vc.positive())
                MCcost += vc * (*sPtr)->curUb();
            else
                MCcost += vc * (*sPtr)->curLb();
        }

        int sense(problemPtr()->objStatus());
        /// Do not let artificial variable have zero cost;
        if (zeroTest(cgSpConfPtr()->fixedCost() + MCcost))
            costrhs(1 * sense);
        else {
            /**
             *  Otherwise it cannot be excluded from lp solution
             * by way of multiplying its cast by a constant factor.
             */
            costrhs((cgSpConfPtr()->fixedCost() + MCcost) * sense);
        }

        setCost(costrhs());
        name(std::string() + "MI" + mcref() + "sp" + cgSpConfPtr()->ref());
        type('C');

        if (cgSpConfPtr()->upperBoundPtr() != NULL) {
            ub(*cgSpConfPtr()->upperBoundPtr());
        }
    }
}

MissingColumn::~MissingColumn()
{
  return;
}

bool MissingColumn::computeCount(ConstVarConstrConstPtr vcPtr)
{
  if (printL(6))
    std::cout << "MissingColumn::computeCount this " << name() << " that " << vcPtr->name() << std::endl;

 if (vcPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    {
        InstMastConvexityConstr * imccPtr = static_cast<InstMastConvexityConstr *> (vcPtr);
        return ((imccPtr->cgSpConfPtr() == cgSpConfPtr()));
    }


  if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask))
  {
      InstanciatedConstr * icPtr = static_cast<InstanciatedConstr *> (vcPtr);

      const Base4NonLinearGenericConstr * const genNLMCPtr =
              dynamic_cast<const Base4NonLinearGenericConstr * const> (icPtr->genVarConstrPtr());
      if (genNLMCPtr != NULL)
          return ((vcPtr->sense() == 'G'));
  }

  if (vcPtr->isTypeOf(VcId::InstMasterConstrMask))
    {
        InstMasterConstr * imcPtr = static_cast<InstMasterConstr *> (vcPtr);

        for (MapSubProbVariablePtr2Double::const_iterator spvPt = imcPtr->subProbVarMember2coefMap().begin();
             spvPt != imcPtr->subProbVarMember2coefMap().end(); ++spvPt)
        {
            if (printL(7))
                std::cout << "MissingColumn::computeCount() test var " << spvPt->first->name() << std::endl;

            if (spvPt->first->cgSpConfPtr() == cgSpConfPtr())
                return (true);
        }

        if (printL(7))
            std::cout << "MissingColumn::computeCount() " << name() << " count = false " << std::endl;
    }

  return (false);
}

const LpCoef MissingColumn::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  if (vcPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    {
        InstMastConvexityConstr * imccPtr = static_cast<InstMastConvexityConstr *> (vcPtr);

        if (imccPtr->cgSpConfPtr() != cgSpConfPtr())
            return LpCoef::ZeroCoef;

        switch (imccPtr->sense())
        {
            case 'G':
                return LpCoef(cgSpConfPtr()->upperBoundPtr() == NULL ? imccPtr->curRhs()
                                                                        : *cgSpConfPtr()->upperBoundPtr());
            case 'L':
                return LpCoef(cgSpConfPtr()->lowerBoundPtr() == NULL ? imccPtr->curRhs()
                                                                        : *cgSpConfPtr()->lowerBoundPtr());
            default:
                return LpCoef(imccPtr->curRhs());
        }
    }

    if (vcPtr->sense() == 'E')
        return LpCoef(vcPtr->curRhs());

    if (vcPtr->isTypeOf(VcId::InstanciatedConstrMask))
    {
        InstanciatedConstr * icPtr = static_cast<InstanciatedConstr *> (vcPtr);
        const Base4NonLinearGenericConstr * const genNLMCPtr =
                dynamic_cast<const Base4NonLinearGenericConstr * const> (icPtr->genVarConstrPtr());

        if (genNLMCPtr != NULL)
        {
            switch (vcPtr->sense())
            {
                case 'G':
                    return LpCoef(vcPtr->curRhs());
                case 'L':
                    return LpCoef::ZeroCoef;
                default:
                    return LpCoef(vcPtr->curRhs());
            }

        }
    }
  if (vcPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
    {
      Base4NonLinearConstraint * nliConstrPtr =
	dynamic_cast<Base4NonLinearConstraint *> (vcPtr);

      switch (vcPtr->sense())
	{
	case 'G':
	  return LpCoef(vcPtr->curRhs());
	case 'L':
	  return LpCoef::ZeroCoef;
	default:
	  return LpCoef(vcPtr->curRhs());
	}
    }

  Double c(0);
  Double vc(0);

  if (vcPtr->isTypeOf(VcId::InstMasterConstrMask))
    {
      InstMasterConstr * imcPtr = static_cast<InstMasterConstr *> (vcPtr);

      for (MapSubProbVariablePtr2Double::const_iterator spvPt =
	     imcPtr->subProbVarMember2coefMap().begin();
	   spvPt != imcPtr->subProbVarMember2coefMap().end(); ++spvPt)
	{
	  if (printL(7))
	    std::cout << "MissingColumn::computeCoef() test var "
		      << spvPt->first->name() << std::endl;

	  if (spvPt->first->cgSpConfPtr() == cgSpConfPtr())
	    {
	      switch (vcPtr->sense())
		{
		case 'L':
		  {
		    vc = spvPt->first->lhsMinContrib(vcPtr);
		    break;
		  }
		case 'G':
		  {
		    vc = spvPt->first->lhsMaxContrib(vcPtr);
		    break;
		  }
		default:
		  vc = 0;
		  break;
		}

	      c += vc;
	      if (printL(7))
		std::cout << "MissingColumn::computeCoef() cur c =  " << c
			  << std::endl;
	    }
	}

      if (printL(7))
	std::cout << "MissingColumn::computeCoef() c = " << c << std::endl;

      switch (vcPtr->sense())
	{
	case 'L':
	  {
	    c = Dmax(c, 0);
	    break;
	  }
	case 'G':
	  {
	    c = Dmin(c, vcPtr->curRhs());
	    break;
	  }
	default:
	  {
	    c = Dmax(c, 0);
	    c = Dmin(c, vcPtr->curRhs());
	    break;
	  }
	}
    }

  return LpCoef(c);
}

void MissingColumn::addMember(VarConstr * vcPtr)
{
  if (printL(7))
    std::cout << "MissingColumn::addMember() " << vcPtr->name() << std::endl;

  MastColumn::addMember(vcPtr);

  return;
}

bool MissingColumn::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::MissingColumnMask, vcIdentifier);
}
