/**
 *
 * This file bcAlg4PreprocessingOfNode.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcColGenSpConfC.hpp"

#define PreprocessingPrintLevel 3
#define PreprocessingPrintLevelPlus 5

bool Alg4PreprocessingOfNode::updateMaxSlack(VarConstr * varConstrPtr, Double const delta)
{
  Constraint * constrPtr = static_cast<Constraint *>(varConstrPtr);
  if (!constrPtr->toBeUsedInPreprocessing())
    return false;

  if (printL(PreprocessingPrintLevel) && !delta.isZero())
    std::cout << "PreprocessingBase::updateMaxSlack() Constraint " << constrPtr->name()
	          << " changed max slack from " << constrPtr->curMaxSlack()
	          << " to " << constrPtr->curMaxSlack() + delta << std::endl;

    constrPtr->incrCurMaxSlack(delta);
    if ((constrPtr->curMaxSlack() < 0) && lessOrEqualConstraint(constrPtr)) {
        if (printL(0)) {
            std::cout << "Constraint " << constrPtr->name() << " induces infeasibility" << std::endl;
        }
        return true; /// infeasible
    }
    if ((constrPtr->curMaxSlack() <= 0) && greaterConstraint(constrPtr)) {
        if (printL(PreprocessingPrintLevel)) {
            std::cout << "PreprocessingBase::updateMaxSlack() Constraint "
                    << constrPtr->name() << "  is redundant" << std::endl;
        }
        /// list of preprocessed constraints contains constraints to desactivate
        /// for now we don t desactivate because they may contain spVars of sub-problems with multiplicity > 1
        /// for such an spVar, if we assume the sp cardinality is 2 the constraint Xi >= 1
        /// (together with the sp cardinality bounds) induces an implicit
        /// constraint that is Xi1+xi2 >= 1. if we say that Xi >= 0 is redundant and we delete it,
        /// we lose in the same time the constraint Xi1+xi2 >= 1.
        //constrPtr->addToPreprocessedList();
    } //else if (delta < 0) /// this condition is needed for the case when called from computeInitialConstrsSlacks()
    //else
    if ((delta < 0) && lessOrEqualConstraint(constrPtr)) {
        constrPtr->addToPropagateList(_constrsListToPropagate);
    }
    return false;
}

bool Alg4PreprocessingOfNode::updateMinSlack(VarConstr * varConstrPtr, Double const delta)
{
    Constraint * constrPtr = static_cast<Constraint *> (varConstrPtr);
    if (!constrPtr->toBeUsedInPreprocessing())
      return false;

    if (printL(PreprocessingPrintLevel))
      std::cout << "PreprocessingBase::updateMinSlack() Constraint " << constrPtr->name()
		<< " changed min slack from " << constrPtr->curMinSlack()
		<< " to " << constrPtr->curMinSlack() + delta << std::endl;

    constrPtr->incrCurMinSlack(delta);
    if ((constrPtr->curMinSlack() > 0) && greaterOrEqualConstraint(constrPtr)) {
      if (printL(0)) {
            std::cout << "Constraint " << constrPtr->name() << " induces infeasibility" << std::endl;
        }
        return true; /// infeasible
    }
    if ((constrPtr->curMinSlack() >= 0) && lessConstraint(constrPtr)) {
        if (printL(PreprocessingPrintLevel)) {
            std::cout << "PreprocessingBase::updateMinSlack() Constraint "
                    << constrPtr->name() << " is redundant" << std::endl;
        }
        /// list of preprocessed constraints contains constraints to desactivate
        /// for now we don t desactivate because they may contain spVars of sub-problems with multiplicity > 1
        /// for such an spVar, if we assume the sp cardinality is 2 the constraint Xi >= 1
        /// (together with the sp cardinality bounds) induces an implicit
        /// constraint that is Xi1+xi2 >= 1. if we say that Xi >= 0 is redundant and we delete it,
        /// we lose in the same time the constraint Xi1+xi2 >= 1.
        //constrPtr->addToPreprocessedList();
    } //else if (delta > 0) /// this condition is needed for the case when called from computeInitialConstrsSlacks()
    //else
    if (delta > 0 && greaterOrEqualConstraint(constrPtr)) {
        constrPtr->addToPropagateList(_constrsListToPropagate);
    }
    return false;
}

ConstVarConstrPtr2Double & masterConstrMap(Variable * varPtr, bool const isSpVar) {
    if (isSpVar) {
        SubProbVariable* spVarPtr = dynamic_cast<SubProbVariable *> (varPtr);
        return spVarPtr->masterConstrMember2coefMap();
    } else {
        return varPtr->member2coefMap();
    }
}

bool Alg4PreprocessingOfNode::updateLowerBound(Variable * varPtr,
					      Double newBound,
					      Constraint * modifyingConsPtr,
					      bool const isSpVar)
{
  Double curLb = varPtr->globalCurLb();
  Double curUb = varPtr->globalCurUb();

  if (newBound > curLb) {
        if (printL(PreprocessingPrintLevel)) {
            if (modifyingConsPtr != NULL)
                std::cout << "PreprocessingBase::updateLowerBound() Constraint " << modifyingConsPtr->name()
                          << " induces a better LB for variable " << varPtr->name() << " FROM " << curLb << " TO "
                          << newBound << std::endl;
            else
                std::cout << "PreprocessingBase::updateLowerBound() The other Bounds "
                          << " induces a better LB for variable " << varPtr->name()
                          << " FROM " << curLb << " TO " << newBound << std::endl;
        }

        if (varPtr->type() != 'C') /// integer Variable
            newBound.Cceil();
        if (newBound > curUb)
        {
            if (printL(0))
              std::cout << "Variable " << varPtr->name() << " new lower bound " << newBound
                        << " induces infeasibility (ub = " << curUb << ")" << std::endl;
            return true; /// infeasible
        }
        Double diff(curLb - newBound);

        ConstVarConstrPtr2Double & constrMap = masterConstrMap(varPtr, isSpVar);

        for (ConstVarConstrPtr2Double::iterator it = constrMap.begin(); it != constrMap.end(); ++it)
            if ((it->first->vcIndexStatus() == VcIndexStatus::Active) && (!it->first->inPreprocessedList())) {
                if ((it->second < 0)) {
                    if (updateMinSlack(it->first, diff * it->second))
                        return true;
                }
                if ((it->second > 0)) {
                    if (updateMaxSlack(it->first, diff * it->second))
                        return true;
                }
            }

        varPtr->globalCurLb(newBound);
        /// list of preprocessed variables contains variables with changed bounds
        varPtr->addToPreprocessedList();
    }
  return false;
}

bool Alg4PreprocessingOfNode::updateUpperBound(Variable * varPtr,
					      Double newBound,
					      Constraint * modifyingConsPtr,
					      bool const isSpVar)
{
  Double curLb = varPtr->globalCurLb();
  Double curUb = varPtr->globalCurUb();

    if (newBound < curUb) {
        if (printL(PreprocessingPrintLevel)) {
            if (modifyingConsPtr != NULL)
                std::cout << "PreprocessingBase::updateUpperBound() Constraint " << modifyingConsPtr->name()
                          << " induces a better UB for variable " << varPtr->name() << " FROM " << curUb << " TO "
                          << newBound << std::endl;
            else
                std::cout << "PreprocessingBase::updateUpperBound() The other Bounds "
                          << " induces a better UB for variable " << varPtr->name() << " FROM " << curUb << " TO "
                          << newBound << std::endl;
        }

        if (varPtr->type() != 'C') /// integer Variable
            newBound.Cfloor();
        if (newBound < curLb) {
            if (printL(0))
              std::cout << "PreprocessingBase::updateUpperBound(): variable " << varPtr->name() << " new upper bound "
                        << newBound << " induces infeasibility (ub = " << curLb << ")" << std::endl;
            return true; /// infeasible
        }
        Double diff(curUb - newBound);

        ConstVarConstrPtr2Double & constrMap = masterConstrMap(varPtr, isSpVar);
        for (ConstVarConstrPtr2Double::iterator it = constrMap.begin(); it != constrMap.end(); ++it)
        {
            if ((it->first->vcIndexStatus() == VcIndexStatus::Active)
                    && (!it->first->inPreprocessedList())) {
                if ((it->second > 0)) {
                    if (updateMinSlack(it->first, diff * it->second))
                        return true;
                }
                if ((it->second < 0)) {
                    if (updateMaxSlack(it->first, diff * it->second))
                        return true;
                }
            }
        }

        varPtr->globalCurUb(newBound);
        /// list of preprocessed variables contains variables with changed bounds
        varPtr->addToPreprocessedList();
    }
    return false;
}

bool Alg4PreprocessingOfNode::updateLocalLowerBound(SubProbVariable * varPtr)
{
    if (varPtr->probConfPtr()->upperBoundPtr() != NULL)
      {
        Double spUb = *(varPtr->probConfPtr()->upperBoundPtr());

        if ((spUb > 0) && updateLocalLowerBound(varPtr, varPtr->globalCurLb() - (spUb - 1) * varPtr->localCurUb()))
          return true;
      }
    return false;
}

bool Alg4PreprocessingOfNode::updateGlobalLowerBound(SubProbVariable * varPtr, const Double & newBound,
						                             Constraint * modifyingConsPtr)
{
  if (updateLowerBound(varPtr, newBound, modifyingConsPtr, true))
        return true;
  return updateLocalLowerBound(varPtr);
}

bool Alg4PreprocessingOfNode::updateLocalUpperBound(SubProbVariable * varPtr)
{
    if (varPtr->probConfPtr()->lowerBoundPtr() != NULL)
      {
        Double spUb = *(varPtr->probConfPtr()->upperBoundPtr());
        Double spLb = *(varPtr->probConfPtr()->lowerBoundPtr());
        if (spLb == 0)
            spLb = 1;

        if ((spUb >= 1)
            && updateLocalUpperBound(varPtr, varPtr->globalCurUb() - (spLb - 1) * varPtr->localCurLb()))
          return true;
      }
    return false;
}

bool Alg4PreprocessingOfNode::updateGlobalUpperBound(SubProbVariable * varPtr, const Double & newBound,
						                             Constraint * modifyingConsPtr)
{
  if (updateUpperBound(varPtr, newBound, modifyingConsPtr, true))
        return true;
  return updateLocalUpperBound(varPtr);
}


bool Alg4PreprocessingOfNode::updateLocalUpperBound(SubProbVariable * varPtr,
						   Double newBound,
						   Constraint * modifyingConsPtr)
{
    if (!(param().PreprocessVariablesLocalBounds())) {
        return false;
    }
    if (newBound < varPtr->localCurUb()) {
        if (printL(PreprocessingPrintLevel)) {
            if (modifyingConsPtr != NULL)
                std::cout << "PreprocessingBase::updateLocalUpperBound() Constraint " << modifyingConsPtr->name()
                          << " induces a better local UB for variable " << varPtr->name() << " FROM "
                          << varPtr->localCurUb() << " TO " << newBound << std::endl;
            else
                std::cout << "PreprocessingBase::updateLocalUpperBound() The other Bounds "
                          << " induces a better local UB for variable " << varPtr->name() << " FROM "
                          << varPtr->localCurUb() << " TO " << newBound << std::endl;
        }

        if (varPtr->type() != 'C') /// integer Variable
            newBound.Cfloor();
        if (newBound < varPtr->localCurLb())
        {
            if (printL(0))
              std::cout << "PreprocessingBase::updateLocalUpperBound(): variable " << varPtr->name()
                        << " new local upper bound " << newBound << " induces infeasibility (lb = "
                        << varPtr->localCurLb() << ")" << std::endl;
            return true; /// infeasible
        }
        Double diff(varPtr->localCurUb() - newBound);

        for (ConstVarConstrPtr2Double::iterator it = varPtr->member2coefMap().begin();
             it != varPtr->member2coefMap().end(); ++it)
            /// only subproblem constraints are considered here
            if ((!it->first->isTypeOf(VcId::MasterConstrMask)) && (it->first->vcIndexStatus() == VcIndexStatus::Active)
                    && (!it->first->inPreprocessedList())) {
                if ((it->second > 0)) {
                    if (updateMinSlack(it->first, diff * it->second))
                        return true;
                }
                if ((it->second < 0)) {
                    if (updateMaxSlack(it->first, diff * it->second))
                        return true;
                }
            }

        varPtr->localCurUb(newBound);
        varPtr->addToPreprocessedList();

        if (varPtr->probConfPtr()->upperBoundPtr() != NULL) {
            Double spUb = *(varPtr->probConfPtr()->upperBoundPtr());

            if (updateGlobalUpperBound(varPtr, varPtr->localCurUb() * spUb))
                return true;

            if (updateLocalLowerBound(varPtr, varPtr->globalCurLb() - (spUb - 1) * varPtr->localCurUb()))
                return true;
        }
    }
    return false;
}

bool Alg4PreprocessingOfNode::updateLocalLowerBound(SubProbVariable * varPtr,
						   Double newBound,
						   Constraint * modifyingConsPtr)
{
    if (!(param().PreprocessVariablesLocalBounds())) {
        return false;
    }
    if (newBound > varPtr->localCurLb()) {
        if (printL(PreprocessingPrintLevel)) {
            if (modifyingConsPtr != NULL)
                std::cout << "PreprocessingBase::updateLocalLowerBound() Constraint " << modifyingConsPtr->name()
                          << " induces a better local LB for variable " << varPtr->name() << " FROM "
                          << varPtr->localCurLb() << " TO " << newBound << std::endl;
            else
                std::cout << "PreprocessingBase::updateLocalLowerBound() The other Bounds induce a better local LB"
                          << " for variable " << varPtr->name() << " FROM " << varPtr->localCurLb()
                          << " TO " << newBound << std::endl;
        }

        if (varPtr->type() != 'C') /// integer Variable
            newBound.Cceil();
        if (newBound > varPtr->localCurUb())
        {
            if (printL(0))
              std::cout << "PreprocessingBase::updateLocalLowerBound(): variable " << varPtr->name()
                        << " new local lower bound " << newBound << " induces infeasibility (ub = "
                        << varPtr->localCurUb() << ")" << std::endl;
            return true; /// infeasible
        }
        Double diff(varPtr->localCurLb() - newBound);

        for (ConstVarConstrPtr2Double::iterator it = varPtr->member2coefMap().begin();
             it != varPtr->member2coefMap().end(); ++it)
            /// only subproblem constraints are considered here
            if ((!it->first->isTypeOf(VcId::MasterConstrMask)) && (it->first->vcIndexStatus() == VcIndexStatus::Active)
                    && (!it->first->inPreprocessedList())) {
                if ((it->second < 0)) {
                    if (updateMinSlack(it->first, diff * it->second))
                        return true;
                }
                if ((it->second > 0)) {
                    if (updateMaxSlack(it->first, diff * it->second))
                        return true;
                }
            }

        varPtr->localCurLb(newBound);
        varPtr->addToPreprocessedList();

        if (varPtr->probConfPtr()->lowerBoundPtr() != NULL) {
            Double spLb = *(varPtr->probConfPtr()->lowerBoundPtr());

            if (updateGlobalLowerBound(varPtr, varPtr->localCurLb() * spLb))
                return true;

            if (updateLocalUpperBound(varPtr, varPtr->globalCurUb() - (spLb - 1) * varPtr->localCurLb()))
                return true;
        }
    }
    return false;
}

bool Alg4PreprocessingOfNode::propagate()
{
    while (!_constrsListToPropagate.empty())
    {
        Constraint * constrPtr = _constrsListToPropagate.front();
        _constrsListToPropagate.pop_front();
        constrPtr->removeFromPropagateList();
        if (constrPtr->isTypeOf(VcId::MasterConstrMask))
        {
            for (ConstVarConstrPtr2Double::iterator it = constrPtr->member2coefMap().begin();
                 it != constrPtr->member2coefMap().end(); ++it) {
                if (it->first->isTypeOf(VcId::VariableMask) && it->first->vcIndexStatus() == VcIndexStatus::Active)
                {
                    Variable* varPtr = static_cast<Variable *> (it->first);

                    /// we consider only pure master variables,
                    /// master columns and artificial variables are excluded
                    if (varPtr->isTypeOf(VcId::InstMasterVarMask))
                    {
                        if (it->second > 0) {
                            if (lessOrEqualConstraint(constrPtr)) {
                                if (updateUpperBound(varPtr,
                                                     (constrPtr->curMaxSlack() + it->second * varPtr->globalCurLb())
                                                     / it->second, constrPtr))
                                    return true;
                            }
                            if (greaterOrEqualConstraint(constrPtr))
                            {
                                if (updateLowerBound(varPtr,
                                                     (constrPtr->curMinSlack() + it->second * varPtr->globalCurUb())
                                                     / it->second, constrPtr))

                                    return true;
                            }
                        } else /// it->second < 0
                        {
                            if (lessOrEqualConstraint(constrPtr)) {
                                if (updateLowerBound(varPtr,
                                                     (constrPtr->curMaxSlack() + it->second * varPtr->globalCurUb())
                                                     / it->second, constrPtr))
                                    return true;

                            }
                            if (greaterOrEqualConstraint(constrPtr)) {
                                if (updateUpperBound(varPtr,
                                                     (constrPtr->curMinSlack() + it->second * varPtr->globalCurLb())
                                                     / it->second, constrPtr))
                                    return true;
                            }
                        }
                    }
                }
            }
            InstMasterConstr* mastConstrPtr = static_cast<InstMasterConstr*> (constrPtr);
            for (MapSubProbVariablePtr2Double::iterator it = mastConstrPtr->subProbVarMember2coefMap().begin();
                 it != mastConstrPtr->subProbVarMember2coefMap().end(); ++it)
            {
                if (it->first->isTypeOf(VcId::SubProbVariableMask)
                    && (it->first->vcIndexStatus() == VcIndexStatus::Active))
                {
                    SubProbVariable* varPtr = static_cast<SubProbVariable *> (it->first);
                    if (it->second > 0) {
                        if (lessOrEqualConstraint(constrPtr)) {
                            if (updateGlobalUpperBound(varPtr,
                                                       (constrPtr->curMaxSlack() + it->second * varPtr->globalCurLb())
                                                       / it->second , constrPtr))
                                return true;
                        }
                        if (greaterOrEqualConstraint(constrPtr)) {
                            if (updateGlobalLowerBound(varPtr,
                                                       (constrPtr->curMinSlack() + it->second * varPtr->globalCurUb())
                                                       / it->second, constrPtr))
                                return true;
                        }
                    } else /// it->second < 0
                    {
                        if (lessOrEqualConstraint(constrPtr)) {
                            if (updateGlobalLowerBound(varPtr,
                                                       (constrPtr->curMaxSlack() + it->second * varPtr->globalCurUb())
                                                       / it->second, constrPtr))
                                return true;
                        }
                        if (greaterOrEqualConstraint(constrPtr)) {
                            if (updateGlobalUpperBound(varPtr,
                                                       (constrPtr->curMinSlack() + it->second * varPtr->globalCurLb())
                                                       / it->second, constrPtr))
                                return true;
                        }
                    }

                }
            }
        } else /// constrPtr is a subproblem constraint
        {
            if( *(constrPtr->probConfPtr()->upperBoundPtr()) <= 0)
            {
              continue;
            }

            for (ConstVarConstrPtr2Double::iterator it = constrPtr->member2coefMap().begin();
                  it != constrPtr->member2coefMap().end(); ++it) {
                if (it->first->isTypeOf(VcId::SubProbVariableMask)
                    && it->first->vcIndexStatus() == VcIndexStatus::Active)
                {
                    SubProbVariable* varPtr = static_cast<SubProbVariable *> (it->first);

                    if (it->second > 0) {
                        if (lessOrEqualConstraint(constrPtr))
                        {
                            if (updateLocalUpperBound(varPtr,
                                                      (constrPtr->curMaxSlack() + it->second * varPtr->localCurLb())
                                                      / it->second, constrPtr))
                                return true;
                        }
                        if (greaterOrEqualConstraint(constrPtr)) {

                            if (updateLocalLowerBound(varPtr,
                                                      (constrPtr->curMinSlack() + it->second * varPtr->localCurUb())
                                                      / it->second, constrPtr))
                                return true;
                        }
                    } else /// it->second < 0
                    {
                        if (lessOrEqualConstraint(constrPtr)) {
                            if (updateLocalLowerBound(varPtr,
                                                      (constrPtr->curMaxSlack()
                                                       + it->second * varPtr->localCurUb())
                                                      / it->second, constrPtr))
                                return true;
                        }
                        if (greaterOrEqualConstraint(constrPtr)) {
                            if (updateLocalUpperBound(varPtr,
                                                      (constrPtr->curMinSlack() + it->second * varPtr->localCurLb())
                                                      / it->second, constrPtr))
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool Alg4PreprocessingOfNode::columnBecameUnsuitable_SolValNotWithinBounds(MastColumn * colPtr)
{
    Solution * solPtr = colPtr->spSol();

    if (solPtr != NULL)
      {
        for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin();
             it != solPtr->solVarValMap().end(); ++it)
          {
            SubProbVariable * spVarPtr = dynamic_cast<SubProbVariable *> (it->first);
            bool withinBounds = true;
            if (param().GenerateProperColumns()
                && ((it->second < spVarPtr->localCurLb()) || (it->second > spVarPtr->localCurUb())))
              withinBounds = false;
            /// if subproblem solver cannot take into account all bounds, still it can take into account
            /// zero upper bound, so we filter columns which do not respect zero upper bound of an sp. variable
            if (!param().GenerateProperColumns() && spVarPtr->localCurUb().isZero() && !it->second.isZero())
              withinBounds = false;

            if (!withinBounds)
              {
                if (printL(PreprocessingPrintLevel))
                  {
                    std::cout << "PreprocessingBase::columnBecameUnsuitable_SolValNotWithinBounds() : var "
                              << spVarPtr->name() << "[" << spVarPtr->localCurLb() << ", " << spVarPtr->localCurUb()
                              << "] ( " "[" << spVarPtr->globalCurLb() << ", " << spVarPtr->globalCurUb()
                              << "]) in column " << colPtr->name() << std::endl;
                  }
                return true;
              }
          }
      }
    return false;
}

bool Alg4PreprocessingOfNode::columnBecameUnsuitable_InexistantNonZeroVar(MastColumn * colPtr)
{
    Solution * solPtr = colPtr->spSol();
    if (_pcPtrToNonZeroVarPtrMap.find(solPtr->probConfPtr()) != _pcPtrToNonZeroVarPtrMap.end())
    {
        for (std::list<Variable*>::iterator it = _pcPtrToNonZeroVarPtrMap[solPtr->probConfPtr()].begin();
             it != _pcPtrToNonZeroVarPtrMap[solPtr->probConfPtr()].end(); it++)
        {
            if (printL(PreprocessingPrintLevel))
            {
                std::cout << "PreprocessingBase::columnBecameUnsuitable_InexistantNonZeroVar() : check var "
                          << (*it)->name() << " in column " << colPtr->name() << std::endl;
            }
            if (solPtr->solVarValMap().find(*it) == solPtr->solVarValMap().end())
            {
                if (printL(PreprocessingPrintLevel))
                {
                    std::cout << "PreprocessingBase::columnBecameUnsuitable() : inexsitant var "
                              << (*it)->name() << " in column " << colPtr->name() << std::endl;
                }
                return true;
            }
        }
    }

    return false;
}

void Alg4PreprocessingOfNode::applyPreprocessingListsInProbAndForm(bool const initialPreprocessing) {

    Problem * mastProbPtr = _problemPts.front();
    mastProbPtr->updateConstrRhsInForm(_mastConstrsToChangeRhs);
    _mastConstrsToChangeRhs.clear();

    for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
         ++probPtrIt)
    {
        /// bounds of preprocessed variables are changed: we desactivate them only if they are fixed to zero
        /// for each varaible with changed bounds, we verify whether columns which contain this variable are
        /// still suitable columns to verify are inserted to the same list probPtr->preprocessedVarsList() at the end
        VarPtrList::iterator varPtrIt;
        VarPtrList varsToChangeBounds;
        VarPtrList varsToRemoveFromForm;

        Problem * probPtr = *probPtrIt;
        /// preprocessed constraints should be desactivated
        ConstrPtrList::iterator constrPtrIt;
        for (constrPtrIt = probPtr->preprocessedConstrsList().begin();
             constrPtrIt != probPtr->preprocessedConstrsList().end(); ++constrPtrIt)
        {
            if (((*constrPtrIt)->flag() == 'd') || initialPreprocessing)
            {
                probPtr->probConstrSet().insert(*constrPtrIt, VcIndexStatus::Unsuitable);
                deactivateLocalArtVarsOfConstr(probPtr, *constrPtrIt, VcIndexStatus::Unsuitable, varsToRemoveFromForm);
            } else {
                probPtr->probConstrSet().insert(*constrPtrIt, VcIndexStatus::Inactive);
                deactivateLocalArtVarsOfConstr(probPtr, *constrPtrIt, VcIndexStatus::Inactive, varsToRemoveFromForm);
            }
            (*constrPtrIt)->desactivate();
        }
        probPtr->delConstrsSimplyInForm(probPtr->preprocessedConstrsList());

        for (varPtrIt = probPtr->preprocessedVarsList().begin(); varPtrIt != probPtr->preprocessedVarsList().end();
             ++varPtrIt)
        {
            if (initialPreprocessing)
            {
                /// update "permanent" bounds
                (*varPtrIt)->globalLb((*varPtrIt)->globalCurLb());
                (*varPtrIt)->globalUb((*varPtrIt)->globalCurUb());
                if ((*varPtrIt)->isTypeOf(VcId::SubProbVariableMask)) {
                    SubProbVariable * spVarPtr = static_cast<SubProbVariable *> (*varPtrIt);
                    spVarPtr->localLb(spVarPtr->localCurLb());
                    spVarPtr->localUb(spVarPtr->localCurUb());
                }
            }
            varsToChangeBounds.push_back(*varPtrIt);

            if ((*varPtrIt)->isTypeOf(VcId::SubProbVariableMask))
            {
                /// add all suitable columns containing this varaible to the preprocessed list
                for (VarPtrSet::const_iterator aggVarIt = (*varPtrIt)->setOfAggregateVarToWhichItBelongs().begin();
                     aggVarIt != (*varPtrIt)->setOfAggregateVarToWhichItBelongs().end(); ++aggVarIt)
                    if ((*aggVarIt)->vcIndexStatus() != VcIndexStatus::Unsuitable)
                        (*aggVarIt)->addToPreprocessedList();

                SubProbVariable * spVarPtr = static_cast<SubProbVariable *> (*varPtrIt);

                //The follwing if bloc added by Issam. TO DO (by Ruslan): check it.
                if (spVarPtr->localCurLb() > 0 || spVarPtr->localCurUb() < 0)
                {
                    _pcPtrToNonZeroVarPtrMap[(*varPtrIt)->probConfPtr()].push_back(*varPtrIt);
                }
            }
        }

        probPtr->delVarsSimplyInForm(varsToRemoveFromForm);
        probPtr->updateBoundsInForm(varsToChangeBounds);
    }

    /// For columns which have been added to to preprocessedVarList
    /// We need to remove them from the formulation.
    VarPtrList varsToRemoveFromForm;

    for (VarPtrList::iterator varPtrIt = mastProbPtr->preprocessedVarsList().begin();
         varPtrIt != mastProbPtr->preprocessedVarsList().end(); ++varPtrIt)
    {
        if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask)
            && columnBecameUnsuitable_SolValNotWithinBounds(static_cast<MastColumn *> (*varPtrIt)))
            {
                /// column can be inactive, so we verify whether we need to desactivate it
                if ((*varPtrIt)->vcIndexStatus() == VcIndexStatus::Active)
                {
                    (*varPtrIt)->desactivate();
                    varsToRemoveFromForm.push_back(*varPtrIt);
                }
                mastProbPtr->probVarSet().insert(*varPtrIt, VcIndexStatus::Unsuitable);
        }
    }

    if (!_pcPtrToNonZeroVarPtrMap.empty() && param().GenerateProperColumns())
    {
        for (VarIndexManager::iterator varPtrIt = mastProbPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
             varPtrIt != mastProbPtr->probVarSet().end(VcIndexStatus::Active, 'd');)
        {
            if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask)
                && columnBecameUnsuitable_InexistantNonZeroVar(static_cast<MastColumn *> (*varPtrIt)))
            {
                Variable* varPtr = *varPtrIt;
                ++varPtrIt;
                varPtr->desactivate();
                varsToRemoveFromForm.push_back(varPtr);
                mastProbPtr->probVarSet().insert(varPtr, VcIndexStatus::Unsuitable);
            } else {
                ++varPtrIt;
            }

        }

        for (VarIndexManager::iterator varPtrIt = mastProbPtr->probVarSet().begin(VcIndexStatus::Inactive, 'd');
             varPtrIt != mastProbPtr->probVarSet().end(VcIndexStatus::Inactive, 'd');)
        {
            if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask)
                 && columnBecameUnsuitable_InexistantNonZeroVar(static_cast<MastColumn *> (*varPtrIt)))
            {
                Variable* varPtr = *varPtrIt;
                ++varPtrIt;
                mastProbPtr->probVarSet().insert(varPtr, VcIndexStatus::Unsuitable);
            } else {
                ++varPtrIt;
            }
        }
    }
    _pcPtrToNonZeroVarPtrMap.clear();

    if (!_colGenSpPtsWithZeroUb.empty()) {
        for (VarIndexManager::iterator varPtrIt = mastProbPtr->probVarSet().begin(VcIndexStatus::Active, 'd');
             varPtrIt != mastProbPtr->probVarSet().end(VcIndexStatus::Active, 'd');) {
            if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            {
                MastColumn * colPtr = static_cast<MastColumn *> (*varPtrIt);
                if (_colGenSpPtsWithZeroUb.find(colPtr->cgSpConfPtr()) != _colGenSpPtsWithZeroUb.end()) {
                    ++varPtrIt;
                    colPtr->desactivate();
                    varsToRemoveFromForm.push_back(colPtr);
                    mastProbPtr->probVarSet().insert(colPtr, VcIndexStatus::Unsuitable);
                    continue;
                }
            }
            ++varPtrIt;
        }

        for (VarIndexManager::iterator varPtrIt = mastProbPtr->probVarSet().begin(VcIndexStatus::Inactive, 'd');
             varPtrIt != mastProbPtr->probVarSet().end(VcIndexStatus::Inactive, 'd');)
        {
            if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            {
                MastColumn * colPtr = static_cast<MastColumn *> (*varPtrIt);
                if (_colGenSpPtsWithZeroUb.find(colPtr->cgSpConfPtr()) != _colGenSpPtsWithZeroUb.end()) {
                    ++varPtrIt;
                    mastProbPtr->probVarSet().insert(colPtr, VcIndexStatus::Unsuitable);
                    continue;
                }
            }
            ++varPtrIt;
        }

    }
    _colGenSpPtsWithZeroUb.clear();

    mastProbPtr->delVarsSimplyInForm(varsToRemoveFromForm);
}

void Alg4PreprocessingOfNode::clearPreprocessingLists()
{
    for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
         ++probPtrIt)
        (*probPtrIt)->clearPreprocessingLists();
}

void Alg4PreprocessingOfNode::deactivateLocalArtVarsOfConstr(Problem * probPtr, Constraint * constrPtr,
                                                             const VcStatus & status, VarPtrList & varsToRemoveFromForm)
{
  deactivateLocalArtVar(probPtr, constrPtr->posLocalArtVarPtr(), status, varsToRemoveFromForm);
  deactivateLocalArtVar(probPtr, constrPtr->negLocalArtVarPtr(), status, varsToRemoveFromForm);

  if (constrPtr->stabInfoPtr() != NULL)
    {
      deactivateLocalArtVar(probPtr, constrPtr->stabInfoPtr()->negInnerArtVarPtr(), status, varsToRemoveFromForm);
      deactivateLocalArtVar(probPtr, constrPtr->stabInfoPtr()->posInnerArtVarPtr(), status, varsToRemoveFromForm);
      deactivateLocalArtVar(probPtr, constrPtr->stabInfoPtr()->negOuterArtVarPtr(), status, varsToRemoveFromForm);
      deactivateLocalArtVar(probPtr, constrPtr->stabInfoPtr()->posOuterArtVarPtr(), status, varsToRemoveFromForm);
    }
}

void Alg4PreprocessingOfNode::deactivateLocalArtVar(Problem * probPtr, LocalArtificialVar * artVarPtr,
                                                    const VcStatus & status, VarPtrList & varsToRemoveFromForm)
{
  if (artVarPtr != NULL)
    {
      artVarPtr->desactivate();
      probPtr->probVarSet().insert(artVarPtr, status);
      if (printL(PreprocessingPrintLevel))
        std::cout << "Local artificial variable " << artVarPtr->name() << " is deactivated" << std::endl;
      varsToRemoveFromForm.push_back(artVarPtr);
    }
}

bool Alg4PreprocessingOfNode::initialUpdateOfSpVarBounds()
{
    bool isMasterProblem = true;
    for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
         ++probPtrIt)
    {
        //the first problem is the not considered here as it is the master
        if (isMasterProblem) {
            isMasterProblem = false;
            continue;
        }

        Problem* subProblemPtr = *probPtrIt;

        // we only consider static variables.
        for (VarIndexManager::iterator varPtrIt = subProblemPtr->probVarSet().begin(Active, 's');
             varPtrIt != subProblemPtr->probVarSet().end(Active, 's'); ++varPtrIt)
        {
            SubProbVariable* spVarPtr = static_cast<SubProbVariable*> (*varPtrIt);

            Double spLb = *(spVarPtr->probConfPtr()->lowerBoundPtr());
            Double spUb = *(spVarPtr->probConfPtr()->upperBoundPtr());

            if (updateGlobalLowerBound(spVarPtr, spVarPtr->localCurLb() * spLb, NULL))
                return true;
            if (updateGlobalUpperBound(spVarPtr, spVarPtr->localCurUb() * spUb, NULL))
                return true;
        }
    }

    return false;
}

bool Alg4PreprocessingOfNode::computeInitialConstrsSlacks()
{
    for (ConstrPtrList::iterator constrPtrIt = _constrsListToPropagate.begin();
         constrPtrIt != _constrsListToPropagate.end(); ++constrPtrIt)
      {
        Double minSlack = (*constrPtrIt)->costrhs();
        Double maxSlack = (*constrPtrIt)->costrhs();
        if (printL(PreprocessingPrintLevelPlus))
            std::cout  << "PreprocessingBase::computeInitialConstrsSlacks() treating constr " << (*constrPtrIt)->name()
                       << std::endl;
        if ((*constrPtrIt)->isTypeOf(VcId::MasterConstrMask))
        {
            for (ConstVarConstrPtr2Double::iterator it = (*constrPtrIt)->member2coefMap().begin();
                  it != (*constrPtrIt)->member2coefMap().end(); ++it)
            {
                if (printL(PreprocessingPrintLevelPlus)) {
                    std::cout << "PreprocessingBase::computeInitialConstrsSlacks() checking mast var "
                              << it->first->name() << std::endl;
                }
                if (it->first->isTypeOf(VcId::InstMasterVarMask)) {
                    if (printL(PreprocessingPrintLevelPlus)) {
                        std::cout << "PreprocessingBase::computeInitialConstrsSlacks() LB = "
                                  << it->first->lb() << std::endl;
                        std::cout << "PreprocessingBase::computeInitialConstrsSlacks() UB = "
                                  << it->first->ub() << std::endl;
                    }
                    if (it->second > 0) {
                        maxSlack -= it->second * it->first->lb();
                        minSlack -= it->second * it->first->ub();
                    } else {
                        maxSlack -= it->second * it->first->ub();
                        minSlack -= it->second * it->first->lb();
                    }
                }
            }
            InstMasterConstr* mastConstrPtr = dynamic_cast<InstMasterConstr*> (*constrPtrIt);
            for (MapSubProbVariablePtr2Double::iterator it = mastConstrPtr->subProbVarMember2coefMap().begin();
                 it != mastConstrPtr->subProbVarMember2coefMap().end(); ++it) {
                if (printL(PreprocessingPrintLevelPlus))
                    std::cout << "PreprocessingBase::computeInitialConstrsSlacks() checking sp var "
                              << it->first->name() << std::endl;
                if (it->first->isTypeOf(VcId::SubProbVariableMask)) //Issam: This if TEST might be removed i guess
                {
                    if (printL(PreprocessingPrintLevelPlus)) {
                        std::cout << "PreprocessingBase::computeInitialConstrsSlacks() global LB = "
                                  << it->first->globalLb() << std::endl;
                        std::cout << "PreprocessingBase::computeInitialConstrsSlacks() global UB = "
                                  << it->first->globalUb() << std::endl;
                    }
                    if (it->second > 0) {
                        maxSlack -= it->second * it->first->globalLb();
                        minSlack -= it->second * it->first->globalUb();
                    } else {
                        maxSlack -= it->second * it->first->globalUb();
                        minSlack -= it->second * it->first->globalLb();
                    }
                }
            }
        } else /// subproblem constraint
        {
            for (ConstVarConstrPtr2Double::iterator it = (*constrPtrIt)->member2coefMap().begin();
                 it != (*constrPtrIt)->member2coefMap().end(); ++it)
                if (it->second > 0) {
                    maxSlack -= it->second * it->first->lb();
                    minSlack -= it->second * it->first->ub();
                } else {
                    maxSlack -= it->second * it->first->ub();
                    minSlack -= it->second * it->first->lb();
                }
        }
        (*constrPtrIt)->minSlack(minSlack);
        (*constrPtrIt)->maxSlack(maxSlack);

        if (printL(PreprocessingPrintLevel)) {
            std::cout << "PreprocessingBase::computeInitialConstrsSlacks() minSlack = " << minSlack << std::endl;
            std::cout << "PreprocessingBase::computeInitialConstrsSlacks() maxSlack = " << maxSlack << std::endl;
        }

        if (updateMaxSlack(*constrPtrIt, 0) || updateMinSlack(*constrPtrIt, 0))
            return true;
    }
    return (false);
}

bool Alg4PreprocessingOfNode::exitWhenInfeasible()
{
    Problem * mastProbPtr = _problemPts.front();
    mastProbPtr->updateConstrRhsInForm(_mastConstrsToChangeRhs);
    _mastConstrsToChangeRhs.clear();

    for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
         ++probPtrIt)
    {
        (*probPtrIt)->updateBoundsInForm((*probPtrIt)->preprocessedVarsList());
    }
    clearPreprocessingLists();

    for (ConstrPtrList::iterator constrPtrIt = _constrsListToPropagate.begin();
         constrPtrIt != _constrsListToPropagate.end(); ++constrPtrIt)
        (*constrPtrIt)->removeFromPropagateList();
    _constrsListToPropagate.clear();

    _colGenSpPtsWithZeroUb.clear();

    return true;
}

void Alg4PreprocessingOfNode::changeSubProblemBounds(ColGenSpConf * spConfPtr, const Double & value)
{
    Double spUb = *(spConfPtr->upperBoundPtr());
    if (spUb < InstMastConvexityConstr::upperBoundWhenInactive)
    {
        spUb -= value;
        if (spUb == 0)
          _colGenSpPtsWithZeroUb.insert(spConfPtr);
        Constraint * ubConvConstr = spConfPtr->upperBoundMastConstrPtr();
        ubConvConstr->curRhs(spUb);
        spConfPtr->upperBoundPtr(new Double(spUb));
        _mastConstrsToChangeRhs.push_back(ubConvConstr);
    }

    Double spLb = *(spConfPtr->lowerBoundPtr());
    if (spLb > InstMastConvexityConstr::lowerBoundWhenInactive)
    {
        spLb -= value;
        if (spLb < InstMastConvexityConstr::lowerBoundWhenInactive)
            spLb = InstMastConvexityConstr::lowerBoundWhenInactive;
        Constraint * lbConvConstr = spConfPtr->lowerBoundMastConstrPtr();
        lbConvConstr->curRhs(spLb);
        spConfPtr->lowerBoundPtr(new Double(spLb));
        if (spLb == InstMastConvexityConstr::lowerBoundWhenInactive)
          lbConvConstr->addToPreprocessedList();
        else
          _mastConstrsToChangeRhs.push_back(lbConvConstr);
    }
}

bool Alg4PreprocessingOfNode::fixVariableValue(Variable * varPtr, const Double & value)
{
    if (printL(PreprocessingPrintLevel))
        std::cout << "PreprocessingBase::fixVariableValue() var "
                  << varPtr->name() << " with bounds [" << varPtr->globalCurLb() << ", "
                  << varPtr->globalCurUb() << "] fixed value " << value << std::endl;

    varPtr->globalCurUb(varPtr->globalCurUb() - value);
    varPtr->globalCurLb(varPtr->globalCurLb() - value);
    varPtr->addToPreprocessedList();

    if (printL(PreprocessingPrintLevel))
        std::cout << "PreprocessingBase::fixVariableValue() var "
                  << varPtr->name() << " has new bounds [" << varPtr->globalCurLb() << ", "
                  << varPtr->globalCurUb() << "]" << std::endl;

    /// TO DO : review this part, it is not clear what does mean "fix master variable to a value"
    /// whether it is var = value, or var >= value
    if (varPtr->isTypeOf(VcId::InstMasterVarMask))
    {
        for (ConstVarConstrPtr2Double::const_iterator mapIt = varPtr->member2coefMap().begin();
             mapIt != varPtr->member2coefMap().end(); ++mapIt)
            if (mapIt->first->inCurProb())
            {
                Constraint * constrPtr = static_cast<Constraint *> (mapIt->first);
                constrPtr->curRhs(constrPtr->curRhs() - mapIt->second * value);
                _mastConstrsToChangeRhs.push_back(constrPtr);
                constrPtr->addToPropagateList(_constrsListToPropagate);
            }

        if ((varPtr->sense() == 'P') && (varPtr->globalCurLb() < 0) && updateLowerBound(varPtr, 0))
            return true;
        if ((varPtr->sense() == 'N') && (varPtr->globalCurUb() > 0) && updateUpperBound(varPtr, 0))
            return true;
    }

    if (varPtr->isTypeOf(VcId::SubProbVariableMask))
    {
        SubProbVariable * spVarPtr = static_cast<SubProbVariable *> (varPtr);
        for (ConstVarConstrPtr2Double::const_iterator mapIt =
                spVarPtr->masterConstrMember2coefMap().begin();
                mapIt != spVarPtr->masterConstrMember2coefMap().end(); ++mapIt)
            if (mapIt->first->inCurProb())
            {
                Constraint * constrPtr = static_cast<Constraint *> (mapIt->first);
                constrPtr->curRhs(constrPtr->curRhs() - mapIt->second * value);
                _mastConstrsToChangeRhs.push_back(constrPtr);
                constrPtr->addToPropagateList(_constrsListToPropagate);
                if (printL(PreprocessingPrintLevel))
                    std::cout << "PreprocessingBase::fixVariableValue() change rhs of constr "
                              << constrPtr->name() << " to " <<  constrPtr->curRhs() << std::endl;
            }

        Double spLowerBound(InstMastConvexityConstr::lowerBoundWhenInactive);
        Double spUpperBound(InstMastConvexityConstr::upperBoundWhenInactive);
        if (spVarPtr->probConfPtr()->lowerBoundPtr() != NULL)
          spLowerBound = *(spVarPtr->probConfPtr()->lowerBoundPtr());
        if (spVarPtr->probConfPtr()->upperBoundPtr() != NULL)
          spUpperBound = *(spVarPtr->probConfPtr()->upperBoundPtr());

        if (spVarPtr->globalCurUb() > spVarPtr->localCurUb() * spUpperBound)
          if (updateGlobalUpperBound(spVarPtr, spVarPtr->localCurUb() * spUpperBound))
            return true;
        if (spVarPtr->globalCurLb() < spVarPtr->localCurLb() * spLowerBound)
          if (updateGlobalLowerBound(spVarPtr, spVarPtr->localCurLb() * spLowerBound))
            return true;

        /// we update now local bounds based on the current values of global bounds,
        /// as the local bounds may be "not synchronized" at this point
        /// order is important (for the case with zero subproblem upper bound)!
        if (updateLocalLowerBound(spVarPtr) || updateLocalUpperBound(spVarPtr))
          return true;
    }

    return false;
}

bool Alg4PreprocessingOfNode::propagateNonLinearMasterConstraints(const MastColumn * const colPtr, const Double & value)
{
    for (ConstVarConstrPtr2Double::const_iterator mapIt = colPtr->member2coefMap().begin();
         mapIt != colPtr->member2coefMap().end(); ++mapIt)
    {
        VarConstr * vcPtr = mapIt->first;
        if (vcPtr->isTypeOf(VcId::SoftConflictsCutConstrMask))
        {
           Constraint * constrPtr = static_cast<Constraint *> (vcPtr);
            constrPtr->curRhs(constrPtr->curRhs() - mapIt->second * value);
            _mastConstrsToChangeRhs.push_back(constrPtr);
            constrPtr->addToPropagateList(_constrsListToPropagate);
            if (printL(PreprocessingPrintLevel))
            {
                std::cout << "PreprocessingBase::propagateNonLinearMasterConstraints() change rhs of constr "
                          << constrPtr->name() << " to " <<  constrPtr->curRhs() << std::endl;
            }
        }
    }
    return false;
}

bool Alg4PreprocessingOfNode::fixPartialSolution(Solution const * solPtr)
{
    if (solPtr != NULL)
      {
        std::set<ColGenSpConf *> subProblemsWithChangedBounds;
        VarPtr2DoubleMap projectedSolution;
        for (VarPtr2DoubleMap::const_iterator mapIt = solPtr->solVarValMap().begin();
                mapIt != solPtr->solVarValMap().end(); mapIt++)
          {
            if (mapIt->first->isTypeOf(VcId::InstMasterVarMask)) /// pure master variable
                projectedSolution[mapIt->first] = mapIt->second;
            if (mapIt->first->isTypeOf(VcId::MastColumnMask)) /// master column
              {
                MastColumn * colPtr = static_cast<MastColumn *> (mapIt->first);
                colPtr->fillAggregateSol(projectedSolution, mapIt->second);
                changeSubProblemBounds(colPtr->cgSpConfPtr(), mapIt->second);
                subProblemsWithChangedBounds.insert(colPtr->cgSpConfPtr());
                if (propagateNonLinearMasterConstraints(colPtr, mapIt->second))
                    return false;
              }
          }

        for (VarPtr2DoubleMap::const_iterator mapIt = projectedSolution.begin();
                mapIt != projectedSolution.end(); mapIt++)
          {
            if (fixVariableValue(mapIt->first, mapIt->second))
                return true;
          }

        /// we need to change the bounds of all variables from subproblems with changed bounds
        for (std::set<ColGenSpConf *>::iterator cgspPtrIt = subProblemsWithChangedBounds.begin();
             cgspPtrIt != subProblemsWithChangedBounds.end(); ++cgspPtrIt)
          {
            const VarIndexManager & probVarSet = (*cgspPtrIt)->probPtr()->probVarSet();
            for (VarIndexManager::const_iterator varPtrIt = probVarSet.begin(VcIndexStatus::Active, 's');
                 varPtrIt != probVarSet.end(VcIndexStatus::Active, 's'); varPtrIt++)
              if (fixVariableValue(*varPtrIt, 0))
                return true;
          }
      }
    return false;
}

bool Alg4PreprocessingOfNode::preprocess(Solution const * solPtr, bool const initialPreprocessing)
{
    if (computeInitialConstrsSlacks())
      {
        if (printL(-1))
	      std::cout << "Preprocessing determines infeasibility (init. constraint slacks)" << std::endl;
        return exitWhenInfeasible();
      }

    if (fixPartialSolution(solPtr))
     {
         if (printL(-1))
             std::cout << "Preprocessing determines infeasibility (after fix of partial solution)" << std::endl;
         return exitWhenInfeasible();
     }

    if (initialPreprocessing)
    {
        if (initialUpdateOfSpVarBounds())
        {
            if (printL(-1))
                std::cout << "Preprocessing determines infeasibility (init. update of sp.var. bounds)" << std::endl;
            return exitWhenInfeasible();
        }
    }

    if(_needToPreprocessCompSetBrConstr && computeCompSetBrConstrInducedGlobalBdOnSpVar())
      {
          if (printL(-1))
              std::cout << "Preprocessing determines infeasibility (comp. set. branching)" << std::endl;
        return exitWhenInfeasible();
      }

    if (propagate())
      {
          if (printL(-1))
              std::cout << "Preprocessing determines infeasibility" << std::endl;
        return exitWhenInfeasible();
      }

    applyPreprocessingListsInProbAndForm(initialPreprocessing);
    clearPreprocessingLists();

    return false;
}

Alg4PreprocessingOfNode::Alg4PreprocessingOfNode(Alg4PreprocessingOfNode & that) :
  _problemPts(that._problemPts), _constrsListToPropagate(), _colGenSpPtsWithZeroUb(), _pcPtrToNonZeroVarPtrMap(),
  _mastConstrsToChangeRhs(), _needToPreprocessCompSetBrConstr(false)
{
}

Alg4PreprocessingOfNode::Alg4PreprocessingOfNode(std::list<Problem *> const & problemPts) :
  _problemPts(problemPts), _constrsListToPropagate(), _colGenSpPtsWithZeroUb(), _pcPtrToNonZeroVarPtrMap(),
  _mastConstrsToChangeRhs(), _needToPreprocessCompSetBrConstr(false)
{
}

bool Alg4PreprocessingOfNode::preprocessRoot()
{
  /// in the beginning, we preprocess all static constraints
  ConstrIndexManager::iterator constrPtrIt;
  for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
       ++probPtrIt)
    for (constrPtrIt = (*probPtrIt)->probConstrSet().begin(VcIndexStatus::Active, 's');
         constrPtrIt != (*probPtrIt)->probConstrSet().end(VcIndexStatus::Active, 's'); ++constrPtrIt)
      {
        if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask) && (*constrPtrIt)->toBeUsedInPreprocessing())
          {
            (*constrPtrIt)->addToPropagateList(_constrsListToPropagate);
          }
      }
  return preprocess(NULL, true);
}

bool Alg4PreprocessingOfNode::preprocessConstraints(ConstrPtrList constrsToPropagate)
{
  for(ConstrPtrList::iterator it = constrsToPropagate.begin(); it != constrsToPropagate.end(); it++)
    {
      if ((*it)->toBeUsedInPreprocessing()
          && !(*it)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
        (*it)->addToPropagateList(_constrsListToPropagate);
      if ((*it)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
        {
          /// deactivated by Ruslan as Issam mentionned that there were bugs
          //_needToPreprocessCompSetBrConstr = true;
        }
    }
  return preprocess(NULL, false);
}

bool Alg4PreprocessingOfNode::preprocessAfterFixingPartialSolution(Solution const * solPtr)
{
    return preprocess(solPtr, false);
}

Algorithm4PreprocessingAtRoot::Algorithm4PreprocessingAtRoot(std::list<Problem *> const & problemPts) :
  Alg4PreprocessingOfNode(problemPts)
{
}

Algorithm4PreprocessingAtNodeOtherThanRoot::Algorithm4PreprocessingAtNodeOtherThanRoot(
  std::list<Problem *> const & problemPts) :
  Alg4PreprocessingOfNode(problemPts)
{
}


bool Alg4PreprocessingOfNode::computeCompSetBrConstrInducedGlobalBdOnSpVar()
{
    bool isMasterProblem = true;
    for (std::list<Problem *>::const_iterator probPtrIt = _problemPts.begin(); probPtrIt != _problemPts.end();
         ++probPtrIt)
    {
        //the first problem is the not considered here as it is the master
        if (isMasterProblem) {
            isMasterProblem = false;
            continue;
        }

        Problem* subProblemPtr = *probPtrIt;

        // we only consider static variables.
        for (VarIndexManager::iterator varPtrIt = subProblemPtr->probVarSet().begin(Active, 's');
             varPtrIt != subProblemPtr->probVarSet().end(Active, 's'); ++varPtrIt)
            ///@todo:iterate only on SP var that are invoved in compSetBrConstr
        {
            SubProbVariable* spVarPtr = static_cast<SubProbVariable*> (*varPtrIt);

            Double globLbInducedByBCS(0);
            Double globUbInducedByBCS(0);
            Double lowerCardinality = *(spVarPtr->probConfPtr()->lowerBoundPtr());
            Double upperCardinality = *(spVarPtr->probConfPtr()->upperBoundPtr());


            for (CompSetConstrPtr2DoubleMap::const_iterator itlbPt = spVarPtr->mapCompSetBrConstr2lowerBd().begin();
                 itlbPt != spVarPtr->mapCompSetBrConstr2lowerBd().end(); itlbPt++)
                if (itlbPt->first->inCurProb())
                {
                    globLbInducedByBCS += itlbPt->first->marginLvalue() * itlbPt->second;
                    lowerCardinality -= itlbPt->first->marginLvalue();
                }
            if (lowerCardinality.positive())
                globLbInducedByBCS += lowerCardinality * spVarPtr->localCurLb(); // newCurLb(); //modified by Issam.

            if (globLbInducedByBCS > spVarPtr->globalCurLb() && printL(PreprocessingPrintLevel))
                std::cout << "PreprocessingNonRoot::computeCompSetBrConstrInducedGlobalBdOnSpVar()"
                          << " induces a better LB for variable " << spVarPtr->name()
                          << " FROM " << spVarPtr->globalCurLb() << " TO " << globLbInducedByBCS << std::endl;

            if (updateGlobalLowerBound(spVarPtr, globLbInducedByBCS, NULL))
                return true;

            for (CompSetConstrPtr2DoubleMap::const_iterator itubPt = spVarPtr->mapCompSetBrConstr2upperBd().begin();
                 itubPt != spVarPtr->mapCompSetBrConstr2upperBd().end(); itubPt++)
                if (itubPt->first->inCurProb())
                {
                    globUbInducedByBCS += itubPt->first->marginLvalue() * itubPt->second;
                    upperCardinality -= itubPt->first->marginLvalue();
                }
            if (upperCardinality.positive())
                globUbInducedByBCS += upperCardinality * spVarPtr->localCurUb(); //newCurUb(); //modified by Issam.

            if(globUbInducedByBCS < spVarPtr->globalCurUb() && printL(PreprocessingPrintLevel))
                std::cout << "PreprocessingNonRoot::computeCompSetBrConstrInducedGlobalBdOnSpVar()"
                          << " induces a better UB for variable " << spVarPtr->name()
                          << " FROM " << spVarPtr->globalCurUb() << " TO " << globUbInducedByBCS << std::endl;

            if (updateGlobalUpperBound(spVarPtr, globUbInducedByBCS, NULL))
                return true;
        }
    }

    return false;
}
