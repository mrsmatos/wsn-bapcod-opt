/**
 *
 * This file bcLimMemStrongKPathCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

//
//  bcLimMemRankOneCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 27/10/2015.
//
//

#ifdef BCP_RCSP_IS_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcNonPublicCuts.hpp"
#include "bcProblemC.hpp"
#include "bcMastColumnC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "rcsp_interface.hpp"

#define LimMemStrongKPathCutsPrintLevel 3


LimMemKPathCut::LimMemKPathCut(const IndexCell & id,
                               GenericLimMemStrongKPathCutConstr * genConstrPtr,
                               ProbConfig * probConfigPtr,
                               const std::string & name,
                               const bcp_rcsp::StrongKPathCut * rcspCutPtr):
  InstMasterConstr(id, genConstrPtr, probConfigPtr, name, rcspCutPtr->rightHandSide, 'G',
                   genConstrPtr->defaultType(), genConstrPtr->defaultKind(),
                   genConstrPtr->defaultFlag()),
  Base4NonLinearConstraint(), _rcspCutPtr(rcspCutPtr), _genLimMemStrongKPathCutConstr(genConstrPtr)
{
}

LimMemKPathCut::~LimMemKPathCut()
{
    delete _rcspCutPtr;
}

bool LimMemKPathCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::LimMemoryKPathCutConstrMask, vcIdentifier);
}

void LimMemKPathCut::nicePrint(std::ostream & os) const
{
    _genLimMemStrongKPathCutConstr->nicePrintCut(_rcspCutPtr, os);
}

const std::vector<int> & LimMemKPathCut::setIds() const
{
    return _rcspCutPtr->setIds;
}

void LimMemKPathCut::setMembership()
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
        LpCoef lpCoeff = _genLimMemStrongKPathCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
        if (lpCoeff.first)
          includeMember(*it, lpCoeff.second, cumulativeCoef);
      }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
      {
        LpCoef lpCoeff = _genLimMemStrongKPathCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
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
          LpCoef lpCoeff = _genLimMemStrongKPathCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }

  Constraint::setMembership();

  return;
}

GenericLimMemStrongKPathCutConstr::GenericLimMemStrongKPathCutConstr(Model * modelPtr,
                                                                     ProbConfig * probConfPtr,
                                                                     const std::string & name,
                                                                     const Double & nonRootPriorityLevel,
                                                                     const Double & rootPriorityLevel,
                                                                     const bool & isFacultative,
                                                                     const bool & equalityCase,
                                                                     const int & maxCapacity,
                                                                     const std::vector<int> & demands,
                                                                     const int & twoPathCutsResId) :
  GenericCutConstr(modelPtr, probConfPtr, name, isFacultative ? 'F' : 'C', SelectionStrategy::MostFractional,
                   nonRootPriorityLevel, rootPriorityLevel, false),
  Base4NonLinearGenericConstr(NULL), _equalityCase(equalityCase), _maxCapacity(maxCapacity),
  _demands(demands), _twoPathCutsResId(twoPathCutsResId), _interfacePtr(nullptr)
{
}

GenericLimMemStrongKPathCutConstr::~GenericLimMemStrongKPathCutConstr()
{
    delete _interfacePtr;
}

bool GenericLimMemStrongKPathCutConstr::prepareSeparation()
{
    bcp_rcsp::StrongKPathSeparatorParameters cutSepParams;
    cutSepParams.maxNumPerRound = param().RCSPrankOneCutsMaxNumPerRound();
    cutSepParams.equalityCase = _equalityCase;
    if (printL(0))
        cutSepParams.printLevel = 1;
    else
        cutSepParams.printLevel = 0;
    std::vector<const bcp_rcsp::GraphData *> graphs;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() != NULL)
            graphs.push_back((*cgSpConfPtrIt)->rcspGraphPtr());
    }

    std::vector<double> demands(_demands.size());
    for (int index = 0; index < (int)_demands.size(); ++index)
        demands[index] = (double)_demands[index];
    _interfacePtr = bcp_rcsp::createAndPrepareStrongKPathCutSeparation(graphs, demands, _maxCapacity, _twoPathCutsResId,
                                                                       cutSepParams);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "BaPCod error : could not prepare strong k-path cuts separation" << std::endl;
        return false;
    }

    return true;
}

void GenericLimMemStrongKPathCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

void GenericLimMemStrongKPathCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset < InstanciatedConstr *, CutSeparationPriorityComp > & generatedCutConstrSet)

{
    if (_interfacePtr == NULL)
        return;

    if (printL(LimMemStrongKPathCutsPrintLevel))
        std::cout << "Separation of strong k-path cuts started" << std::endl;

    bcp_rcsp::FractionalMasterSolution rcspFracSolution;
    rcspFracSolution.solPts.reserve(curSol.size());
    rcspFracSolution.values.reserve(curSol.size());

    for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
        if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;

        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);

        auto rcspSolPtr = colPtr->spSol()->rcspSolPtr();
        if (rcspSolPtr != nullptr)
        {
            rcspFracSolution.solPts.push_back(rcspSolPtr);
            rcspFracSolution.values.push_back(colPtr->val());
        }

        if (printL(LimMemStrongKPathCutsPrintLevel))
        {
            std::cout << "Route from sp." << colPtr->cgSpConfPtr()->name()
                      << " with val = " << colPtr->val() << " and cost = " << colPtr->costrhs() << " : ";
            colPtr->spSol()->printOrderedSolution();
        }
    }

    std::vector<const bcp_rcsp::StrongKPathCut *> rcspCutPts;
    if (_interfacePtr->separate(rcspFracSolution, rcspCutPts))
    {
        for (const auto rcspCutPtr : rcspCutPts)
        {
            std::string name("SKP");
            MultiIndex newCutId(rcspCutPtr->id);
            newCutId.appendRef2name(name, multiIndexNames());
            generatedCutConstrSet.insert(new LimMemKPathCut(newCutId, this, probConfPtr(), name, rcspCutPtr));
        }
    }
    else
    {
        for (auto & rcspCutPtr : rcspCutPts)
            delete rcspCutPtr;
        rcspCutPts.clear();
    }
}


const LpCoef GenericLimMemStrongKPathCutConstr
             ::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
    if (!icPtr->isTypeOf(VcId::LimMemoryKPathCutConstrMask))
        return LpCoef(false, 0.0);

    const bcp_rcsp::StrongKPathCut * rcspCutPtr = static_cast<LimMemKPathCut *>(icPtr)->rcspCutPtr();
    double coeff = _interfacePtr->coefficient(colPtr->spSol()->rcspSolPtr(), rcspCutPtr);
    if (coeff > 0.0)
        return LpCoef(true, coeff);
    return LpCoef(false, 0.0);
}

void GenericLimMemStrongKPathCutConstr::nicePrintCut(const bcp_rcsp::StrongKPathCut * rcspCutPtr, std::ostream & os)
{
    if (_interfacePtr != nullptr)
        _interfacePtr->niceCutPrint(rcspCutPtr, os);
}

#endif /* BCP_RCSP_IS_FOUND */
