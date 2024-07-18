/**
 *
 * This file bcLimMemRankOneCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcNetworkBasedCuts.hpp"
#include "bcProblemC.hpp"
#include "bcMastColumnC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "rcsp_interface.hpp"

#define LimMemRankOneCutsPrintLevel 3


LimMemRankOneCut::LimMemRankOneCut(const IndexCell & id,
                                   GenericLimMemRankOneCutConstr * genConstrPtr,
                                   ProbConfig * probConfigPtr,
                                   const std::string & name,
                                   const bcp_rcsp::RankOneCut * rcspCutPtr):
  InstMasterConstr(id, genConstrPtr, probConfigPtr, name, rcspCutPtr->rightHandSide,
                   (rcspCutPtr->cutClass == SET_COVERING_CUT) ? 'G' : 'L',
                   genConstrPtr->defaultType(), genConstrPtr->defaultKind(),
                   genConstrPtr->defaultFlag()),
  Base4NonLinearConstraint(), _rcspCutPtr(rcspCutPtr), _genLimMemRankOneCutConstr(genConstrPtr)
{
}

LimMemRankOneCut::~LimMemRankOneCut()
{
    delete _rcspCutPtr;
}

bool LimMemRankOneCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::LimMemoryRankOneCutConstrMask, vcIdentifier);
}

void LimMemRankOneCut::nicePrint(std::ostream& os) const
{
   _genLimMemRankOneCutConstr->nicePrintCut(_rcspCutPtr, os);
}

const std::vector<int> & LimMemRankOneCut::setIds() const
{
    return _rcspCutPtr->setIds;
}

const std::vector<int> & LimMemRankOneCut::coeffs() const
{
    return _rcspCutPtr->coeffs;
}

int LimMemRankOneCut::numRows() const
{
    return _rcspCutPtr->numRows;
}

int LimMemRankOneCut::denominator() const
{
    return _rcspCutPtr->denominator;
}

void LimMemRankOneCut::setMembership()
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
        LpCoef lpCoeff = _genLimMemRankOneCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
        if (lpCoeff.first)
          includeMember(*it, lpCoeff.second, cumulativeCoef);
      }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
      {
        LpCoef lpCoeff = _genLimMemRankOneCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
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
          LpCoef lpCoeff = _genLimMemRankOneCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }

  Constraint::setMembership();

  return;
}

GenericLimMemRankOneCutConstr::GenericLimMemRankOneCutConstr(Model * modelPtr,
                                                             ProbConfig * probConfPtr,
                                                             const std::string & name,
                                                             const Double & nonRootPriorityLevel,
                                                             const Double & rootPriorityLevel,
                                                             const std::string & spVarName,
                                                             const int & memoryType,
                                                             const bool & isFacultative) :
  GenericCutConstr(modelPtr, probConfPtr, name, isFacultative? 'F' : 'C', SelectionStrategy::MostFractional,
                   nonRootPriorityLevel, rootPriorityLevel, false),
  Base4NonLinearGenericConstr(NULL), _interfacePtr(nullptr), _nbPackingSets(0), _memoryType(memoryType),
  _currentPhase(1), _maxNumberOfRows(param().RCSPrankOneCutsMaxNumRows()), _standaloneFileName(""),
  _spVarName(spVarName), _genIndicVarPts(), _genIndicConstrPts()
{
}

GenericLimMemRankOneCutConstr::~GenericLimMemRankOneCutConstr()
{
    delete _interfacePtr;
}

bool GenericLimMemRankOneCutConstr::prepareSeparation()
{
    if (_spVarName != "")
    {
        _nbPackingSets = 0;
        std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt;
        for (cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
             cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
        {
            GenericVar * genVarPtr = (*cgSpConfPtrIt)->getGenericVar(_spVarName);
            const IndexCell2InstancVarPtrMap & varPtrMap = genVarPtr->indexCell2InstancVarPtrMap();
            for (IndexCell2InstancVarPtrMap::const_iterator mapIt = varPtrMap.begin(); mapIt != varPtrMap.end(); ++mapIt)
            {
                int firstId = mapIt->first.first();
                if (_nbPackingSets < firstId + 1)
                    _nbPackingSets = firstId + 1;
            }
        }

        for (cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
             cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
        {
            _genIndicVarPts[*cgSpConfPtrIt] = modelPtr()->createGenericVar(*cgSpConfPtrIt, BcVarConstrType::local2Formulation,
                                                                           defaultName() + "V", MultiIndexNames('k'),
                                                                           'I', 0);
            _genIndicConstrPts[*cgSpConfPtrIt] = modelPtr()->createGenericConstr(*cgSpConfPtrIt,
                                                                                 BcVarConstrType::local2Formulation,
                                                                                 defaultName() + "C",
                                                                                 MultiIndexNames('k'), 'L');
        }

        if ((param().RCSPrankOneCutsMaxNumRows() >= 4) && (param().RCSPrankOneCutsMaxNumTriplets() <= 0))
        {
            std::cerr << "lm-1Rank cuts separator error: for separation of 4 and 5-rows cuts, tuple based separation "
                      << "should be used" << std::endl;
            return false;
        }

        if (param().RCSPrankOneCutsMaxNumRows() >= 6)
        {
            std::cerr << "lm-1Rank cuts separator error: separation of 6-row and more cuts only allowed if network(s) "
                      << "are defined" << std::endl;
            return false;
        }

        return true;
    }

    bcp_rcsp::LimMemRankOneCutsSeparatorParameters cutSepParams;
    cutSepParams.maxNumRows = param().RCSPrankOneCutsMaxNumRows();
    cutSepParams.maxNumPerRound = param().RCSPrankOneCutsMaxNumPerRound();
    cutSepParams.memoryType = param().RCSPrankOneCutsMemoryType();
    cutSepParams.arcMemoryType = param().RCSPrankOneCutsArcMemoryType();
    cutSepParams.spDependentMemory = param().RCSPrankOneCutsSpDependentMemory();
    if (param().RCSPrankOneCutsLSnumIterations() >= 0) {
        cutSepParams.numLocSearchIterations = param().RCSPrankOneCutsLSnumIterations();
        cutSepParams.useLocSearchOnlyForSixRowsAndMore = false;
    } else {
        cutSepParams.numLocSearchIterations = - param().RCSPrankOneCutsLSnumIterations();
        cutSepParams.useLocSearchOnlyForSixRowsAndMore = true;
    }
    cutSepParams.neighbourhoodSize = param().RCSPrankOneCutsNeighbourhoodSize();
    cutSepParams.maxNumTriplets = param().RCSPrankOneCutsMaxNumTriplets();
    cutSepParams.cutsTypeToSeparate = param().RCSPrankOneCutsTypeToSeparate();
    cutSepParams.cutViolationTolerance = param().BapCodCutViolationTolerance();
    cutSepParams.useCoveringSetsForSeparatingOneRowPackingCuts = param().RCSPuseCovSetsForSepOneRowPackingCuts();
    cutSepParams.printLevel = param().DEFAULTPRINTLEVEL();
    std::vector<const bcp_rcsp::GraphData *> graphs;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() == NULL)
            (*cgSpConfPtrIt)->fillRCSPGraph();

        if ((*cgSpConfPtrIt)->rcspGraphPtr() != NULL)
            graphs.push_back((*cgSpConfPtrIt)->rcspGraphPtr());
    }

    _interfacePtr = bcp_rcsp::createAndPrepareRankOneCutSeparation(graphs, cutSepParams);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "BaPCod error : could not prepare rank-1 cuts separation" << std::endl;
        return false;
    }

    return true;
}

void GenericLimMemRankOneCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

void GenericLimMemRankOneCutConstr::updateSubprobemsWithIndicatorVarAndConstr(const bcp_rcsp::RankOneCut * cutPtr)
{
    if (param().colGenSubProbSolMode().status() == SolutionMethod::customSolver)
        return;

    int coeffSum = 0;
    for (int rowId = 0; rowId < cutPtr->numRows; ++rowId)
        coeffSum += cutPtr->coeffs[rowId];
    int indicVarUpperBound = coeffSum / cutPtr->denominator;


    /// subproblems are solved by MIP, thus we add conflict indicator variables and linking constraints in subproblems
    for (std::vector< ColGenSpConf * >::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        std::list<SubProbVariable *> newIndicVarPts;
        std::list<Constraint *> newIndicConstrPts;
        MultiIndex indicId(cutPtr->id);

        GenericVar * indicGenVarPtr = _genIndicVarPts[*cgSpConfPtrIt];
        InstanciatedVar * ivPtr = modelPtr()->createVariable(*cgSpConfPtrIt, indicGenVarPtr, indicId, 0, 'I',
                                                             indicVarUpperBound);
        ivPtr = (*cgSpConfPtrIt)->castAndAddVariable(ivPtr, false);
        SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(ivPtr);
        newIndicVarPts.push_back(spVarPtr);

        GenericConstr * indicGenConstrPtr = _genIndicConstrPts[*cgSpConfPtrIt];
        InstanciatedConstr * icPtr = modelPtr()->createConstraint(*cgSpConfPtrIt, indicGenConstrPtr, indicId,
                                                                  cutPtr->denominator - 1, 'L');
        icPtr->includeMember(ivPtr, -cutPtr->denominator, false);

        GenericVar * setIdGenVarPtr = (*cgSpConfPtrIt)->getGenericVar(_spVarName);
        for (int rowId = 0; rowId < cutPtr->numRows; ++rowId)
        {
            MultiIndex setId(cutPtr->setIds[rowId]);
            InstanciatedVar * setIdVarPtr = setIdGenVarPtr->checkIfInstanciationAlreadyExist(setId);
            if (setIdVarPtr != NULL)
            {
                icPtr->includeMember(setIdVarPtr, cutPtr->coeffs[rowId], false);
            }
        }
        newIndicConstrPts.push_back(icPtr);

        (*cgSpConfPtrIt)->probPtr()->addVarSet(newIndicVarPts, 1, 2);
        (*cgSpConfPtrIt)->probPtr()->addConstrSet(newIndicConstrPts, 1, 2);
    }
}

void GenericLimMemRankOneCutConstr::cutSeparationRoutine(const std::string & fileName)
{
    VarPtrSet emptySol;
    std::multiset < InstanciatedConstr *, CutSeparationPriorityComp > cutConstrSet;
    _standaloneFileName = fileName;
    cutSeparationRoutine(emptySol, cutConstrSet);
}

void GenericLimMemRankOneCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset < InstanciatedConstr *, CutSeparationPriorityComp > & generatedCutConstrSet)

{
    if (_interfacePtr == NULL)
        return;

    if (printL(LimMemRankOneCutsPrintLevel))
        std::cout << "Separation of Rank-1 cuts started" << std::endl;

    bcp_rcsp::RankOneCutSeparationInput input;
    input.memoryTypeForThisRound = _memoryType;
    input.arcIdsInPathsArePackSetIds = (_spVarName != "");
    input.phase = _currentPhase;

    /// we retrieve the set of active rank one cuts, they are needed to join the arc memory set
    if (input.memoryTypeForThisRound == ARC_MEMORY_TYPE)
    {
        ConstrIndexManager::const_iterator constrPtrIt;
        for (constrPtrIt = probConfPtr()->probPtr()->probConstrSet().begin(VcIndexStatus::Active, 'd');
             constrPtrIt != probConfPtr()->probPtr()->probConstrSet().end(VcIndexStatus::Active, 'd');
             ++constrPtrIt)
            if ((*constrPtrIt)->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
            {
                LimMemRankOneCut *cutPtr = static_cast<LimMemRankOneCut *>(*constrPtrIt);
                input.rankOneCuts.push_back(std::make_pair(cutPtr->_rcspCutPtr, cutPtr->valOrSepPointVal()));
            }
    }

    input.fracMasterSol.solPts.reserve(curSol.size());
    input.fracMasterSol.values.reserve(curSol.size());
    input.solutionFile = _standaloneFileName;
    for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
        if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;

        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);

        if (_spVarName == "")
        {
            auto rcspSolPtr = colPtr->spSol()->rcspSolPtr();
            if (rcspSolPtr != nullptr)
            {
                input.fracMasterSol.solPts.push_back(rcspSolPtr);
                input.fracMasterSol.values.push_back(colPtr->val());
            }
        }
        else
        {
            ProbConfig * probConfPtr = colPtr->cgSpConfPtr();
            auto * rcspSolPtr = new bcp_rcsp::Solution(probConfPtr->id().first());
            rcspSolPtr->enumeratedFlag = colPtr->spSol()->enumeratedFlag();
            const VarPtr2DoubleMap & varValMap = colPtr->spSol()->solVarValMap();
            for (VarPtr2DoubleMap::const_iterator mapIt = varValMap.begin(); mapIt != varValMap.end(); ++mapIt)
            {
                if (mapIt->first->isTypeOf(VcId::InstanciatedVarMask))
                {
                    const InstanciatedVar * iVarPtr = static_cast<const InstanciatedVar *>(mapIt->first);
                    if ((iVarPtr->genVarPtr()->defaultName() == _spVarName) && (iVarPtr->id().endPosition == 1))
                    {
                        int multiplicity = (int)Dfloor(mapIt->second);
                        int setId = iVarPtr->id().first();
                        for (int mult = 0; mult < multiplicity; ++mult)
                            rcspSolPtr->arcIds.push_back(setId);
                    }
                }
            }
        }
        if (printL(LimMemRankOneCutsPrintLevel))
        {
            std::cout << "Route from sp." << colPtr->cgSpConfPtr()->name()
                      << " with val = " << colPtr->val() << " and cost = " << colPtr->costrhs() << " : ";
            colPtr->spSol()->printOrderedSolution();
        }
    }

    std::vector<const bcp_rcsp::RankOneCut *> rcspCutPts;
    if (_interfacePtr->separate(input, rcspCutPts))
    {
        for (const auto rcspCutPtr : rcspCutPts)
        {
            _currentPhase = (std::max)(_currentPhase, rcspCutPtr->numRows);

            std::string name("R1C");
            MultiIndex newCutId(rcspCutPtr->id);
            newCutId.appendRef2name(name, multiIndexNames());
            if (_spVarName != "")
                updateSubprobemsWithIndicatorVarAndConstr(rcspCutPtr);
            generatedCutConstrSet.insert(new LimMemRankOneCut(newCutId, this, probConfPtr(), name, rcspCutPtr));
        }
    }
    else
    {
        for (auto & rcspCutPtr : rcspCutPts)
            delete rcspCutPtr;
        rcspCutPts.clear();
    }
}


const LpCoef GenericLimMemRankOneCutConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
    if (!icPtr->isTypeOf(VcId::LimMemoryRankOneCutConstrMask))
        return LpCoef(false, 0.0);

    const bcp_rcsp::RankOneCut * rcspCutPtr = static_cast<LimMemRankOneCut *>(icPtr)->rcspCutPtr();
    if (_spVarName != "")
    {
        int state = 0;
        const VarPtr2DoubleMap & varValMap = colPtr->spSol()->solVarValMap();
        for (VarPtr2DoubleMap::const_iterator mapIt = varValMap.begin(); mapIt != varValMap.end(); ++mapIt)
        {
            if (mapIt->first->isTypeOf(VcId::InstanciatedVarMask))
            {
                const InstanciatedVar * iVarPtr = static_cast<const InstanciatedVar *>(mapIt->first);
                if ((iVarPtr->genVarPtr()->defaultName() == _spVarName) && (iVarPtr->id().endPosition == 1))
                {
                    int multiplicity = (int)Dfloor(mapIt->second);
                    int setId = iVarPtr->id().first();
                    int rowId = 0;
                    while ((rowId < rcspCutPtr->numRows) && (rcspCutPtr->setIds[rowId] != setId))
                        ++rowId;
                    if (rowId < rcspCutPtr->numRows)
                        state += rcspCutPtr->coeffs[rowId] * multiplicity;
                }
            }
        }
        int coeff = state / rcspCutPtr->denominator;
        if (coeff > 0)
            return LpCoef(true, coeff);
        return LpCoef(false, 0.0);
    }

    double coeff = _interfacePtr->coefficient(colPtr->spSol()->rcspSolPtr(), rcspCutPtr);
    if (coeff > 0.0)
        return LpCoef(true, coeff);
    return LpCoef(false, 0.0);
}

void GenericLimMemRankOneCutConstr::nicePrintCut(const bcp_rcsp::RankOneCut * rcspCutPtr, std::ostream & os)
{
    if (_interfacePtr != nullptr)
        _interfacePtr->niceCutPrint(rcspCutPtr, os);
}

void GenericLimMemRankOneCutConstr::resetSeparationPhase()
{
  _currentPhase = 1;
}

void GenericLimMemRankOneCutConstr::increaseSeparationPhase()
{
  if (_currentPhase == 1)
    _currentPhase = 3;
  else if (_currentPhase < _maxNumberOfRows)
    _currentPhase += 1;
}

bool GenericLimMemRankOneCutConstr::separationPhaseIsMaximum() const
{
  return (_currentPhase == _maxNumberOfRows);
}

void GenericLimMemRankOneCutConstr::setArcMemory()
{
    _memoryType = ARC_MEMORY_TYPE;
}

void GenericLimMemRankOneCutConstr::setVertexMemory()
{
    _memoryType = VERTEX_MEMORY_TYPE;
}

#endif /* BCP_RCSP_IS_FOUND */
