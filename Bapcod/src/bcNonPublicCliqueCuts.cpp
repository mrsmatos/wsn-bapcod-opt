/**
 *
 * This file bcNonPublicCliqueCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcGenNonRobustCutsC.hpp
//  Project
//
//  Created by Ruslan Sadykov on 16/02/2015.
//
//

#ifdef BCP_RCSP_IS_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcNonPublicCuts.hpp"

#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

#include "rcsp_interface.hpp"


#ifdef CLIQUE_SEP_IS_FOUND

#include "CglEClique.hpp"

void getActiveCliqueCuts(const BcFormulation & spPtr, std::vector<const BcCliqueCut *> & cutPts)
{
  if (spPtr.probConfPtr() == NULL)
  {
    std::cerr << "ERROR Model BcFormulation == NULL in getLimMemRankOneActiveMasterCutsList" << std::endl;
    exit(1);
  }

  MasterConf * mastConfPtr = spPtr.probConfPtr()->modelPtr()->master();
  GenericCutConstr * _genR1CutConstrPtr = mastConfPtr->getGenericCutConstr("CLQ");
  if (_genR1CutConstrPtr == NULL)
    return;

  const IndexCell2InstancConstrPtrMap & constrPtrMap = _genR1CutConstrPtr->indexCell2InstancConstrPtrMap();
  for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
    if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
        && it->second->isTypeOf(VcId::CliqueCutConstrMask))
    {
      CliqueCut * cliqueCutPtr = static_cast<CliqueCut *>(it->second);
      cutPts.push_back(new BcCliqueCut(cliqueCutPtr));
    }
}

/***************************************************************************
 *******************   TODO: Methods for CliqueCut   ***********************
 ***************************************************************************/


CliqueCut::CliqueCut(const IndexCell& id, GenericCliqueCutConstr * genConstrPtr, ProbConfig * probConfigPtr,
                     const std::string & name, const std::set<std::vector<int> > & setIds):
        InstMasterConstr(id, genConstrPtr, probConfigPtr, name, 1.0, 'L', genConstrPtr->defaultType(),
                         genConstrPtr->defaultKind(), genConstrPtr->defaultFlag()),
        Base4NonLinearConstraint(), _setIds(), _genCliqueCutConstr(genConstrPtr)
{
   for (std::set<std::vector<int> >::const_iterator vectIt = setIds.begin(); vectIt != setIds.end(); ++vectIt)
       _setIds.push_back(*vectIt);
}

CliqueCut::~CliqueCut()
{
}

void CliqueCut::setMembership()
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
            LpCoef lpCoeff = _genCliqueCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
            if (lpCoeff.first)
                includeMember(*it, lpCoeff.second, cumulativeCoef);
        }
    for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
         it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
        if ((*it)->isTypeOf(VcId::MastColumnMask))
        {
            LpCoef lpCoeff = _genCliqueCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
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
                LpCoef lpCoeff = _genCliqueCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
                if (lpCoeff.first)
                    includeMember(*it, lpCoeff.second, cumulativeCoef);
            }

    Constraint::setMembership();

    return;
}

bool CliqueCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
    return compareIdentifier(VcId::CliqueCutConstrMask, vcIdentifier);
}

void CliqueCut::nicePrint(std::ostream & os) const
{
    os << "Clique cut " << name() << ": packSetIds = (";
    for (std::vector<std::vector<int> >::const_iterator intVectIt = _setIds.begin(); intVectIt != _setIds.end();
         ++intVectIt)
    {
       os << "{";
       for (std::vector<int>::const_iterator setIdIt = intVectIt->begin(); setIdIt != intVectIt->end(); ++setIdIt)
       {
           if (setIdIt != intVectIt->begin())
               os << ", ";
           os << *setIdIt;
       }
       os << "}";
    }
    os << ")" << std::endl;
}

/***************************************************************************
 *************   TODO: Methods for GenericCliqueCutConstr   ****************
 ***************************************************************************/

GenericCliqueCutConstr::GenericCliqueCutConstr(Model * modelPtr,
                                               ProbConfig * probConfPtr,
                                               const std::string & name,
                                               const Double & nonRootPriorityLevel,
                                               const Double & rootPriorityLevel,
                                               const int & maxNumCutsPerRound,
                                               const double & separationMaxTime):
  GenericCutConstr(modelPtr, probConfPtr, name, 'F', SelectionStrategy::MostFractional,
                   nonRootPriorityLevel, rootPriorityLevel, false),
  Base4NonLinearGenericConstr(NULL), _numGeneratedCuts(0), _maxNumCutsPerRound(maxNumCutsPerRound),
  _separationMaxTime(separationMaxTime), _nbPackingSets(0)
{
}

GenericCliqueCutConstr::~GenericCliqueCutConstr()
{
}

bool GenericCliqueCutConstr::prepareSeparation()
{
    bool networkDefined = false;
    std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt;
    for (cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        const NetworkFlow * networkPtr = (*cgSpConfPtrIt)->networkFlowPtr();
        if (networkPtr != NULL)
        {
            networkDefined = true;
            if (_nbPackingSets < networkPtr->packingSetPts().size())
                _nbPackingSets = networkPtr->packingSetPts().size();
        }
    }

    if (!networkDefined)
    {
        std::cerr << "Clique cuts separator error: a network should be defined for at least one CG subproblem "
                  << std::endl;
        return false;
    }
    if (_nbPackingSets >= MAX_NUMBER_OF_ELEMENTARITY_SETS)
    {
        std::cerr << "Clique cuts separator error: number of packing sets cannot be more than "
                  << MAX_NUMBER_OF_ELEMENTARITY_SETS - 1 << std::endl;
        return false;
    }
    return true;
}

const LpCoef GenericCliqueCutConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
    if (!icPtr->isTypeOf(VcId::CliqueCutConstrMask))
        return LpCoef(false, 0.0);

    return getMastColumnCoeff(static_cast<CliqueCut *>(icPtr), colPtr);
}

void GenericCliqueCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
    iconstrPtr->presetMembership(true);
    return;
}

void GenericCliqueCutConstr::printColsToConstrsMatrix(int spConfNum, std::vector<std::vector<int> > & colsToConstrsMatrix)
{
    int numVarsInClique = colsToConstrsMatrix.size();
    int mastConstrNum = colsToConstrsMatrix[0].size();
    
    for (int constrIndex = 0; constrIndex < mastConstrNum; constrIndex++)
    {
        if (constrIndex == spConfNum)
        {
            for (int varIndex = 0; varIndex < numVarsInClique; varIndex++)
                std::cout << "-\t";
            std::cout << std::endl;
        }
        
        for (int varIndex = 0; varIndex < numVarsInClique; varIndex++)
            std::cout << colsToConstrsMatrix[varIndex][constrIndex] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

const LpCoef GenericCliqueCutConstr::getMastColumnCoeff(CliqueCut * cutPtr, MastColumn * colPtr) const
{
    const std::vector<int> & orderedArcIds = colPtr->spSol()->orderedIds();
    const ProbConfig * probConfPtr = colPtr->cgSpConfPtr();
    const NetworkFlow * netFlowPtr = probConfPtr->networkFlowPtr();

    LpCoef lpCoeff(false, 0.0);
    /// only an enumerated column may have a non-zero coefficient in a clique cut
    if ((netFlowPtr == NULL) || !colPtr->spSol()->enumeratedFlag())
        return lpCoeff;

    bool setIdIsPresent[_nbPackingSets];
    std::memset(setIdIsPresent, 0, _nbPackingSets * sizeof(bool));

    std::vector<int>::const_iterator arcIdIt = orderedArcIds.begin();
    const BcArcInfo * arcInfoPtr = netFlowPtr->getArcInfoPtr(*arcIdIt);
    if ((arcInfoPtr->tailPackSetId >= 0) && (arcInfoPtr->tailPackSetId < _nbPackingSets))
        setIdIsPresent[arcInfoPtr->tailPackSetId] = true;
    while (arcIdIt != orderedArcIds.end())
    {
        arcInfoPtr = probConfPtr->networkFlowPtr()->getArcInfoPtr(*arcIdIt);
        int setId = (arcInfoPtr->headPackSetId >= 0) ? arcInfoPtr->headPackSetId : arcInfoPtr->arcPackSetId;
        if ((setId >= 0) && (setId < _nbPackingSets))
            setIdIsPresent[setId] = true;
        ++arcIdIt;
    }

    for (std::vector<std::vector<int> >::const_iterator intVectIt = cutPtr->_setIds.begin();
         intVectIt != cutPtr->_setIds.end(); ++intVectIt)
    {
        bool allSetIdsArePresent = true;
        std::vector<int>::const_iterator setIdIt = intVectIt->begin();
        while ((setIdIt != intVectIt->end()) && allSetIdsArePresent)
        {
           allSetIdsArePresent = allSetIdsArePresent && setIdIsPresent[*setIdIt];
           ++setIdIt;
        }
        if (allSetIdsArePresent)
        {
            LpCoef lpCoeff(true, 1.0);
            return lpCoeff;
        }
    }
    return lpCoeff;
}

void GenericCliqueCutConstr::minimizeColumnToSetIdMatrix(std::vector<std::vector<bool> > & columnToSetIdMatrix,
                                                         std::set<std::vector<int> > & cutSetIds)
{
    if (columnToSetIdMatrix.empty())
        return;
    
    //std::cout << "columnToSetIdMatrix before" << std::endl;
    //printColsToConstrsMatrix(spConfNum, columnToSetIdMatrix);
    
    int numColsInClique = columnToSetIdMatrix.size();

    /// count the number of columns containing each setId
    int maxNbColumns = 0;
    std::vector<int> nbColumns(_nbPackingSets, 0);
    for (int setId = 0; setId < _nbPackingSets; setId++)
    {
        for (int columnId = 0; columnId < numColsInClique; columnId++)
            if (columnToSetIdMatrix[columnId][setId])
                nbColumns[setId]++;
        if (nbColumns[setId] > maxNbColumns)
            maxNbColumns = nbColumns[setId];
    }
    
    // bucket sort the constraints by number of vars with non-zero coeff
    std::vector< std::vector<int> > bucket(maxNbColumns + 1);
    for (int setId = _nbPackingSets - 1; setId >= 0; setId--)
        bucket[nbColumns[setId]].push_back(setId);
    
    // try to remove the vertices from the cut starting from the ones with fewer routes
    for (int bucketId = 1; bucketId <= maxNbColumns; bucketId++)
    {
        for (int setOrder = 0; setOrder < (int)bucket[bucketId].size(); setOrder++)
        {
            int setId = bucket[bucketId][setOrder];
            for (int columnId = 0; columnId < numColsInClique; columnId++)
                if (columnToSetIdMatrix[columnId][setId])
                {
                    /// try to remove the assignment of the set to the column
                    bool canBeRemoved = true;
                    if (nbColumns[setId] > 1)
                    {
                        /// check every other column containing this packing set
                        for (int otherColumnId = 0; (otherColumnId < numColsInClique) && canBeRemoved; otherColumnId++)
                        {
                            if ((columnId == otherColumnId) || !columnToSetIdMatrix[otherColumnId][setId])
                                continue;
                            bool otherConflictingSetFound = false;
                            // try to find other packing sets that ensure the conflict
                            for (int otherSetId = 0; (otherSetId < _nbPackingSets) && !otherConflictingSetFound;
                                 otherSetId++)
                            {
                                if ((setId != otherSetId) && columnToSetIdMatrix[columnId][otherSetId]
                                    && columnToSetIdMatrix[otherColumnId][otherSetId])
                                    otherConflictingSetFound = true;
                            }
                            if (!otherConflictingSetFound)
                                canBeRemoved = false;
                        }
                    }
                    
                    if (canBeRemoved)
                    {
                        columnToSetIdMatrix[columnId][setId] = false;
                        nbColumns[setId]--;
                    }
                }
        }
    }

    /// we now remove the same columns from the matrix
    for (int colId = 0; colId < (int)columnToSetIdMatrix.size(); ++colId)
    {
        std::vector<int> localSetIds;
        for (int setId = 0; setId < _nbPackingSets; ++setId)
        {
            if (columnToSetIdMatrix[colId][setId])
                localSetIds.push_back(setId);
        }
        cutSetIds.insert(localSetIds);
    }

    //    std::cout << "columnToSetIdMatrix after" << std::endl;
    //    printColsToConstrsMatrix(spConfNum, columnToSetIdMatrix);

//    if (printL(0)) {
//        std::cout << " with packing sets {";
//        for (int setId = 0; setId < _nbPackingSets; setId++)
//            if (nbColumns[setId] > 0)
//                std::cout << setId << ", ";
//        std::cout << "}" << std::endl;
//    }
    
}

void GenericCliqueCutConstr::cutSeparationRoutine(const VarPtrSet & curSol,
                                                  std::multiset < InstanciatedConstr * ,
                                                                  CutSeparationPriorityComp > & generatedCutConstrSet)
{
    /// we first check whether there is at least one enumerated subproblem, if not, we quit immediately
    bool atLeastOneEnumeratedSubproblem = false;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->enumeratedStatus())
            atLeastOneEnumeratedSubproblem = true;
    }

    if (!atLeastOneEnumeratedSubproblem)
        return;

    /// we create data structures needed for the clique cut separation object
    std::vector<std::vector<int> > cGraphVertices(_nbPackingSets);
    std::vector<std::vector<int> > setIdsOfColumn;
    std::vector<double> varValues;

    /// we fill the these data structures
    int columnId = 0;
    for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
        if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;
        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);

        if (!colPtr->spSol()->enumeratedFlag())
            continue;

        setIdsOfColumn.push_back(std::vector<int>());

        const std::vector<int> & orderedArcIds = colPtr->spSol()->orderedIds();
        const ProbConfig * probConfPtr = colPtr->cgSpConfPtr();
        const NetworkFlow * netFlowPtr = probConfPtr->networkFlowPtr();
        std::vector<int>::const_iterator arcIdIt = orderedArcIds.begin();
        const BcArcInfo * arcInfoPtr = netFlowPtr->getArcInfoPtr(*arcIdIt);
        if ((arcInfoPtr->tailPackSetId >= 0) && (arcInfoPtr->tailPackSetId < _nbPackingSets))
        {
            cGraphVertices[arcInfoPtr->tailPackSetId].push_back(columnId);
            setIdsOfColumn[columnId].push_back(arcInfoPtr->tailPackSetId);
        }
        while (arcIdIt != orderedArcIds.end())
        {
            arcInfoPtr = probConfPtr->networkFlowPtr()->getArcInfoPtr(*arcIdIt);
            int setId = (arcInfoPtr->headPackSetId >= 0) ? arcInfoPtr->headPackSetId : arcInfoPtr->arcPackSetId;
            if ((setId >= 0) && (setId < _nbPackingSets))
            {
                cGraphVertices[setId].push_back(columnId);
                setIdsOfColumn[columnId].push_back(setId);
            }
            ++arcIdIt;
        }

        varValues.push_back((*varPtrIt)->val());
        columnId++;
    }


    /// we create the clique cut separation object
    int numEnumColumnsInSolution = columnId;
    CglEClique cliqueSep;
    cliqueSep.resetCGraph(numEnumColumnsInSolution);
    for (int setId = 0; setId < _nbPackingSets; setId++) {
        if (!cGraphVertices[setId].empty())
            cliqueSep.addCliqueToCGraph(cGraphVertices[setId]);
    }
    
    /// we call the cut separation
    std::vector< std::vector<double> > coeff;
    std::vector<double> RHS;

    cliqueSep.setGenOddHoles(false);
    cliqueSep.generateCuts(varValues, coeff, RHS, _separationMaxTime, printL(0));

    std::set<std::pair<double, int> > bestCutIndices;
    for (int cutIndex = 0; cutIndex < (int)coeff.size(); cutIndex++)
    {
        /// for each generated cut
        if (RHS[cutIndex] <= 1.0)
        {
            double lhs = 0.0;
            for (int colId = 0; colId < numEnumColumnsInSolution; ++colId)
            {
                if (coeff[cutIndex][colId] > 0)
                {
                    lhs += varValues[colId];
                }
            }
            bestCutIndices.insert(std::make_pair(1.0 - lhs, cutIndex));
        }
    }

    int numAddedCuts = 0;
    for (std::set<std::pair<double, int> >::iterator pairIt = bestCutIndices.begin(); pairIt != bestCutIndices.end();
            ++pairIt)
    {
        int cutIndex = pairIt->second;
        std::vector<std::vector<bool> > columnToSetIdMatrix;
        for (int colId = 0; colId < numEnumColumnsInSolution; ++colId)
        {
            if (coeff[cutIndex][colId] > 0)
            {
                columnToSetIdMatrix.push_back(std::vector<bool>(_nbPackingSets, false));
                for (std::vector<int>::iterator setIdIt = setIdsOfColumn[colId].begin();
                     setIdIt != setIdsOfColumn[colId].end(); ++setIdIt)
                    columnToSetIdMatrix.back()[*setIdIt] = true;
            }
        }

//            if (printL(0))
//              std::cout << "Found clique cut " << _numGeneratedCuts + 1 << " with lhs " << lhs << " ";

        /// we try to reduce the number of involved constraints and the number of sub-columns coefficients
        std::set<std::vector<int> > cutSetIds;
        minimizeColumnToSetIdMatrix(columnToSetIdMatrix, cutSetIds);

        MultiIndex newCutId(_numGeneratedCuts++);
        std::string cutName(defaultName());
        newCutId.appendRef2name(cutName, multiIndexNames());
        CliqueCut * newCutPtr = new CliqueCut(newCutId, this, probConfPtr(), cutName, cutSetIds);

        generatedCutConstrSet.insert(newCutPtr);

        if (++numAddedCuts >= _maxNumCutsPerRound)
            break;
    }
}

#endif /* CLIQUE_SEP_IS_FOUND */

#endif /* BCP_RCSP_IS_FOUND */
