/**
 *
 * This file bcNonPublicCapacityCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcCapacityCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 11/09/2017.
//
//

#include "bcNonPublicCuts.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "bcModelNetworkFlow.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif

#ifdef CVRPSEP_IS_FOUND
#include "basegrph.h"
#include "cnstrmgr.h"
#include "capsep.h"

const int MaxOldCutsMemory = 1000;

GenericCapacityCutConstr::GenericCapacityCutConstr(Model * modelPtr,
                                                   ProbConfig * probConfPtr,
                                                   const std::string & name,
                                                   const Double & nonRootPriorityLevel,
                                                   const Double & rootPriorityLevel,
                                                   const bool & isFacultative,
                                                   const int & maxNumCutsPerRound,
                                                   const int & maxCapacity,
                                                   const std::vector<int> & demands):
  GenericCutConstr(modelPtr, probConfPtr, name, isFacultative ? 'F' : 'C',
                   SelectionStrategy::MostViolated, nonRootPriorityLevel, rootPriorityLevel),
  _maxCapacity(maxCapacity), _demands(demands.size() + 1, 0), _myCutsCMP(NULL), _myOldCutsCMP(NULL), _cutCount(0),
  _cutRound(0), _maxNumCutsPerRound(maxNumCutsPerRound), _numElemSets(demands.size()),
  _elemSetPairVarPts(_numElemSets, std::vector<std::vector<InstanciatedVar *> >(_numElemSets,
                                                                                std::vector<InstanciatedVar *>()))
{
  if (_maxNumCutsPerRound > param().MaxNbOfCutGeneratedAtEachIter())
    _maxNumCutsPerRound = param().MaxNbOfCutGeneratedAtEachIter();
  CMGR_CreateCMgr(&_myCutsCMP, _maxNumCutsPerRound);
  CMGR_CreateCMgr(&_myOldCutsCMP, MaxOldCutsMemory);
  /// CVRPSEP uses 1 as the first index, index 0 is not used, so _demands are "shifted" by one index
  for (int elemSetId = 0; elemSetId < demands.size(); ++elemSetId)
    _demands[elemSetId + 1] = demands[elemSetId];
}

GenericCapacityCutConstr::~GenericCapacityCutConstr()
{
  CMGR_FreeMemCMgr(&_myOldCutsCMP);
  CMGR_FreeMemCMgr(&_myCutsCMP);
}

/// we say that a vertex belong to an elementarity set only if the demand of this set is positive
/// we use the following definitions
/// 1) "connecting" arc - arc incident to vertices both belonging to elementary sets
/// 2) "appropriate" variable - all connecting arcs associated with this variable are associated with coefficient 1 and
///    are incident to vertices belonging to the same pair of elementarity sets
/// here we verify conditions which should be satisfied to apply capacity cuts
/// - elementarity sets should be defined only on vertices
/// - every connecting arc should be associated to exactly one appropriate variable
bool GenericCapacityCutConstr::prepareSeparation()
{
  /// first, we find all appropriate variables and put them to _elemSetPairVarPts
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
      NetworkFlow * networkPtr = (*cgSpConfPtrIt)->networkFlowPtr();
      if (networkPtr == NULL)
        continue;

      const VarIndexManager & probVarPts = (*cgSpConfPtrIt)->probPtr()->probVarSet();
      for (VarIndexManager::const_iterator varPtrIt = probVarPts.begin(VcIndexStatus::Active, 's');
           varPtrIt != probVarPts.end(VcIndexStatus::Active, 's'); ++varPtrIt)
        if ((*varPtrIt)->isTypeOf(VcId::InstanciatedVarMask))
          {
            int smallerElemSetId = -1;
            int largerElemSetId = -1;
            bool appropriateVar = true;
            InstanciatedVar * instVarPtr = static_cast<InstanciatedVar *>(*varPtrIt);
            for (std::map<int, double>::const_iterator mapIt = instVarPtr->arcIdToCoeff().begin();
                 mapIt != instVarPtr->arcIdToCoeff().end(); ++mapIt)
              {
                const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(mapIt->first);
                if ((arcInfoPtr != NULL)
                    && (arcInfoPtr->tailElemSetId >= 0) && (arcInfoPtr->tailElemSetId < _numElemSets)
                    && (_demands[arcInfoPtr->tailElemSetId + 1] > 0)
                    && (arcInfoPtr->headElemSetId >= 0) && (arcInfoPtr->headElemSetId < _numElemSets)
                    && (_demands[arcInfoPtr->headElemSetId + 1] > 0)
						)
                  {
                    if (mapIt->second != 1.0)
                      {
                        appropriateVar = false;
                        if (printL(2))
						  std::cout << "Variable " << instVarPtr->name() << " is not appropriate due to coefficient "
						  		<< mapIt->second << std::endl;
                        break;
                      }
                    int thisSmallerElemSetId = arcInfoPtr->tailElemSetId;
                    int thisLargerElemSetId = arcInfoPtr->headElemSetId;
                    if (thisSmallerElemSetId > thisLargerElemSetId)
                      {
                        thisSmallerElemSetId = arcInfoPtr->headElemSetId;
                        thisLargerElemSetId = arcInfoPtr->tailElemSetId;
                      }
                   if (smallerElemSetId < 0)
                      {
                        smallerElemSetId = thisSmallerElemSetId;
                        largerElemSetId = thisLargerElemSetId;
                      }
                    else if ((smallerElemSetId != thisSmallerElemSetId) || (largerElemSetId != thisLargerElemSetId))
                      {
                        appropriateVar = false;
                        if (printL(2))
						  {
							std::cout << "Variable " << instVarPtr->name() << " is not appropriate since";
							std::cout << " mapped to pairs (" << smallerElemSetId << ", " << largerElemSetId;
							std::cout << ") and (" << thisSmallerElemSetId << ", "  << thisLargerElemSetId << ")" << std::endl;
						  }
                        break;
                      }
                  }
              }
            if (appropriateVar && (smallerElemSetId >= 0))
			  {
				_elemSetPairVarPts[smallerElemSetId][largerElemSetId].push_back(instVarPtr);
				if (printL(2))
				  std::cout << "Variable " << instVarPtr->name() << " is appropriate" << std::endl;
			  }
            else if (appropriateVar && printL(2))
			  std::cout << "Variable " << instVarPtr->name() << " connects " << smallerElemSetId << " to "
			  	<< largerElemSetId << std::endl;
          }
    }

  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
      NetworkFlow * networkPtr = (*cgSpConfPtrIt)->networkFlowPtr();
      if (networkPtr == NULL)
        continue;
      for (lemon::ListDigraph::ArcIt lemonArc(networkPtr->digraph()); lemonArc != lemon::INVALID; ++lemonArc)
        {
          NetworkArc * netArcPtr = networkPtr->netArcPtr(lemonArc);
          const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(netArcPtr->id());
          
          if (arcInfoPtr->arcElemSetId >= 0)
            {
              std::cerr << "BaPCod error : capacity cuts are not compatible with elementarity sets on arcs"
                        << std::endl;
              return false;
            }

           bool connectingArc = (arcInfoPtr->tailElemSetId >= 0) && (arcInfoPtr->tailElemSetId < _numElemSets)
                                && (_demands[arcInfoPtr->tailElemSetId + 1] > 0)
                                && (arcInfoPtr->headElemSetId >= 0) && (arcInfoPtr->headElemSetId < _numElemSets)
                                && (_demands[arcInfoPtr->headElemSetId + 1] > 0);
           if (!connectingArc)
             continue;

          if (netArcPtr->varToCoeffMap().empty())
            {
              std::cerr << "BaPCod error : capacity cuts cannot be used as arc " << arcInfoPtr->tailVertId
                        << "->" << arcInfoPtr->headVertId << " is not associated to a variable " << std::endl;
              return false;
            }

           int smallerElemSetId = arcInfoPtr->tailElemSetId;
           int largerElemSetId = arcInfoPtr->headElemSetId;
           if (smallerElemSetId > largerElemSetId)
             {
               smallerElemSetId = arcInfoPtr->headElemSetId;
               largerElemSetId = arcInfoPtr->tailElemSetId;
             }
           std::vector<InstanciatedVar *> & associatedVarPts = _elemSetPairVarPts[smallerElemSetId][largerElemSetId];
           std::vector<InstanciatedVar *>::iterator iVarPtrIt;
           int numAssociateVars = 0;
           for (iVarPtrIt = associatedVarPts.begin(); iVarPtrIt != associatedVarPts.end(); ++iVarPtrIt)
             for (std::map<InstanciatedVar *, double>::const_iterator mapIt = netArcPtr->varToCoeffMap().begin();
                  mapIt != netArcPtr->varToCoeffMap().end(); ++mapIt)
               {
                 if (*iVarPtrIt == mapIt->first)
                   numAssociateVars += 1;
               }

          if (numAssociateVars != 1)
            {
              std::cerr << "BaPCod error : capacity cuts cannot be used as the number of appropriate variables "
                        << "associated to arc " << arcInfoPtr->tailVertId << "->" << arcInfoPtr->headVertId
                        << " is different from one " << std::endl;
              return false;
            }
        }
    }
  return true;
}

inline double roundP(double d)
{
    return (floor(d*100000000 + 0.5) / 100000000);
}

void GenericCapacityCutConstr::cutSeparationRoutine(const VarPtrSet & curSol,
                                                    std::multiset < InstanciatedConstr * ,
                                                                    CutSeparationPriorityComp > & generatedCutConstrSet)
{
#ifdef BCP_RCSP_IS_FOUND

  _cutRound++;

  /// we need to aggregate edge variables values for different subproblems
  std::vector<std::vector<double> > edgeSol(_numElemSets + 1, std::vector<double>(_numElemSets + 1, 0.0));
//  std::vector<double> adjacentArcsToElement(numElements, 0.0);

  for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
      if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
        continue;

      MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
      if ((colPtr->spSol()->rcspPathPtr() == nullptr) || colPtr->spSol()->rcspPathPtr()->arcIds.empty())
          continue;

      ProbConfig * probConfPtr = colPtr->cgSpConfPtr();
      NetworkFlow * networkPtr = probConfPtr->networkFlowPtr();
      if (networkPtr == NULL)
        continue;

      for (const auto & arcId : colPtr->spSol()->rcspPathPtr()->arcIds)
        {
          const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(arcId);
          if (arcInfoPtr != NULL)
            {
              int tailElemId = arcInfoPtr->tailElemSetId + 1; /// remember that indices in CVRPSEP are shifted by 1
              if ((tailElemId <= 0) || (tailElemId > _numElemSets) || (_demands[tailElemId] == 0))
                tailElemId = 0;
              int headElemId = arcInfoPtr->headElemSetId + 1; /// remember that indices in CVRPSEP are shifted by 1
              if ((headElemId <= 0) || (headElemId > _numElemSets) || (_demands[headElemId] == 0))
                headElemId = 0;
              if ((tailElemId > 0) || (headElemId > 0))
                {
                  if (tailElemId <= headElemId)
                    edgeSol[tailElemId][headElemId] += roundP((*varPtrIt)->val());
                  else
                    edgeSol[headElemId][tailElemId] += roundP((*varPtrIt)->val());
                }
            }
        }
    }

  /// define the data structure needed for the cut separation routine
  std::vector<int> edgeTail(1, 0);
  std::vector<int> edgeHead(1, 0);
  std::vector<double> edgeX(1, 0.0);

  /// Store the information on the current LP solution in EdgeTail, EdgeHead, EdgeX.
  for (int firstElemId = 0; firstElemId <= _numElemSets; ++firstElemId)
    for (int secondElemId = 0; secondElemId <= _numElemSets; ++secondElemId)
      if (edgeSol[firstElemId][secondElemId] > 0.0)
        {
          edgeTail.push_back((firstElemId == 0) ? _numElemSets + 1 : firstElemId); /// _numElemSets + 1 means depot
          edgeHead.push_back((secondElemId == 0) ? _numElemSets + 1 : secondElemId);
          edgeX.push_back(edgeSol[firstElemId][secondElemId]);
        }

  /// Call separation. Pass the previously found cuts in MyOldCutsCMP:
  char integerAndFeasible = 0;
  double maxViolation = 0;
  CAPSEP_SeparateCapCuts(_numElemSets, &_demands[0], _maxCapacity, edgeX.size() - 1, &edgeTail[0], &edgeHead[0],
                         &edgeX[0], _myOldCutsCMP, _maxNumCutsPerRound, param().BapCodIntegralityTolerance(),
                         &integerAndFeasible, &maxViolation, _myCutsCMP);

  if (_myCutsCMP->Size == 0)
    return; /// No cuts found

  /// Read the cuts from MyCutsCMP, and add them to the LP
  for (int k = 0; k < _myCutsCMP->Size; k++)
    {
      if (printL(2))
      {
        std::cout << "New capacity cut " << _cutCount << " for set";
        for (int i = 1; i <= _myCutsCMP->CPL[k]->IntListSize; i++)
          std::cout << _myCutsCMP->CPL[k]->IntList[i] - 1 << " ";
        std::cout << ": ";
      }

      std::stringstream ss;
      ss << "CAP_" << _cutCount;
      MultiIndex newCutId(_cutCount++);
      InstMasterConstr * newCutPtr = new InstMasterConstr(newCutId, this, probConfPtr(), ss.str(),
                                                          _myCutsCMP->CPL[k]->RHS, 'L', defaultType(), defaultKind(),
                                                          defaultFlag());
      std::set<Variable *> varsInCut;
      for (int ii = 1; ii < _myCutsCMP->CPL[k]->IntListSize; ii++)
        for (int jj = ii+1; jj <= _myCutsCMP->CPL[k]->IntListSize; jj++)
          {
            int smallerElemSetId = _myCutsCMP->CPL[k]->IntList[ii] - 1; /// remember that indices in CVRPSEP
            int largerElemSetId = _myCutsCMP->CPL[k]->IntList[jj] - 1;  /// are shifted by 1
            std::vector<InstanciatedVar *> & associatedVarPts = _elemSetPairVarPts[smallerElemSetId][largerElemSetId];
            std::vector<InstanciatedVar *>::iterator iVarPtrIt;
            for (iVarPtrIt = associatedVarPts.begin(); iVarPtrIt != associatedVarPts.end(); ++iVarPtrIt)
              {
                varsInCut.insert(*iVarPtrIt);
                newCutPtr->includeMember(*iVarPtrIt, 1.0, true);
                if (printL(2))
                  std::cout << "+" << (*iVarPtrIt)->name();
              }
          }

      if (printL(2))
      {
        double lhs = 0.0;
        for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
        {
          if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;

          MastColumn *colPtr = static_cast<MastColumn *>(*varPtrIt);
          if (colPtr->spSol()->rcspPathPtr()->arcIds.empty())
            continue;

          for (VarPtr2DoubleMap::const_iterator mapIt = colPtr->spSol()->solVarValMap().begin();
               mapIt != colPtr->spSol()->solVarValMap().end(); ++mapIt)
            if (varsInCut.count(mapIt->first))
              lhs += (*varPtrIt)->val() * mapIt->second;
        }
        std::cout << " <= " << _myCutsCMP->CPL[k]->RHS << "(lhs is " << lhs << ")" << std::endl;
      }

      generatedCutConstrSet.insert(newCutPtr);
    }


  // Move the new cuts to the list of old cuts:
  for (int index = 0; index < _myCutsCMP->Size; index++)
    CMGR_MoveCnstr(_myCutsCMP, _myOldCutsCMP, index, 0);
  _myCutsCMP->Size = 0;
#endif
}

#endif
