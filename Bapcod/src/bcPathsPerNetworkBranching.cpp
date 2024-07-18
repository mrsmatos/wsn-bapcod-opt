/**
 *
 * This file bcPathsPerNetworkBranching.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcPathsPerNetworkBranching.cpp
//  Project
//
//  Created by Ruslan Sadykov on 18/09/2017.
//
//

#include "bcNetworkBasedBranching.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcModelNetworkFlow.hpp"

GenPathsPerNetworkBranchingConstr
::GenPathsPerNetworkBranchingConstr(Model * modelPtr,
			                        ProbConfig * probConfPtr,
			                        const std::string & name,
                                    const SelectionStrategy & priorityRule,
                                    const Double & nonRootPriorityLevel,
                                    const Double & rootPriorityLevel,
			                        const bool & toBeUsedInPreprocessing) :
  GenericBranchingConstr(modelPtr, probConfPtr, name, priorityRule,
                         nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing),
  _numGeneratedBrConstrs(0), _branchingVars()
{
}

GenPathsPerNetworkBranchingConstr::~GenPathsPerNetworkBranchingConstr()
{
}

/// we use the following definitions
/// 1) "appropriate" variable - all arcs associated with this variable are associated with coefficient 1 and
///    incident to the source or to the sink, but not to the both
/// he we verify the condition which should be satisfied to apply this type of branching
/// - every arc incident to the source or to the sink should be associated to exactly one appropriate variable
bool GenPathsPerNetworkBranchingConstr::prepareSeparation()
{
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
      NetworkFlow * networkPtr = (*cgSpConfPtrIt)->networkFlowPtr();
      if (networkPtr == NULL)
        continue;

      if (networkPtr->sourceList().empty())
        {
          std::cerr << "BaPCod PathsPerNetwork branching error : source vertex is not defined for the network "
                       "associated with subproblem " << (*cgSpConfPtrIt)->name() << std::endl;
          return false;
        }

      if (networkPtr->sinkList().empty())
        {
          std::cerr << "BaPCod PathsPerNetwork branching error : sink vertex is not defined for the network "
                       "associated with subproblem " << (*cgSpConfPtrIt)->name() << std::endl;
          return false;
        }

      int sourceVertId = networkPtr->netVertexPtr(networkPtr->sourceList().front())->id();
      int sinkVertId = networkPtr->netVertexPtr(networkPtr->sinkList().front())->id();

      std::pair<std::map<const ProbConfig *, std::set<InstanciatedVar *> >::iterator, bool> insertResult;
      insertResult = _branchingVars.insert(std::make_pair(*cgSpConfPtrIt, std::set<InstanciatedVar *>()));

      const VarIndexManager & probVarPts = (*cgSpConfPtrIt)->probPtr()->probVarSet();
      for (VarIndexManager::const_iterator varPtrIt = probVarPts.begin(VcIndexStatus::Active, 's');
           varPtrIt != probVarPts.end(VcIndexStatus::Active, 's'); ++varPtrIt)
        if ((*varPtrIt)->isTypeOf(VcId::InstanciatedVarMask))
          {
            bool associatedToAnArc = false;
            bool appropriateVar = true;
            InstanciatedVar * instVarPtr = static_cast<InstanciatedVar *>(*varPtrIt);
            for (std::map<int, double>::const_iterator mapIt = instVarPtr->arcIdToCoeff().begin();
                 mapIt != instVarPtr->arcIdToCoeff().end(); ++mapIt)
              {
                const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(mapIt->first);
                if (arcInfoPtr == NULL)
                  continue;
                associatedToAnArc = true;
                if (mapIt->second != 1.0)
                  appropriateVar = false;
                if ((arcInfoPtr->headVertId != sinkVertId) && (arcInfoPtr->tailVertId != sourceVertId))
                  appropriateVar = false;
                if ((arcInfoPtr->headVertId == sinkVertId) && (arcInfoPtr->tailVertId == sourceVertId))
                  appropriateVar = false;
              }
            if (associatedToAnArc && appropriateVar)
              insertResult.first->second.insert(instVarPtr);
          }

      for (lemon::ListDigraph::ArcIt lemonArc(networkPtr->digraph()); lemonArc != lemon::INVALID; ++lemonArc)
        {
          NetworkArc * netArcPtr = networkPtr->netArcPtr(lemonArc);
          const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(netArcPtr->id());
          if ((arcInfoPtr->tailVertId == sourceVertId) || (arcInfoPtr->headVertId == sinkVertId))
            {
                auto & mapping = netArcPtr->varToCoeffMaps().front();
                if (mapping.empty())
                {
                  std::cerr << "BaPCod PathsPerNetwork branching error : arc " << arcInfoPtr->tailVertId
                            << arcInfoPtr->tailVertId << "->" << arcInfoPtr->headVertId
                            << " is not associated to a variable" << std::endl;
                  return false;
                }

              int numAssociatedVars = 0;
              std::set<InstanciatedVar *>::iterator setIt;
              for (setIt = insertResult.first->second.begin(); setIt != insertResult.first->second.end(); ++setIt)
                for (auto & pair : mapping)
                  {
                    if (*setIt == pair.first)
                      numAssociatedVars += 1;
                  }
              
              if (numAssociatedVars != 1)
                {
                  std::cerr << "BaPCod PathsPerNetwork branching error : number of appropriate variables "
                            << " associated to arc " << arcInfoPtr->tailVertId << "->" << arcInfoPtr->headVertId
                            << " is different from one" << std::endl;
                  return false;
                }
            }
        }
    }
  return true;
}

void GenPathsPerNetworkBranchingConstr
     ::branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                         const MasterColSolution & curListOfMasterCol,
                                         const int & maxNumOfCandidates,
		  			                     BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  if (probConfPtr() == NULL)
    return;

  std::map<const ProbConfig *, double> lhsMap;
  std::pair<std::map<const ProbConfig *, double>::iterator, bool> mapInsertResult;

  MasterVarSolution::const_iterator solIt;
  for (solIt = curListOfMastAndSubprobVar.begin(); solIt != curListOfMastAndSubprobVar.end(); ++solIt)
    if (solIt->first->isTypeOf(VcId::SubProbVariableMask))
      {
        SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(solIt->first);
        std::map<const ProbConfig *, std::set<InstanciatedVar *> >::iterator mapIt;
        mapIt = _branchingVars.find(spVarPtr->cgSpConfPtr());
        if ((mapIt != _branchingVars.end()) && (mapIt->second.find(spVarPtr) != mapIt->second.end()))
          {
            double varValue = solIt->second._value;
            mapInsertResult = lhsMap.insert(std::make_pair(spVarPtr->cgSpConfPtr(), 0.5 * varValue));
            if (!mapInsertResult.second)
              (mapInsertResult.first)->second += 0.5 * varValue;
          }
      }

  /// this code is non-deterministic : sorting by pointers to ProbConfig
  /// TO DO : make it deterministic
  std::set<std::pair<double, const ProbConfig * > > sortedIdsByPriority;
  std::map<const ProbConfig *, double>::iterator ivMapIt;

  for (ivMapIt = lhsMap.begin(); ivMapIt != lhsMap.end(); ++ivMapIt)
    {
      double intPart;
      double fracPart = std::modf(ivMapIt->second, &intPart);
      if ((fracPart < param().BapCodIntegralityTolerance()) || (fracPart > 1 - param().BapCodIntegralityTolerance()))
        continue;
      
      double priority = std::abs(fracPart - 0.5);
      sortedIdsByPriority.insert(std::make_pair(priority, ivMapIt->first));
    }
    
  std::set<std::pair<double, const ProbConfig *> >::iterator setIt;
  int candCounter = 0;
  for (setIt = sortedIdsByPriority.begin(); setIt != sortedIdsByPriority.end(); ++setIt)
    {
      if (++candCounter > maxNumOfCandidates)
        break;
      
      MultiIndex constrId(_numGeneratedBrConstrs++);
      InstMasterConstr * prototypeConstrPtr = new InstMasterConstr(constrId, this, probConfPtr(), defaultName(),
                                                                   0, defaultSense());
      
      std::map<const ProbConfig *, std::set<InstanciatedVar *> >::iterator varMapIt;
      std::set<InstanciatedVar *>::iterator varSetIt;

      varMapIt = _branchingVars.find(setIt->second);
      for (varSetIt = varMapIt->second.begin(); varSetIt != varMapIt->second.end(); ++varSetIt)
        prototypeConstrPtr->includeMember(*varSetIt, 0.5, true);
      
      std::string description(defaultName());
      MultiIndex spId = setIt->second->id();
      BranchingConstrGenerator * gentorPtr = new BranchingConstrGenerator(this, 'U', 0.0, prototypeConstrPtr,
                                                                          spId.appendRef2name(description));
      generatedBrConstrGeneratorSet.insert(gentorPtr);
    }
    
  return;
}
