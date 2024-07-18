/**
 *
 * This file bcPackSetAssignBranching.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcPackSetAssignBranching.cpp
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

GenPackSetAssignBranchingConstr
::GenPackSetAssignBranchingConstr(Model * modelPtr,
			                      ProbConfig * probConfPtr,
			                      const std::string & name,
                                  const SelectionStrategy & priorityRule,
                                  const Double & nonRootPriorityLevel,
                                  const Double & rootPriorityLevel,
			                      const bool & toBeUsedInPreprocessing) :
  GenericBranchingConstr(modelPtr, probConfPtr, name, priorityRule,
                         nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing),
  _numGeneratedBrConstrs(0), _numElemSets(0), _incidBranchingVars(), _belongBranchingVars()
{
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
      NetworkFlow * networkPtr = (*cgSpConfPtrIt)->networkFlowPtr();
      if (networkPtr == NULL)
        continue;
      
      int thisNetworkNumElemSets = networkPtr->elementaritySetPts().size();
      if (_numElemSets < thisNetworkNumElemSets)
        _numElemSets = thisNetworkNumElemSets;
    }
}

GenPackSetAssignBranchingConstr::~GenPackSetAssignBranchingConstr()
{
}

/// we use the following definitions
/// 1) "appropriate" variable - all arcs associated with this variable are associated with coefficient 1 and
///    -- either belong to the same elem. set
///    -- or incident to the vertices belonging to the same pair of elem. sets
/// he we verify the condition which should be satisfied to apply this type of branching
/// - every arc incident to a vertex belonging to an elem. set every arc belonging to an elem. set is
///   associated to exactly one appropriate variable
bool GenPackSetAssignBranchingConstr::prepareSeparation()
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

      std::pair<std::map<const ProbConfig *, std::vector<std::set<InstanciatedVar *> > >::iterator, bool>
        incidInsertResult = _incidBranchingVars.insert(
          std::make_pair(*cgSpConfPtrIt,
                         std::vector<std::set<InstanciatedVar *> >(_numElemSets, std::set<InstanciatedVar *>())));

      std::pair<std::map<const ProbConfig *, std::vector<std::set<InstanciatedVar *> > >::iterator, bool>
        belongInsertResult = _belongBranchingVars.insert(
          std::make_pair(*cgSpConfPtrIt,
                         std::vector<std::set<InstanciatedVar *> >(_numElemSets, std::set<InstanciatedVar *>())));

      const VarIndexManager & probVarPts = (*cgSpConfPtrIt)->probPtr()->probVarSet();
      for (VarIndexManager::const_iterator varPtrIt = probVarPts.begin(VcIndexStatus::Active, 's');
           varPtrIt != probVarPts.end(VcIndexStatus::Active, 's'); ++varPtrIt)
        if ((*varPtrIt)->isTypeOf(VcId::InstanciatedVarMask))
          {
            int smallerElemSetId = -1;
            int largerElemSetId = -1;
            int arcElemSetId = -1;
            bool appropriateVar = true;
            InstanciatedVar * instVarPtr = static_cast<InstanciatedVar *>(*varPtrIt);
            for (std::map<int, double>::const_iterator mapIt = instVarPtr->arcIdToCoeff().begin();
                 mapIt != instVarPtr->arcIdToCoeff().end(); ++mapIt)
              {
                const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(mapIt->first);
                if (arcInfoPtr == NULL)
                  continue;
                
                if ((arcInfoPtr->arcElemSetId >= 0) && (arcInfoPtr->arcElemSetId < _numElemSets))
                  {
                    if (largerElemSetId >= 0)
                      appropriateVar = false;
                    else if ((arcElemSetId >= 0) && (arcElemSetId != arcInfoPtr->arcElemSetId))
                      appropriateVar = false;
                    else
                      arcElemSetId = arcInfoPtr->arcElemSetId;
                  }

                int thisSmallerElemSetId = arcInfoPtr->tailElemSetId;
                int thisLargerElemSetId = arcInfoPtr->headElemSetId;
                if (thisSmallerElemSetId > thisLargerElemSetId)
                  {
                    thisSmallerElemSetId = arcInfoPtr->headElemSetId;
                    thisLargerElemSetId = arcInfoPtr->tailElemSetId;
                  }
                
                if ((thisLargerElemSetId >= 0) && (thisLargerElemSetId < _numElemSets))
                  {
                    if ((largerElemSetId >= 0)
                        && ((largerElemSetId != thisLargerElemSetId) || (smallerElemSetId != thisSmallerElemSetId)))
                      appropriateVar = false;
                    else if (arcElemSetId >= 0)
                      appropriateVar = false;
                    else
                      {
                        largerElemSetId = thisLargerElemSetId;
                        smallerElemSetId = thisSmallerElemSetId;
                      }
                  }
              }
            if (appropriateVar && (arcElemSetId >= 0))
              belongInsertResult.first->second[arcElemSetId].insert(instVarPtr);
            if (appropriateVar && (smallerElemSetId >= 0))
              incidInsertResult.first->second[smallerElemSetId].insert(instVarPtr);
            if (appropriateVar && (largerElemSetId >= 0))
              incidInsertResult.first->second[largerElemSetId].insert(instVarPtr);
          }

      for (lemon::ListDigraph::ArcIt lemonArc(networkPtr->digraph()); lemonArc != lemon::INVALID; ++lemonArc)
        {
          NetworkArc * netArcPtr = networkPtr->netArcPtr(lemonArc);
          const BcArcInfo * arcInfoPtr = networkPtr->getArcInfoPtr(netArcPtr->id());
          if (arcInfoPtr->arcElemSetId >= 0 || arcInfoPtr->headElemSetId >= 0 || arcInfoPtr->headElemSetId >= 0)
          {
              auto & mapping = netArcPtr->varToCoeffMaps().front();
              if (mapping.empty())
                {
                  std::cerr << "BaPCod PackSetAssign branching error : arc " << arcInfoPtr->tailVertId
                            << arcInfoPtr->tailVertId << "->" << arcInfoPtr->headVertId
                            << " is not associated to a variable" << std::endl;
                  return false;
                }

              int numAssociatedArcVars = 1;
              int numAssociatedTailVars = 1;
              int numAssociatedHeadVars = 1;
              std::set<InstanciatedVar *>::iterator setIt;
              if (arcInfoPtr->arcElemSetId >= 0)
                {
                  numAssociatedArcVars = 0;
                  for (setIt = belongInsertResult.first->second[arcInfoPtr->arcElemSetId].begin();
                       setIt != belongInsertResult.first->second[arcInfoPtr->arcElemSetId].end(); ++setIt)
                    for (auto & pair : mapping)
                      {
                        if (*setIt == pair.first)
                          numAssociatedArcVars += 1;
                      }
                }
              if (arcInfoPtr->tailElemSetId >= 0)
                {
                  numAssociatedTailVars = 0;
                  for (setIt = incidInsertResult.first->second[arcInfoPtr->tailElemSetId].begin();
                       setIt != incidInsertResult.first->second[arcInfoPtr->tailElemSetId].end(); ++setIt)
                    for (auto & pair : mapping)
                      {
                        if (*setIt == pair.first)
                          numAssociatedTailVars += 1;
                      }
                }
              if (arcInfoPtr->headElemSetId >= 0)
                {
                  numAssociatedHeadVars = 0;
                  for (setIt = incidInsertResult.first->second[arcInfoPtr->headElemSetId].begin();
                       setIt != incidInsertResult.first->second[arcInfoPtr->headElemSetId].end(); ++setIt)
                    for (auto & pair : mapping)
                      {
                        if (*setIt == pair.first)
                          numAssociatedHeadVars += 1;
                      }
                }

              if ((numAssociatedArcVars != 1) || (numAssociatedTailVars != 1) || (numAssociatedHeadVars != 1))
                {
                  std::cerr << "BaPCod PackSetAssign branching error : number of appropriate variables "
                            << " associated to arc " << arcInfoPtr->tailVertId << "->" << arcInfoPtr->headVertId
                            << " is different from one" << std::endl;
                  return false;
                }

            }
        }
    }
  return true;
}

void GenPackSetAssignBranchingConstr::augmentLhs(const ProbConfig * probConfPtr,
                                                 std::map<InstanciatedVar *, double> & varValueMap,
                                                 std::vector<double> & lhsVector)
{
  std::map<InstanciatedVar *, double>::iterator varValMapIt;
  for (int elemSetId = 0; elemSetId < _numElemSets; ++elemSetId)
    {
      std::set<InstanciatedVar *> incidSet = _incidBranchingVars.find(probConfPtr)->second[elemSetId];
      std::set<InstanciatedVar *>::iterator incidSetIt = incidSet.begin();
      varValMapIt = varValueMap.begin();
      while ((incidSetIt != incidSet.end()) && (varValMapIt != varValueMap.end()))
        {
          if (*incidSetIt == varValMapIt->first)
            {
              lhsVector[elemSetId] += 0.5 * varValMapIt->second;
              ++incidSetIt;
              ++varValMapIt;
            }
          else if (*incidSetIt < varValMapIt->first)
            {
              ++incidSetIt;
            }
          else if (*incidSetIt > varValMapIt->first)
            {
              ++varValMapIt;
            }
        }
      std::set<InstanciatedVar *> belongSet = _belongBranchingVars.find(probConfPtr)->second[elemSetId];
      std::set<InstanciatedVar *>::iterator belongSetIt = belongSet.begin();
      varValMapIt = varValueMap.begin();
      while ((belongSetIt != belongSet.end()) && (varValMapIt != varValueMap.end()))
        {
          if (*belongSetIt == varValMapIt->first)
            {
              lhsVector[elemSetId] += varValMapIt->second;
              ++belongSetIt;
              ++varValMapIt;
            }
          else if (*belongSetIt < varValMapIt->first)
            {
              ++belongSetIt;
            }
          else if (*belongSetIt > varValMapIt->first)
            {
              ++varValMapIt;
            }
        }
   }
}

void GenPackSetAssignBranchingConstr
     ::branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                         const MasterColSolution & curListOfMasterCol,
                                         const int & maxNumOfCandidates,
		  			                     BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  if (probConfPtr() == NULL)
    return;

  std::map<const ProbConfig *, std::map<InstanciatedVar *, double> > probConfigSolMap;
  std::map<const ProbConfig *, std::map<InstanciatedVar *, double> >::iterator pcVarMapIt;
  std::map<const ProbConfig *, std::vector<double> > lhsVectMap;
  std::map<const ProbConfig *, std::vector<double> >::iterator lhsVectMapIt;

  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    if ((*cgSpConfPtrIt)->networkFlowPtr() != NULL)
      {
        probConfigSolMap.insert(std::make_pair(*cgSpConfPtrIt, std::map<InstanciatedVar *, double>()));
        lhsVectMap.insert(std::make_pair(*cgSpConfPtrIt, std::vector<double>(_numElemSets)));
      }
    
  /// first we distribute the variables in the solution to the subproblem maps
  MasterVarSolution::const_iterator solIt;
  for (solIt = curListOfMastAndSubprobVar.begin(); solIt != curListOfMastAndSubprobVar.end(); ++solIt)
    if (solIt->first->isTypeOf(VcId::SubProbVariableMask))
      {
        SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(solIt->first);
        pcVarMapIt = probConfigSolMap.find(spVarPtr->cgSpConfPtr());
        std::pair<std::map<InstanciatedVar *, double>::iterator, bool> solMapInsert;
        double varValue = solIt->second._value;
        solMapInsert = pcVarMapIt->second.insert(std::make_pair(spVarPtr, varValue));
        if (!solMapInsert.second)
          (solMapInsert.first)->second += varValue;
      }

  /// then, for each subproblem, we pass over its variables in the solution and fill lhsVectMap,
  /// which contains, for each pair (subproblem, elem. set), the cumulative solution value
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    if ((*cgSpConfPtrIt)->networkFlowPtr() != NULL)
      {
        pcVarMapIt = probConfigSolMap.find(*cgSpConfPtrIt);
        lhsVectMapIt = lhsVectMap.find(*cgSpConfPtrIt);
        augmentLhs(*cgSpConfPtrIt, pcVarMapIt->second, lhsVectMapIt->second);
      }

  /// now we sort solution values for pairs (subproblem, elem. set) in the "most fractional value first" order
  /// this code is non-deterministic : sorting by pointers to ProbConfig
  /// TO DO : make it deterministic
  std::set<std::pair<double, std::pair<const ProbConfig *, int> > > sortedPairsByPriority;
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
      lhsVectMapIt = lhsVectMap.find(*cgSpConfPtrIt);
      for (int elemSetId = 0; elemSetId < _numElemSets; ++elemSetId)
        {
          double thisPairValue = lhsVectMapIt->second[elemSetId];
          if ((thisPairValue > param().BapCodIntegralityTolerance())
              && (thisPairValue < 1.0 - param().BapCodIntegralityTolerance()))
            {
              sortedPairsByPriority.insert(std::make_pair(std::abs(thisPairValue - 0.5),
                                                          std::make_pair(*cgSpConfPtrIt, elemSetId)));
            }
        }
    }

  /// we loop over solution values for pairs (subproblem, elem. set) in the "most fractional value first" order
  /// for each pair, until the limit is reached, we create the branching constraint generator
  /// only one branching generator per elem. est is allowed
  std::vector<bool> elemSetPushed(_numElemSets, false);

  std::set<std::pair<double, std::pair<const ProbConfig *, int> > >::iterator setIt;
  int candCounter = 0;
  for (setIt = sortedPairsByPriority.begin(); setIt != sortedPairsByPriority.end(); ++setIt)
    {
      const ProbConfig * probConfigPtr = setIt->second.first;
      int elemSetId = setIt->second.second;
      if (elemSetPushed[elemSetId])
        continue;
      else
        elemSetPushed[elemSetId] = true;
      
      if (++candCounter > maxNumOfCandidates)
        break;
      
      MultiIndex constrId(_numGeneratedBrConstrs++);
      InstMasterConstr * prototypeConstrPtr = new InstMasterConstr(constrId, this, probConfPtr(), defaultName(),
                                                                   0, defaultSense());
      
      std::set<InstanciatedVar *> incidSet = _incidBranchingVars.find(probConfigPtr)->second[elemSetId];
      for (std::set<InstanciatedVar *>::iterator incSetIt = incidSet.begin(); incSetIt != incidSet.end(); ++incSetIt)
        prototypeConstrPtr->includeMember(*incSetIt, 0.5, true);
      
      std::set<InstanciatedVar *> belongSet = _belongBranchingVars.find(probConfigPtr)->second[elemSetId];
      for (std::set<InstanciatedVar *>::iterator belSetIt = belongSet.begin(); belSetIt != belongSet.end(); ++belSetIt)
        prototypeConstrPtr->includeMember(*belSetIt, 1.0, true);
      
      std::string description(defaultName());
      MultiIndex spId = probConfigPtr->id();
      spId += elemSetId;
      BranchingConstrGenerator * gentorPtr = new BranchingConstrGenerator(this, 'U', 0.0, prototypeConstrPtr,
                                                                          spId.appendRef2name(description));
      generatedBrConstrGeneratorSet.insert(gentorPtr);
    }
    
  return;
}
