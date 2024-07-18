/**
 *
 * This file bcMasterCommons.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

// bcMasterCommons.cpp
//

#include "bcMasterCommons.hpp"
#include "bcMasterConfC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

MasterCommons4ProblemSetup::MasterCommons4ProblemSetup(MasterConf & masterConf) :
    _masterConf(masterConf)
{
}

const std::list<Problem *> & MasterCommons4ProblemSetup::problemList() const
{
  return _masterConf.problemList();
}

const std::vector<ColGenSpConf *> & MasterCommons4ProblemSetup::colGenSubProbConfPts() const
{
  return _masterConf.colGenSubProbConfPts();
}

Solution * MasterCommons4ProblemSetup::debugSolution() const
{
    return _masterConf.getDebugSolution();
}

MasterCommons4PrimalHeuristic::MasterCommons4PrimalHeuristic(MasterConf & masterConf) :
    _masterConf(masterConf)
{
}

const std::list<Problem *> & MasterCommons4PrimalHeuristic::problemList() const
{
  return _masterConf.problemList();
}

const std::vector<ColGenSpConf *> & MasterCommons4PrimalHeuristic::colGenSubProbConfPts() const
{
  return _masterConf.colGenSubProbConfPts();
}


MasterCommons4ProblemSetup & MasterCommons4PrimalHeuristic::masterCommons4ProblemSetup()
{
  return _masterConf._masterCommons4ProblemSetup;
}

MasterCommons4EvalAlg & MasterCommons4PrimalHeuristic::masterCommons4EvalAlg()
{
  return _masterConf._masterCommons4EvalAlg;
}

MasterCommons4GenChildNodesAlgorithm & MasterCommons4PrimalHeuristic::masterCommons4GenChildNodes()
{
  return _masterConf._masterCommons4GenChildNodes;
}

MasterCommons4GenChildNodesAlgorithm::MasterCommons4GenChildNodesAlgorithm(MasterConf & masterConf) :
    _masterConf(masterConf)
{
}

const std::vector<ColGenSpConf *> & MasterCommons4GenChildNodesAlgorithm::colGenSubProbConfPts() const
{
  return _masterConf.colGenSubProbConfPts();
}

const std::set<GenericBranchingConstr *, DynamicGenConstrSort> & MasterCommons4GenChildNodesAlgorithm
                                                                 ::candidateBranchingGenericConstr() const
{
  return _masterConf.candidateBranchingGenericConstr();
}

const std::set<GenericCutConstr *, DynamicGenConstrSort> & MasterCommons4GenChildNodesAlgorithm
                                                           ::candidateCutGenericConstr() const
{
  return _masterConf.candidateCutGenericConstr();
}

bool MasterCommons4GenChildNodesAlgorithm::getAverageSubtreeSize(const int & depth, Double & averSubtreeSize) const
{
  return _masterConf.getAverageSubtreeSize(depth, averSubtreeSize);
}

const std::list<Problem *> & MasterCommons4GenChildNodesAlgorithm::problemList() const
{
  return _masterConf.problemList();
}

MasterCommons4ProblemSetup & MasterCommons4GenChildNodesAlgorithm::masterCommons4ProblemSetup()
{
  return _masterConf._masterCommons4ProblemSetup;
}

MasterCommons4EvalAlg & MasterCommons4GenChildNodesAlgorithm::masterCommons4EvalAlg()
{
  return _masterConf._masterCommons4EvalAlg;
}

int MasterCommons4GenChildNodesAlgorithm::getNodeCountAndIncreaseIt()
{
   int nodeCount = _masterConf.pcNodeCount();
   _masterConf.increasePCNodeCount();
   return nodeCount;
}

BcObjStatus::MinMaxIntFloat MasterCommons4GenChildNodesAlgorithm::objStatus() const
{
    return _masterConf.modelPtr()->objectiveSense();
}


MasterCommons4EvalAlg::MasterCommons4EvalAlg(MasterConf & masterConf) :
    _masterConf(masterConf)
{
}

const std::vector<ColGenSpConf *> & MasterCommons4EvalAlg::colGenSubProbConfPts() const
{
  return _masterConf._colGenSubProbConfPts;
}

long int MasterCommons4EvalAlg::totalNumberOfEnumeratedSubprobSolutions() const
{
  long int numEnumSolutions = 0;
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterConf._colGenSubProbConfPts.begin();
       spcPt != _masterConf._colGenSubProbConfPts.end(); ++spcPt)
  {
    if (!(*spcPt)->enumeratedStatus())
      return -1;
    numEnumSolutions += (*spcPt)->probPtr()->getNumberOfEnumeratedSolutions();
  }
  return numEnumSolutions;
}

const std::set<GenericBranchingConstr *, DynamicGenConstrSort> & MasterCommons4EvalAlg
                                                                 ::candidateBranchingGenericConstr() const
{
  return _masterConf.candidateBranchingGenericConstr();
}

BcObjStatus::MinMaxIntFloat MasterCommons4EvalAlg::objStatus() const
{
    return _masterConf.modelPtr()->objectiveSense();
}

InstanciatedConstr * MasterCommons4EvalAlg::castAndAddConstraint(InstanciatedConstr * iconstrPtr,
                                                                 bool const & insertImmediately)
{
  return _masterConf.castAndAddConstraint(iconstrPtr, insertImmediately);
}

std::list<Constraint *> & MasterCommons4EvalAlg::pcDelayedConstrPtrList()
{
  return _masterConf._pcDelayedConstrPtrList;
}

std::list<Constraint *> & MasterCommons4EvalAlg::pcConstrPtrList()
{
  return _masterConf._pcConstrPtrList;
}

const std::set<GenericCutConstr *, DynamicGenConstrSort> & MasterCommons4EvalAlg::candidateCutGenericConstr()
{
  return _masterConf.candidateCutGenericConstr();
}

Solution * MasterCommons4EvalAlg::debugSolution() const
{
    return _masterConf.getDebugSolution();
}

