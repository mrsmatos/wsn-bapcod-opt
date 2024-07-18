/**
 *
 * This file bcGreedyHeuristic.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcGreedyHeuristic.hpp"
#include "bcColGenSpConfC.hpp"

#define GreedyHeuristicPrintLevel 3

bool Algorithm4GreedyEval::setupAlgo(Node * nodePtr)
{
  if (_thereArePureMasterVariables)
    {
      if (Alg4EvalByMip::setupAlgo(nodePtr))
        return true;
    }
  else
    {
      if (Alg4EvalOfNode::setupAlgo(nodePtr))
        return true;
    }

  GreedyEvalInfo * greedyInfoPtr = dynamic_cast<GreedyEvalInfo *>(nodePtr->nodeEvalInfoPtr());

  bapcodInit().require(greedyInfoPtr != NULL,
                       "BaPCod error: NodeEvalInfo for GreedyEvalALg is not of type GreedyEvalInfo.");

  _subProblemIndex = greedyInfoPtr->subProblemIndex;

  if (nodePtr->depth() == 0)
    updateSubProbVarsCosts();

  _needBeConqueredIfSolIsInteger = false;

  return false;
}

NodeEvalInfo * Algorithm4GreedyEval::recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr)
{
  GreedyEvalInfo * greedyEvalInfoPtr = NULL;
  if (nodeEvalInfoPtr != NULL)
    {
      greedyEvalInfoPtr = dynamic_cast<GreedyEvalInfo *>(nodeEvalInfoPtr);

      bapcodInit().require(greedyEvalInfoPtr != NULL,
                           "BaPCod error: nodeEvalInfoPtr passed to GreedyEvalAlg::recordNodeEvalInfo "
                           "is not of type GreedyEvalInfo");

      greedyEvalInfoPtr->subProblemIndex = _subProblemIndex;
    }
  else
    greedyEvalInfoPtr = new GreedyEvalInfo(_subProblemIndex);

  return Alg4EvalOfNode::recordNodeEvalInfo(globalTreatOrder, greedyEvalInfoPtr);
}

void Algorithm4GreedyEval::setDownAlgo()
{
  if (_thereArePureMasterVariables)
      Alg4EvalByMip::setDownAlgo();
  else
      Alg4EvalOfNode::setDownAlgo();
}

bool Algorithm4GreedyEval::partialSolutionIsFeasible()
{
  for (std::vector<ColGenSpConf *>::const_iterator spcIt = _masterCommons.colGenSubProbConfPts().begin();
       spcIt != _masterCommons.colGenSubProbConfPts().end(); ++spcIt)
    if (*((*spcIt)->lowerBoundPtr()) > 0)
      return false;

  for (ConstrIndexManager::const_iterator constrPtrIt = _probPtr->probConstrSet().begin(Active, 's');
       constrPtrIt != _probPtr->probConstrSet().end(Active, 's'); ++constrPtrIt)
    if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask))
      {
        if (((*constrPtrIt)->sense() == 'E') && ((*constrPtrIt)->curRhs() != 0))
          return false;
        if (((*constrPtrIt)->sense() == 'G') && ((*constrPtrIt)->curRhs() > 0))
          return false;
        if (((*constrPtrIt)->sense() == 'L') && ((*constrPtrIt)->curRhs() < 0))
          return false;
      }

  /// TO DO: it remains to check whether all subproblem variables contain zero
  /// in their feasibility interval

  return true;
}

/// this work only for minimization problem for now
void Algorithm4GreedyEval::updateSubProbVarsCosts()
{
  /// first, we change the cost of the subproblem variables
  std::vector<ColGenSpConf *>::const_iterator spcIt;
  for (spcIt = _masterCommons.colGenSubProbConfPts().begin(); spcIt != _masterCommons.colGenSubProbConfPts().end();
       ++spcIt)
    {
      Double maxCost(-BapcodInfinity);
      Double minCost(BapcodInfinity);
      Problem * probPtr = (*spcIt)->probPtr();

      for (VarIndexManager::const_iterator varPtrIt =
          probPtr->probVarSet().begin(Active, 's');
          varPtrIt != probPtr->probVarSet().end(Active, 's'); ++varPtrIt)
        {
          const Double & cost = (*varPtrIt)->costrhs();
          if (maxCost < cost)
            maxCost = cost;
          if (minCost > cost)
            minCost = cost;
        }

      if (maxCost >= 0)
        {
          if (maxCost > 0)
            maxCost *= 1.1;
          else
            maxCost = 1;
          for (VarIndexManager::const_iterator varPtrIt =
              probPtr->probVarSet().begin(Active, 's');
              varPtrIt != probPtr->probVarSet().end(Active, 's'); ++varPtrIt)
            (*varPtrIt)->resetCurCostByValue((*varPtrIt)->costrhs() - maxCost);
        }
    }
}

void Algorithm4GreedyEval::recordSolution(ColGenSpConf * cgSpConfPtr, Solution * solPtr)
{
  MastColumn * colPtr = cgSpConfPtr->recordSubproblemSolution(solPtr, false, 2);
  int multiplicity = colPtr->maxValueInCurrentMasterProblem();
  if (multiplicity > 0)
    _currentNodePtr->localFixedSolution()->includeVar(colPtr,
        multiplicity, true);
}

bool Algorithm4GreedyEval::eval()
{
  if (_thereArePureMasterVariables)
    Alg4EvalByMip::eval();
  else
    {
      if (partialSolutionIsFeasible())
        {
          VarPtrSet emptySol;
          _algCurLpPrimalBound = Bound(_probPtr->partialSolutionValue(), _masterCommons.objStatus());
          updatePrimalIpSolAndBnds(emptySol, _probPtr->partialSolution());
        }
    }

  _currentNodePtr->clearLocalFixedSolution();

  /// we now determine which columns to fix and store them in _currentNodePtr->localFixedSolution()
  /// for the GreedyDiveAlgorithm

  if ((_subProblemIndex < 0) || (_subProblemIndex >= _masterCommons.colGenSubProbConfPts().size()))
    return false;

  std::vector<ColGenSpConf *>::const_iterator cgspIt = _masterCommons.colGenSubProbConfPts().begin();
  cgspIt += _subProblemIndex;

  int maxLevelOfSubProbRestriction(0);

  Problem * spProbPtr = (*cgspIt)->probPtr();
  spProbPtr->setPrimalLpBound(Bound::infPrimalBound(_masterCommons.objStatus()));
  spProbPtr->setDualBound(Bound::infDualBound(_masterCommons.objStatus()));
  if (spProbPtr->solveProb(maxLevelOfSubProbRestriction) > 0)
    {
        recordSolution(*cgspIt, spProbPtr->retrieveCurPrimalLpSol());
        for (std::list<Solution *>::const_iterator spSolPt = spProbPtr->recordedSolList().begin();
             spSolPt != spProbPtr->recordedSolList().end(); spSolPt++)
            recordSolution(*cgspIt, *spSolPt);
    }

  if (printL(GreedyHeuristicPrintLevel))
    _currentNodePtr->printFixedSolution(std::cout);

  int cgspCounter(0);
  do
    {
      ++cgspIt;
      _subProblemIndex += 1;
      cgspCounter += 1;
      if (cgspIt == _masterCommons.colGenSubProbConfPts().end())
        {
          cgspIt = _masterCommons.colGenSubProbConfPts().begin();
          _subProblemIndex = 0;
        }
    }
  while ((*((*cgspIt)->upperBoundPtr()) == 0)
      && (cgspCounter <= _masterCommons.colGenSubProbConfPts().size()));

  if (cgspCounter > _masterCommons.colGenSubProbConfPts().size())
    _subProblemIndex = -1;

  return false;
}

bool Algorithm4GreedyDive::setupAlgo(Node * nodePtr)
{
  return Alg4GenChildrenOfNode::setupAlgo(nodePtr);
}


void Algorithm4GreedyDive::run(int & globalTreatOrder)
{
  if (_currentNodePtr->localFixedSolution()->solVarValMap().empty())
    {
      /// Nothing is fixed, we stop the greedy heuristic
      return;
    }

  std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
  Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                    tmpLocalNodeBrConstrList,_currentNodePtr->localFixedSolution()->clone());
  _currentNodePtr->sons().push_back(newChildNodePtr);
}

void Algorithm4GreedyDive::setDownAlgo()
{
  Alg4GenChildrenOfNode::setDownAlgo();
}

bool GreedyHeuristic::prepareNodeForTreatment(Node * nodePtr, const int globalTreatOrder)
{
  nodePtr->setEvalAlg(new Algorithm4GreedyEval(_problemPtr,
                                               _masterCommons.masterCommons4EvalAlg(),
                                               _thereArePureMasterVariables));

  if (globalTreatOrder == 0)
    nodePtr->setPreprocessor(new Algorithm4PreprocessingAtRoot(_masterCommons.problemList()));
  else
    nodePtr->setPreprocessor(new Algorithm4PreprocessingInDive(_masterCommons.problemList()));

  if (globalTreatOrder == 0)
    nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupRootNode(
        _masterCommons.masterCommons4ProblemSetup()));
  else
    nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupOfNode(
        _masterCommons.masterCommons4ProblemSetup()));

  nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(
        _masterCommons.masterCommons4ProblemSetup()));

  nodePtr->setGenChildNodesAlgorithm(new Algorithm4GreedyDive(
        _masterCommons.masterCommons4GenChildNodes()));

  return true;
}

void GreedyHeuristic::runBody(int & globalTreatOrder)
{
  /// Specificity of this heuristic is that it's root node is passed as parameter, it is not created

  std::cout << "Initial greedy heuristic is started " << std::endl;

  Node * curNodePtr(_currentBaPNodePtr);

  while (curNodePtr != NULL)
    {
      prepareNodeForTreatment(curNodePtr, globalTreatOrder);

      if (curNodePtr->treat(globalTreatOrder,
                            _currentBaPNodePtr->nodeIncIpPrimalBound()) == false)
        {
          std::cout << "ERROR: Initial Greedy Heuristic is interrupted" << std::endl;
          break;
        }

      if (curNodePtr->primalBoundIsUpdated())
        {
          _currentBaPNodePtr->updateNodeIncPrimalSolution(curNodePtr->nodeIncIpPrimalSolPtr());
        }

      /// maximum one child node can be generated
      if (!curNodePtr->sons().empty())
        {
          curNodePtr = curNodePtr->sons().front();
          _generatedNodes.push_back(curNodePtr);
        }
      else
        curNodePtr = NULL;
    }

}

GreedyHeuristic::GreedyHeuristic(Problem * probPtr, MasterCommons4PrimalHeuristic & masterCommons,
                                 const bool pureMasterVars) :
    Alg4PrimalHeuristicOfNode(probPtr, masterCommons), _thereArePureMasterVariables(pureMasterVars)
{
}

GreedyHeuristic::~GreedyHeuristic()
{
}

