/**
 *
 * This file bcRestrictedMasterIpHeuristicC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcRestrictedMasterIpHeuristicC.hpp"

RestrictedMasterIpHeuristic::RestrictedMasterIpHeuristic(Problem * probPtr,
    MasterCommons4PrimalHeuristic & masterCommons) :
    Alg4PrimalHeuristicOfNode(probPtr, masterCommons), _activateAllColumns(false),
    _localNodeBrConstrList()
{
}

RestrictedMasterIpHeuristic::~RestrictedMasterIpHeuristic()
{
}

Node * RestrictedMasterIpHeuristic::createRootNode()
{
  /// we make a "clone" of the node which launched the heuristic
  int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
  Node * nodePtr = new Node(newNodeRef, _currentBaPNodePtr, _localNodeBrConstrList, NULL);

  return nodePtr;
}

bool RestrictedMasterIpHeuristic::prepareNodeForTreatment(Node * nodePtr,
                                                     const int globalNodesTreatOrder)
{
  bool nodeShouldBeTreated = true;

  if (nodePtr->isToBePruned(_currentBaPNodePtr->nodeIncIpPrimalBound()))
  {
    nodePtr->prunedAtBeginningOfTreatNode(true);
    return nodeShouldBeTreated = false;
  }

  if (!nodePtr->solved())
  {
    Alg4EvalByMip * evalAlgPtr = new Alg4EvalByMip(_problemPtr, _masterCommons.masterCommons4EvalAlg());
    evalAlgPtr->setOptionNeedBeConqueredIfSolIsInteger(false);
    evalAlgPtr->setParamMaxTime(param().MaxTimeForRestrictedMasterIpHeur());
    evalAlgPtr->setParamExactSolution(false);
    evalAlgPtr->setParamRepeatIfCoreCutAdded(true);
    nodePtr->setEvalAlg(evalAlgPtr);

    if (param().ApplyPreprocessing())
      {
        nodePtr->setPreprocessor(new Algorithm4PreprocessingAtNodeOtherThanRoot(_masterCommons.problemList()));
      }

    /// set problem setup algorithm
    if ((nodePtr->probSetupInfoPtr()->treatOrderId == globalNodesTreatOrder) && !_activateAllColumns
        && !nodePtr->probSetupInfoPtr()->fullSetupIsObligatory)
      {
        /// the problem was not altered after solving the parent node, so we do not do problem setup
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupOfNode(_masterCommons.masterCommons4ProblemSetup()));
      }
    else
      {
        Alg4ProblemSetupFull * fullProblemSetupAlgPtr
          = new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup());
        fullProblemSetupAlgPtr->setOptionMakeAllColumnsActive(_activateAllColumns);
        nodePtr->setProblemSetupAlgorithm(fullProblemSetupAlgPtr);
      }

    /// set problem set down algorithm
    nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(_masterCommons.masterCommons4ProblemSetup()));
  }

  /// for now, no heuristics inside diving heuristic

  /// set the child nodes generation algorithm
  nodePtr->setGenChildNodesAlgorithm(NULL);

  return nodeShouldBeTreated;
}

/// Implementation of the diving heuristic is similar to the implementation of
/// the branch-and-price algorithm (which is in MasterConf for now),
/// so we might think of creation a generic SearchTreeAlgorithm class

void RestrictedMasterIpHeuristic::runBody(int & globalTreatOrder)
{
  if (printL(0))
    std::cout << "------------------------------------------------" << std::endl;
  if (printL(-1))
    std::cout << "---- Restricted Master IP Heuristic started ----" << std::endl;
  if (printL(0))
    std::cout << "------------------------------------------------" << std::endl;

  Node * curNodePtr = createRootNode(); //Only one node is needed for this heuristic

  if (prepareNodeForTreatment(curNodePtr, globalTreatOrder))
  {
    if (curNodePtr->treat(globalTreatOrder, _currentBaPNodePtr->nodeIncIpPrimalBound()) == false)
    {
      if (printL(0))
        std::cout << "ERROR: RestrictedMasterIpHeuristic is interrupted" << std::endl;
    }

    /// nodePtr->nodeIncIpDualBound() is not used here,
    /// as the dual bound of a heuristic node is not a valid global bound

    if (curNodePtr->primalBoundIsUpdated())
      _currentBaPNodePtr->updateNodeIncPrimalSolution(curNodePtr->nodeIncIpPrimalSolPtr());
  }
}

