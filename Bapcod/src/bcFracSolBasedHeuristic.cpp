/**
 *
 * This file bcFracSolBasedHeuristic.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcFracSolBasedHeuristic.hpp"
#include "bcAlg4ProblemSetup.hpp"

#define FracSolBasedHeuristicPrintLevel 3

bool Algorithm4FracSolBasedHeuristicEval::eval()
{
  if (printL(0))
    std::cout << "Started frac. solution based heuristic." << std::endl;
  
  /// we first try to complete the current partial solution with a heuristic solution for the residual problem
  if (_probPtr->runFracSolBasedHeuristicFunctor(_heurPtr->_partialSol, _heurPtr->_masterSol)
      && checkIfCurSolIsMasterLpFeasible())
    {
      _algCurLpPrimalBound = Bound(totalObjVal(), _masterCommons.objStatus());
        updateAlgPrimalLpBounds();
      if (checkIfCurSolIsInteger() && !addCutToMaster('C'))
        {
          std::cout << "Frac. solution based heuristic returned solution of value "
                    << _algCurLpPrimalBound << std::endl;
          updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
        }
    }
  return false;
}

bool FracSolBasedHeuristic::prepareNodeForTreatment(Node * nodePtr, const int globalTreatOrder)
{
  nodePtr->setEvalAlg(new Algorithm4FracSolBasedHeuristicEval(_problemPtr, _masterCommons.masterCommons4EvalAlg(),
                                                              this));

  nodePtr->setPreprocessor(NULL);

  if ((nodePtr->probSetupInfoPtr()->treatOrderId == globalTreatOrder)
      && !nodePtr->probSetupInfoPtr()->fullSetupIsObligatory)
    nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupBranchingOnly(_masterCommons.masterCommons4ProblemSetup()));
  else
    nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup()));

  nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(_masterCommons.masterCommons4ProblemSetup()));

  return true;
}

void FracSolBasedHeuristic::runBody(int & globalTreatOrder)
{
  ProblemSetupInfo * probSetupInfoPtr = _currentBaPNodePtr->probSetupInfoPtr();
  std::list<VariableSolInfo> & solInfoList = probSetupInfoPtr->masterPartialSolutionInfo;
  const SolutionVarInfoPtrList & varList = _currentBaPNodePtr->primalSol();

  _partialSol = std::list<VariableSolInfo>(_currentBaPNodePtr->probSetupInfoPtr()->masterPartialSolutionInfo);
  _masterSol = SolutionVarInfoPtrList(_currentBaPNodePtr->primalSol());

  std::list<BranchingConstrBaseType *> localNodeBrConstrList; /// this should be empty
  int newNodeRef = _masterCommons.masterCommons4GenChildNodes().getNodeCountAndIncreaseIt();
  Node * curNodePtr = new Node(newNodeRef, _currentBaPNodePtr, localNodeBrConstrList, NULL);

  prepareNodeForTreatment(curNodePtr, globalTreatOrder);

  Bound primalBound(_currentBaPNodePtr->nodeIncIpPrimalBound());
  curNodePtr->treat(globalTreatOrder, primalBound);
  
  if (curNodePtr->primalBoundIsUpdated())
    _currentBaPNodePtr->updateNodeIncPrimalSolution(curNodePtr->nodeIncIpPrimalSolPtr());
  
  _generatedNodes.push_back(curNodePtr); /// curNodePtr will be deleted in ~Alg4PrimalHeuristicOfNode()
}

FracSolBasedHeuristic::FracSolBasedHeuristic(Problem * probPtr, MasterCommons4PrimalHeuristic & masterCommons) :
    Alg4PrimalHeuristicOfNode(probPtr, masterCommons)
{
}

FracSolBasedHeuristic::~FracSolBasedHeuristic()
{
}

