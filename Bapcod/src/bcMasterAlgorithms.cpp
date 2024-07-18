/**
 *
 * This file bcMasterAlgorithms.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMasterConfC.hpp"
#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"

#include "bcOvfVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"

#include "bcAlg4EvalByLagrangianDuality.hpp"
#include "bcAlg4EvalBySimplexBasedColGen.hpp"
#include "bcAlg4EvalByColAndCutGen.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4GenChildrenInBranching.hpp"
#include "bcAlg4DivingHeuristic.hpp"
#include "bcFracSolBasedHeuristic.hpp"
#include "bcRestrictedMasterIpHeuristicC.hpp"
#include "bcGreedyHeuristic.hpp"

using namespace std;
using namespace VcIndexStatus;

void MasterConf::runGreedyHeuristic(int & globalNodesTreatOrder)
{

  if (!param().UseGreedyHeur())
    return;

  std::vector<Variable *> activeColumnsToInitializeMaster;
  std::vector<Variable *> inactiveColumnsToInitializeMaster;
  ProblemSetupInfo * problemSetupInfoPtr = new ProblemSetupInfo(0, activeColumnsToInitializeMaster,
                                                                inactiveColumnsToInitializeMaster);

  Node * greedyNodePtr = new Node(this, dualIncBound(), problemSetupInfoPtr, new GreedyEvalInfo());
  greedyNodePtr->nodeIncIpPrimalBound(primalIncBound());
  GreedyHeuristic * greedyHeuristicPtr = new GreedyHeuristic(_probPtr, _masterCommons4PrimalHeuristic,
                                                             !_setOfPureMasterVar.empty());
  greedyHeuristicPtr->run(greedyNodePtr, globalNodesTreatOrder);
  if (greedyNodePtr->primalBoundIsUpdated())
    updatePrimalIncSolution(greedyNodePtr->nodeIncIpPrimalBound(), greedyNodePtr->nodeIncIpPrimalSolPtr());

  delete greedyHeuristicPtr;
  delete greedyNodePtr;
}

Node * MasterConf::createRootNode()
{
  std::vector<Variable *> activeColumnsToInitializeMaster;
  std::vector<Variable *> inactiveColumnsToInitializeMaster;
  if ((param().mastInitMode().status() == MasterInitMode::incSolCol)
      || (param().mastInitMode().status() == MasterInitMode::incSolColAndGac)
      || (param().mastInitMode().status() == MasterInitMode::incSolColAndLac))
    {
      if (_primalSolPtr != NULL)
        for (VarPtr2DoubleMap::const_iterator varPt = _primalSolPtr->solVarValMap().begin();
             varPt != _primalSolPtr->solVarValMap().end(); varPt++)
          activeColumnsToInitializeMaster.push_back(varPt->first);
    }

  if (_initialSetOfActiveColumns != NULL)
    {
      for (VarPtr2DoubleMap::const_iterator varPt = _initialSetOfActiveColumns->solVarValMap().begin();
           varPt != _initialSetOfActiveColumns->solVarValMap().end(); varPt++)
        activeColumnsToInitializeMaster.push_back(varPt->first);
      delete _initialSetOfActiveColumns;
      _initialSetOfActiveColumns = NULL;
    }
  if (_initialSetOfInactiveColumns != NULL)
    {
      for (VarPtr2DoubleMap::const_iterator varPt = _initialSetOfInactiveColumns->solVarValMap().begin();
           varPt != _initialSetOfInactiveColumns->solVarValMap().end(); varPt++)
        inactiveColumnsToInitializeMaster.push_back(varPt->first);
      delete _initialSetOfInactiveColumns;
      _initialSetOfInactiveColumns = NULL;
    }

  ProblemSetupInfo * problemSetupInfoPtr = new ProblemSetupInfo(0, activeColumnsToInitializeMaster,
                                                                inactiveColumnsToInitializeMaster);
  StabilizationInfo * stabInfoPtr = new StabilizationInfo(probPtr(), param());
  LpBasisRecord * masterLpBasisPtr = new LpBasisRecord("Basis0");

  NodeEvalInfo * nodeEvalInfoPtr;

  switch (_probPtr->solMode().status()) // method to solve the restricted master
    {
    case SolutionMethod::lpSolver:
    case SolutionMethod::mipSolver:
    {
      nodeEvalInfoPtr = new ColGenEvalInfo(stabInfoPtr, masterLpBasisPtr, BapcodInfinity);
      break;
    }
    case SolutionMethod::customSolver:
    case SolutionMethod::custom2mipSolver:
    case SolutionMethod::undefined:
    case SolutionMethod::none:
    default:
    {
      nodeEvalInfoPtr = NULL;
      bapcodInit().check(true, "MasterConf::createRootNode(): ERROR unsupported solution method");
    }
    }
  return new Node(this, dualIncBound(), problemSetupInfoPtr, nodeEvalInfoPtr,
                  _debugSolutionPtr != NULL);
}

bool MasterConf::getAverageSubtreeSize(const int & depth, Double & averSubtreeSize) const
{
  if ((_subTreeSizeByDepth.size() <= depth) || (_subTreeSizeByDepth[depth].first == 0))
    return false;
  averSubtreeSize = _subTreeSizeByDepth[depth].second;
  return true;
}

/// this function sets the algorithms for the node treatment
/// return false in does not need to be treated
bool MasterConf::prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder,
                                         const int thisSearchTreeTreatedNodesNumber)
{
  bool nodeShouldBeTreated = true;

  //Issam: this if bloc should be out of this method.
  if (nodePtr->isToBePruned(primalIncBound()) || !progStatus().doRun())
    {
      nodePtr->prunedAtBeginningOfTreatNode(true);
      return nodeShouldBeTreated = false;
    }

  if (!nodePtr->solved())
    {
        /// set solver
        int colGenLogFrequency = param().ColGenLogFrequency();
        if (!printL(0) && (colGenLogFrequency < 10))
            colGenLogFrequency = 10;

        switch (_probPtr->solMode().status())
        {
            case SolutionMethod::lpSolver:
            case SolutionMethod::mipSolver:
            {
                /// set solver
                Alg4EvalOfNode * evalAlgPtr = new Alg4EvalByColAndCutGen(_probPtr, _masterCommons4EvalAlg);

                Alg4EvalByLagrangianDuality * colAndCutGenSolverPtr
                        = dynamic_cast<Alg4EvalByLagrangianDuality *>(evalAlgPtr);

                if (colAndCutGenSolverPtr != NULL)
                {
                    if (param().ReducedCostFixingThreshold() > 0)
                        colAndCutGenSolverPtr->setOptionDoRedCostFixingAndEnumeration(1);
                    else
                        colAndCutGenSolverPtr->setOptionDoRedCostFixingAndEnumeration(0);
                    colAndCutGenSolverPtr->setOptionMaxNbOfCgIterations(param().MaxNbOfCgIterations());
                    colAndCutGenSolverPtr->setOptionMaxNbOfPenaltyUpdates(param().ArtVarMaxNbOfPenaltyUpdates());
                    colAndCutGenSolverPtr->setOptionLogPrintFrequency(colGenLogFrequency);
                    colAndCutGenSolverPtr->setOptionMinNbOfCutRounds(param().MinNumOfCutRoundsBeforeStopBySp());
                }
                nodePtr->setEvalAlg(evalAlgPtr);

                if (nodePtr->isRoot() && (param().RCSPrankOneCutsMemoryType() == 0))
                {
                    AutoRankOneCutsMemoryInfo * infoPtr = new AutoRankOneCutsMemoryInfo(_modelPtr->objectiveSense());
                    Alg4EvalByColAndCutGen * evalAlgPtr = new Alg4EvalByColAndCutGen(_probPtr, _masterCommons4EvalAlg);
                    if (param().ReducedCostFixingThreshold() > 0)
                        evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(1);
                    else
                        evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(0);
                    evalAlgPtr->setOptionMaxNbOfCgIterations(param().MaxNbOfCgIterations());
                    evalAlgPtr->setOptionMaxNbOfPenaltyUpdates(param().ArtVarMaxNbOfPenaltyUpdates());
                    evalAlgPtr->setOptionLogPrintFrequency(colGenLogFrequency);
                    evalAlgPtr->setOptionMinNbOfCutRounds(param().MinNumOfCutRoundsBeforeStopBySp());
                    infoPtr->evalAlgPtr = evalAlgPtr;
                    infoPtr->problemSetupAlgPtr = new Alg4ProblemSetupFull(_masterCommons4ProblemSetup);
                    infoPtr->problemSetDownAlgPtr = new ProblemFullSetDownAlgorithm(_masterCommons4ProblemSetup);
                    nodePtr->applyAutoRankOneCutsMemorySearch(infoPtr);
                }

                break;
            }
            case SolutionMethod::customSolver:
            case SolutionMethod::custom2mipSolver:
            case SolutionMethod::undefined:
            case SolutionMethod::none:
            default:
                bapcodInit().check(true,"MasterConf::prepareNodeForTreatment(): ERROR undefined solution method");
        }


      /// set preprocessing algorithm
      if (param().ApplyPreprocessing())
        if (globalNodesTreatOrder == 0)
          nodePtr->setPreprocessor(new Algorithm4PreprocessingAtRoot(_problemList));
        else
          nodePtr->setPreprocessor(new Algorithm4PreprocessingAtNodeOtherThanRoot(_problemList));
      else
        nodePtr->setPreprocessor(NULL);

      /// set problem setup algorithm
      if (globalNodesTreatOrder == 0)
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupRootNode(_masterCommons4ProblemSetup));

      else if ((nodePtr->probSetupInfoPtr()->treatOrderId == globalNodesTreatOrder)
               && !nodePtr->probSetupInfoPtr()->fullSetupIsObligatory)
        /// the problem was not altered after solving the parent node,
        /// so we do setup of only branching constraints
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupBranchingOnly(_masterCommons4ProblemSetup));
      else
        nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons4ProblemSetup));

    }

  /// set the list of heuristic algorithms
  if (nodePtr->depth() <= param().DivingHeurUseDepthLimit())
    {
      if ((param().CallFrequencyOfDivingHeur() <= 0 && (nodePtr->depth() == 0)) ||
          ((param().CallFrequencyOfDivingHeur() > 0)
           && (thisSearchTreeTreatedNodesNumber % param().CallFrequencyOfDivingHeur() == 0)))
        {

          DivingHeuristic * divingHeuristicPtr = new DivingHeuristic(_probPtr, _masterCommons4PrimalHeuristic,
                                                                     param().RoundingColSelectionCriteria());
          divingHeuristicPtr->setOptionPriority(2);
          divingHeuristicPtr->setOptionMaxDepth(param().MaxLDSdepth());
          divingHeuristicPtr->setOptionMaxDiscrepancy(param().MaxLDSbreadth());
          divingHeuristicPtr->setOptionStopAfterFirstSolutionFound(param().DivingHeurStopsWithFirstFeasSol());
          nodePtr->addPrimalHeuristic(divingHeuristicPtr);
        }
    }

  if (nodePtr->depth() <= param().LocalSearchHeurUseDepthLimit())
    {
      LocalSearchHeuristic *
        localSearchHeuristicPtr = new LocalSearchHeuristic(_probPtr, _masterCommons4PrimalHeuristic,
                                                           param().LocalSearchColSelectionCriteria());
      localSearchHeuristicPtr->setOptionFixedVarsRatio(param().MaxFactorOfColFixedByLocalSearchHeur());
      localSearchHeuristicPtr->setOptionMaxIterationsNumber(param().MaxLocalSearchIterationCounter());
      localSearchHeuristicPtr->setOptionPriority(3);
      nodePtr->addPrimalHeuristic(localSearchHeuristicPtr);
    }

  if (param().UseCustomFracSolBasedHeur())
    {
      FracSolBasedHeuristic * fracSolBasedHeurPtr = new FracSolBasedHeuristic(_probPtr,
                                                                              _masterCommons4PrimalHeuristic);
      fracSolBasedHeurPtr->setOptionPriority(5);
      nodePtr->addPrimalHeuristic(fracSolBasedHeurPtr);
    }

  if (param().MaxTimeForRestrictedMasterIpHeur() > 0)
    {
      if ((param().CallFrequencyOfRestrictedMasterIpHeur() <= 0 && nodePtr->depth() == 0) ||
          ((param().CallFrequencyOfRestrictedMasterIpHeur() > 0)
           && (globalNodesTreatOrder % param().CallFrequencyOfRestrictedMasterIpHeur() == 0)))
        {
          if (param().RCSPmaxNumOfLabelsInHeurEnumeration() > 0)
          {
            /// diving heuristic is used here because the enumeration with false gap is implemented
            /// inside Algorithm4DivingEval
            DivingHeuristic * divingHeuristicPtr = new DivingHeuristic(_probPtr, _masterCommons4PrimalHeuristic,
                                                                       param().RoundingColSelectionCriteria());
            divingHeuristicPtr->setOptionPriority(6);
            divingHeuristicPtr->setOptionRunRestrMasterAfterFalseGapEnumeration(true);
            nodePtr->addPrimalHeuristic(divingHeuristicPtr);
          }
          else
          {
            RestrictedMasterIpHeuristic *restrictedMasterIpHeuristicPtr
                    = new RestrictedMasterIpHeuristic(_probPtr, _masterCommons4PrimalHeuristic);
            restrictedMasterIpHeuristicPtr->setOptionPriority(6);
            restrictedMasterIpHeuristicPtr->setOptionActivateAllColumns
                    (param().ActivateAllColumnsForRestrictedMasterIpHeur());
            nodePtr->addPrimalHeuristic(restrictedMasterIpHeuristicPtr);
          }
        }
    }

  /// set the child nodes generation algorithm
  if (nodePtr->depth() < param().MaxDepthInBBtree()
      && (thisSearchTreeTreatedNodesNumber + 1 < param().MaxNbOfBBtreeNodeTreated()))
    nodePtr->setGenChildNodesAlgorithm(new Alg4GenChildrenInBranching(_masterCommons4GenChildNodes));
  else
    nodePtr->setGenChildNodesAlgorithm(NULL);

  /// set problem set down algorithm
  /// this should be done after setting heuristic and gen. child. nodes algorithms
  if (!nodePtr->solved())
  {
    if (nodePtr->needProblemFullSetDownAlgorithm())
      nodePtr->setProblemSetDownAlgorithm(new ProblemFullSetDownAlgorithm(_masterCommons4ProblemSetup));
    else
      nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(_masterCommons4ProblemSetup));
  }

  return nodeShouldBeTreated;
}

/// solvePC function reimplemented by Ruslan and Issam
Solution * MasterConf::solvePC()
{
    if (param().UseInitialPrimalHeur())
    {
        Solution * solPtr = _modelPtr->initPrimalHeur();
        if (solPtr != nullptr)
        {
            solPtr->resetCost();
            Bound primalBound(solPtr->cost(), modelPtr()->objectiveSense());
            updatePrimalIncSolution(primalBound, solPtr);
            delete solPtr;
        }
    }

  int globalNodesTreatOrder = 0;
  runGreedyHeuristic(globalNodesTreatOrder);

  std::priority_queue<Node *, std::vector<Node *>, smallerNodePt> searchTree;

  /// Once the number of nodes in searchTree reaches openNodesLimit we put additional nodes in secondary_queue
  std::priority_queue<Node *, std::vector<Node *>, secondarySmallerNodePt> secondarySearchTree;

  Node * curNodePtr = createRootNode();
  searchTree.push(curNodePtr);
  _bapTreeNodes.push_back(curNodePtr);

  int bapTreatOrder = 1; //usefull only for printing. (Issam)
  int thisSearchTreeTreatedNodesNumber = 0;
  initializeBaPTreeDotFile();

  Time bcTimeBaP;

  while (!searchTree.empty() && (thisSearchTreeTreatedNodesNumber++ < param().MaxNbOfBBtreeNodeTreated()))
    {
      if (!progStatus().doRun())
        break;

      bool thisIsPrimaryTreeNode = secondarySearchTree.empty();
      if (secondarySearchTree.empty())
        {
          curNodePtr = searchTree.top();
          searchTree.pop();
        }
      else
        {
          curNodePtr = secondarySearchTree.top();
          secondarySearchTree.pop();
        }

      bapcodInit().statistics().incrCounter("bcCountNodeGen");

      bool curNodeSolvedBefore = curNodePtr->solved();
      if (prepareNodeForTreatment(curNodePtr, globalNodesTreatOrder,
                                  thisSearchTreeTreatedNodesNumber - 1))
        {
          if (!curNodeSolvedBefore)
            curNodePtr->branchAndPriceOrder(bapTreatOrder++);
          if (printL(-1))
            printInfoBeforeSolvingNode(curNodePtr, searchTree.size() + ((thisIsPrimaryTreeNode) ? 1 : 0),
                                       secondarySearchTree.size() + ((thisIsPrimaryTreeNode) ? 0 : 1),
                                       calculateTreeSizeEstimation());

          if (!curNodePtr->treat(globalNodesTreatOrder, primalIncBound()))
            {
              if (printL(-1))
                std::cout << "BaPCod WARNING: Branch-and-Price is interrupted" << std::endl;
              std::cerr << "BaPCod WARNING: Branch-and-Price is interrupted" << std::endl;
              break;
            }

          if (curNodePtr->depth() == 0)
            {
              bapcodInit().statistics().recTime("bcTimeRoot", bapcodInit().startTime().getElapsedTime_dbl());
              bapcodInit().statistics().recValue("bcRecRootInc", curNodePtr->nodeIncIpPrimalBound());
              bapcodInit().statistics().recValue("bcRecRootDb", curNodePtr->nodeIncIpDualBound());
              bapcodInit().statistics().recValue("bcRecRootLpVal", curNodePtr->nodeIncLpPrimalBound());
            }

          /// the output of the treated node are the generated child nodes and possibly the updated bounds
          /// and the updated solution, we should update primal bound before dual one
          /// as the dual bound will be limited by the primal one
          if (curNodePtr->primalBoundIsUpdated())
            updatePrimalIncSolution(curNodePtr->nodeIncIpPrimalBound(), curNodePtr->nodeIncIpPrimalSolPtr());
          if (curNodePtr->dualBoundIsUpdated())
            updateCurValidDualBound(curNodePtr);
          for (std::list<Node *>::const_iterator childPtrIt = curNodePtr->sons().begin();
               childPtrIt != curNodePtr->sons().end(); ++childPtrIt)
            {
              _bapTreeNodes.push_back(*childPtrIt);
              if ((*childPtrIt)->dualBoundIsUpdated())
                updateCurValidDualBound(*childPtrIt);
              if (searchTree.size() < param().OpenNodesLimit)
                searchTree.push(*childPtrIt);
              else
                secondarySearchTree.push(*childPtrIt);
            }
        }
      if (!curNodeSolvedBefore)
          addNodeToBaPTreeDotFile(curNodePtr);

      for (std::list<Node*>::iterator childPtrIt = curNodePtr->sons().begin();
           childPtrIt != curNodePtr->sons().end(); childPtrIt++)
        {
          if ((*childPtrIt)->solved())
            {
              (*childPtrIt)->branchAndPriceOrder(bapTreatOrder++);
              addNodeToBaPTreeDotFile(*childPtrIt);
            }
        }
      /// if the node $n$ has been conquered, we update the subtree size for every node $n'$
      /// such that $n$ was the last conquered in the subtree rooted at $n'$
      if (curNodePtr->sons().empty())
        curNodePtr->calculateSubtreeSize(_subTreeSizeByDepth);
    }

  if (printL(-1))
  {
    if (!searchTree.empty() && (thisSearchTreeTreatedNodesNumber >= param().MaxNbOfBBtreeNodeTreated()))
      std::cout << "SEACH IS INTERRUPTED as the limit on the number of branch-and-bound nodes is reached!" << std::endl;
    std::cout << "************************************************************************************************"
              << std::endl;
    std::cout << "Search is finished, global bounds : [ " << std::setprecision(6)
              << dualIncBound() << " , " << primalIncBound() << " ], ";
    printTime(bapcodInit().startTime().getElapsedTime(), std::cout);
    std::cout << "************************************************************************************************"
              << std::endl;
  }

  bapcodInit().statistics().incrCounter("bcCountNodeTreat", globalNodesTreatOrder);
  bapcodInit().statistics().setCounter("bcEstimTreeSize", calculateTreeSizeEstimation());

  bapcodInit().statistics().incrTimer("bcTimeBaP", bcTimeBaP.getElapsedTime_dbl());

  Solution * temp = _primalSolPtr;
  if (_primalSolPtr != NULL)
    {
      _primalSolPtr = getDissagregatedSolution(_primalSolPtr);
      delete temp;
    }

  bapcodInit().statistics().incrTimer("bcTimeMain", bapcodInit().startTime().getElapsedTime_dbl());

  if (printL(2) && (_primalSolPtr != NULL))
    _primalSolPtr->print();

  return _primalSolPtr;
}

int MasterConf::calculateTreeSizeEstimation()
{
    int nbSolvedChildNodes = 0;
    double totalDbImprovement = 0.0;
    for (auto * nodePtr : _bapTreeNodes)
    {
        if (nodePtr->solved() && !nodePtr->isRoot())
        {
            nbSolvedChildNodes += 1;
            if (nodePtr->nodeIncLpDualBound() < _cutOffValue)
                totalDbImprovement += nodePtr->nodeIncLpDualBound() - nodePtr->father()->nodeIncLpDualBound();
            else
                totalDbImprovement += _cutOffValue - nodePtr->father()->nodeIncLpDualBound();
        }
    }

    if (nbSolvedChildNodes > 0)
    {
        double remainingTreeSizeEstimation = 0.0;
        double averDbImprovement = totalDbImprovement / nbSolvedChildNodes;
//        std::cout << "Aver. DB improvement : " << averDbImprovement;
        for (auto * nodePtr: _bapTreeNodes)
            if (nodePtr->sons().empty())
            {
                double absoluteGap = _cutOffValue - nodePtr->nodeIncLpDualBound();
                if (absoluteGap > Double::precision)
                {
//                    std::cout << " " << absoluteGap;
                    remainingTreeSizeEstimation += std::pow(2, absoluteGap / averDbImprovement);
                }
            }
        if (remainingTreeSizeEstimation > INT_MAX)
            return 0;
        auto currTreeSize = bapcodInit().statistics().getCounter("bcCountNodeProc");
        return (int)currTreeSize + (int)remainingTreeSizeEstimation;
    }
    return 0;
}


Solution * MasterConf::getDebugSolution()
{
    return _debugSolutionPtr;
}

Solution * MasterConf::enumerateAllColumns(int & nbEnumColumns)
{
    if (_enumSolutionPtr != NULL)
    {
        _enumSolutionPtr->deleteSolutionsChain();
        delete _enumSolutionPtr;
    }
    _enumSolutionPtr = new Solution();


  nbEnumColumns = 0;
  for (std::vector<ColGenSpConf *>::iterator cgSpConfPtrIt = _colGenSubProbConfPts.begin();
       cgSpConfPtrIt != _colGenSubProbConfPts.end(); ++cgSpConfPtrIt)
  {
    int nbEnumColumnsThisCgSpConf = (*cgSpConfPtrIt)->probPtr()->enumerateAllColumns(_enumSolutionPtr);
    if (nbEnumColumnsThisCgSpConf < 0)
    {
      nbEnumColumns = -1;
      break;
    }
    else
    {
      nbEnumColumns += nbEnumColumnsThisCgSpConf;
    }
  }

  return _enumSolutionPtr;
}


void MasterConf::initializeBaPTreeDotFile() const
{
    std::string fileName = param().baPTreeDot_file().c_str();
    if (fileName == "")
        return;

    std::ofstream ofile(param().baPTreeDot_file().c_str(), std::ios::out);
    ofile << "##Command to get a nice layout: dot -Tpdf thisfile > thisfile.pdf" << std::endl;
    ofile << std::endl;
    ofile << "digraph " << name() << "_BaP_Tree {" << std::endl;
    ofile << "edge[fontname = \"Courier\", fontsize = 10];" << std::endl;
    ofile << "}";
    ofile.close();
}

void MasterConf::addNodeToBaPTreeDotFile(Node * nodePtr) const
{
  bapcodInit().statistics().incrCounter("bcCountNodeProc");

  std::string fileName = param().baPTreeDot_file().c_str();
  if (fileName == "")
    {
      nodePtr->clearLocalNodeBrConstrList(false);
      return;
    }
  std::ifstream infile(param().baPTreeDot_file().c_str());
  stringstream os;
  std::string line;
  while (std::getline(infile, line))
    {
      if (line != "}")
        os << line << std::endl;
    }
  infile.close();

  os << "n" << nodePtr->branchAndPriceOrder() << " [label="
     << "\"N_" << nodePtr->branchAndPriceOrder() << " ("
     << nodePtr->evalEndTimeString() << ")";
  if (nodePtr->debugSolutionAtThisNode())
      os << "(D)";
  os << " \\n";
  if (nodePtr->infeasible())
    os << "INFEAS\"];" << std::endl;
  else if (nodePtr->isConquered() || nodePtr->prunedAtBeginningOfTreatNode())
    os << "BOUND [" << fixed << setprecision(2)
       << nodePtr->nodeIncIpPrimalBound() <<"]\"];" << std::endl;
  else
    {
      os.precision(2);
      os << "[" << fixed << setprecision(2)
         << nodePtr->nodeIncLpDualBound() << ","
         << nodePtr->nodeIncIpPrimalBound() << "]\"];" << std::endl;
    }

  if (nodePtr->depth() > 0)
    {
      os << "n" << nodePtr->fatherBranchAndPriceOrder()
         << " -> n" << nodePtr->branchAndPriceOrder();
      os << " [ label = \"";
      std::vector<std::string>::size_type maxSize = 0;
      std::list<BranchingConstrBaseType*>::const_iterator constrPtrIt;
      for (constrPtrIt = nodePtr->localNodeBrConstrList().begin();
           constrPtrIt != nodePtr->localNodeBrConstrList().end(); constrPtrIt++)
        maxSize = (std::max)((*constrPtrIt)->forDotPrint().size(), maxSize);

      for (constrPtrIt = nodePtr->localNodeBrConstrList().begin();
           constrPtrIt != nodePtr->localNodeBrConstrList().end(); constrPtrIt++)
        {
          std::vector<std::string> constrLines = (*constrPtrIt)->forDotPrint();
          for (int i = 0; i < maxSize; i++)
            {
              if (constrPtrIt != nodePtr->localNodeBrConstrList().begin())
                {
                  if (i == 0)
                    os << " & ";
                  else
                    os << "   ";
                }
              else
                {
                  if ((i == 0) && (nodePtr->localNodeBrConstrList().size() > 1))
                    os << "[";
                  else
                    os << " ";
                }
              if (constrLines.size() > i)
                os << constrLines[i];
              else
                os << string(constrLines[0].size(), ' ');

              if (constrPtrIt != nodePtr->localNodeBrConstrList().end())
                {
                  if ((i == 0) && (nodePtr->localNodeBrConstrList().size() > 1))
                    os << "]";
                  else
                    os << " ";
                  if (i != maxSize - 1)
                    os << "\\n";
                }
            }
        }
      os << "\" ];";
      os << std::endl;
    }

  os << "}";

  std::ofstream outfile(param().baPTreeDot_file().c_str());
  outfile << os.str();
  outfile.close();

  /// we now can clear local branching constraints of the node
  nodePtr->clearLocalNodeBrConstrList(false);
}
