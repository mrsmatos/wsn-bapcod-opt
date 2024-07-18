/**
 *
 * This file bcNodeC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcNodeC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMasterConfC.hpp"
#include "bcMastColumnC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcProblemC.hpp"
#include "bcProbConfigC.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcNetworkBasedCuts.hpp"
#include "bcAlg4EvalOfNode.hpp"
#include "bcAlg4GenChildrenOfNode.hpp"
#include "bcAlg4EvalByColAndCutGen.hpp"
#include "bcRestrictedMasterIpHeuristicC.hpp"
#include "bcColGenSpConfC.hpp"

AutoRankOneCutsMemoryInfo::AutoRankOneCutsMemoryInfo(BcObjStatus::MinMaxIntFloat objStatus_) :
    objStatus(objStatus_), probSetupInfoPtr(NULL), nodeEvalInfoPtr(NULL), evalAlgPtr(NULL), problemSetupAlgPtr(NULL),
    problemSetDownAlgPtr(NULL), savedLpDualBound(Bound::infDualBound(objStatus)),
    startTime(0.0), elapsedTime(0.0), quitByTailingOff(false)
{
}

AutoRankOneCutsMemoryInfo::~AutoRankOneCutsMemoryInfo()
{
  deleteInfoAndBounds();
  if (evalAlgPtr != NULL)
    delete evalAlgPtr;
  evalAlgPtr = NULL;
  if (problemSetupAlgPtr != NULL)
    delete problemSetupAlgPtr;
  problemSetupAlgPtr = NULL;
  if (problemSetDownAlgPtr != NULL)
    delete problemSetDownAlgPtr;
  problemSetDownAlgPtr = NULL;
}

void AutoRankOneCutsMemoryInfo::deleteInfoAndBounds()
{
  if (probSetupInfoPtr != NULL)
  {
    probSetupInfoPtr->numberOfNodes -= 1;
    if (probSetupInfoPtr->numberOfNodes == 0)
      delete probSetupInfoPtr;
  }
  probSetupInfoPtr = NULL;
  if (nodeEvalInfoPtr != NULL)
  {
    nodeEvalInfoPtr->numberOfNodes -= 1;
    if (nodeEvalInfoPtr->numberOfNodes == 0)
      delete nodeEvalInfoPtr;
  }
  nodeEvalInfoPtr = NULL;
  savedLpDualBound = Bound::infDualBound(objStatus);
}


Node::Node(ProbConfig * probConfigPtr,
           const Bound & dualBound,
           ProblemSetupInfo * problemSetupInfoPtr,
           NodeEvalInfo * nodeSolverInfoPtr,
           bool debugSolutionIsDefined) :
  _probConfigPtr(probConfigPtr), _objStatus(probConfigPtr->probPtr()->objStatus()),
  _ref(probConfigPtr->pcNodeCount()), _param(probConfigPtr->param()), _father(NULL), _sons(), _depth(0),
  _prunedAtBeginningOfTreatNode(false), _debugSolutionAtThisNode(debugSolutionIsDefined),
  _estimatedSubtreeSize(BapcodInfinity), _branchingPriorityLevel(BapcodInfinity), _subtreeSize(-1),
  _nodeIncLpDualBound(dualBound), _nodeIncIpDualBound(dualBound),
  _nodeIncLpPrimalBound(probConfigPtr->primalIncBound()), _nodeIncIpPrimalBound(probConfigPtr->primalIncBound()),
  _subtreeDualBound(dualBound), _dualBoundIsUpdated(false), _ipPrimalBoundIsUpdated(false),
  _nodeIncIpPrimalSolPtr(NULL), _localNodeBrConstrList(), _localFixedSolution(NULL),
  _evalEndTime(-1), _treatOrder(-1), _branchAndPriceOrder(-1), _fatherBranchAndPriceOrder(-1),
  _infeasible(false), _evaluated(false), _treated(false),
  _problemAndEvalAlgInfoSaved(false), _primalSol(), _probSetupInfoPtr(NULL), _nodeEvalInfoPtr(NULL),
  _genChildNodesInfoPtr(NULL), _branchEvalInfoPtr(NULL), _autoRankOneCutsMemoryInfoPtr(NULL),
  _applyAutoRankOneCutsMemorySearch(false), _strongBranchPhaseNumber(0), _strongBranchNodeNumber(-1),
  _preprocessorPtr(NULL), _evalAlgPtr(NULL), _problemSetupAlgPtr(NULL), _problemSetDownAlgPtr(NULL),
  _genChildNodesAlgPtr(NULL), _PrimalHeuristicPts(), _treeOfColClasses(), _cgSpConfTreeOfColClassesMap()
{
  probConfigPtr->increasePCNodeCount();

  _nodeIncIpDualBound.round();

  associateProblemSetupInfoPtr(problemSetupInfoPtr);
  associateNodeEvalInfoPtr(nodeSolverInfoPtr);
  associateGenChildNodesInfoPtr(new GenChildNodesInfo());

  for (std::list<BranchingConstrBaseType *>::iterator it = _localNodeBrConstrList.begin();
          it != _localNodeBrConstrList.end(); it++)
    {
      Constraint * brConstrPtr = dynamic_cast<Constraint *> (*it);
      if (brConstrPtr != NULL)
        brConstrPtr->incrParticipation(0);
      if (printL(7))
        std::cout << "Node::node() participation of brConstr " << brConstrPtr->name()
                  << " at " << static_cast<void*> (brConstrPtr)
                  << " was incremented to " << brConstrPtr->participation() << std::endl;
    }
}

Node::Node(int ref,
           Node * father,
           const std::list<BranchingConstrBaseType *> & localNodeBrConstrList,
           Solution * localFixedSolution,
           bool inheritDualBound):
  _probConfigPtr(father->_probConfigPtr), _objStatus(father->_objStatus), _ref(ref), _param(_probConfigPtr->param()),
  _father(father), _sons(), _depth(father->depth() + 1), _prunedAtBeginningOfTreatNode(false),
  _debugSolutionAtThisNode(father->debugSolutionAtThisNode()),
  _estimatedSubtreeSize(BapcodInfinity), _branchingPriorityLevel(BapcodInfinity), _subtreeSize(-1),
  _nodeIncLpDualBound(Bound::infDualBound(_objStatus)), _nodeIncIpDualBound(Bound::infDualBound(_objStatus)),
  _nodeIncLpPrimalBound(Bound::infPrimalBound(_objStatus)), _nodeIncIpPrimalBound(father->nodeIncIpPrimalBound()),
  _subtreeDualBound(father->nodeIncIpDualBound()), _dualBoundIsUpdated(false), _ipPrimalBoundIsUpdated(false),
  _nodeIncIpPrimalSolPtr(NULL), _localNodeBrConstrList(localNodeBrConstrList),
  _localFixedSolution(localFixedSolution), _evalEndTime(-1), _treatOrder(-1), _branchAndPriceOrder(-1),
  _fatherBranchAndPriceOrder(father->branchAndPriceOrder()),
  _infeasible(false), _evaluated(false), _treated(false), _problemAndEvalAlgInfoSaved(false), _primalSol(),
  _probSetupInfoPtr(NULL), _nodeEvalInfoPtr(NULL), _genChildNodesInfoPtr(NULL), _branchEvalInfoPtr(NULL),
  _autoRankOneCutsMemoryInfoPtr(NULL), _applyAutoRankOneCutsMemorySearch(false), _strongBranchPhaseNumber(0),
  _strongBranchNodeNumber(-1), _preprocessorPtr(NULL), _evalAlgPtr(NULL), _problemSetupAlgPtr(NULL),
  _problemSetDownAlgPtr(NULL), _genChildNodesAlgPtr(NULL), _PrimalHeuristicPts(), _treeOfColClasses(),
  _cgSpConfTreeOfColClassesMap()
{
  if (inheritDualBound)
    {
      _nodeIncIpDualBound = father->nodeIncIpDualBound();
      _nodeIncLpDualBound = father->nodeIncLpDualBound();
    }

  associateProblemSetupInfoPtr(father->_probSetupInfoPtr);
  associateNodeEvalInfoPtr(father->_nodeEvalInfoPtr);

  for (std::list<BranchingConstrBaseType *>::iterator it = _localNodeBrConstrList.begin();
       it != _localNodeBrConstrList.end(); it++)
    {
      Constraint * brConstrPtr = dynamic_cast<Constraint *> (*it);
      if (brConstrPtr != NULL)
        brConstrPtr->incrParticipation(0);
      if (printL(7))
        {
          std::cout << "Node::node() participation of brConstr " << brConstrPtr->name()
                    << " at " << static_cast<void*> (brConstrPtr)
                    << " was incremented to " << brConstrPtr->participation() << std::endl;
        }
    }
}

Problem * Node::probPtr() const
{
  return _probConfigPtr->probPtr();
}

Node::~Node()
{
  /// Does the memory clean-up normally done after node treatment
  /// (enumeration can be interrupted before all nodes have been treated)
  if (!_treated)
    exitTreatment();

  clearLocalNodeBrConstrList(true);

  if (_nodeIncIpPrimalSolPtr != NULL)
    delete _nodeIncIpPrimalSolPtr;
  _nodeIncIpPrimalSolPtr = NULL;

  for (std::set<Alg4PrimalHeuristicOfNode *, PrimalHeurSort>::iterator heurPtrIt = _PrimalHeuristicPts.begin();
       heurPtrIt != _PrimalHeuristicPts.end(); heurPtrIt++)
    delete (*heurPtrIt);
  _PrimalHeuristicPts.clear();
}

void Node::clearLocalNodeBrConstrList(const bool & calledFromDestructor)
{
  std::list<BranchingConstrBaseType *>::iterator brConstrPtrIt = _localNodeBrConstrList.begin();
  while (brConstrPtrIt != _localNodeBrConstrList.end())
    {
      /// we delete component set branching constraints only in the destructor
      /// as they are needed for component set branching of node's descendants
      if (!(*brConstrPtrIt)->isTypeOf(VcId::CompSetInstMastBranchConstrMask) || calledFromDestructor)
        {
          Constraint * brConstrPtr = dynamic_cast<Constraint*>(*brConstrPtrIt);
          if (brConstrPtr != NULL)
            brConstrPtr->decrParticipation(0);
          if (printL(7))
            std::cout << "Node::clearLocalNodeBrConstrList() participation of brConstr " << brConstrPtr->name()
                      <<  " at " << static_cast<void*>(brConstrPtr)
                      << " was decremented to " << brConstrPtr->participation() << std::endl;
          if ((brConstrPtr != NULL) && (brConstrPtr->problemPtr() == NULL) && (brConstrPtr->participation() == 0))
          {
            /// when constraint does not belong to a problem, it should be deleted, otherwise there is a memory leak
            delete brConstrPtr;
          }
          brConstrPtrIt = _localNodeBrConstrList.erase(brConstrPtrIt);
        }
      else
        ++brConstrPtrIt;
    }
}

void Node::setProblemSetupAlgorithm(Alg4ProblemSetupOfNode * problemSetupAlgPtr)
{
  _problemSetupAlgPtr = problemSetupAlgPtr;

}
void Node::setProblemSetDownAlgorithm(Alg4ProblemSetDownOfNode * problemSetDownAlgPtr)
{
  _problemSetDownAlgPtr = problemSetDownAlgPtr;
}

void Node::setPreprocessor(Alg4PreprocessingOfNode * preprocessorPtr)
{
  _preprocessorPtr = preprocessorPtr;
}

void Node::setEvalAlg(Alg4EvalOfNode * nodeSolverPtr)
{
  _evalAlgPtr = nodeSolverPtr;
}

void Node::setGenChildNodesAlgorithm(Alg4GenChildrenOfNode * genChildNodesAlgPtr)
{
  _genChildNodesAlgPtr = genChildNodesAlgPtr;
}

void Node::addPrimalHeuristic(Alg4PrimalHeuristicOfNode * PrimalHeuristicPtr)
{
  _PrimalHeuristicPts.insert(PrimalHeuristicPtr);
}

bool Node::needProblemFullSetDownAlgorithm()
{
  return !_PrimalHeuristicPts.empty() || (_genChildNodesAlgPtr != NULL);
}

void Node::calculateSubtreeSize(std::vector<std::pair<int, double> > & subTreeSizeByDepth)
{
  bool allSonsSubreesWereCalculated = true;
  int totalSonsSubreesSize = 0;
  for (std::list<Node *>::const_iterator childPtrIt = _sons.begin();
       childPtrIt != _sons.end(); ++childPtrIt)
    if ((*childPtrIt)->subtreeSize() > 0)
      totalSonsSubreesSize += (*childPtrIt)->subtreeSize();
    else
      allSonsSubreesWereCalculated = false;
  if (allSonsSubreesWereCalculated)
    {
      _subtreeSize = totalSonsSubreesSize + 1;

      if (subTreeSizeByDepth.size() <= _depth)
        subTreeSizeByDepth.resize(_depth + 1, std::make_pair(0, 0.0));
      int & curNumEvals = subTreeSizeByDepth[_depth].first;
      double & curAverTreeSize = subTreeSizeByDepth[_depth].second;
      curAverTreeSize = (curAverTreeSize * curNumEvals + _subtreeSize ) / (curNumEvals + 1);
      curNumEvals += 1;

      if (!isRoot())
        _father->calculateSubtreeSize(subTreeSizeByDepth);
    }
}

Node * Node::father() const
{
  return _father;
}

void Node::clearLocalFixedSolution()
{
  if (_localFixedSolution == NULL)
    _localFixedSolution = new Solution();
  else
    _localFixedSolution->clear();
}

void Node::printFixedSolution(std::ostream& os, const bool printMasterSolution) const
{
  if (_localFixedSolution == NULL)
    return;

  int solSize = _localFixedSolution->solVarValMap().size();
  if (solSize == 0)
    {
      os << "Fixed master solution is empty" << std::endl;
      return;
    }

  VarPtr2DoubleMap projectedSolution;
  VarPtr2DoubleMap::const_iterator mapIt;
  for (mapIt = _localFixedSolution->solVarValMap().begin();
      mapIt != _localFixedSolution->solVarValMap().end(); mapIt++)
    {
      if (mapIt->first->isTypeOf(VcId::InstMasterVarMask)) /// pure master variable
        projectedSolution[mapIt->first] = mapIt->second;
      if (mapIt->first->isTypeOf(VcId::MastColumnMask)) /// master column
        mapIt->first->fillAggregateSol(projectedSolution, mapIt->second);
    }
  if (printMasterSolution)
  {
    if (solSize <= 20)
    {
      for (mapIt = _localFixedSolution->solVarValMap().begin();
           mapIt != _localFixedSolution->solVarValMap().end(); mapIt++)
      {
        if (mapIt == _localFixedSolution->solVarValMap().begin())
          os << "Fixed master solution = [";
        else
          os << ", ";
        os << mapIt->first->name() << "=" << mapIt->second;
        if (mapIt->first->isTypeOf(VcId::MastColumnMask))
          os << " from " << mapIt->first->spSol()->probConfPtr()->name();
        else
          os << " (pure master var) ";
      }
      os << "]" << std::endl;
    }
    else
    {
      os << "Fixed " << solSize << " master variables " << std::endl;
    }
  }
  if (projectedSolution.size() <= 100)
    {
      for (mapIt = projectedSolution.begin(); mapIt != projectedSolution.end(); mapIt++)
        {
          if (mapIt == projectedSolution.begin())
            os << "Fixed projected solution = [";

          if (mapIt->first->genVarPtr()->defaultName() == "TLCCV")
            continue;

          if (mapIt != projectedSolution.begin())
            os << ", ";
          os << mapIt->first->name() << "=" << mapIt->second;
        }
      os << "]" << std::endl;
    }
  else
    os << "Fixed " << solSize << " projected variables " << std::endl;
}

ConstrPtrList Node::upCastedBranchingConstrList()
{
  std::list<Constraint *> branchingConstrList;
  if (_depth == 0)
    return branchingConstrList;
  else
    {
      for (std::list<BranchingConstrBaseType *>::iterator it = _localNodeBrConstrList.begin();
          it != _localNodeBrConstrList.end(); ++it)
        {
          if ((*it)->isTypeOf(VcId::ConstraintMask))
            {
              Constraint * constrPtr = dynamic_cast<Constraint *>(*it);
              if (constrPtr != NULL)
                  branchingConstrList.push_back(constrPtr);
            }
        }
      return branchingConstrList;
    }
}

bool Node::runEnumeratedMIP()
{
  MasterConf * mastConfPtr = dynamic_cast<MasterConf *>(_probConfigPtr);
  if (mastConfPtr == NULL)
    {
      std::cerr << "Node::runEnumeratedMIP() error: node does not belong to the master" << std::endl;
      exit(1);
    }

  Time startTime;

  if (printL(-1))
    std::cout << "----- Terminating the node by MIP -----" << std::endl;

  /// enumerated MIP could increase the Lp bound, therefore we should reset it now
  _nodeIncLpPrimalBound = Bound::infPrimalBound(_objStatus);

  MasterCommons4EvalAlg masterCommons(*mastConfPtr);
  Alg4EvalByMip alg4EvalByMip(_probConfigPtr->probPtr(), masterCommons);
  alg4EvalByMip.setParamMaxTime(param().MipSolverMaxTime());
  alg4EvalByMip.setParamExactSolution(true);
  alg4EvalByMip.setParamRepeatIfCoreCutAdded(true);

  bool returnValue = alg4EvalByMip.setupAlgo(this);
  if (!returnValue)
    {
      alg4EvalByMip.eval();

      if (alg4EvalByMip.algIncIpPrimalBoundUpdated())
        recordIpPrimalSolAndUpdateIpPrimalBound(alg4EvalByMip.algIncIpPrimalBound(),
                                                alg4EvalByMip.algIncIpPrimalSolMap());
      nodeIncLpPrimalBound(alg4EvalByMip.algIncLpPrimalBound());
      nodeIncLpDualBound(alg4EvalByMip.algIncLpDualBound());
      nodeIncIpDualBound(alg4EvalByMip.algIncIpDualBound());

      alg4EvalByMip.setDownAlgo();
    }

  bapcodInit().statistics().incrTimer("bcTimeEnumMPsol", startTime.getElapsedTime_dbl());
  return returnValue;
}

/// assumes that _branchEvalInfoPtr is not NULL
void Node::storeBranchingEvaluationInfo(double time)
{
  double lpValue = ((nodeIncLpPrimalBound() > nodeIncIpPrimalBound()) ? nodeIncIpPrimalBound()
                                                                      : nodeIncLpPrimalBound());
  _branchEvalInfoPtr->addEvalEntry(_strongBranchPhaseNumber, _strongBranchNodeNumber, lpValue, time);
}

void Node::setDualBoundEqualToIncPrimalBound()
{
  updateDualBounds(_nodeIncIpPrimalBound, _nodeIncIpPrimalBound);
}

#ifdef BCP_RCSP_IS_FOUND
bool Node::autoRankOneCutMemoryEvaluation(int & globalTreatOrder, const Bound & incPrimalBound)
{
  Alg4EvalByColAndCutGen * cutGenEvalAlg = dynamic_cast<Alg4EvalByColAndCutGen *>(_evalAlgPtr);
  if (cutGenEvalAlg == NULL)
    return evaluation(globalTreatOrder, incPrimalBound);

  GenericLimMemRankOneCutConstr * rankOneGenCutConstrPtr = dynamic_cast<GenericLimMemRankOneCutConstr *>
                                                           (probConfPtr()->getGenericCutConstr("R1C"));
  if (rankOneGenCutConstrPtr == NULL)
    return evaluation(globalTreatOrder, incPrimalBound);

  cutGenEvalAlg->setOptionSaveInfoForAutoRankOneCutsMemory(true);
  _autoRankOneCutsMemoryInfoPtr->startTime = bapcodInit().startTime().getElapsedTime_dbl();

  if (printL(-1))
    std::cout << "--------------------------------------------------------------------------" << std::endl
              << "BaPCod info : first we run the root node with rank-1 cuts with node memory" << std::endl
              << "--------------------------------------------------------------------------" << std::endl;


  rankOneGenCutConstrPtr->setVertexMemory();
  bool returnValue = evaluation(globalTreatOrder, incPrimalBound);

  /// if the node is not conquered and we did not quit by tailing off, we repeat the root with arc memory for cuts
  if (returnValue && !_treated && (_nodeEvalInfoPtr != NULL) &&  !_autoRankOneCutsMemoryInfoPtr->quitByTailingOff)
  {
    MasterConf * masterConfPtr = static_cast<MasterConf *>(probConfPtr());

    Node * arcMemNode = new Node(masterConfPtr, _autoRankOneCutsMemoryInfoPtr->savedLpDualBound,
                                 _autoRankOneCutsMemoryInfoPtr->probSetupInfoPtr,
                                 _autoRankOneCutsMemoryInfoPtr->nodeEvalInfoPtr,
                                 masterConfPtr->getDebugSolution() != NULL);
    arcMemNode->setProblemSetupAlgorithm(_autoRankOneCutsMemoryInfoPtr->problemSetupAlgPtr);
    arcMemNode->setProblemSetDownAlgorithm(_autoRankOneCutsMemoryInfoPtr->problemSetDownAlgPtr);
    arcMemNode->setEvalAlg(_autoRankOneCutsMemoryInfoPtr->evalAlgPtr);

    if (printL(-1))
    {
      std::cout << "-----------------------------------------------------------------------" << std::endl
                << "BaPCod info : now we run the root node with rank-1 cuts with arc memory" << std::endl
                << "-----------------------------------------------------------------------" << std::endl;
    }

    double startArcTime = bapcodInit().startTime().getElapsedTime_dbl();

    rankOneGenCutConstrPtr->setArcMemory();
    returnValue = arcMemNode->evaluation(globalTreatOrder, _nodeIncIpPrimalBound);

    if (arcMemNode->primalBoundIsUpdated())
      recordIpPrimalSolAndUpdateIpPrimalBound(arcMemNode->nodeIncIpPrimalBound(), arcMemNode->nodeIncIpPrimalSolPtr());


    if (!returnValue || arcMemNode->treated())
    {
      /// arcMemNode is conquered
      deleteProblemAndEvalAlgInfo(); /// also deletes _primalSol;
      _nodeIncLpDualBound = arcMemNode->nodeIncLpDualBound();
      _nodeIncLpPrimalBound = arcMemNode->nodeIncLpPrimalBound();
      _nodeIncIpDualBound = arcMemNode->nodeIncIpDualBound();
    }
    else
    {
      double arcElapsedTime = bapcodInit().startTime().getElapsedTime_dbl();

      double gapRatio = (_nodeIncIpPrimalBound - arcMemNode->nodeIncLpDualBound())
                        / (_nodeIncIpPrimalBound - _nodeIncLpDualBound);
      double timeRatio = arcElapsedTime / _autoRankOneCutsMemoryInfoPtr->elapsedTime;
      double arcMemoryPrefScore = pow(2, 10) / pow(2, 10 * gapRatio) / timeRatio;
      bool arcMemoryIsPreferred = (arcMemNode->nodeIncLpDualBound() > _nodeIncLpDualBound) && (arcMemoryPrefScore > 1);

      if (printL(0))
        std::cout << "-----------------------------------------------------------------------" << std::endl
                  << "Run with node memory : bound is " << _nodeIncLpDualBound
                  << ", elapsed time is " <<  _autoRankOneCutsMemoryInfoPtr->elapsedTime / 100.0 << std::endl
                  << "Run with arc memory : bound is " << arcMemNode->nodeIncLpDualBound()
                  << ", elapsed time is " << arcElapsedTime / 100.0 << std::endl
                  << (arcMemoryIsPreferred ? "Arc" : "Node") << " memory is preferred (score = "
                  << arcMemoryPrefScore << ")" << std::endl
                  << "-----------------------------------------------------------------------" << std::endl;
      else if (printL(-1))
        std::cout << "-----------------------------------------------------------------------" << std::endl
                  << (arcMemoryIsPreferred ? "Arc" : "Node") << " memory is preferred (score = "
                  << arcMemoryPrefScore << ")" << std::endl
                  << "-----------------------------------------------------------------------" << std::endl;

      if (arcMemoryIsPreferred)
      {
        deleteProblemAndEvalAlgInfo(); /// also deletes _primalSol;
        _nodeIncLpDualBound = arcMemNode->nodeIncLpDualBound();
        _nodeIncLpPrimalBound = arcMemNode->nodeIncLpPrimalBound();
        _nodeIncIpDualBound = arcMemNode->nodeIncIpDualBound();
        associateNodeEvalInfoPtr(arcMemNode->_nodeEvalInfoPtr);
        associateProblemSetupInfoPtr(arcMemNode->_probSetupInfoPtr);
        /// we copy the solution
        for (SolutionVarInfoPtrList::const_iterator solIt = arcMemNode->primalSol().begin();
             solIt != arcMemNode->primalSol().end(); ++solIt)
          _primalSol.push_back(new SolutionVarInfo(**solIt));
      }
      else
      {
        rankOneGenCutConstrPtr->setVertexMemory();
      }
    }

    /// algorithms will be deleted together with the node
    _autoRankOneCutsMemoryInfoPtr->problemSetupAlgPtr = NULL;
    _autoRankOneCutsMemoryInfoPtr->problemSetDownAlgPtr = NULL;
    _autoRankOneCutsMemoryInfoPtr->evalAlgPtr = NULL;
    delete arcMemNode;
  }

  if (_autoRankOneCutsMemoryInfoPtr != NULL)
    delete _autoRankOneCutsMemoryInfoPtr;
  _autoRankOneCutsMemoryInfoPtr = NULL;
  return returnValue;
}
#endif //BCP_RCSP_IS_FOUND

bool Node::evaluation(int & globalTreatOrder, const Bound & incPrimalBound)
{
  _treatOrder = ++globalTreatOrder;
  _nodeIncIpPrimalBound = incPrimalBound;
  _ipPrimalBoundIsUpdated = false;
  _dualBoundIsUpdated = false;

  bapcodInit().require(_problemSetupAlgPtr != NULL, "Problem setup algorithm should be defined for a node!");
  bapcodInit().require(_problemSetDownAlgPtr != NULL, "Problem set down algorithm should be defined for a node!");
  bapcodInit().require(_evalAlgPtr != NULL, "Node solver should be defined for a node!");

  Time start;

  /// Problem setup
  if (_problemSetupAlgPtr->run(this))
    {
      _problemSetDownAlgPtr->run();
      return markInfeasibleAndExitTreatment();
    }
  removeProblemSetupInfoAssociation();

  if ((_preprocessorPtr != NULL) && _preprocessorPtr->run(this))
    {
      _problemSetDownAlgPtr->run();
      return markInfeasibleAndExitTreatment();
    }

  if (_evalAlgPtr->setupAlgo(this))
    {
      _evalAlgPtr->setDownAlgo();
      _problemSetDownAlgPtr->run();
      return markInfeasibleAndExitTreatment();
    }
  removeNodeEvalInfoAssociation();

  bapcodInit().statistics().incrTimer("bcTimeSetMast", start.getElapsedTime_dbl());

  if (_evalAlgPtr->eval())
    {
      _evalAlgPtr->setDownAlgo();
      _problemSetDownAlgPtr->run();
      return markInfeasibleAndExitTreatment();
    }
  _evaluated = true;

  //Issam : this should be also called after the heuristics.
  if (_evalAlgPtr->algIncIpPrimalBoundUpdated())
    recordIpPrimalSolAndUpdateIpPrimalBound(_evalAlgPtr->algIncIpPrimalBound(), _evalAlgPtr->algIncIpPrimalSolMap());
  nodeIncLpPrimalBound(_evalAlgPtr->algIncLpPrimalBound());
  updateDualBounds(_evalAlgPtr->algIncLpDualBound(), _evalAlgPtr->algIncIpDualBound());

  if (printL(1))
  {
      std::cout << "Alg IP bounds =  [" << std::setprecision(12) << _evalAlgPtr->algIncIpDualBound() << ","
                << _evalAlgPtr->algIncIpPrimalBound() << "], node IP bounds = ["
                << _nodeIncIpDualBound << "," << _nodeIncIpPrimalBound << "]"
                << std::setprecision(6) << std::endl;
  }


  if (!isConquered() && _evalAlgPtr->checkIfSubProbSolutionsEnumeratedToMIP() && runEnumeratedMIP())
    {
      _evalAlgPtr->setDownAlgo();
      _problemSetDownAlgPtr->run();

      if (_branchEvalInfoPtr != NULL)
        storeBranchingEvaluationInfo(start.getElapsedTime_dbl());
      return markInfeasibleAndExitTreatment();
    }

  if (isConquered())
    {
      _evalAlgPtr->setDownAlgo();
      _problemSetDownAlgPtr->run();

      if (_branchEvalInfoPtr != NULL)
        storeBranchingEvaluationInfo(start.getElapsedTime_dbl());

      return exitTreatment(true);
    }

  Time start2;

  if (!_problemAndEvalAlgInfoSaved)
    saveProblemAndEvalAlgInfo();

  _evalAlgPtr->setDownAlgo();
  _problemSetDownAlgPtr->run();

  bapcodInit().statistics().incrTimer("bcTimeSetMast", start2.getElapsedTime_dbl());

  _evalEndTime = bapcodInit().startTime().getElapsedTime();
  if (_branchEvalInfoPtr != NULL)
      storeBranchingEvaluationInfo(start.getElapsedTime_dbl());

  return true;
}

void Node::updateDualBounds(const Bound & lpDualBound, const Bound & ipDualBound)
{
  if (_nodeIncLpDualBound < lpDualBound)
    {
      _nodeIncLpDualBound = lpDualBound;
      _dualBoundIsUpdated = true;
    }
  if (_nodeIncIpDualBound < ipDualBound)
    {
      _nodeIncIpDualBound = ipDualBound;
      _dualBoundIsUpdated = true;
    }
}

void Node::recordIpPrimalSolAndUpdateIpPrimalBound(const Bound & primalBound, Solution * solPtr)
{
    if (solPtr == NULL)
        return;

    solPtr->resetCost();
    //Bound primalBound(solPtr->cost(), _objStatus);
    if (printL(5))
        std::cout << "Node::recordIpPrimalSolAndUpdateIpPrimalBound node " << _treatOrder
                  << "  primalBound " << primalBound << " _nodeIncIpPrimalBound " << _nodeIncIpPrimalBound << std::endl;

    if (_nodeIncIpPrimalBound > primalBound)
    {
        _nodeIncIpPrimalBound = primalBound;
        if (_nodeIncIpPrimalSolPtr != NULL)
            delete _nodeIncIpPrimalSolPtr;
        _nodeIncIpPrimalSolPtr = solPtr->clone();
        _ipPrimalBoundIsUpdated = true;
    }

}

void Node::recordIpPrimalSolAndUpdateIpPrimalBound(const Bound & primalBound, const VarPtr2DoubleMap & primalSolMap)
{
  Solution * intermSolPtr = new Solution(_probConfigPtr, primalSolMap);
  recordIpPrimalSolAndUpdateIpPrimalBound(primalBound, intermSolPtr);
  delete intermSolPtr;
}

/// the main function in which the node is treated
/// globalTreatOrder count the overall number of treated nodes
/// (additional nodes can be created and treated inside heuristics and/or branching)
/// returns false if a error occured
bool Node::treat(int & globalTreatOrder, const Bound & incPrimalBound)
{
  /// In strong branching, part I of treat (setup, preprocessing and solve) is separated
  /// from part II (heuristics and children generation).
  /// Therefore, treat() can be called two times, one inside strong branching, second inside the branch-and-price tree.
  /// Thus, variables _solved is used to know whether part I has already been done or not.
  if (!_evaluated)
    {
#ifdef BCP_RCSP_IS_FOUND
      if (_applyAutoRankOneCutsMemorySearch)
      {
        if (!autoRankOneCutMemoryEvaluation(globalTreatOrder, incPrimalBound))
          return false;
      }
      else
#endif //BCP_RCSP_IS_FOUND
      {
        if (!evaluation(globalTreatOrder, incPrimalBound))
          return false;
      }
    }
  else
    {
      if (incPrimalBound <= _nodeIncIpPrimalBound)
        {
          _nodeIncIpPrimalBound = incPrimalBound;
          _ipPrimalBoundIsUpdated = false;
        }
    }

  if (depth() == 0)
    bapcodInit().statistics().recTime("bcTimeRootEval", bapcodInit().startTime().getElapsedTime_dbl());

  if (_treated)
    return true;

  /// check time limit
  if (bapcodInit().startTime().getElapsedTime() > param().GlobalTimeLimitInTick())
    return exitTreatment();

  if (isConquered())
    return exitTreatment();

  for (std::set<Alg4PrimalHeuristicOfNode *, PrimalHeurSort>::iterator heurIt = _PrimalHeuristicPts.begin();
       heurIt != _PrimalHeuristicPts.end(); ++heurIt)
    {

        Time startTime;

        (*heurIt)->run(this, globalTreatOrder);

        bapcodInit().statistics().incrTimer("bcTimePrimalHeur", startTime.getElapsedTime_dbl());

        //TODO remove node bound updates from inside heuristics and put it here.

       if (isConquered())
         return exitTreatment();
    }

  /// the generation child nodes algorithm fills the sons
  if (_genChildNodesAlgPtr == NULL)
    return exitTreatment();

  if (_genChildNodesAlgPtr->setupAlgo(this))
    {
      _genChildNodesAlgPtr->setDownAlgo();
      return exitTreatment();
    }
  removeGenChildNodesInfoAssociation();

  _genChildNodesAlgPtr->run(globalTreatOrder);
  _genChildNodesAlgPtr->setDownAlgo();

  return exitTreatment();
}

bool Node::markInfeasibleAndExitTreatment()
{
  _infeasible = true;
  _subtreeDualBound = _nodeIncLpDualBound = _nodeIncIpDualBound = Bound::infPrimalBound(_objStatus);

  if (printL(1))
    std::cout << " Node:: EARLY TERMINATION of node treatment : infeasibility is detected" << std::endl;

  return exitTreatment();
}

bool Node::exitTreatment(const bool exitValue)
{
  if (_evalEndTime == -1)
    _evalEndTime = bapcodInit().startTime().getElapsedTime();

  if (_preprocessorPtr != NULL)
    delete _preprocessorPtr;
  _preprocessorPtr = NULL;

  if (_evalAlgPtr != NULL)
    delete _evalAlgPtr;
  _evalAlgPtr = NULL;

  if (_problemSetDownAlgPtr != NULL)
    delete _problemSetDownAlgPtr;
  _problemSetDownAlgPtr = NULL;

  if (_problemSetupAlgPtr != NULL)
    delete _problemSetupAlgPtr;
  _problemSetupAlgPtr = NULL;

  if (_genChildNodesAlgPtr != NULL)
    delete _genChildNodesAlgPtr;
  _genChildNodesAlgPtr = NULL;

  if (_autoRankOneCutsMemoryInfoPtr != NULL)
    delete _autoRankOneCutsMemoryInfoPtr;
  _autoRankOneCutsMemoryInfoPtr = NULL;

  removeNodeEvalInfoAssociation();
  removeProblemSetupInfoAssociation();
  removeGenChildNodesInfoAssociation();

  _treeOfColClasses.clear();
  _cgSpConfTreeOfColClassesMap.clear();

  while (!_primalSol.empty())
  {
     delete _primalSol.back();
     _primalSol.pop_back();
  }

  if (_localFixedSolution != NULL)
    delete _localFixedSolution;
  _localFixedSolution = NULL;

  /// cannot clear, as used in generateBaPTreeDotFile
  //_localNodeBrConstrList.clear();

  _evaluated = true;
  _treated = true;
  return exitValue;
}

void Node::deleteProblemAndEvalAlgInfo()
{
  if (_problemAndEvalAlgInfoSaved)
    {
      while (!_primalSol.empty())
        {
          delete _primalSol.back();
          _primalSol.pop_back();
        }
      removeNodeEvalInfoAssociation();
      removeProblemSetupInfoAssociation();
    }
  _problemAndEvalAlgInfoSaved = false;
}

void Node::saveProblemAndEvalAlgInfo()
{
  Time start;

  deleteProblemAndEvalAlgInfo();
  if (_evalAlgPtr != NULL)
    associateNodeEvalInfoPtr(_evalAlgPtr->recordNodeEvalInfo(_treatOrder));
  if (_problemSetDownAlgPtr)
    associateProblemSetupInfoPtr(_problemSetDownAlgPtr->recordProblemInfo(_treatOrder));
  _problemAndEvalAlgInfoSaved = true;

  bapcodInit().statistics().incrTimer("bcTimeSetMast", start.getElapsedTime_dbl());
}

void Node::saveAutoRankOneCutsMemoryInfo(const bool quitByTailingOff, const bool rank1CutsArePresent)
{
  if (_autoRankOneCutsMemoryInfoPtr == NULL)
    return;

  _autoRankOneCutsMemoryInfoPtr->quitByTailingOff = quitByTailingOff;

  if (rank1CutsArePresent)
  {
    _autoRankOneCutsMemoryInfoPtr->elapsedTime = bapcodInit().startTime().getElapsedTime_dbl()
                                                 - _autoRankOneCutsMemoryInfoPtr->startTime;

  } else
  {
    _autoRankOneCutsMemoryInfoPtr->startTime = bapcodInit().startTime().getElapsedTime_dbl();
    _autoRankOneCutsMemoryInfoPtr->deleteInfoAndBounds();
    _autoRankOneCutsMemoryInfoPtr->nodeEvalInfoPtr = _evalAlgPtr->recordNodeEvalInfo(_treatOrder);
    _autoRankOneCutsMemoryInfoPtr->nodeEvalInfoPtr->numberOfNodes += 1;
    _autoRankOneCutsMemoryInfoPtr->probSetupInfoPtr = _problemSetDownAlgPtr->recordProblemInfo(_treatOrder);
    _autoRankOneCutsMemoryInfoPtr->probSetupInfoPtr->numberOfNodes += 1;
    _autoRankOneCutsMemoryInfoPtr->savedLpDualBound = _evalAlgPtr->algIncLpDualBound();
  }

}

void Node::recordPrimalSol(const VarPtrSet & varPtrSet, const bool printSol)
{
  /// verify whether primal solution has been already recorded
  if (!_primalSol.empty())
    return;

  if (printSol)
    std::cout << "Node solution is : ";

  for (VarPtrSet::const_iterator varPtrIt = varPtrSet.begin(); varPtrIt != varPtrSet.end(); ++varPtrIt)
    {
      _primalSol.push_back(new SolutionVarInfo(*varPtrIt));
      if (printSol)
        std::cout << (*varPtrIt)->name() << "(" << (*varPtrIt)->curCost() << ","
                  << _primalSol.back()->canRoundUp << " ) = " << (*varPtrIt)->val() << "  ";
    }

  if (printSol)
    std::cout << std::endl;
}

void Node::addToPrimalSol(Variable * varPtr)
{
  _primalSol.push_back(new SolutionVarInfo(varPtr));
}

bool Node::isConquered()
{
  return gapSmallerThanTol(_nodeIncIpDualBound, _nodeIncIpPrimalBound, _param);
}

/// Minimisation problem assumed by convention
bool Node::isToBePruned(const Bound & primalBound)
{
  bool nodeShouldBePruned = gapSmallerThanTol(_nodeIncIpDualBound, primalBound, _param);

  if (printL(2))
    std::cout << "try to prune node ref " << ref() << " _nodeIncIpDualBound  = " << _nodeIncIpDualBound
              << " optimalityGapTolerance = " << _param.optimalityGapTolerance()
              << " primalBound = " << primalBound << "  nodeShouldBePruned = " << nodeShouldBePruned << std::endl;

  return nodeShouldBePruned;
}

std::string Node::evalEndTimeString() const
{
  stringstream ss;
  if (_evalEndTime == -1)
    ss << "NE";
  else if (_evalEndTime < 6000)
    {
      ss << _evalEndTime / 100 << "."
         << _evalEndTime % 100 << "s";
    }
  else if (_evalEndTime < 360000)
    {
      ss << _evalEndTime / 6000 << "m"
         << (_evalEndTime % 6000) / 100 << "s";
    }
  else
    {
      ss << _evalEndTime / 360000 << "h"
         << (_evalEndTime % 360000) / 6000 << "m";
    }
  return ss.str();
}

Solution * Node::nodeIncIpPrimalSolPtr()
{
  return _nodeIncIpPrimalSolPtr;
}

bool Node::updateNodeIncPrimalSolution(Solution * solPtr)
{
  if (solPtr == NULL)
    return false;

  Bound primalBound(solPtr->cost(), _objStatus);
  if (printL(5))
    std::cout << "Node::updateNodePrimalBound() node " << "  dualBound " << primalBound
    << " _nodeIncIpPrimalBound "
    << _nodeIncIpPrimalBound << std::endl;

  if (_nodeIncIpPrimalBound > primalBound)
    {
      _nodeIncIpPrimalBound = primalBound;
      if (_nodeIncIpPrimalSolPtr != NULL)
        delete _nodeIncIpPrimalSolPtr;
      _nodeIncIpPrimalSolPtr = solPtr->clone();
      _ipPrimalBoundIsUpdated = true;
    }

  return false;
}


bool Node::updateSubtreeDualBound(const Bound & dualBound)
{
  if (printL(5))
    std::cout << "Node::updateSubtreeDualBound() node " << _treatOrder << "  dualBound " << dualBound << std::endl;

  if (_subtreeDualBound < dualBound)
    {
      if (printL(5))
        std::cout << "Node::updateSubtreeDualBound() node " << _treatOrder << " dualBound = " << dualBound
                  << " better than  _subtreeDualBound = " << _subtreeDualBound << std::endl;

      _subtreeDualBound = dualBound;
      return (true);
    }

  return (false);
}

std::ostream& Node::print(std::ostream& os) const
{
  os << std::endl << "Node  ref = " << ref() << std::endl;

  os << "   depth = " << _depth << std::endl;
  os << "   _BaPOrder = " << _branchAndPriceOrder << std::endl;
  for (std::list<Node *>::const_iterator it =  _sons.begin();
       it != _sons.end(); ++it)
    if ((*it) != NULL)
      os << "   son ref = " << (*it)->ref() << std::endl;
  // os << "   _nodeIncLpPrimalBound = " << _nodeIncLpPrimalBound << std::endl;
  os << "   _nodeIncIpDualBound = " << _nodeIncIpDualBound << std::endl;
  os << "   _subtreeDualBound = " << _subtreeDualBound << std::endl;
  if (_father != NULL)
    os << "   father ref = " << _father->ref() << std::endl;
  else
    os << "   this node is root of the tree " << std::endl;

  if (!_localNodeBrConstrList.empty())
    for (std::list<BranchingConstrBaseType *>::const_iterator it =
         _localNodeBrConstrList.begin(); it != _localNodeBrConstrList.end(); ++it)
      {
        os << "   localNodeBrConstrList = ";
        (*it)->print(os) << std::endl;
      }
  else
    os << "   localNodeBrConstrList undefined " << std::endl;

  return (os);
}

void Node::nicePrint(std::ostream& os) const
{
  if (isRoot())
    {
      os << "BaB tree root node";
      if (printL(3))
        os << " " << ref();
      os << std::endl;
      os << "**** Local DB = " << _nodeIncLpDualBound;
    }
  else
    {
      os << "BaB tree node ";
      if (printL(3))
        os << ref();
      os << "(N° " << _branchAndPriceOrder << ", parent N° " << _fatherBranchAndPriceOrder << ", depth " << _depth;
      if (printL(0))
        os << ", treatOrderId " << _treatOrder;
      os << ")";

      if (_debugSolutionAtThisNode)
          os << "(with debug solution)";

      os << std::endl;
      os << "**** Local DB = " << _nodeIncLpDualBound;

      if (!_localNodeBrConstrList.empty())
        {
          os << ", branch: ";
          for (std::list<BranchingConstrBaseType *>::const_iterator it =
               _localNodeBrConstrList.begin(); it != _localNodeBrConstrList.end(); ++it)
            {
              if (printL(7))
                os << "Printing constr at 0x" << hex << (long)(*it) << dec << endl;
              (*it)->shortPrint(os);
            }
        }
    }
}


bool Node::operator<(const Node & that) const
{
  /// Returns true if this has Less priority than that and zero otherwise
  switch (param().treeSearchStrategy())
  {
    case Node::BestDualBoundThanDF:
    {
      /// Priority to the best dual bound
      if (_nodeIncLpDualBound < that.nodeIncLpDualBound())
        return (false);

      if (_nodeIncLpDualBound > that.nodeIncLpDualBound())
        return (true);

      /// Priority to the deepest node
      if (_depth > that.depth())
        return (false);
      if (_depth < that.depth())
        return (true);

      /// Then give priority to according to the Directives
      if (ref() < that.ref())
        /// Priority to the last child node that was created
        return (false);
      else
        return (true);

      break;
    }
    case Node::DepthFirstWithWorseBound:
    {
      /// Priority to the deepest node
      if (_depth > that.depth())
        return (false);
      if (_depth < that.depth())
        return (true);

      /// nodes with the same depth are solved in the non-increasing order of their dual bounds
      if (_nodeIncLpDualBound < that.nodeIncLpDualBound())
        return (false);

      if (_nodeIncLpDualBound > that.nodeIncLpDualBound())
        return (true);


      if (ref() < that.ref())
        /// Priority to the last child node that was created
        return (false);
      else
        return (true);

      break;
    }
    case Node::BestLpBound:
    {
      /// Priority to the best dual bound
      if (_subtreeDualBound < that.subtreeDualBound())
        return (false);
      if (_subtreeDualBound > that.subtreeDualBound())
        return (true);

      if (_nodeIncLpPrimalBound < that.nodeIncLpPrimalBound())
        return (false);

      if (_nodeIncLpPrimalBound > that.nodeIncLpPrimalBound())
        return (true);

      if (_depth > that.depth())
        return (false); // priority to the deepest node
      if (_depth < that.depth())
        return (true);

      /// Priority to the oldest father node
      if (_fatherBranchAndPriceOrder < that.fatherBranchAndPriceOrder())
        return (false);
      if (_fatherBranchAndPriceOrder > that.fatherBranchAndPriceOrder())
        return (true);
      // else they hang from the same father node

      /// Then give priority to according to the DIRECTIVES
      if (ref() < that.ref())
        /// Priority to the last child node that was created
        return (false);
      else
        return (true);

      break;
    }
    case Node::DepthFirstWithBetterBound:
    {
      /// Priority to the deepest node
      if (_depth > that.depth())
        return (false);
      if (_depth < that.depth())
        return (true);

      /// nodes with the same depth are solved in the non-decreasing order of their dual bounds
      if (_nodeIncLpDualBound > that.nodeIncLpDualBound())
        return (false);
      if (_nodeIncLpDualBound < that.nodeIncLpDualBound())
        return (true);


      if (ref() < that.ref())
        /// Priority to the last child node that was created
        return (false);
      else
        return (true);

      break;
    }
  }

  return (true);
}

bool Node::secondaryLessThan(const Node & that) const
{
    /// Secondary search strategy is fixed to Node::DepthFirstWithWorseBound

    /// Priority to the deepest node
    if (_depth > that.depth())
      return (false);
    if (_depth < that.depth())
      return (true);

    /// nodes with the same depth are solved in the non-increasing order of their dual bounds
    if (_nodeIncLpDualBound < that.nodeIncLpDualBound())
      return (false);
    if (_nodeIncLpDualBound > that.nodeIncLpDualBound())
      return (true);

    /// Priority to the last child node that was created
    return ref() > that.ref();
}

bool smallerNodePt::operator()(Node * a, Node * b) const
{
  return (*a < *b);
}

bool secondarySmallerNodePt::operator()(Node * a, Node * b) const
{
  return (*a).secondaryLessThan(*b);
}

BapcodInit & Node::bapcodInit() const
{
  return probConfPtr()->bapcodInit();
}

BapcodInit * Node::bapcodInitPtr() const
{
  return probConfPtr()->bapcodInitPtr();
}

const ControlParameters& Node::param() const
{
  return probConfPtr()->param();
}

ControlParameters & Node::param()
{
  return probConfPtr()->param();
}

ProblemSetupInfo * Node::probSetupInfoPtr() const
{
  return _probSetupInfoPtr;
}

void Node::makeProbSetupInfoObligatoryForFullSetup()
{
  _probSetupInfoPtr->setFullSetupIsObligatory(true);
}

void Node::removeProbSetupObligationForFullSetup()
{
  _probSetupInfoPtr->setFullSetupIsObligatory(false);
}

NodeEvalInfo * Node::nodeEvalInfoPtr() const
{
  return _nodeEvalInfoPtr;
}

GenChildNodesInfo * Node::genChildNodesInfoPtr() const
{
  return _genChildNodesInfoPtr;
}

void Node::associateProblemSetupInfoPtr(ProblemSetupInfo * const problemSetupInfoPtr)
{
  _probSetupInfoPtr = problemSetupInfoPtr;
  if (_probSetupInfoPtr != NULL)
    _probSetupInfoPtr->numberOfNodes += 1;
}

void Node::associateNodeEvalInfoPtr(NodeEvalInfo * const nodeSolverInfoPtr)
{
  _nodeEvalInfoPtr = nodeSolverInfoPtr;
  if (_nodeEvalInfoPtr != NULL)
    _nodeEvalInfoPtr->numberOfNodes += 1;

}

void Node::associateGenChildNodesInfoPtr(GenChildNodesInfo * const genChildNodesInfoPtr)
{
  _genChildNodesInfoPtr = genChildNodesInfoPtr;
  if (_genChildNodesInfoPtr != NULL)
    _genChildNodesInfoPtr->numberOfNodes += 1;
}

void Node::removeProblemSetupInfoAssociation()
{
  if (_probSetupInfoPtr != NULL)
    {
      _probSetupInfoPtr->numberOfNodes -= 1;
      if (_probSetupInfoPtr->numberOfNodes == 0)
        delete _probSetupInfoPtr;
      _probSetupInfoPtr = NULL;
    }
}

void Node::removeNodeEvalInfoAssociation()
{
  if (_nodeEvalInfoPtr != NULL)
    {
      _nodeEvalInfoPtr->numberOfNodes -= 1;
      if (_nodeEvalInfoPtr->numberOfNodes == 0)
        delete _nodeEvalInfoPtr;
      _nodeEvalInfoPtr = NULL;
    }
}

void Node::removeGenChildNodesInfoAssociation()
{
  if (_genChildNodesInfoPtr != NULL)
    {
      _genChildNodesInfoPtr->numberOfNodes -= 1;
      if (_genChildNodesInfoPtr->numberOfNodes == 0)
        delete _genChildNodesInfoPtr;
      _genChildNodesInfoPtr = NULL;
    }
}

void Node::setBranchEvaluationInfo(BranchingEvaluationInfo * branchEvalInfoPtr,
                                   const int & strongBranchNodeNumber)
{
  _branchEvalInfoPtr = branchEvalInfoPtr;
  _strongBranchPhaseNumber = 1;
  _strongBranchNodeNumber = strongBranchNodeNumber;
}

void Node::setBranchEvaluationInfoFromNode(const Node * prevPhaseNodePtr)
{
  _branchEvalInfoPtr = prevPhaseNodePtr->_branchEvalInfoPtr;
  _strongBranchPhaseNumber = prevPhaseNodePtr->_strongBranchPhaseNumber + 1;
  _strongBranchNodeNumber = prevPhaseNodePtr->_strongBranchNodeNumber;
}

void Node::applyAutoRankOneCutsMemorySearch(AutoRankOneCutsMemoryInfo * infoPtr)
{
  _applyAutoRankOneCutsMemorySearch = true;
  if (_autoRankOneCutsMemoryInfoPtr != NULL)
    delete _autoRankOneCutsMemoryInfoPtr;
  _autoRankOneCutsMemoryInfoPtr = infoPtr;
}
