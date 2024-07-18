/**
 *
 * This file bcNodeC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BRANCHaBOUNDCLASSES_H
#define BRANCHaBOUNDCLASSES_H
#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcProblemC.hpp"
#include "bcProbConfigC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMasterCommons.hpp"
#include "bcAlg4PrimalHeurOfNode.hpp"
#include "bcGenBranchingConstrC.hpp"

class Node;
class BranchingConstrBaseType;
class BranchingConstrGenerator;
class Solution;
struct LpBasisRecord;
struct NodeEvalInfo;
struct ProblemSetupInfo;
struct GenChildNodesInfo;
struct BranchingEvaluationInfo;
class Alg4PreprocessingOfNode;
class Alg4EvalOfNode;
class Alg4ProblemSetupOfNode;
class Alg4ProblemSetDownOfNode;
class Alg4GenChildrenOfNode;
struct Alg4MasterSolAndBounds;

struct AutoRankOneCutsMemoryInfo
{
  BcObjStatus::MinMaxIntFloat objStatus;
  ProblemSetupInfo * probSetupInfoPtr;
  NodeEvalInfo * nodeEvalInfoPtr;
  Alg4EvalOfNode * evalAlgPtr;
  Alg4ProblemSetupOfNode * problemSetupAlgPtr;
  Alg4ProblemSetDownOfNode * problemSetDownAlgPtr;
  Bound savedLpDualBound;
  double startTime;
  double elapsedTime;
  bool quitByTailingOff;

  void deleteInfoAndBounds();
  AutoRankOneCutsMemoryInfo(BcObjStatus::MinMaxIntFloat objStatus_);
  ~AutoRankOneCutsMemoryInfo();
};


class Node
{
public:

  enum SearchStrategy
  {
    BestDualBoundThanDF = 0,
    DepthFirstWithWorseBound = 1,
    BestLpBound = 2,
    DepthFirstWithBetterBound = 3
  } ;
private:

  ProbConfig * _probConfigPtr;
  BcObjStatus::MinMaxIntFloat _objStatus;
  int _ref;
  const ControlParameters & _param;
  Node * _father;
  std::list<Node *> _sons;
  int _depth;
  bool _prunedAtBeginningOfTreatNode;
  bool _debugSolutionAtThisNode;
    
  Double _estimatedSubtreeSize;
  Double _branchingPriorityLevel;
  int _subtreeSize;

  /// members related to bounds where algorith store their results
  Bound _nodeIncLpDualBound;  /// incumbent valid dual bound for the master LP // _nodeFractionalDualBound;
  Bound _nodeIncIpDualBound; /// incumbent valid dual bound for the master IP //_curIpDualVal;
  Bound _nodeIncLpPrimalBound; /// incumbent valid dual bound for the master LP
  Bound _nodeIncIpPrimalBound; /// incumbent valid dual bound for the master IP
  Bound _subtreeDualBound;
  
  bool _dualBoundIsUpdated;
  bool _ipPrimalBoundIsUpdated;

  Solution * _nodeIncIpPrimalSolPtr;  /// used in MasterConf::extractNodeLpSol()

  std::list<BranchingConstrBaseType *> _localNodeBrConstrList;
  Solution * _localFixedSolution;

  //This will be the time and order in which the node is treated in the BaP tree
  long _evalEndTime;
  int _treatOrder;
  int _branchAndPriceOrder; // only used for printing BaP nodes (represents evaluation order in the BaP tree)
  int _fatherBranchAndPriceOrder; //  = _father->branchAndPriceOrder()
  bool _infeasible; //Issam: added mainly for dot file generation
  bool _evaluated;
  bool _treated;

  bool _problemAndEvalAlgInfoSaved;
  SolutionVarInfoPtrList _primalSol;
  ProblemSetupInfo * _probSetupInfoPtr;
  NodeEvalInfo * _nodeEvalInfoPtr;
  GenChildNodesInfo * _genChildNodesInfoPtr;

  AutoRankOneCutsMemoryInfo * _autoRankOneCutsMemoryInfoPtr;
  bool _applyAutoRankOneCutsMemorySearch;

  /// pointer to the branching evaluation information
  /// (used to update the branching history after this node being evaluated)
  BranchingEvaluationInfo * _branchEvalInfoPtr;
  int _strongBranchPhaseNumber;
  int _strongBranchNodeNumber;

  Alg4PreprocessingOfNode * _preprocessorPtr;
  Alg4EvalOfNode * _evalAlgPtr;
  Alg4ProblemSetupOfNode * _problemSetupAlgPtr;
  Alg4ProblemSetDownOfNode * _problemSetDownAlgPtr;
  Alg4GenChildrenOfNode * _genChildNodesAlgPtr;
  std::set<Alg4PrimalHeuristicOfNode *, PrimalHeurSort> _PrimalHeuristicPts;

  ColClassesVector _treeOfColClasses;
  ColClassesVectorsMap _cgSpConfTreeOfColClassesMap;
  
  Node() = delete;
  
  //This function stores a solution in the node
  //This should be called by node after Eval and after running a heuristic.
  void recordIpPrimalSolAndUpdateIpPrimalBound(const Bound & primalBound, Solution *solPtr);
  void recordIpPrimalSolAndUpdateIpPrimalBound(const Bound & primalBound, const VarPtr2DoubleMap & primalSolMap);
  void updateDualBounds(const Bound & lpDualBound, const Bound & ipDualBound);

public:

  /**
   * To be used to create the root node
   */
  Node(ProbConfig * probConfigPtr,
       const Bound & dualBound,
       ProblemSetupInfo * varConstrResetInfoPtr,
       NodeEvalInfo * nodeEvalInfoPtr,
       bool debugSolutionIsDefined = false);

  /**
   * To be use to create child nodes
   */
  Node(int ref,
       Node * father,
       const std::list<BranchingConstrBaseType *> & localNodeBrConstrList,
       Solution * localFixedSolution = NULL,
       bool inheritDualBound = true);

  virtual ~Node();

  Problem * probPtr() const;
  int ref() const {return _ref;}

  /// start added by Ruslan

  ConstrPtrList upCastedBranchingConstrList();

  bool markInfeasibleAndExitTreatment();

  bool exitTreatment(const bool exitValue = true);

  bool autoRankOneCutMemoryEvaluation(int & globalTreatOrder, const Bound & incPrimalBound);
  bool evaluation(int & globalTreatOrder, const Bound & incPrimalBound);
  bool treat(int & globalTreatOrder, const Bound & incPrimalBound);    

  void storeBranchingEvaluationInfo(double time);
  bool runEnumeratedMIP();
    
  void clearLocalNodeBrConstrList(const bool & calledFromDestructor);

  void setProblemSetupAlgorithm(Alg4ProblemSetupOfNode * problemSetupAlgPtr);
  void setProblemSetDownAlgorithm(Alg4ProblemSetDownOfNode * problemSetDownAlgPtr);
  void setPreprocessor(Alg4PreprocessingOfNode * preprocessorPtr);
  void setEvalAlg(Alg4EvalOfNode * nodeEvalAlgPtr);
  void setGenChildNodesAlgorithm(Alg4GenChildrenOfNode * genChildNodesAlgPtr);
  void setDualBoundEqualToIncPrimalBound();
  void addPrimalHeuristic(Alg4PrimalHeuristicOfNode * PrimalHeuristicPtr);
  bool needProblemFullSetDownAlgorithm();
  /// end added by Ruslan

  virtual const int & depth() const
  {
    return (_depth);
  }

  const bool debugSolutionAtThisNode() const
  {
      return _debugSolutionAtThisNode;
  }

  void removeDebugSolutionFromThisNode()
  {
      _debugSolutionAtThisNode = false;
  }

  const bool isRoot() const
  {
    return (_depth == 0);
  }

  const Double & estimatedSubtreeSize() const
  {
    return _estimatedSubtreeSize;
  }

  const Double & branchingPriorityLevel() const
  {
    return _branchingPriorityLevel;
  }

  void calculateSubtreeSize(std::vector<std::pair<int, double> > & subTreeSizeByDepth);
  
  const int & subtreeSize() const
  {
    return _subtreeSize;
  }

  void setEstimatedSubtreeSize(const Double & value)
  {
    _estimatedSubtreeSize = value;
  }

  void setBranchingPriorityLevel(const Double & value)
  {
    _branchingPriorityLevel = value;
  }

    const bool solved() const
  {
    return _evaluated;
  }
  
  void setSolved(bool evaluated)
  {
    _evaluated = evaluated;
  }
  
  const bool treated() const
  {
    return _treated;
  }

  const SolutionVarInfoPtrList & primalSol() const
  {
    return _primalSol;
  }
    
  void saveProblemAndEvalAlgInfo();

  void saveAutoRankOneCutsMemoryInfo(const bool quitByTailingOff, const bool rank1CutsArePresent = true);

  void deleteProblemAndEvalAlgInfo();
  
  void recordPrimalSol(const VarPtrSet & varPtrSet, const bool printSol = false);

  void addToPrimalSol(Variable * varPtr);
  
  virtual bool isConquered();

  ColClassesVector & treeOfColClasses()
  {
    return _treeOfColClasses;
  }

  ColClassesVectorsMap & cgSpConfTreeOfColClassesMap()
  {
    return _cgSpConfTreeOfColClassesMap;
  }

  /**
   * Minimisation problem assumed by convention
   */
  virtual bool isToBePruned(const Bound & primalBound);

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  virtual void nicePrint(std::ostream& os = std::cout) const;

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;
  const ControlParameters& param() const;
  ControlParameters& param();

  virtual bool operator<(const Node & that) const;
  
  virtual bool secondaryLessThan(const Node & that) const;

  virtual ProbConfig * probConfPtr() const
  {
    return _probConfigPtr;
  }

  virtual std::list<Node *> & sons()
  {
    return (_sons);
  }

  const Bound & nodeIncIpDualBound() const
  {
    return _nodeIncIpDualBound;
  }

  const Bound & nodeIncLpDualBound() const
  {
    return _nodeIncLpDualBound;
  }
 
 
  const Bound & nodeIncIpPrimalBound() const
  {
    return _nodeIncIpPrimalBound;
  }

  const Bound & nodeIncLpPrimalBound() const
  {
    return _nodeIncLpPrimalBound;
  }

  
  void nodeIncLpDualBound(const Bound & newVal)
  {
    _nodeIncLpDualBound = newVal;
  } 
  
  void nodeIncIpDualBound(const Bound & newVal)
  {
    _nodeIncIpDualBound = newVal;
  } 
  
   void nodeIncIpPrimalBound(const Bound & newVal)
  {
    _nodeIncIpPrimalBound = newVal;
  } 
  void nodeIncLpPrimalBound(const Bound & newVal)
  {
    _nodeIncLpPrimalBound = newVal;
  } 

  const bool dualBoundIsUpdated() const
  {
    return _dualBoundIsUpdated;
  }

  const bool primalBoundIsUpdated() const
  {
    return _ipPrimalBoundIsUpdated;
  }
    
  virtual bool prunedAtBeginningOfTreatNode() const {return _prunedAtBeginningOfTreatNode;}
  virtual void prunedAtBeginningOfTreatNode(bool newVal) {_prunedAtBeginningOfTreatNode = newVal;}
  virtual bool infeasible() const {return _infeasible;}

  std::string evalEndTimeString() const;

  virtual void treatOrder(int newVal) {
    _treatOrder = newVal;
  }

  virtual int treatOrder() const
  {
    return _treatOrder;
  }

  inline const int & branchAndPriceOrder() const
  {
    return _branchAndPriceOrder;
  }

  inline const int & fatherBranchAndPriceOrder() const
  {
    return _fatherBranchAndPriceOrder;
  }
  
  
  virtual void branchAndPriceOrder(int newVal)
  {
    _branchAndPriceOrder = newVal;
  }
    
  virtual Node * father() const;

  virtual Solution * nodeIncIpPrimalSolPtr();
  
  virtual bool updateNodeIncPrimalSolution(Solution * solPtr);

  virtual bool updateSubtreeDualBound(const Bound & dualBound);

  const Bound & subtreeDualBound() const
  {
    return _subtreeDualBound;
  }

  virtual const std::list<BranchingConstrBaseType *> & localNodeBrConstrList() const
  {
    return (_localNodeBrConstrList);
  }

  virtual Solution * localFixedSolution()
  {
    return _localFixedSolution;
  }

  virtual const Solution * localFixedSolution() const
  {
    return _localFixedSolution;
  }
  void clearLocalFixedSolution();
  void printFixedSolution(std::ostream& os = std::cout, const bool printMasterSolution = true) const;

  ProblemSetupInfo * probSetupInfoPtr() const;
  void makeProbSetupInfoObligatoryForFullSetup();
  void removeProbSetupObligationForFullSetup();
  NodeEvalInfo * nodeEvalInfoPtr() const;
  GenChildNodesInfo * genChildNodesInfoPtr() const;

  void associateProblemSetupInfoPtr(ProblemSetupInfo * const problemSetupInfoPtr);
  void associateNodeEvalInfoPtr(NodeEvalInfo * const nodeEvalInfoPtr);
  void associateGenChildNodesInfoPtr(GenChildNodesInfo * const genChildNodesInfoPtr);
  void removeProblemSetupInfoAssociation();
  void removeNodeEvalInfoAssociation();
  void removeGenChildNodesInfoAssociation();
  
  /// this function is for a node of phase one
  void setBranchEvaluationInfo(BranchingEvaluationInfo * branchEvalInfoPtr, const int & strongBranchNodeNumber);
  /// this function is for a node of other phases
  void setBranchEvaluationInfoFromNode(const Node * prevPhaseNodePtr);
  
  void setBranchPhaseNumber(const int & strongBranchPhaseNumber)
  {
    _strongBranchPhaseNumber = strongBranchPhaseNumber;
  }

  void applyAutoRankOneCutsMemorySearch(AutoRankOneCutsMemoryInfo * infoPtr);
};

inline std::ostream& operator<<(std::ostream& os, const Node & that)
{
  return that.print(os);
}

struct smallerNodePt
{
  bool operator()(Node *  a, Node * b) const;
};

struct secondarySmallerNodePt
{
  bool operator()(Node *  a, Node * b) const;
};

struct BestNodeLpValue
{
  bool operator()(Node *  nodeA, Node  * nodeB) const
  {
    if (nodeA->nodeIncLpPrimalBound() == nodeB->nodeIncLpPrimalBound())
        return (nodeA->ref() < nodeB->ref());
    return (nodeA->nodeIncLpPrimalBound() > nodeB->nodeIncLpPrimalBound());
  }
};

#endif // BRANCHaBOUNDCLASSES_H
