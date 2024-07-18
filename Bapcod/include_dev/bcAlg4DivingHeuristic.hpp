/**
 *
 * This file bcAlg4DivingHeuristic.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCDIVINGHEURISTIC_HPP__
#define BCDIVINGHEURISTIC_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcNodeC.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4GenChildrenOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4EvalByColAndCutGen.hpp"

struct DiveInfo : public GenChildNodesInfo
{
private:
  VarPtrSet _tabuVariables;

public:
  int remainingDepth;
  int maxTabuSize;

  /// This following param was introduced for the following use-case:
  /// 1) When enabling dive Heur with usingDiveOnMastColOnly, we wanted to fix pure master varsat the end of each
  ///    dive by using a rest. mip heur while keeping only the columns that were fixed during the dive.
  /// 2) To implement the restricted master heuristic with the false gap enumeration
  ///    (enumeration is performed in Algorithm4DivingEval, and then the restricted master is called without diving)
  /// = 0 if restricted master is not called
  /// = 1 if both evaluation algorithm and restricted master are called
  /// = 2 if only restricted master is called
  int performRestrictedMaster;

  DiveInfo(const int depth, const int tabuSize, const int performRestrictedMaster_ = 0);
  DiveInfo(const DiveInfo & that);
  DiveInfo(const VarPtrSet & tabuVars, const int depth, const int tabuSize, const int performRestrictedMaster_ = 0);
  ~DiveInfo();

  const VarPtrSet & tabuVariables() const;
  void tabuVariables(const VarPtrSet & tabuVars);

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os << "DiveInfo with number of nodes = " << numberOfNodes << std::endl;
    return os;
  }

};

struct SmallestNodeLpValue
{
  bool operator()(Node *  nodeA, Node  * nodeB) const
  {
    return (nodeA->nodeIncLpPrimalBound() < nodeB->nodeIncLpPrimalBound());
  }
} ;

struct LargestNodeDepth
{
  bool operator()(Node *  nodeA, Node  * nodeB) const
  {
    if (nodeA->depth() != nodeB->depth())
      return (nodeA->depth() < nodeB->depth());
    if ((nodeA->genChildNodesInfoPtr() != NULL) && (nodeB->genChildNodesInfoPtr() != NULL))
      {
        DiveInfo * diveInfoA = static_cast<DiveInfo *>(nodeA->genChildNodesInfoPtr());
        DiveInfo * diveInfoB = static_cast<DiveInfo *>(nodeB->genChildNodesInfoPtr());
        if (diveInfoA->tabuVariables().size() != diveInfoB->tabuVariables().size())
          return (diveInfoA->tabuVariables().size() > diveInfoB->tabuVariables().size());
      }

    return (nodeA->ref() > nodeB->ref());
  }
};

struct DivingEvalInfo : public ColGenEvalInfo
{
  int nbNeededProperColumns;
  bool runCutGeneration;
  int enumerationWithFalseGapMode;

  DivingEvalInfo() :
    ColGenEvalInfo(), nbNeededProperColumns(0), runCutGeneration(true), enumerationWithFalseGapMode(0)
  {
  }

  DivingEvalInfo(const ColGenEvalInfo & conGenEvalInfo, const int nbNeededProperColumns_,
                 const bool runCutGeneration_ = true, const int enumerationWithFalseGapMode_ = 0) :
    ColGenEvalInfo(conGenEvalInfo), nbNeededProperColumns(nbNeededProperColumns_),
    runCutGeneration(runCutGeneration_), enumerationWithFalseGapMode(enumerationWithFalseGapMode_)
  {
  }

  DivingEvalInfo(const int nbNeededProperColumns_, const bool runCutGeneration_ = true,
                 const int enumerationWithFalseGapMode_ = 0) :
    ColGenEvalInfo(), nbNeededProperColumns(nbNeededProperColumns_), runCutGeneration(runCutGeneration_),
    enumerationWithFalseGapMode(enumerationWithFalseGapMode_)
  {
    if (printL(5))
      std::cout << "ColGenEvalInfo with " << *masterLpBasisPtr << " is created " << std::endl;
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os << "DivingEvalInfo with number of nodes = " << numberOfNodes
       << " and nbNeededProperColumns = " << nbNeededProperColumns << std::endl;
    return os;
  }

};

class Algorithm4DivingEval : public Alg4EvalByColAndCutGen
{
    int _nbNeededProperColumns;
    bool _runCutGeneration;
    int _enumerationWithFalseGapMode;

    void addFixedSolutionBasedCutsToMaster();

public:

    Algorithm4DivingEval(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons):
            Alg4EvalByColAndCutGen(probPtr, masterCommons), _nbNeededProperColumns(0), _runCutGeneration(true),
            _enumerationWithFalseGapMode(0)
    {
    }

    virtual bool setupAlgo(Node * nodePtr);
    virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr);
    virtual bool eval();
};

class DiveAlgorithm : public Alg4GenChildrenOfNode
{
protected:
  VarPtrSet _tabuVariables;
  int _remainingDepth;
  int _maxTabuSize;
  MultitokenSelectionStrategyVector & _colSelectionCriteria;

  bool fixVariables();
  bool solutionCausesInfeasibility(int & globalTreatOrder, Solution * solPtr);
  Solution * roundVariable(int & globalTreatOrder);

  virtual void prepareCandidateNodeInStrongDive(Node * nodePtr, int & globalTreatOrder);
public:

  DiveAlgorithm(MasterCommons4GenChildNodesAlgorithm & masterCommons, MultitokenSelectionStrategyVector & colSelectionCriteria) :
      Alg4GenChildrenOfNode(masterCommons), _tabuVariables(), _colSelectionCriteria(colSelectionCriteria)
  {
  }

  virtual ~DiveAlgorithm()
  {
  }

  virtual bool setupAlgo(Node * nodePtr);

  void runStrongDive(int & globalTreatOrder, const int maxNumberOfCandidates,
                     const int maxNumberOfDives, const bool inheritParentDualBound);

  virtual void run(int & globalTreatOrder);

  virtual void setDownAlgo();
};

class DivingHeuristic : public Alg4PrimalHeuristicOfNode
{

  void printDivingNodeInformation(const Node * nodePtr, const int diveNumber);
  virtual Alg4EvalOfNode * createEvaluationAlgorithm(const bool runColGen);
  void replaceNodeInNonProperStrongDiving(Node **nodePtr);
  bool prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder);
  virtual void runBody(int & globalTreatOrder);

protected:
  bool _stopAfterFirstSolutionFound;
  int _maxDepth;
  int _maxDiscrepancy;
  bool _runRestrMasterAfterFalseGapEnumeration;

  MultitokenSelectionStrategyVector & _colSelectionCriteria;
  VarPtr2DoubleMap _referenceSolution; /// reference solution is needed for the local search heuristic
                                       /// even if infeasibility is reached, a partial solution is kept here
  bool runDiving(int & globalTreatOrder, Node * nodePtr);

public:
  void setOptionMaxDepth(const int value);
  void setOptionMaxDiscrepancy(const int value);
  void setOptionStopAfterFirstSolutionFound(const bool value);
  void setOptionRunRestrMasterAfterFalseGapEnumeration(const bool value);

  DivingHeuristic(Problem * probPtr,
                  MasterCommons4PrimalHeuristic & masterCommons,
                  MultitokenSelectionStrategyVector & colSelectionCriteria);

  virtual ~DivingHeuristic();
};

class LocalSearchHeuristic : public DivingHeuristic
{
  int _maxIterationsNumber;
  double _fixedVarsRatio;

  Solution * fixPartialSolution();
  virtual void runBody(int & globalTreatOrder);

public:
  void setOptionMaxIterationsNumber(const int value);
  void setOptionFixedVarsRatio(const double value);

  LocalSearchHeuristic(Problem * probPtr,
                       MasterCommons4PrimalHeuristic & masterCommons,
                       MultitokenSelectionStrategyVector & colSelectionCriteria);

  virtual ~LocalSearchHeuristic();
};

#endif /* BCDIVINGHEURISTIC_HPP_ */
