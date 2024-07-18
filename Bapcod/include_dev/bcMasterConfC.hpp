/**
 *
 * This file bcMasterConfC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMASTERCONFC_H_
#define BCMASTERCONFC_H_

#include "bcProbConfigC.hpp"
#include "bcSolutionC.hpp"
#include "bcNodeC.hpp"
#include "bcAlg4EvalOfNode.hpp"

class ConvexityGenConstr;

class MasterConf: public ProbConfig
{
  friend class MasterCommons4ProblemSetup;
  friend class MasterCommons4EvalAlg;
  friend class MasterCommons4GenChildNodesAlgorithm;  
  friend class MasterCommons4PrimalHeuristic;

protected:
  
  MasterCommons4ProblemSetup _masterCommons4ProblemSetup;
  MasterCommons4EvalAlg _masterCommons4EvalAlg;
  MasterCommons4GenChildNodesAlgorithm _masterCommons4GenChildNodes;
  MasterCommons4PrimalHeuristic _masterCommons4PrimalHeuristic;

  std::list<ConvexityGenConstr *> _convexityGenConstrPtrList;

  /**
   * To hold variables other than dynamically generated columns
   */
  InstMasterVarPtrSet _setOfPureMasterVar;

  /**
   * To hold constraints other than dynamically generated convexity constraints
   */
  InstMasterConstrPtrSet _setOfPureMasterConstr;

  Solution * _enumSolutionPtr;
  Solution * _debugSolutionPtr;

  std::list<Variable * > _nonStabStaticArtVarPtrList;

  Solution * _initialSetOfActiveColumns; /// these columns do not necessarily form a solution
  Solution * _initialSetOfInactiveColumns; /// these columns do not necessarily form a solution

  Bound _bestPrimalBound;
  Bound _bestDualBound;
  
  std::vector< Node *> _bapTreeNodes;
  std::vector<std::pair<int, double> > _subTreeSizeByDepth;

  void initializeBaPTreeDotFile() const;
  void addNodeToBaPTreeDotFile(Node * nodePtr) const;
 
  void bestPrimalBound(const Bound & val) {_bestPrimalBound = val;}
  void bestDualBound(const Bound & val) {_bestDualBound = val;}

  void updateCurValidDualBound(Node * nodePtr);
  std::ostream& printInfoBeforeSolvingNode(const Node * nodePtr, int numPrimOpenNodes, int numSecOpenNodes,
                                           int estimTreeSize, std::ostream& os = std::cout) const;

  friend class ColGenSpConf;

  void insertColGenSpConf(ColGenSpConf * cgSpConfPtr);

  void readDebugSolutionFromFile(const std::string & fileName);

  int calculateTreeSizeEstimation();

public:
    
  const Bound & bestPrimalBound() const {return _bestPrimalBound;}
  const Bound & bestDualBound() const {return _bestDualBound;}

  virtual bool updateDualIncBound(const Bound & newDualBound);
  
  virtual bool updatePrimalIncBound(const Bound & newIncVal);

  MasterConf(Model * modelPtr, 
	         Problem * problemPtr,
	         const Bound & initialPrimalBound,
	         const Bound & initialDualBound);

  virtual ~MasterConf();

  void resetColGenSpListOfFractMastCol(const MasterColSolution & curListOfFractMastCol);

  virtual Constraint * castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately = false);
  virtual InstanciatedConstr * castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately = false);
  virtual Variable * castAndAddVariable(Variable * varPtr, const bool & insertImmediately = false);
  virtual InstanciatedVar * castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately = false);

  virtual void insertPureVar(InstanciatedVar * ivarPtr);
  virtual void insertPureConstr(InstanciatedConstr * iconstrPtr);

  /**
   * Set initial variables (including artificial and missing columns)
   * and constraints (including convexity constraints) and build formulation
   */
  virtual void prepareProbConfig();

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  std::ostream& printSol(std::ostream& os = std::cout) const;

  void runGreedyHeuristic(int & globalNodesTreatOrder);
  virtual Node * createRootNode();
  virtual bool prepareNodeForTreatment(Node * nodePtr,
                                       const int globalNodesTreatOrder,
                                       const int thisSearchTreeTreatedNodesNumber);
  bool getAverageSubtreeSize(const int & depth, Double & averSubtreeSize) const;

  virtual Solution * enumerateAllColumns(int & nbEnumColumns);
  virtual Solution * getDebugSolution();
  virtual Solution * solvePC();
  virtual void recordInitialActiveSetOfColumns(Solution * solPtr);
  virtual void recordInitialInactiveSetOfColumns(Solution * solPtr);
  virtual void updatePrimalIncSolution(const Bound & primalBound, Solution * solPtr);
  virtual Solution * getDissagregatedSolution(Solution * curSolPtr); /// project solution in original space if need be
  virtual Solution * getAggregatedSolution(Solution * curSolPtr);

  const InstMasterVarPtrSet & setOfPureMasterVar() const
  {
    return _setOfPureMasterVar;
  }

  const InstMasterConstrPtrSet & setOfPureMasterConstr() const
  {
    return _setOfPureMasterConstr;
  }

  void addArtVar(MissingColumn * artVarPtr);
  void addGlobalArtVar(GlobalArtificialVar * artVarPtr);
  void addNonStabilizationLocalArtVar(LocalArtificialVar * artVarPtr);

  // To add instanciated variables to the master problem
  virtual void addVariablesToForm();

  virtual bool isTypeOf(const PcId::PcIdentifier & pcIdentifier) const
  {
      return compareIdentifier(PcId::MasterMask, pcIdentifier);
  }

};

inline std::ostream& operator<<(std::ostream& os, const MasterConf & that)
{
  return that.print(os);
}


#endif /* BCMASTERCONFC_H_ */
