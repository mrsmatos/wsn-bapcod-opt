/**
 *
 * This file bcMasterCommons.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

// bcMasterCommons.h
//

#ifndef LZZ_bcMasterCommons_hpp
#define LZZ_bcMasterCommons_hpp

#include "bcUsefulHeadFil.hpp"

#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"

class MasterConf;
class Bound;
class Double;
class Node;
class ControlParameters;
class Problem;
class BapcodInit;
struct BranchingGeneratorHistory;

class MasterCommons4ProblemSetup
{
  MasterConf & _masterConf;

public:
  MasterCommons4ProblemSetup(MasterConf & masterConf);

  const std::list<Problem *> & problemList() const;
  const std::vector<ColGenSpConf *> & colGenSubProbConfPts() const;
  Solution * debugSolution() const;
};

class MasterCommons4EvalAlg
{
  MasterConf & _masterConf;  
  
public:
  MasterCommons4EvalAlg(MasterConf & masterConf);

  const std::vector<ColGenSpConf *> & colGenSubProbConfPts() const;

  /// returns -1 if at least one subproblem is not enumerated. needed for cut generation
  long int totalNumberOfEnumeratedSubprobSolutions() const;

  BcObjStatus::MinMaxIntFloat objStatus() const;
  InstanciatedConstr * castAndAddConstraint(InstanciatedConstr * iconstrPtr, bool const & insertImmediately);
  std::list<Constraint *> & pcDelayedConstrPtrList();
  std::list<Constraint *> & pcConstrPtrList();
  const std::set <GenericBranchingConstr * , DynamicGenConstrSort> & candidateBranchingGenericConstr() const;
  const std::set <GenericCutConstr *, DynamicGenConstrSort> & candidateCutGenericConstr();
  Solution * debugSolution() const;
} ;

class MasterCommons4GenChildNodesAlgorithm
{
  MasterConf & _masterConf;

public:
  MasterCommons4GenChildNodesAlgorithm(MasterConf & masterConf);

  const std::vector<ColGenSpConf *> & colGenSubProbConfPts() const;
  const std::set<GenericBranchingConstr * , DynamicGenConstrSort> & candidateBranchingGenericConstr() const;
  const std::set<GenericCutConstr *, DynamicGenConstrSort> & candidateCutGenericConstr() const;
  bool getAverageSubtreeSize(const int & depth, Double & averSubtreeSize) const;
  BcObjStatus::MinMaxIntFloat objStatus() const;

  const std::list<Problem *> & problemList() const;
  MasterCommons4ProblemSetup & masterCommons4ProblemSetup();
  MasterCommons4EvalAlg & masterCommons4EvalAlg();
  int getNodeCountAndIncreaseIt();
};

class MasterCommons4PrimalHeuristic
{
  MasterConf & _masterConf;

public:
  MasterCommons4PrimalHeuristic(MasterConf & masterConf);

  const std::list<Problem *> & problemList() const;
  const std::vector<ColGenSpConf *> & colGenSubProbConfPts() const;
  MasterCommons4ProblemSetup & masterCommons4ProblemSetup();
  MasterCommons4EvalAlg & masterCommons4EvalAlg();
  MasterCommons4GenChildNodesAlgorithm & masterCommons4GenChildNodes();
};

#endif
