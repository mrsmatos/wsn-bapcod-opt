/**
 *
 * This file bcGreedyHeuristic.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCGREEDYHEURISTIC_HPP__
#define BCGREEDYHEURISTIC_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcNodeC.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4GenChildrenOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4EvalBySimplexBasedColGen.hpp"
#include "bcAlg4EvalByMip.hpp"

struct GreedyEvalInfo: public NodeEvalInfo
{
  int subProblemIndex;

  GreedyEvalInfo() :
    NodeEvalInfo(), subProblemIndex(0)
  {
  }

  GreedyEvalInfo(const int subProbIndex) :
    NodeEvalInfo(), subProblemIndex(subProbIndex)
  {
  }

  virtual ~GreedyEvalInfo()
  {
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os << "GreedyEvalInfo with number of nodes = " << numberOfNodes
        << " and subProblemIndex = " << subProblemIndex << std::endl;
    return os;
  }

};


class Algorithm4GreedyEval : public Alg4EvalByMip
{
  bool _thereArePureMasterVariables;
  int _subProblemIndex;

  bool partialSolutionIsFeasible();
  void updateSubProbVarsCosts();
  void recordSolution(ColGenSpConf * cgSpConfPtr, Solution * solPtr);
public:

  Algorithm4GreedyEval(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons, const bool pureMasterVars) :
    Alg4EvalByMip(probPtr, masterCommons), _thereArePureMasterVariables(pureMasterVars), _subProblemIndex(0)
  {
  }

  virtual ~Algorithm4GreedyEval()
  {
  }

  bool setupAlgo(Node * nodePtr);

  virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL);

  void setDownAlgo();

  virtual bool eval();

} ;


class Algorithm4GreedyDive : public Alg4GenChildrenOfNode
{
public:

  Algorithm4GreedyDive(MasterCommons4GenChildNodesAlgorithm & masterCommons) :
      Alg4GenChildrenOfNode(masterCommons)
  {
  }

  virtual ~Algorithm4GreedyDive()
  {
  }

  virtual bool setupAlgo(Node * nodePtr);

  virtual void run(int & globalTreatOrder);

  virtual void setDownAlgo();
};

class GreedyHeuristic : public Alg4PrimalHeuristicOfNode
{
  bool _thereArePureMasterVariables;

protected:
  virtual void runBody(int & globalTreatOrder);

  bool prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder);

public:
  GreedyHeuristic(Problem * probPtr, MasterCommons4PrimalHeuristic & masterCommons, const bool pureMasterVars);

  virtual ~GreedyHeuristic();
};


#endif /* BCGREEDYHEURISTIC_HPP_ */
