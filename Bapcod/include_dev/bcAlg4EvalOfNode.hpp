/**
 *
 * This file bcAlg4EvalOfNode.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCEVALUATIONALGORITHMC_HPP
#define BCEVALUATIONALGORITHMC_HPP

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4Master.hpp"

struct NodeEvalInfo
{
  int numberOfNodes;
  int treatOrderId;

  NodeEvalInfo() :
  numberOfNodes(0), treatOrderId(-1)
  {
  }

  NodeEvalInfo(const int treatOrder) :
  numberOfNodes(0), treatOrderId(treatOrder)
  {
  }

  virtual ~NodeEvalInfo()
  {
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os << "NodeEvalInfo with number of Nodes = " << numberOfNodes << std::endl;
    return os;
  }

} ;

inline std::ostream& operator<<(std::ostream& os, const NodeEvalInfo & that)
{
  return that.print(os);
}

class Alg4EvalOfNode: public Alg4Master
{
protected:
  MasterCommons4EvalAlg & _masterCommons;
  Double _lastPriorityLevel;
  int _pricingSolverCutsMessageId;
  bool _addCutToMasterFirstCall;

  bool _solIsInteger;
  bool _needBeConqueredIfSolIsInteger;

  std::list<Variable *> _nonStabArtVarPtrList;

  virtual bool isConquered()
  {
    return gapSmallerThanTol(algIncIpDualBound(),algIncIpPrimalBound(), param());
  }

  bool checkIfCurSolIsInteger();

  virtual bool checkIfCurSolIsMasterLpFeasible();

  void recordCutSeparationStats(const std::string & cutClass,
                                const std::multiset<InstanciatedConstr *,
                                                    CutSeparationPriorityComp> & generatedCutConstrSet,
                                int numCutsBefore, double sepTime);

  void addCutsToProblem(char const & flag, const double & cutSepTime,
                        std::map<std::string, int> & numCutGeneratedPerType,
                        std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> & generatedCutConstrSet);

  bool addCutToMaster(char flag, bool masterConverged = false);

  bool runColGenSpRelaxationImprovement(bool masterConverged,
                                        std::set<ColGenSpConf *> & subprobsWithImprovedRelaxation);

  void removeColumnsNotSatisfyingColGenSpRelaxation(const std::set<ColGenSpConf *> & subprobsWithImprovedRelaxation);


  void setRollbackPointSavedStatus(const bool & status);

  const Node * currentNodePtr()
  {
	  return _currentNodePtr;
  }

public:

  Alg4EvalOfNode(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons) ;

  virtual ~Alg4EvalOfNode()
  {
  }

  virtual bool setupAlgo(Node * nodePtr);
  /**
   * Return true if the evalAlg determines infeasibility
   */
  virtual bool eval()
  {
    return false;
  }

  virtual bool checkIfSubProbSolutionsEnumeratedToMIP()
  {
    return false;
  }

  virtual void setOptionNeedBeConqueredIfSolIsInteger(const bool value);

  virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL);

  virtual void setDownAlgo();
} ;

class Alg4EvalByNone : public Alg4EvalOfNode
{

protected:

public:

  Alg4EvalByNone(MasterCommons4EvalAlg & masterCommons) : Alg4EvalOfNode(NULL, masterCommons)
  {
  }

  virtual ~Alg4EvalByNone()
  {
  }

  virtual bool setupAlgo(Node * nodePtr) {return false;}

  virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL)
  {
    return NULL;
  }

  virtual void setDownAlgo() {}

  virtual bool eval() {return false;}

} ;

#endif  /* BCEVALUATIONALGORITHM_HPP */

