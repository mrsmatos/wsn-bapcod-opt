/**
 *
 * This file bcFracSolBasedHeuristic.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCFRACSOLBASEDHEURISTIC_HPP__
#define BCFRACSOLBASEDHEURISTIC_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalOfNode.hpp"
#include "bcAlg4PrimalHeurOfNode.hpp"

class FracSolBasedHeuristic;

class Algorithm4FracSolBasedHeuristicEval : public Alg4EvalOfNode
{
  FracSolBasedHeuristic * _heurPtr;
public:

  Algorithm4FracSolBasedHeuristicEval(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons,
                                      FracSolBasedHeuristic * heurPtr) :
    Alg4EvalOfNode(probPtr, masterCommons), _heurPtr(heurPtr)
  {
  }

  virtual ~Algorithm4FracSolBasedHeuristicEval()
  {
  }

  virtual bool eval();
} ;


class FracSolBasedHeuristic : public Alg4PrimalHeuristicOfNode
{
  friend class Algorithm4FracSolBasedHeuristicEval;
  std::list<VariableSolInfo> _partialSol;
  SolutionVarInfoPtrList _masterSol;

protected:
  Algorithm4FracSolBasedHeuristicEval * _fracSolBasedEvalAlgPtr;

  virtual void runBody(int & globalTreatOrder);

  bool prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder);

public:
  FracSolBasedHeuristic(Problem * probPtr, MasterCommons4PrimalHeuristic & masterCommons);

  virtual ~FracSolBasedHeuristic();
};


#endif /* BCGREEDYHEURISTIC_HPP_ */
