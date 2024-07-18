/**
 *
 * This file bcAlg4EvalByLp.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCLPEVALALG_HPP_
#define BCLPEVALALG_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcVarConstrC.hpp"
#include "bcProblemC.hpp"
#include "bcAlg4EvalOfNode.hpp"
#include "bcAlg4EvalBySimplexBasedColGen.hpp"
#include "bcMathProgSolverInterfaceC.hpp"

class Alg4EvalByLp : public Alg4EvalOfNode
{
  bool _updateIncDualBound;

protected:

  // We need to keep this information to provide it possibly for the next node using it.
  StabilizationInfo * _keptStabInfoPtr;
  Double _keptLatestRedCostFixingGap;

  int solveRestrictedMastLP();

public:

  Alg4EvalByLp(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons);

  virtual bool setupAlgo(Node * nodePtr);

  virtual bool eval();

  NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL);

  virtual void setDownAlgo();

  void setOptionUpdateIncDualBound(const bool & value);
} ;

#endif /* BCLPEVALALG_HPP_ */
