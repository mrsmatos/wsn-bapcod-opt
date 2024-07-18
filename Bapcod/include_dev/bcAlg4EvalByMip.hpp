/**
 *
 * This file bcAlg4EvalByMip.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMIPEVALALG_HPP_
#define BCMIPEVALALG_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcVarConstrC.hpp"
#include "bcProblemC.hpp"
#include "bcAlg4EvalOfNode.hpp"

//MipNodeEval does not need a specific nodeEvalInfo.

class Alg4EvalByMip : public Alg4EvalOfNode
{
  double _paramMaxTime;
  bool _paramExactSolution;
  bool _paramRepeatIfCoreCutAddedAfterMIPHeur;

protected:

  // Added by Boris
  SolutionStatus _problemSolutionRequiredStatus;

  int solveRestrictedMastIP(); //This name is probably not the most adequate at it is not only for Restricted Mast.

public:

  Alg4EvalByMip(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons) :
          Alg4EvalOfNode(probPtr, masterCommons), _paramMaxTime(1e+12), _paramExactSolution(true),
          _paramRepeatIfCoreCutAddedAfterMIPHeur(false)
  {
  }

  virtual ~Alg4EvalByMip()
  {
  }

  bool setupAlgo(Node * nodePtr);

  void setDownAlgo();

  virtual bool eval();

  inline void setParamMaxTime(const double maxTime) {_paramMaxTime = maxTime;}
  inline void setParamExactSolution(const bool exactSolution) {_paramExactSolution = exactSolution;}
  inline void setParamRepeatIfCoreCutAdded(const bool repeat) { _paramRepeatIfCoreCutAddedAfterMIPHeur = repeat; }

} ;

#endif /* BCMIPEVALALG_HPP_ */
