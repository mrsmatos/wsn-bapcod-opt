/**
 *
 * This file bcAlg4EvalByColAndCutGen.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef COLANDCUTGENEVALALGC_HPP_
#define COLANDCUTGENEVALALGC_HPP_

#include "bcAlg4EvalBySimplexBasedColGen.hpp"

class Alg4EvalByColAndCutGen : public Alg4EvalBySimplexBasedColGen
{
protected:

  bool _saveInfoForAutoRankOneCutsMemory;

  void reintroduceArtVarInMast();

  void drawPrimalSolutionToDotFile(std::string & filename);

  void printAndRecordActiveCutsStats(bool recordStats);

public:

  Alg4EvalByColAndCutGen(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons);

  virtual bool setupAlgo(Node * nodePtr);

  virtual bool eval();

  virtual void setDownAlgo();

  virtual void setOptionSaveInfoForAutoRankOneCutsMemory(const bool value)
  {
    _saveInfoForAutoRankOneCutsMemory = true;
  }

};

#endif /* COLANDCUTGENEVALALGC_HPP_ */
