/**
 *
 * This file bcAlg4EvalBySimplexBasedColGen.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCCOLGENEVALALG_HPP_
#define BCCOLGENEVALALG_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalOfNode.hpp"
#include "bcAlg4EvalByLagrangianDuality.hpp"

struct StabilizationInfo;

struct ColGenEvalInfo : public NodeEvalInfo
{
  StabilizationInfo * stabilizationInfoPtr;
  LpBasisRecord * masterLpBasisPtr;
  Double latestReducedCostFixingGap;

  ColGenEvalInfo() :
    stabilizationInfoPtr(NULL), masterLpBasisPtr(NULL), latestReducedCostFixingGap(BapcodInfinity)
  {
  }

  ColGenEvalInfo(const ColGenEvalInfo & that);

  ColGenEvalInfo(StabilizationInfo * stabInfoPtr, LpBasisRecord * lpBasisPtr,
                 const Double & latestReducedCostFixingGap_) :
    NodeEvalInfo(), stabilizationInfoPtr(stabInfoPtr), masterLpBasisPtr(lpBasisPtr),
    latestReducedCostFixingGap(latestReducedCostFixingGap_)
  {
    if (printL(5))
      std::cout << "ColGenEvalInfo with " << *masterLpBasisPtr << " is created " << std::endl;
  }

  ~ColGenEvalInfo();

  virtual std::ostream & print(std::ostream& os = std::cout) const;

};

inline std::ostream& operator<<(std::ostream& os, const ColGenEvalInfo & that)
{
  return that.print(os);
}

class Alg4EvalBySimplexBasedColGen : public Alg4EvalByLagrangianDuality
{
  Double _latestReducedCostFixingGap;
  bool _subProbSolutionsEnumeratedToMIP;
  bool _masterConverged; /// whether dual solution is value is close to the primal solution value
                         /// (optimalityGapTolerance is used)
  bool _canRunReducedCostFixing;
  int _cutSepRound;

protected:

  inline const bool & masterConverged()
  {
    return _masterConverged;
  }

  inline const bool & subProbSolutionsEnumeratedToMIP()
  {
    return _subProbSolutionsEnumeratedToMIP;
  }

  bool setPurePhaseI();

  bool unsetPurePhaseI();

  bool solveMastLpPrimPh1A2(int & nbOfPenaltiesUpdates, int & nbCgIterations);

  bool solveMastLpPrimPh2(int & nbCgIterations);

  int solveRestrictedMastLP();

  bool updatePenalties(const Double & factor);

  bool runColGenSpRelaxationLightening(const int callMode);

  bool addEnumColumnsToMaster(const int maxNumOfColumns);
  bool runEnumSolBasedUserHeuristicFunctor();
  bool addEnumColumnsWithSmallestRedCostInMaster(const int maxNumOfColumns);
  bool runReducedCostFixingAndEnumeration(const int enumerationMode = 0, const double falseGap = -1.0);

  void setCutSepRound(int cutSepRound) {_cutSepRound = cutSepRound;}

public:

  Alg4EvalBySimplexBasedColGen(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons);

  virtual ~Alg4EvalBySimplexBasedColGen();

  virtual bool setupAlgo(Node * nodePtr);

  virtual bool eval();

  virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL);

  virtual void setDownAlgo();

  virtual bool checkIfSubProbSolutionsEnumeratedToMIP();
};

#endif /* BCCOLGENEVALALG_HPP_ */
