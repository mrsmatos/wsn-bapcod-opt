/**
 *
 * This file bcStabilizationColgen.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

// bcStabilization.hpp
//


#ifndef BCSTABILIZATIONCOLGEN_HPP__
#define BCSTABILIZATIONCOLGEN_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcProbConfigC.hpp"
#include "bcStabilizationInfo.hpp"

class ColGenStabilization
{
  Problem * _probPtr;
  const std::vector<ColGenSpConf *> & _colGenSubProbConfPts;
  const ControlParameters & _param;
  bool _incumbentSoValWasUpdated;
  Bound _bestIntermediateBoundRecord;
  int _currentColGenStage;
  bool _dualPriceSmoothingIsActive;
  Double _stabFunctionCurvature;
  Double _baseDualPriceSmoothingAlpha;
  Double _mispriceDualPriceSmoothingAlpha;
  int _nbOfDualPriceSmoothingMisprices;
  Double _incumbentAngle;
  Double _pricingPointAngle;
  bool _incumbentAngleIsDefined;
  bool _pricingPointAngleIsDefined;
  bool _subgradientAtIncumbentIsSet;
  bool _outputHeaderIsPrinted;
  Double _savedOuterInterval;
  Double _savedInnerInterval;
  Double _outerHalfInterval;
  Double _innerHalfInterval;
  Double _dynamicVarsHalfInterval;
  Double _maxHalfInterval;
  Double _averDifference;
  Double _dynamicAverDifference;
  Double _maxInnerArtVarValue;
  Double _maxOuterArtVarValue;
  Double _innerAngle;
  Double _outerAngle;
  Double _dynamicVarsAngle;
  bool _inRootNode;
  std::vector<StabilizationConstraint *> _multiPointStabConstraints;
  Variable * _multiPointBetaVarPtr;
  int _multiPointStabConstraintsCounter;
  int _maxMultiPointStabConstraintsNumber;
  std::list<Variable *> _stabilizationCandArtVarPtrList;
  std::list<VarConstrStabInfo *> _stabilizationCandConstrList;
  void createMultiPointBetaVarPtr();
  void addMultiPointStabConstraint();
  void checkMultiPointStabConstraints();
  void setArtCostAndBound(VarConstrStabInfo * constrStabInfoPtr, LocalArtificialVar * artVarPtr);
  void reset();
  void setStabArtVarsCostsAndBounds();
  void addConstrAndAssociatedArtVarsToStabCandList(VarConstrStabInfo * constrStabInfoPtr);
  void updatePenaltyFunctionBasedOnCurvature();
  void saveAverageHalfIntervals(int const numbMasterIterations);
  void saveNormalizedIncumbentToKelleyDirection();
  void calculateAngleAtIncumbent();
  bool computeDirectionalOutPointValues();
  void getSubgradientVectContrib(Variable * varPtr);
  void getSubgradientVectContrib(MastColumn * colPtr);
  void calculateAngleAtPricingPoint(bool const forAutoSmoothing);
  void dynamicUpdateSmoothingDualValFactor(const bool misprice);
  void setStabLocalArtVarMaxValues();
  void saveNormalizedSubgradientAtIncumbent();

public:
  ColGenStabilization(Problem * probPtr, const std::vector<ColGenSpConf *> & colGenSubProbConfPts,
                      const ControlParameters & param);
  void setupStab(StabilizationInfo * stabInfoPtr, const Bound & lagrBound,
                 const int & currentColGenStage, const int & nodeDepth);
  void deactivate();
  void printDetailedStabilizationInformation(std::ostream & os, const int nbCgIterations, const long cpuTime,
                                             const Double & objVal);
  void initializationAfterSolvingRestrictedMaster(Double const & optimalityGap,
                                                  int const numbMasterIterations,
                                                  const int & currentColGenStage);

  bool stabVarsInSolution();
  bool isActive();
  bool solValueSmoothingIsActive();
  void changePureMasterVarsReducedCostUsingSepValues(std::map<VarConstr *, double> & modifiedRedCost);
  void updateOnLagrBoundChange(Bound const & LagrBound, const int & currentColGenStage,
                               const bool & imposeBoundChange = false);
  bool updateAfterPricingProblemSolution(int nbAddedNegRedCostCol);
  void resetOnColGenTermination();
  void updateAfterColGenIteration();
  bool updateOnArtVarsInFinalSolution();
  StabilizationInfo * recordStabilizationInfo();
  void setDownStab();
  double curAlphaValue();
};

#endif
