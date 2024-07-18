/**
 *
 * This file bcMasterConfData.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMASTERCONFDATA_HPP
#define	BCMASTERCONFDATA_HPP



#endif	/* BCMASTERCONFDATA_HPP */

struct LagrangianStabilizationInfo{
  bool  _dualPriceSmoothingIsActive;
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
  Double _averOuterHalfInterval;
  Double _averInnerHalfInterval;
  Double _averDifference;
  Double _innerHIBase;
  Double _maxInnerArtVarValue;
  Double _maxOuterArtVarValue;
  Double _innerAngle;
  Double _outerAngle;
  int _curvatureModifStage;
  Double _prevPrimalVal;
  Double _firstHalfInterval;
  std::list<Variable * > _stabilizationCandArtVarPtrList;
  std::list<Constraint *> _stabilizationCandConstrList;
};
