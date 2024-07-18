/**
 *
 * This file bcStabilizationInfo.hpp is a part of BaPCod - a generic Branch-And-Price Code.
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


#ifndef bcStabilizationinfo_hpp
#define bcStabilizationinfo_hpp

#include "bcUsefulHeadFil.hpp"
#include "bcProbConfigC.hpp"

const int StabilizationPrintLevel = 2;

/// Stabilization information relative to a constraint or Var
/// varconstr has a pointer to its stab. information
class VarConstrStabInfo
{
  friend class ColGenStabilization;
  Constraint * _constrPtr;
  Variable * _varPtr;
  int _stabilizationParticipationFlag;
  LocalArtificialVar * _posOuterArtVarPtr;
  LocalArtificialVar * _negOuterArtVarPtr;
  LocalArtificialVar * _posInnerArtVarPtr;
  LocalArtificialVar * _negInnerArtVarPtr;
  Double _sepPointVal;
  Double _prevPointVal;
  Double _inPointVal;
  Double _dirOutPointVal;
  Double _subgradient;
  Double _incumbentNormalizedSubgradient;
  Double _incToKelleyNormalizedDirection;
  Double const & dirOutPointVal(Double const & beta, Double const & inOutNorm);
  void normalizeDirSmoothDualVal(Double const & inOutNorm, Double const & dirInOutNorm);
  void recomputeSmoothedValue(Double const & alpha, bool const dirSmoothing = false);
  void subgradientInit();
  void subgradientAdd(Double const & subgr);
  inline Double const & valOrSepPointVal() const;

public:
  VarConstrStabInfo(VarConstr * varConstrPtr);
  VarConstrStabInfo(VarConstrStabInfo const & that, VarConstr * varConstrPtr);

  ~VarConstrStabInfo();
  LocalArtificialVar * posOuterArtVarPtr() const;
  LocalArtificialVar * negOuterArtVarPtr() const;
  LocalArtificialVar * posInnerArtVarPtr() const;
  LocalArtificialVar * negInnerArtVarPtr() const;
  void posOuterArtVarPtr(LocalArtificialVar * ptr);
  void negOuterArtVarPtr(LocalArtificialVar * ptr);
  void posInnerArtVarPtr(LocalArtificialVar * ptr);
  void negInnerArtVarPtr(LocalArtificialVar * ptr);
  int const & stabilizationParticipationFlag() const;
  Double const & sepPointVal() const;
  std::ostream & print(std::ostream & os = std::cout) const;
  void setMembership();
  bool computeCount(ConstVarConstrConstPtr vcPtr);
  const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
};

typedef std::pair<Constraint *, Double> StabCenterPair;

typedef std::list<StabCenterPair> StabCenterList;

/// Stabilization information to be saved to the (col.gen.) evaluation algoritm information
/// (to be passed to a child node)
struct StabilizationInfo
{
  StabCenterList stabilityCenter;
  Double stabilizationBySmoothingAutoAlpha;
  Double stabFunctionCurvature;
  Double maxHalfInterval;
  Double averInnerHalfInterval;
  Double averOuterHalfInterval;
  Double innerAngle;
  Double outerAngle;
  StabilizationInfo();
  StabilizationInfo(const StabilizationInfo & that);
  StabilizationInfo(Problem * probPtr, const ControlParameters & param);
  ~StabilizationInfo();
  std::ostream & print(std::ostream & os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const StabilizationInfo & that)
{
  return that.print(os);
}

/// Stabilization constraint relative to the multi point stabilization
class StabilizationConstraint : public Constraint
{
  StabCenterList _dualSolution;
  Variable * _betaVarPtr;

public:
  int nonActiveCounter;
  int genCounter;
  bool stabCenter;
  StabilizationConstraint(ProbConfig * probConfPtr, const std::string & name,
                          std::list<VarConstrStabInfo *> & stabilizationCandConstrList,
                          Variable * varPtr, const int counter);
  virtual void setMembership();
};

/// Structure for comparison of multi point stabilization constraints
struct StabilizationConstraintSort
{
  /// return true iff constraint 'a' has more priority than 'b' to remain in the master
  bool operator()(StabilizationConstraint * a, StabilizationConstraint  * b) const
  {
    if (a->stabCenter != b->stabCenter)
      return a->stabCenter;

    if (a->val() != b->val())
      return (a->val() > b->val());

    if (a->nonActiveCounter != b->nonActiveCounter)
      return (a->nonActiveCounter < b->nonActiveCounter);

    return (a->genCounter > b->genCounter);
  }
};

inline Double const & VarConstrStabInfo::valOrSepPointVal() const
{
  if (_stabilizationParticipationFlag == 2)
    return _sepPointVal;
  else if ( _varPtr != NULL )
    return _varPtr->val();
  else
    return _constrPtr->val();
}


#endif
