/**
 *
 * This file bcStabilizationInfo.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

// bcStabilization.cpp
//

#include "bcStabilizationInfo.hpp"
#include "bcColGenSpConfC.hpp"

/// computes the current "directional" OUT point for a constraint
/// ("real" OUT point pushed towards the subgradient direction at the dual incumbent solution)
Double const & VarConstrStabInfo::dirOutPointVal(Double const & beta, Double const & inOutNorm)
{
  if (beta > 0)
    _dirOutPointVal = (-_incumbentNormalizedSubgradient) * inOutNorm * beta
                       + (1 - beta) * (_constrPtr->val() - _inPointVal) + _inPointVal;
  else
    _dirOutPointVal = _constrPtr->val();
  return _dirOutPointVal;
}

/// computes the current "normalized" directional OUT point for a constraint
void VarConstrStabInfo::normalizeDirSmoothDualVal(Double const & inOutNorm, Double const & dirInOutNorm)
{
  _dirOutPointVal = (_dirOutPointVal - _inPointVal) / dirInOutNorm * inOutNorm + _inPointVal;
}

/// calculates the current SEP point for a constraint

void VarConstrStabInfo::recomputeSmoothedValue(Double const & alpha, bool const dirSmoothing)
{
  if (dirSmoothing)
  {
    _sepPointVal = alpha * _inPointVal + (1 - alpha) * _dirOutPointVal;
    if (_constrPtr != NULL)
    {
      if ((_constrPtr->sense() == 'L') && (_sepPointVal < 0))
        _sepPointVal = 0;
      else if ((_constrPtr->sense() == 'G') && (_sepPointVal > 0))
        _sepPointVal = 0;
    }
    else if (_varPtr != NULL)
    {
      if ((_varPtr->sense() == 'P') && (_sepPointVal < 0))
        _sepPointVal = 0;
      else if ((_varPtr->sense() == 'N') && (_sepPointVal > 0))
        _sepPointVal = 0;
    }

  }
  else
  {
    if (_constrPtr != NULL)
    {
      _sepPointVal = alpha * _inPointVal + (1 - alpha) * _constrPtr->val();
    }
    else if (_varPtr != NULL)
    {
      _sepPointVal = alpha * _inPointVal + (1 - alpha) * _varPtr->val();
    }
  }
  _stabilizationParticipationFlag = 2;

  if (printL(StabilizationPrintLevel))
  {
    if (_constrPtr != NULL)
    {
      std::cout << "OUT = kelley val[" << _constrPtr->name() << "] = " << _constrPtr->val() << std::endl;
      if (dirSmoothing)
        std::cout << "DIROUT = kelley val[" << _constrPtr->name() << "] = " << _dirOutPointVal << std::endl;
      std::cout << "IN = val[" << _constrPtr->name() << "] = " << _inPointVal << std::endl;
      std::cout << "SEP = smoothed val[" << _constrPtr->name() << "] = " << sepPointVal() << std::endl;
    }
    else if (_varPtr != NULL)
    {
      std::cout << "IN = val[" << _varPtr->name() << "] = " << _inPointVal << std::endl;
      std::cout << "OUT = kelley val[" << _varPtr->name() << "] = " << _varPtr->val() << std::endl;
      if (dirSmoothing)
        std::cout << "DIROUT = kelley val[" << _varPtr->name() << "] = " << _dirOutPointVal << std::endl;
      std::cout << "SEP = smoothed val[" << _varPtr->name() << "] = " << sepPointVal() << std::endl;
    }
  }

  return;
}

/// initialization for subgradient calculation (at the SEP or IN point)
void VarConstrStabInfo::subgradientInit()
{
  _subgradient = _constrPtr->curRhs();
}

/// subgradient update during its calculation (at the SEP or IN point)
void VarConstrStabInfo::subgradientAdd(Double const & subgr)
{
  _subgradient += subgr;
}

/// constructor for empty stabilization information
VarConstrStabInfo::VarConstrStabInfo(VarConstr * varConstrPtr) :
        _constrPtr(NULL), _varPtr(NULL),
        _stabilizationParticipationFlag(0), _posOuterArtVarPtr(NULL),
        _negOuterArtVarPtr(NULL), _posInnerArtVarPtr(NULL), _negInnerArtVarPtr(NULL),
        _sepPointVal(0), _prevPointVal(0), _inPointVal(0), _dirOutPointVal(0), _subgradient(0),
        _incumbentNormalizedSubgradient(0), _incToKelleyNormalizedDirection(0)
{
  if (varConstrPtr->isTypeOf(VcId::ConstraintMask))
    {
      _constrPtr = static_cast<Constraint *>(varConstrPtr);
    }
  else
    {
      _varPtr = static_cast<Variable *>(varConstrPtr);
    }

}

/// copy constructor stabilization information (called from copy constructor of the constraint)
VarConstrStabInfo::VarConstrStabInfo(VarConstrStabInfo const & that, VarConstr * varConstrPtr) :
  _constrPtr(NULL),
  _varPtr(NULL),
  _stabilizationParticipationFlag(that._stabilizationParticipationFlag),
  _posOuterArtVarPtr(that._posOuterArtVarPtr),
  _negOuterArtVarPtr(that._negOuterArtVarPtr),
  _posInnerArtVarPtr(that._posInnerArtVarPtr),
  _negInnerArtVarPtr(that._negInnerArtVarPtr),
  _sepPointVal(that._sepPointVal),
  _prevPointVal(that._prevPointVal),
  _inPointVal(that._inPointVal),
  _dirOutPointVal(that._dirOutPointVal),
  _subgradient(that._subgradient),
  _incumbentNormalizedSubgradient(that._incumbentNormalizedSubgradient),
  _incToKelleyNormalizedDirection(that._incToKelleyNormalizedDirection)
{
  if (varConstrPtr->isTypeOf(VcId::ConstraintMask))
    {
      _constrPtr =  static_cast<Constraint *>(varConstrPtr);

      if (_posOuterArtVarPtr != NULL)
	    _posOuterArtVarPtr->setConstraintPtr(_constrPtr);
      if (_negOuterArtVarPtr != NULL)
	    _negOuterArtVarPtr->setConstraintPtr(_constrPtr);
      if (_posInnerArtVarPtr != NULL)
	    _posInnerArtVarPtr->setConstraintPtr(_constrPtr);
      if (_negInnerArtVarPtr != NULL)
	    _negInnerArtVarPtr->setConstraintPtr(_constrPtr);
    }
  else
    {
      _varPtr  =  static_cast<Variable *>(varConstrPtr);
    }
}

VarConstrStabInfo::~VarConstrStabInfo()
{
  if ((_constrPtr != NULL) && (_constrPtr->flag() == 'd'))
    {
      /// if associated local artificial variables were not added to the problem
      /// we delete them, otherwise there will be memory leak
      /// (if they were added to the problem, they will be deleted when erased from the
      ///  problem's var. manager)
      if ((_posOuterArtVarPtr != NULL) && (_posOuterArtVarPtr->problemPtr() == NULL))
        delete _posOuterArtVarPtr;
      _posOuterArtVarPtr = NULL;
      if ((_negOuterArtVarPtr != NULL) && (_negOuterArtVarPtr->problemPtr() == NULL))
        delete _negOuterArtVarPtr;
      _negOuterArtVarPtr = NULL;
      if ((_posInnerArtVarPtr != NULL) && (_posInnerArtVarPtr->problemPtr() == NULL))
        delete _posInnerArtVarPtr;
      _posInnerArtVarPtr = NULL;
      if ((_negInnerArtVarPtr != NULL) && (_negInnerArtVarPtr->problemPtr() == NULL))
        delete _negInnerArtVarPtr;
      _negInnerArtVarPtr = NULL;
    }
}

LocalArtificialVar * VarConstrStabInfo::posOuterArtVarPtr() const
{
  return (_posOuterArtVarPtr);
}

LocalArtificialVar * VarConstrStabInfo::negOuterArtVarPtr() const
{
  return (_negOuterArtVarPtr);
}

LocalArtificialVar * VarConstrStabInfo::posInnerArtVarPtr() const
{
  return (_posInnerArtVarPtr);
}

LocalArtificialVar * VarConstrStabInfo::negInnerArtVarPtr() const
{
  return (_negInnerArtVarPtr);
}

void VarConstrStabInfo::posOuterArtVarPtr(LocalArtificialVar * ptr)
{
  _posOuterArtVarPtr = ptr;
}

void VarConstrStabInfo::negOuterArtVarPtr(LocalArtificialVar * ptr)
{
  _negOuterArtVarPtr = ptr;
}

void VarConstrStabInfo::posInnerArtVarPtr(LocalArtificialVar * ptr)
{
  _posInnerArtVarPtr = ptr;
}

void VarConstrStabInfo::negInnerArtVarPtr(LocalArtificialVar * ptr)
{
  _negInnerArtVarPtr = ptr;
}

int const & VarConstrStabInfo::stabilizationParticipationFlag() const
{
  return _stabilizationParticipationFlag;
}

Double const & VarConstrStabInfo::sepPointVal() const
{
  return _sepPointVal;
}

std::ostream & VarConstrStabInfo::print(std::ostream & os) const
{
  if (_posOuterArtVarPtr != NULL)
    os << "  posThetaLocArtVar = " << _posOuterArtVarPtr->name() << std::endl;
  if (_negOuterArtVarPtr != NULL)
    os << "  negThetaLocArtVar = " << _negOuterArtVarPtr->name() << std::endl;
  if (_posInnerArtVarPtr != NULL)
    os << "  posGammaLocArtVar = " << _posInnerArtVarPtr->name() << std::endl;
  if (_negInnerArtVarPtr != NULL)
    os << "  negGammaLocArtVar = " << _negInnerArtVarPtr->name() << std::endl;
  return (os);
}

/// sets membership of stabilization artifical variables
/// (called from setMembership of the stabilized constraint)
void VarConstrStabInfo::setMembership()
{
  if (_posOuterArtVarPtr != NULL)
    _posOuterArtVarPtr->addMember(_constrPtr);
  if (_negOuterArtVarPtr != NULL)
    _negOuterArtVarPtr->addMember(_constrPtr);
  if (_posInnerArtVarPtr != NULL)
    _posInnerArtVarPtr->addMember(_constrPtr);
  if (_negInnerArtVarPtr != NULL)
    _negInnerArtVarPtr->addMember(_constrPtr);
}

/// sets membership of stabilization artifical variables
bool VarConstrStabInfo::computeCount(ConstVarConstrConstPtr vcPtr)
{
  return ((vcPtr == _posOuterArtVarPtr) || (vcPtr == _negOuterArtVarPtr)
          || (vcPtr == _posInnerArtVarPtr) || (vcPtr == _negInnerArtVarPtr));
}

/// sets membership of stabilization artifical variables
const LpCoef VarConstrStabInfo::computeCoef(ConstVarConstrConstPtr vcPtr)
{
  Double coeff(_constrPtr->costrhs());
  if (zeroTest(coeff))
    coeff = 1;
  if (vcPtr == _posOuterArtVarPtr)
    return LpCoef(coeff);
  else if (vcPtr == _negOuterArtVarPtr)
    return LpCoef(-coeff);
  else if (vcPtr == _posInnerArtVarPtr)
    return LpCoef(coeff);
  else if (vcPtr == _negInnerArtVarPtr)
    return LpCoef(-coeff);
  return LpCoef::ZeroCoef;
}

/// copy constructor of stabilization info
StabilizationInfo::StabilizationInfo(const StabilizationInfo & that):
    stabilityCenter(that.stabilityCenter), stabilizationBySmoothingAutoAlpha(that.stabilizationBySmoothingAutoAlpha),
    stabFunctionCurvature(that.stabFunctionCurvature), maxHalfInterval(that.maxHalfInterval),
    averInnerHalfInterval(that.averInnerHalfInterval), averOuterHalfInterval(that.averOuterHalfInterval),
    innerAngle(that.innerAngle), outerAngle(that.outerAngle)
{
  for (StabCenterList::iterator listIt = stabilityCenter.begin();
       listIt != stabilityCenter.end(); ++listIt)
    {
      listIt->first->incrParticipation(3);
      if(printL(7))
        std::cout << "StabilizationInfo::StabilizationInfo() participation of constr "
                  << listIt->first->name()
                  << " was incremented to " << listIt->first->participation() << std::endl;
    }
}

/// constructor of empty stabilization info
StabilizationInfo::StabilizationInfo() :
    stabilityCenter(), stabilizationBySmoothingAutoAlpha(0.5), stabFunctionCurvature(BapcodInfinity),
    maxHalfInterval(BapcodInfinity),
    averInnerHalfInterval(BapcodInfinity), averOuterHalfInterval(BapcodInfinity),
    innerAngle(0), outerAngle(0)
{
}

/// constructor of initial stabilization info (called in the beginning of the root node)
StabilizationInfo::StabilizationInfo(Problem * probPtr, const ControlParameters & param) :
  stabilityCenter(),
  stabilizationBySmoothingAutoAlpha(0.5),
  stabFunctionCurvature(BapcodInfinity),
  maxHalfInterval(BapcodInfinity),
  averInnerHalfInterval(BapcodInfinity),
  averOuterHalfInterval(BapcodInfinity),
  innerAngle(0),
  outerAngle(0)
{
  ConstrIndexManager::iterator constrIndIt;
  for (constrIndIt = probPtr->probConstrSet().begin(VcIndexStatus::Active, 's');
       constrIndIt != probPtr->probConstrSet().end(VcIndexStatus::Active, 's'); constrIndIt++)
    {
      if ((*constrIndIt)->stabInfoPtr() != NULL)
        {
          (*constrIndIt)->incrParticipation(3);
          if (printL(7))
            std::cout << "StabilizationInfo::StabilizationInfo() participation of constr "
                      << (*constrIndIt)->name()
                      << " was incremented to " << (*constrIndIt)->participation() << std::endl;
          stabilityCenter.push_back(StabCenterPair(*constrIndIt, (*constrIndIt)->incumbentVal()));
        }
    }

  if (param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::none)
    return;

  innerAngle = param.StabilFuncInnerAngle();
  outerAngle = param.StabilFuncOuterAngle();
  if (param.StabilFuncKappa() <= 0) // manual setting of curvature or half intervals
    {
      if (param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
        stabFunctionCurvature = param.StabilFuncCurvature();
      if (param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
        {
          maxHalfInterval = averInnerHalfInterval = param.StabilFuncInnerHalfInterval();
          averOuterHalfInterval = param.StabilFuncOuterHalfInterval();
        }
    }

}

StabilizationInfo::~StabilizationInfo()
{
  for (StabCenterList::iterator stabIt = stabilityCenter.begin();
       stabIt != stabilityCenter.end(); ++stabIt)
    {
      stabIt->first->decrParticipation(3);
      if (printL(7))
        std::cout << "StabilizationInfo::~StabilizationInfo participation ofConstr "
                  << stabIt->first->name()
                  << " was decremented to " << stabIt->first->participation() << std::endl;
    }
}

std::ostream & StabilizationInfo::print(std::ostream & os) const
{
  os << "StabilizationInfo:" << std::endl;
  os << "stabilizationBySmoothingAutoAlpha = " << stabilizationBySmoothingAutoAlpha << std::endl;
  os << "stabFunctionCurvature = " << stabFunctionCurvature << std::endl;
  os << "averInnerHalfInterval = " << averInnerHalfInterval << std::endl;
  os << "averOuterHalfInterval = " << averOuterHalfInterval << std::endl;
  os << "innerAngle = " << innerAngle << std::endl;
  os << "outerAngle = " << outerAngle << std::endl;

  StabCenterList::const_iterator stabIt;
  for (stabIt = stabilityCenter.begin(); stabIt != stabilityCenter.end(); ++stabIt)
    os << "Constraint " << stabIt->first->name() << " stab. center value is " << stabIt->second << std::endl;
  return (os);
}

StabilizationConstraint::StabilizationConstraint(ProbConfig * probConfPtr,
    const std::string & name,
    std::list<VarConstrStabInfo *> & stabilizationCandConstrList,
    Variable * varPtr, const int counter) :
    Constraint(probConfPtr->modelPtr(), name, 0, 'L', 'F', 'E', 'd'), _dualSolution(),
               _betaVarPtr(varPtr), nonActiveCounter(0), genCounter(counter), stabCenter(false)
{
  Problem * probPtr = probConfPtr->probPtr();

  if (printL(StabilizationPrintLevel))
    std::cout << "StabilizationConstraint::StabilizationConstraint() : constraint" << name << " is created"
              << std::endl;

  for (ConstrPtrSet::const_iterator cpsIt = probPtr->inDualSol().begin(); cpsIt != probPtr->inDualSol().end(); ++cpsIt)
    if ((*cpsIt)->stabInfoPtr() != NULL)
        _dualSolution.push_back(StabCenterPair(*cpsIt, (*cpsIt)->valOrSepPointVal()));
}

void StabilizationConstraint::setMembership()
{
  bool cumulativeCoef(false);

  includeMember(_betaVarPtr, -1, cumulativeCoef);

  if (printL(StabilizationPrintLevel))
    std::cout << "Included variable " << _betaVarPtr->name() << " with coeff -1 to " << name() << std::endl;

  for (StabCenterList::const_iterator it = _dualSolution.begin(); it != _dualSolution.end(); ++it)
    {
      if (printL(StabilizationPrintLevel))
        std::cout << "Check constraint " << it->first->name() << " with coeff " << it->second << std::endl;
      LocalArtificialVar * posLocArtVarPtr = it->first->stabInfoPtr()->posInnerArtVarPtr();
      if (posLocArtVarPtr != NULL)
        {
          includeMember(posLocArtVarPtr, -it->second, cumulativeCoef);
          if (printL(StabilizationPrintLevel))
            std::cout << "Included variable " << posLocArtVarPtr->name() << " with coeff " << -it->second
                      << " to " << name() << std::endl;
        }
      LocalArtificialVar * negLocArtVarPtr = it->first->stabInfoPtr()->negInnerArtVarPtr();
      if (negLocArtVarPtr != NULL)
        {
          includeMember(negLocArtVarPtr, it->second, cumulativeCoef);
          if (printL(StabilizationPrintLevel))
            std::cout << "Included variable " << negLocArtVarPtr->name() << " with coeff " << it->second
                      << " to " << name() << std::endl;
        }
    }
}
