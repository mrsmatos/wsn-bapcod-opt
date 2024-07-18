/**
 *
 * This file bcStabilizationColgen.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcStabilizationColgen.hpp"
#include "bcColGenSpConfC.hpp"


ColGenStabilization::ColGenStabilization(Problem * probPtr,
                                         const std::vector<ColGenSpConf *> & colGenSubProbConfPts,
                                         const ControlParameters & param) :
  _probPtr(probPtr),
  _colGenSubProbConfPts(colGenSubProbConfPts),
  _param(param),
  _incumbentSoValWasUpdated(false),
  _bestIntermediateBoundRecord(Bound::infDualBound(probPtr->objStatus())),
  _currentColGenStage(-1),
  _dualPriceSmoothingIsActive(false),
  _stabFunctionCurvature(0),
  _baseDualPriceSmoothingAlpha(0),
  _mispriceDualPriceSmoothingAlpha(0),
  _nbOfDualPriceSmoothingMisprices(0),
  _incumbentAngle(0),
  _pricingPointAngle(0),
  _incumbentAngleIsDefined(false),
  _pricingPointAngleIsDefined(false),
  _subgradientAtIncumbentIsSet(false),
  _outputHeaderIsPrinted(false),
  _savedOuterInterval(0),
  _savedInnerInterval(0),
  _outerHalfInterval(0),
  _innerHalfInterval(0),
  _dynamicVarsHalfInterval(0),
  _maxHalfInterval(0),
  _averDifference(0),
  _dynamicAverDifference(0),
  _maxInnerArtVarValue(0),
  _maxOuterArtVarValue(0),
  _innerAngle(0),
  _outerAngle(0),
  _dynamicVarsAngle(0),
  _inRootNode(false),
  _multiPointStabConstraints(),
  _multiPointBetaVarPtr(NULL),
  _multiPointStabConstraintsCounter(0),
  _maxMultiPointStabConstraintsNumber(0),
  _stabilizationCandArtVarPtrList(),
  _stabilizationCandConstrList()
{
}

/// creates the variable beta needed for the multi point stabilization
void ColGenStabilization::createMultiPointBetaVarPtr()
{
  _multiPointBetaVarPtr = new Variable(_probPtr->probConfPtr()->modelPtr(),
                                       "betaArtVar", 1, 'F', 'C', 'E', BapcodInfinity,
                                       -BapcodInfinity,'a');
  _multiPointBetaVarPtr->addToProb(_probPtr);
  _probPtr->probVarSet().insert(_multiPointBetaVarPtr, VcIndexStatus::Active);
  _multiPointBetaVarPtr->activate();
  std::list<Variable *> varsToAddToForm;
  varsToAddToForm.push_back(_multiPointBetaVarPtr);
  _probPtr->addVarsSimplyInForm(varsToAddToForm);
}

/// adds a multi point stabilization constraint relative to the current dual solution (SEP point)
void ColGenStabilization::addMultiPointStabConstraint()
{
  if (_multiPointBetaVarPtr == NULL)
    createMultiPointBetaVarPtr();
  StabilizationConstraint * newStabConstrPtr
    = new StabilizationConstraint(_probPtr->probConfPtr(),
                                  std::string("stabConstr") + _multiPointStabConstraintsCounter,
                                  _stabilizationCandConstrList, _multiPointBetaVarPtr,
                                  _multiPointStabConstraintsCounter);
  _multiPointStabConstraintsCounter += 1;
  _multiPointStabConstraints.push_back(newStabConstrPtr);
  newStabConstrPtr->addToProb(_probPtr);
  _probPtr->probConstrSet().insert(newStabConstrPtr, VcIndexStatus::Active);
  newStabConstrPtr->activate();
  std::list<Constraint *> constrsToAddToForm;
  constrsToAddToForm.push_back(newStabConstrPtr);
  _probPtr->addConstrsSimplyInForm(constrsToAddToForm);
}

/// updates the set of multi point stabilization constraints
/// 1) updates multi point stabilization constraints information
/// 2) sorts the constraint according to a priority
/// 3) removes some least priority constraints
void ColGenStabilization::checkMultiPointStabConstraints()
{
  if (_multiPointStabConstraints.empty())
    return;

  std::vector<StabilizationConstraint *>::iterator it;
  for (it = _multiPointStabConstraints.begin(); it != --_multiPointStabConstraints.end(); ++it)
    {
      if ((*it)->val() == 0)
        (*it)->nonActiveCounter += 1;
      else
        (*it)->nonActiveCounter = 0;
      if (_incumbentSoValWasUpdated)
        (*it)->stabCenter = false;
    }
  _multiPointStabConstraints.back()->stabCenter = _incumbentSoValWasUpdated;

  std::stable_sort(_multiPointStabConstraints.begin(),
      _multiPointStabConstraints.end(), StabilizationConstraintSort());

  it = _multiPointStabConstraints.begin();
  int numConstrsToLeave = 0;
  while ((it != _multiPointStabConstraints.end()) && (numConstrsToLeave < _maxMultiPointStabConstraintsNumber)
         && ((*it)->stabCenter || ((*it)->nonActiveCounter <= 2)))
    {
      ++it;
      numConstrsToLeave += 1;
    }

  std::list<Constraint *> constrsToRemoveFromForm;
  for (; it != _multiPointStabConstraints.end(); ++it)
    {
      (*it)->desactivate();
      _probPtr->probConstrSet().insert(*it, VcIndexStatus::Unsuitable);
      constrsToRemoveFromForm.push_back(*it);
    }
  _probPtr->delConstrsSimplyInForm(constrsToRemoveFromForm);

  _multiPointStabConstraints.resize(numConstrsToLeave);
}

/// changes costs and bounds of a stabilization artificial variable based on current penalty function
void ColGenStabilization::setArtCostAndBound(VarConstrStabInfo * constrStabInfoPtr, LocalArtificialVar * artVarPtr)
{
  /// TO DO: currently dual value is multiplied by -1, verify whether this always works
  /// (regardless min or max objective, and the sense of constraints)
  int sign = -1;
  Double mult(constrStabInfoPtr->_constrPtr->costrhs());
  Double newCost;
  Double innerHalfInterval(_innerHalfInterval);
  Double innerAngle(_innerAngle);

  if (mult == 0)
    mult = 1;
  if (constrStabInfoPtr->_stabilizationParticipationFlag)
  {
    switch (artVarPtr->localClassId())
    {
      case LocalArtificialVar::PosOuterId:
        newCost = (sign * constrStabInfoPtr->_constrPtr->incumbentVal() + _outerHalfInterval)
        * mult;
        artVarPtr->resetCurCostByValue(newCost);
        artVarPtr->curUb(_outerAngle);
        break;
      case LocalArtificialVar::NegOuterId:
        newCost = -(sign * constrStabInfoPtr->_constrPtr->incumbentVal() - _outerHalfInterval)
        * mult;
        artVarPtr->resetCurCostByValue(newCost);
        artVarPtr->curUb(_outerAngle);
        break;
      case LocalArtificialVar::PosInnerId:
        if (_param.colGenProximalStabilizationRule().status() != ColGenProximalStabilizationMode::multiPointMode)
          artVarPtr->resetCurCostByValue((sign * constrStabInfoPtr->_constrPtr->incumbentVal() + innerHalfInterval)
                                         * mult);
        else
          artVarPtr->resetCurCostByValue(0);
        artVarPtr->curUb(_innerAngle);
        break;
      case LocalArtificialVar::NegInnerId:
        if (_param.colGenProximalStabilizationRule().status() != ColGenProximalStabilizationMode::multiPointMode)
          artVarPtr->resetCurCostByValue(-(sign * constrStabInfoPtr->_constrPtr->incumbentVal() - innerHalfInterval)
                                         * mult);
        else
          artVarPtr->resetCurCostByValue(0);
        artVarPtr->curUb(innerAngle);
        break;
      default:
        break;
    };


    if (printL(StabilizationPrintLevel))
      std::cout << "name is " << constrStabInfoPtr->_constrPtr->name() << std::setprecision(12)
                << ", incumbVal = " << constrStabInfoPtr->_constrPtr->incumbentVal().val()
                << ", halfInterval = " << _innerHalfInterval.val() << ", curCost = "
                << artVarPtr->curCost().val() << std::endl;
  }
  else
  {
    artVarPtr->resetCurCostByValue(0);
    artVarPtr->curUb(0);
  }
}

/// called from the stabilization deactivation procedure
void ColGenStabilization::reset()
{
  _dualPriceSmoothingIsActive = false;
  _baseDualPriceSmoothingAlpha = 0;
  _mispriceDualPriceSmoothingAlpha = 0;
  _nbOfDualPriceSmoothingMisprices = 0;
  _incumbentAngleIsDefined = false;
  _pricingPointAngleIsDefined = false;
  _subgradientAtIncumbentIsSet = false;
  _outputHeaderIsPrinted = false;
  return;
}

/// changes costs and bounds of stabilization artificial variables based on current penalty function
void ColGenStabilization::setStabArtVarsCostsAndBounds()
{
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end();
      ++stabIt)
    {
      if ((*stabIt)->_posOuterArtVarPtr != NULL)
        setArtCostAndBound(*stabIt, (*stabIt)->_posOuterArtVarPtr);
      if ((*stabIt)->_negOuterArtVarPtr != NULL)
        setArtCostAndBound(*stabIt, (*stabIt)->_negOuterArtVarPtr);
      if ((*stabIt)->_posInnerArtVarPtr != NULL)
        setArtCostAndBound(*stabIt, (*stabIt)->_posInnerArtVarPtr);
      if ((*stabIt)->_negInnerArtVarPtr != NULL)
        setArtCostAndBound(*stabIt, (*stabIt)->_negInnerArtVarPtr);
    }
}

/// adds a constraint and the associated stabilization artificial variables to the stabilization lists
void ColGenStabilization::addConstrAndAssociatedArtVarsToStabCandList(VarConstrStabInfo * constrStabInfoPtr)
{
  _stabilizationCandConstrList.push_back(constrStabInfoPtr);
  constrStabInfoPtr->_constrPtr->incrParticipation(4);
  if (printL(7))
    std::cout << "ColGenStabilization::addConstrAndAssociatedArtVarsToStabCandList() participation of constr "
              << constrStabInfoPtr->_constrPtr->name()
              << " was incremented to " << constrStabInfoPtr->_constrPtr->participation() << std::endl;
  if (_param.colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
    {
      if (constrStabInfoPtr->_posOuterArtVarPtr != NULL)
        _stabilizationCandArtVarPtrList.push_back(constrStabInfoPtr->_posOuterArtVarPtr);
      if (constrStabInfoPtr->_negOuterArtVarPtr != NULL)
        _stabilizationCandArtVarPtrList.push_back(constrStabInfoPtr->_negOuterArtVarPtr);
      if (constrStabInfoPtr->_posInnerArtVarPtr != NULL)
        _stabilizationCandArtVarPtrList.push_back(constrStabInfoPtr->_posInnerArtVarPtr);
      if (constrStabInfoPtr->_negInnerArtVarPtr != NULL)
        _stabilizationCandArtVarPtrList.push_back(constrStabInfoPtr->_negInnerArtVarPtr);
    }

  if (printL(StabilizationPrintLevel))
    {
      std::cout << "add constraint " << constrStabInfoPtr->_constrPtr->name()
          << " to stabilizationCandConstrList with the following artificial variables:"
          << std::endl;
      constrStabInfoPtr->print();
    }
}

/// updates the current piecewise linear penalty function
/// (mode based on the curvature of the quadratic function being approximated)
void ColGenStabilization::updatePenaltyFunctionBasedOnCurvature()
{
  if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
    {
      if (_averDifference < _savedInnerInterval)
        {
          _savedOuterInterval = _savedOuterInterval - (_savedOuterInterval - _savedInnerInterval)
                                                      * _param.StabilFuncCurvatureAdvanceRate();
          _savedInnerInterval = _savedInnerInterval - (_savedInnerInterval - _averDifference)
                                                      * _param.StabilFuncCurvatureAdvanceRate();
        }
      else if (_averDifference < _savedOuterInterval)
        {
          _savedInnerInterval = _savedInnerInterval + (_averDifference - _savedInnerInterval)
                                                      * _param.StabilFuncCurvatureAdvanceRate();
        }
      else
        {
          _savedInnerInterval = _savedInnerInterval + (_savedOuterInterval - _savedInnerInterval)
                                                      * _param.StabilFuncCurvatureAdvanceRate();
          _savedOuterInterval = _savedOuterInterval + (_averDifference - _savedOuterInterval)
                                                      * _param.StabilFuncCurvatureAdvanceRate();
        }
    }
  else /// three-piece
    {
      _savedInnerInterval = _savedInnerInterval
                            + (_averDifference - _savedInnerInterval) * _param.StabilFuncCurvatureAdvanceRate();
    }

  if (_savedInnerInterval < BapcodInfinity)
    {
      _innerHalfInterval = _savedInnerInterval * 0.5;
      _innerAngle = _savedInnerInterval / _stabFunctionCurvature;
      _dynamicVarsHalfInterval = _savedInnerInterval * 0.5;
      _dynamicVarsAngle = 0;
    }
  else
    {
      _innerHalfInterval = 0;
      _innerAngle = 0;
    }

  if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
    {
      if (_savedOuterInterval < BapcodInfinity)
        {
          _outerHalfInterval = _savedOuterInterval * 0.5;
          _outerAngle = _savedOuterInterval / _stabFunctionCurvature;
        }
      else
        {
          _outerHalfInterval = 0;
          _outerAngle = 0;
        }
    }

  if (printL(StabilizationPrintLevel))
    std::cout << "ColGenStabilization::updatePenaltyFunctionBasedOnCurvature() : outerHalfInterval = "
              << _outerHalfInterval << ", outerAngle = " << _outerAngle << ", innerHalfInterval = "
              << _innerHalfInterval << ", innerAngle = " << _innerAngle << std::endl;

  setStabArtVarsCostsAndBounds();
  _probPtr->updateObjCoeffsInForm(_stabilizationCandArtVarPtrList);
  _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
  return;
}

/// computes average deviation of the current dual solution from the stability center
void ColGenStabilization::saveAverageHalfIntervals(int const numbMasterIterations)
{
    if (printL(StabilizationPrintLevel))
        std::cout << "ColGenStabilization::saveAverageHalfIntervals()" << std::endl;

    int numStabConstraints = 0;
    int numDynStabConstraints = 0;
    Double differenceSum = 0;
    Double dynamicDifferenceSum = 0;

    std::list<VarConstrStabInfo *>::iterator stabIt;
    for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    {
        if ((*stabIt)->_stabilizationParticipationFlag)
        {
            (*stabIt)->_inPointVal = (*stabIt)->_constrPtr->incumbentVal();
            double value = (*stabIt)->_constrPtr->incumbentVal() - (*stabIt)->_constrPtr->val();
            if (value < 0)
                value = -value;
            if ((*stabIt)->_constrPtr->flag() == 'd')
            {
                dynamicDifferenceSum += value;
                numDynStabConstraints += 1;
            }
            else
            {
                differenceSum += value;
                numStabConstraints += 1;
            }
            if (printL(StabilizationPrintLevel))
                std::cout << "diff of constr " << (*stabIt)->_constrPtr->name() << " = " << value
                          << " (" << (*stabIt)->_constrPtr->incumbentVal() << ")" << std::endl;
        }
    }

    if (numStabConstraints > 0)
        _averDifference = differenceSum / (double) numStabConstraints;
    else
        _averDifference = 0;

    if (numDynStabConstraints > 0)
        _dynamicAverDifference = dynamicDifferenceSum / (double) numDynStabConstraints;
    else
        _dynamicAverDifference = 0;


    if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
        setStabLocalArtVarMaxValues();

    if (printL(StabilizationPrintLevel))
    {
        std::cout << "averDifference = " << _averDifference << std::endl;
        std::cout << "dynamicAverDifference = " << _dynamicAverDifference << std::endl;
    }
}

/// saves normalized IN-OUT direction
void ColGenStabilization::saveNormalizedIncumbentToKelleyDirection()
{
  if (printL(StabilizationPrintLevel))
    std::cout << "ColGenStabilization::saveNormalizedIncumbentToKelleyDirection()" << std::endl;

  Double incKelleyNorm = 0;

  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        Double value = (*stabIt)->_constrPtr->incumbentVal() - (*stabIt)->_constrPtr->val();
        incKelleyNorm += value * value;
      }
  incKelleyNorm = sqrt((double) incKelleyNorm);

  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    {
      if ((*stabIt)->_stabilizationParticipationFlag)
        {
          (*stabIt)->_incToKelleyNormalizedDirection = ((*stabIt)->_constrPtr->incumbentVal()
              - (*stabIt)->_constrPtr->val()) / incKelleyNorm;
          if (printL(StabilizationPrintLevel))
            std::cout << "Constraint " << (*stabIt)->_constrPtr->name() << ": incumbentVal = "
                      << (*stabIt)->_constrPtr->incumbentVal() << ", val = "
                      << (*stabIt)->_constrPtr->val() << std::endl;
        }
    }
}

/// calculates the angle between the subgradient at the dual incumbent
/// and the IN-OUT direction (needed for automatic directional smoothing)
void ColGenStabilization::calculateAngleAtIncumbent()
{
  if (printL(StabilizationPrintLevel))
    std::cout << "ColGenStabilization::calculateAngleAtIncumbent()" << std::endl;

  _incumbentAngle = 0;
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        _incumbentAngle += (*stabIt)->_incumbentNormalizedSubgradient * (*stabIt)->_incToKelleyNormalizedDirection;
        if (printL(StabilizationPrintLevel))
          std::cout << "Constraint " << (*stabIt)->_constrPtr->name() << ": incNormSubgrad = "
                    << (*stabIt)->_incumbentNormalizedSubgradient << ", incToKelNormDir = "
                    << (*stabIt)->_incToKelleyNormalizedDirection << ", incumbentAngle = "
                    << _incumbentAngle << std::endl;
      }
  _incumbentAngleIsDefined = true;
  return;
}

/// computes "directional" OUT point by "pushing" the real OUT point towards the
/// subgradient of the incumbent dual solution
/// (used only if directional stabilization is applied)
bool ColGenStabilization::computeDirectionalOutPointValues()
{
  if (_param.colGenDualPriceSmoothingBetaFactor() <= 0)
    return false;

  Double inOutNorm = 0;
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        Double value = (*stabIt)->_inPointVal - (*stabIt)->_constrPtr->val();
        inOutNorm += value * value;
      }
  inOutNorm = sqrt((double) inOutNorm);
  if (printL(StabilizationPrintLevel))
    std::cout << "MasterConf::computeDirectionalOutPointValues(): inOutNorm = " << inOutNorm << std::endl;

  if (inOutNorm.isZero())
    return false;

  Double curBetaFactor = 0;
  if (_param.colGenDualPriceSmoothingBetaFactor() == 1) /// dynamic beta factor rule
    {
      if (_subgradientAtIncumbentIsSet)
        calculateAngleAtIncumbent();
      if (_incumbentAngleIsDefined)
        curBetaFactor = _incumbentAngle;
      if (printL(StabilizationPrintLevel))
        std::cout << "_incumbentAngle = " << _incumbentAngle << std::endl;
    }
  else if (_param.colGenDualPriceSmoothingBetaFactor() > 0)
    curBetaFactor = _param.colGenDualPriceSmoothingBetaFactor();

  Double dirInOutNorm = 0;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        Double value = (*stabIt)->_inPointVal - (*stabIt)->dirOutPointVal(curBetaFactor, inOutNorm);
        dirInOutNorm += value * value;
      }

  dirInOutNorm = sqrt((double) dirInOutNorm);
  if (printL(StabilizationPrintLevel))
    std::cout << "MasterConf::computeDirectionalOutPointValues(): dirInOutNorm = " << dirInOutNorm << std::endl;

  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end();
      ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        (*stabIt)->normalizeDirSmoothDualVal(inOutNorm, dirInOutNorm);
      }

  return true;
}

/// updates the subradient during its calculation for a pure master variable
void ColGenStabilization::getSubgradientVectContrib(Variable * varPtr)
{
  // - it assumes that the subgradient has been initialized with RHSs

  ConstVarConstrPtr2Double::const_iterator itm;
  std::map<Constraint *, Double>::iterator mapIt;
  for (itm = varPtr->member2coefMap().begin(); itm != varPtr->member2coefMap().end(); ++itm)
    if (itm->first->inCurProb())
      {
        Constraint * constrPtr = dynamic_cast<Constraint *>(itm->first);
        if (constrPtr == NULL)
          continue;

        if (printL(7))
          std::cout << "maxcol[" << varPtr->name() << "] in constr[" << constrPtr->name()
                    << "] of sense " << constrPtr->sense() << "    has coef " << itm->second << std::endl;

        if ((constrPtr->stabInfoPtr() != NULL) && constrPtr->stabInfoPtr()->_stabilizationParticipationFlag)
          constrPtr->stabInfoPtr()->subgradientAdd(-itm->second * varPtr->mult());
      }
}

/// updates the subradient during its calculation using the column
/// representing the best solution at the SEP or IN point
void ColGenStabilization::getSubgradientVectContrib(MastColumn * colPtr)
{
  Double mult(0);

  // calculate the number of copies of the current subproblem
  ConstVarConstrPtr2Double::const_iterator itm;
  for (itm = colPtr->member2coefMap().begin(); itm != colPtr->member2coefMap().end(); ++itm)
    {
      InstMastConvexityConstr * convConstrPtr = dynamic_cast<InstMastConvexityConstr *>(itm->first);
      if (convConstrPtr != NULL)
        {
          if (convConstrPtr->curRhs() > mult)
            mult = convConstrPtr->curRhs();
        }
    }

  // update the constraint violation for each constraint
  // - it assumes that the subgradient has been initialized with RHSs
  std::map<Constraint *, Double>::iterator mapIt;
  for (itm = colPtr->member2coefMap().begin(); itm != colPtr->member2coefMap().end(); ++itm)
    if (itm->first->inCurProb())
      {
        Constraint * constrPtr = dynamic_cast<Constraint *>(itm->first);
        if (constrPtr == NULL)
          continue;

        if (printL(7))
          std::cout << "maxcol[" << colPtr->name() << "] in constr[" << constrPtr->name()
              << "] of sense " << constrPtr->sense() << "    has coef " << itm->second << std::endl;

        if ((constrPtr->stabInfoPtr() != NULL) && constrPtr->stabInfoPtr()->stabilizationParticipationFlag())
          constrPtr->stabInfoPtr()->subgradientAdd(-itm->second * mult);
      }
}

/// calculates the angle between the gradient direction at the SEP point
/// and the IN-OUT direction
void ColGenStabilization::calculateAngleAtPricingPoint(bool const forAutoSmoothing)
{
  if (printL(StabilizationPrintLevel))
    std::cout << "ColGenStabilization::calculateAngleAtPricingPoint() " << std::endl;

  Double value, inSepNorm = 0;
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        if (forAutoSmoothing)
          value = (*stabIt)->_inPointVal - (*stabIt)->valOrSepPointVal();
        else
          value = (*stabIt)->_constrPtr->incumbentVal() - (*stabIt)->_constrPtr->val();
        inSepNorm += value * value;
      }
  inSepNorm = sqrt((double) inSepNorm);

  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      (*stabIt)->subgradientInit();

  std::vector<ColGenSpConf *>::const_iterator spcPt;
  for (spcPt = _colGenSubProbConfPts.begin(); spcPt != _colGenSubProbConfPts.end(); spcPt++)
    {
      if (*((*spcPt)->upperBoundPtr()) == 0)
        continue;

      MastColumn* colPtr = (*spcPt)->curBestMastColumnPtr();
      if (colPtr != NULL)
        getSubgradientVectContrib(colPtr);
      else
        {
          if (printL(2))
            std::cout << "BaPCod info: cannot access one of the current best subproblem solutions, "
                      << "thus automatic smoothing cannot be applied" << std::endl;
          return;
        }
    }

  /// sub-gradient contribution of pure master variables
  for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active,'s');
       varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); ++varPtrIt)
    getSubgradientVectContrib(*varPtrIt);

  Double subgradNorm = 0;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        Double value = (*stabIt)->_subgradient;
        subgradNorm += value * value;
      }
  subgradNorm = sqrt((double) subgradNorm);

  /// now we calculate angle
  _pricingPointAngle = 0;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        if (forAutoSmoothing)
          value = (*stabIt)->_inPointVal - (*stabIt)->valOrSepPointVal();
        else
          value = (*stabIt)->_constrPtr->incumbentVal() - (*stabIt)->_constrPtr->val();
        _pricingPointAngle += (*stabIt)->_subgradient * value;
        if (printL(StabilizationPrintLevel))
          {
            std::cout << "Constraint " << (*stabIt)->_constrPtr->name() << ": subgrad = "
                      << (*stabIt)->_subgradient << ", inPointVal = " << (*stabIt)->_inPointVal
                      << ", valOrSepPointVal = " << (*stabIt)->valOrSepPointVal() << ", inToSepDir = "
                      << value << ", pricingPointAngle = " << _pricingPointAngle << std::endl;
          }
      }
  _pricingPointAngle = _pricingPointAngle / (subgradNorm * inSepNorm);
  if (printL(StabilizationPrintLevel))
    std::cout << "Pricing point angle = " << _pricingPointAngle << std::endl;

  _pricingPointAngleIsDefined = true;
  return;
}

/// automatic update of the current alpha parameter based on the current angle
/// between the SEP point subgradient and the IN-OUT direction
void ColGenStabilization::dynamicUpdateSmoothingDualValFactor(const bool misprice)
{
  /// commented by Ruslan : if pricing point angle is not defined, then we do not have
  /// the best column, meaning we are in misprice, that angle should be reduced below
  // if (!_pricingPointAngleIsDefined)
  //   return;

  if (misprice || ((double) _pricingPointAngle > 1e-12)) /// we transform to double for strict comparison
    {
      _baseDualPriceSmoothingAlpha -= 0.1;
    }
  else if (((double) _pricingPointAngle < -1e-12) && (_baseDualPriceSmoothingAlpha < 0.999))
  {
      _baseDualPriceSmoothingAlpha += (1 - _baseDualPriceSmoothingAlpha) * 0.1;
  }
  return;
}

/// computes maximum values of stabilization artificial variables
void ColGenStabilization::setStabLocalArtVarMaxValues()
{
  _maxInnerArtVarValue = 0;
  _maxOuterArtVarValue = 0;
  Variable * varPtr;
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    {
      varPtr = (*stabIt)->_posOuterArtVarPtr;
      if ((varPtr != NULL) && (varPtr->val() > _maxOuterArtVarValue))
        _maxOuterArtVarValue = varPtr->val();
      varPtr = (*stabIt)->negOuterArtVarPtr();
      if ((varPtr != NULL) && (varPtr->val() > _maxOuterArtVarValue))
        _maxOuterArtVarValue = varPtr->val();
      varPtr = (*stabIt)->posInnerArtVarPtr();
      if ((varPtr != NULL) && (varPtr->val() > _maxInnerArtVarValue))
        _maxInnerArtVarValue = varPtr->val();
      varPtr = (*stabIt)->negInnerArtVarPtr();
      if ((varPtr != NULL) && (varPtr->val() > _maxInnerArtVarValue))
        _maxInnerArtVarValue = varPtr->val();
    }
}

/// prints detailed information about current state of stabilization,
/// replaces the standard output of column generation
/// one needs to add '-f <name of file>' in the command line to activate this
void ColGenStabilization::printDetailedStabilizationInformation(std::ostream & os, const int nbCgIterations,
                                                                const long cpuTime, const Double & objVal)
{
  if (!_outputHeaderIsPrinted)
    {
      os << "dbUpd\t" << "iter\t" << "time\t" << "dualB\t" << "primalB\t";
      if (!_stabilizationCandArtVarPtrList.empty())
        {
          if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
          {
            if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
              os << "Curvture\t" << "basInterv\t" << "avInnInt\t" << "avOutInt\t" << "intRatio\t" << "inMxVal\t"
                 << "outMxVal\t";
            else
              os << "Curvture\t" << "basInterv\t" << "avInterv\t" << "intRatio\t" << "maxValue\t";
          }
          if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
          {
            if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
              os << "DeltaIn\t" << "DeltaOut\t" << "avInterv\t" << "boundIn\t" << "boundOut\t" << "inMaxVal\t"
                 << "outMaxVal\t";
            else
              os << "Delta\t" << "avInterv\t" << "bound\t" << "maxValue\t";
          }
          if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::multiPointMode)
            {
              os << "bound\t" << "maxValue\t" << "cnstrNum\t";
            }
        }
      if (_param.colGenDualPriceSmoothingAlphaFactor() > 0)
        os << "alpha\t" << "misprice\t";
      os << "incAngl\t" << "pricAngl" << std::endl;
      _outputHeaderIsPrinted = true;
    }

  std::stringstream ss;
  if (_incumbentSoValWasUpdated)
    ss << "1";
  else
    ss << "0";

  ss << std::setprecision(8) << "\t" << nbCgIterations << "\t" << cpuTime / (double) 100
     << "\t" << _bestIntermediateBoundRecord << "\t" << objVal;

  if (!_stabilizationCandArtVarPtrList.empty())
    {
      if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
        {
          setStabLocalArtVarMaxValues();
          if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
            ss << std::setprecision(8) << "\t" << _stabFunctionCurvature << "\t" << _averDifference << "\t"
               << _innerHalfInterval << '\t' << _outerHalfInterval << '\t'
               << _innerHalfInterval / _averDifference << '\t' << _maxInnerArtVarValue << '\t'
               << _maxOuterArtVarValue;
          else
            ss << std::setprecision(8) << "\t" << _stabFunctionCurvature << "\t" << _averDifference << '\t'
               << _innerHalfInterval << '\t' << _innerHalfInterval / _averDifference << '\t'
               << _maxInnerArtVarValue;
        }
      if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
        {
          //setStabLocalArtVarMaxValues(); /// done in saveAverageHalfIntervals
          if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::fivePiece)
            ss << std::setprecision(8) << "\t" << _innerHalfInterval << "\t" << _outerHalfInterval << '\t'
               << _averDifference << '\t' << _innerAngle << '\t' << _outerAngle << '\t'
               << _maxInnerArtVarValue << '\t' << _maxOuterArtVarValue;
          else
            ss << std::setprecision(8) << '\t' << _innerHalfInterval << '\t' << _averDifference << "\t" << _innerAngle
               << '\t' << _maxInnerArtVarValue;
        }
      if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::multiPointMode)
        {
          setStabLocalArtVarMaxValues();
          ss << std::setprecision(8) << '\t' << _innerAngle << '\t' << _maxInnerArtVarValue
             << '\t' << _multiPointStabConstraints.size();
        }
    }

  if (_param.colGenDualPriceSmoothingAlphaFactor() > 0)
    {
      if (_dualPriceSmoothingIsActive)
        ss << std::setprecision(8) << "\t" << _mispriceDualPriceSmoothingAlpha;
      else
        ss << std::setprecision(8) << "\t" << "na(" << _baseDualPriceSmoothingAlpha << ")";
      ss << std::setprecision(8) << "\t" << _nbOfDualPriceSmoothingMisprices;
    }
  if (_incumbentAngleIsDefined)
    ss << std::setprecision(8) << "\t" << _incumbentAngle;
  else
    ss << std::setprecision(8) << "\t" << "na";
  if (_pricingPointAngleIsDefined)
    ss << std::setprecision(8) << "\t" << _pricingPointAngle;
  else
    ss << std::setprecision(8) << "\t" << "na";

  std::cout << ss.str() << std::endl;
  return;
}

/// saves normalized subgradient at the IN point
/// (needed for directional smoothing, both fixed parameter value and automatic)
void ColGenStabilization::saveNormalizedSubgradientAtIncumbent()
{
  _subgradientAtIncumbentIsSet = false;

  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      (*stabIt)->subgradientInit();

  std::vector<ColGenSpConf *>::const_iterator spcPt;
  for (spcPt = _colGenSubProbConfPts.begin(); spcPt != _colGenSubProbConfPts.end(); spcPt++)
    {
      if (*((*spcPt)->upperBoundPtr()) == 0)
        continue;

      MastColumn* colPtr = (*spcPt)->incBestMastColumnPtr();
      if (colPtr != NULL)
        getSubgradientVectContrib(colPtr);
      else
        {
          if (printL(2))
            std::cout << "BaPCod info: cannot access one of the current best subproblem solutions, "
                      << "thus directional smoothing cannot be applied" << std::endl;
          return;
        }
    }

  /// sub-gradient contribution of pure master constraints
  for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active,'s');
       varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); ++varPtrIt)
    getSubgradientVectContrib(*varPtrIt);

  Double subgradNorm = 0;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      {
        Double value = (*stabIt)->_subgradient;
        subgradNorm += value * value;
        if (printL(StabilizationPrintLevel))
          {
            std::cout << "Constraint " << (*stabIt)->_constrPtr->name() << ": incumbent subgrad = "
                      << value << std::endl;
          }
      }
  subgradNorm = sqrt((double) subgradNorm);

  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      (*stabIt)->_incumbentNormalizedSubgradient = (*stabIt)->_subgradient / subgradNorm;

  _subgradientAtIncumbentIsSet = true;
  return;
}

/// deactivation of the stabilization
/// deactivates stabilization artificial variables by„ putting their upper bound to zero
void ColGenStabilization::deactivate()
{
  reset();

  if (_stabilizationCandConstrList.empty())
    return;

  /// deactivate penalty function by setting bounds to zero
  _outerAngle = 0;
  _innerAngle = 0;

  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    {
      (*stabIt)->_constrPtr->decrParticipation(4);
      if (printL(7))
        std::cout << "ColGenStabilization::deactivate participation ofConstr "
                  << (*stabIt)->_constrPtr->name()
                  << " was decremented to " << (*stabIt)->_constrPtr->participation() << std::endl;
      (*stabIt)->_stabilizationParticipationFlag = 0;
    }

  setStabArtVarsCostsAndBounds();
  _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);

  _stabilizationCandConstrList.clear();
  _stabilizationCandArtVarPtrList.clear();
  return;
}


/// stabilization setup using the stabilization information
/// (in the col. gen. algorithm setup in the beginning of a node or after cut generation)
void ColGenStabilization::setupStab(StabilizationInfo * stabInfoPtr, const Bound & lagrBound,
                                    const int & currentColGenStage, const int & nodeDepth)
{
  _bestIntermediateBoundRecord = lagrBound;
  _currentColGenStage = currentColGenStage;
  _inRootNode = (nodeDepth == 0);

  if (_param.colGenDualPriceSmoothingAlphaFactor() < 1)
    _baseDualPriceSmoothingAlpha = _param.colGenDualPriceSmoothingAlphaFactor();
  else
    _baseDualPriceSmoothingAlpha = stabInfoPtr->stabilizationBySmoothingAutoAlpha;

  ConstrIndexManager::iterator constrIndIt;
  for (constrIndIt = _probPtr->probConstrSet().begin(VcIndexStatus::Active, 's');
      constrIndIt != _probPtr->probConstrSet().end(VcIndexStatus::Active, 's'); constrIndIt++)
    {
      VarConstrStabInfo * constrStabInfoPtr = (*constrIndIt)->stabInfoPtr();
      if (constrStabInfoPtr != NULL)
        {
          constrStabInfoPtr->_stabilizationParticipationFlag = 0;
          addConstrAndAssociatedArtVarsToStabCandList(constrStabInfoPtr);
        }
    }
  for (constrIndIt = _probPtr->probConstrSet().begin(VcIndexStatus::Active, 'd');
      constrIndIt != _probPtr->probConstrSet().end(VcIndexStatus::Active, 'd'); constrIndIt++)
    {
      VarConstrStabInfo * constrStabInfoPtr = (*constrIndIt)->stabInfoPtr();
      if (constrStabInfoPtr != NULL)
        {
            if ((*constrIndIt)->VarConstr::isTypeOf(VcId::InstMasterBranchingConstrMask))
                constrStabInfoPtr->_stabilizationParticipationFlag = 0;
            else
            {
                /// if this is a cut, we stabilize it immediately
                /// otherwise it can cause a very slow convergence
                /// (as many cuts can be added at once )
                constrStabInfoPtr->_stabilizationParticipationFlag = 1;
                constrStabInfoPtr->_prevPointVal = 0;
                (*constrIndIt)->incumbentVal(0);

                if (printL(StabilizationPrintLevel))
                    std::cout << "ColGenStabilization setup : new cut " << (*constrIndIt)->name()
                              <<  " incumbent val is set to 0" << std::endl;
            }
            addConstrAndAssociatedArtVarsToStabCandList(constrStabInfoPtr);
        }
    }

  StabCenterList::iterator stabIt;
  for (stabIt = stabInfoPtr->stabilityCenter.begin(); stabIt != stabInfoPtr->stabilityCenter.end(); ++stabIt)
    {
      stabIt->first->stabInfoPtr()->_stabilizationParticipationFlag = 1;
      stabIt->first->stabInfoPtr()->_prevPointVal = stabIt->second;
      stabIt->first->incumbentVal(stabIt->second);
    }

  if (_param.colGenStabilizationMaxTreeDepth() <= nodeDepth)
    {
      /// we need to do it after computing _stabilizationCandConstrList,
      /// as penalty artificial variables may remain active
      deactivate();
      return;
    }

  if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::none)
    return;

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::multiPointMode)
    {
      _innerAngle = stabInfoPtr->innerAngle;
      _outerAngle = stabInfoPtr->outerAngle;
      setStabArtVarsCostsAndBounds();
      _probPtr->updateObjCoeffsInForm(_stabilizationCandArtVarPtrList);
      _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
      _maxMultiPointStabConstraintsNumber = _param.StabilFuncKappa();
    }

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
    {
      _maxHalfInterval = stabInfoPtr->maxHalfInterval;
      _innerHalfInterval = stabInfoPtr->averInnerHalfInterval;
      _outerHalfInterval = stabInfoPtr->averOuterHalfInterval;
      _innerAngle = stabInfoPtr->innerAngle;
      _outerAngle = stabInfoPtr->outerAngle;
      setStabArtVarsCostsAndBounds();
      _probPtr->updateObjCoeffsInForm(_stabilizationCandArtVarPtrList);
      _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
    }

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
    {
      _stabFunctionCurvature = stabInfoPtr->stabFunctionCurvature;
      _savedInnerInterval = stabInfoPtr->averInnerHalfInterval;
      _averDifference = _savedInnerInterval;
      _savedOuterInterval = stabInfoPtr->averOuterHalfInterval;
      updatePenaltyFunctionBasedOnCurvature();
    }
}

/// initializes stabilization at every col.gen. iteration after solving the restricted master
/// 1) decides whether stabilization should be applied at this particular moment
/// 2) computes IN and SEP points
void ColGenStabilization::initializationAfterSolvingRestrictedMaster(Double const & optimalityGap,
                                                                     int const numbMasterIterations,
                                                                     const int & currentColGenStage)
{
  if (_stabilizationCandConstrList.empty())
    return;

  _dualPriceSmoothingIsActive = false;
  _incumbentSoValWasUpdated = false;
  _incumbentAngleIsDefined = false;
  _pricingPointAngleIsDefined = false;

  if (_param.colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
    saveAverageHalfIntervals(numbMasterIterations);

  /// this is needed for penalty functions
  _nbOfDualPriceSmoothingMisprices = 0;
  saveNormalizedIncumbentToKelleyDirection();

  /// switching smoothing off on the first iteration works better
  /// I do not completely understand why: to review it later
  if (_inRootNode && (numbMasterIterations == 0))
    return;

  /// smoothing is not used if the current phase number is less than the minimum permitted by the parameter
  if (currentColGenStage < _param.StabilizationMinPhaseOfStage())
    {
      /// _stabilizationParticipationFlag may remain equal to 2 after the previous phase
      /// (for example, in the case when col. gen. is interrupted by bound)
      std::list<VarConstrStabInfo *>::iterator stabIt;
      for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
        if ((*stabIt)->_stabilizationParticipationFlag)
          (*stabIt)->_stabilizationParticipationFlag = 1;
      return;
    }

  _mispriceDualPriceSmoothingAlpha = _baseDualPriceSmoothingAlpha;

    /// we use the Wentges rule to update the stability center
    for (std::list<VarConstrStabInfo *>::iterator stabIt = _stabilizationCandConstrList.begin();
         stabIt != _stabilizationCandConstrList.end(); ++stabIt)
        if ((*stabIt)->_stabilizationParticipationFlag)
            (*stabIt)->_inPointVal = (*stabIt)->_constrPtr->incumbentVal();

    /// Neame rule is commented as it is never used in practice
//    for (std::list<VarConstrStabInfo *>::iterator stabIt = _stabilizationCandConstrList.begin();
//         stabIt != _stabilizationCandConstrList.end(); ++stabIt)
//        /// for the moment participatesInStabilization status is set to true only on
//        /// the dual incumbent change, which is valid for the penalty functions and Wentges smoothing
//        /// participation in the Neame smoothing can be activated earlier than that: to review this
//        if ((*stabIt)->_stabilizationParticipationFlag)
//            (*stabIt)->_inPointVal = (*stabIt)->_prevPointVal;

  /// this should be verified after setting inPointVal
  if (((optimalityGap > _param.colGenDualPriceSmoothingMaxGap()) && (_param.colGenDualPriceSmoothingMaxGap() > 0))
      || (optimalityGap < _param.colGenDualPriceSmoothingMinGap()) || (_baseDualPriceSmoothingAlpha <= 0))
    return;

  bool directionalSmoothing = computeDirectionalOutPointValues();
  for (std::list<VarConstrStabInfo *>::iterator stabIt = _stabilizationCandConstrList.begin();
       stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
    {
      (*stabIt)->recomputeSmoothedValue(_mispriceDualPriceSmoothingAlpha, directionalSmoothing);
      (*stabIt)->_prevPointVal = (*stabIt)->_sepPointVal;
    }
  _probPtr->updateInDualSol(); /// val values are changed, we need to update _inDualSol
  _dualPriceSmoothingIsActive = true;
  return;
}

/// detemines whether there are stabilization related artificial variables present in the current master solution
bool ColGenStabilization::stabVarsInSolution()
{
  if (_stabilizationCandArtVarPtrList.empty())
    return (false);

  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (std::list<Variable *>::iterator varIt = _stabilizationCandArtVarPtrList.begin();
       varIt != _stabilizationCandArtVarPtrList.end(); ++varIt)
    {
      if (_probPtr->inPrimalLpSol().count(*varIt) && !(*varIt)->val().isZero())
        return (true);
    }
  return (false);
}

/// returns a flag whether the stabilization is active at this moment
bool ColGenStabilization::isActive()
{
  return (!_stabilizationCandArtVarPtrList.empty() || _dualPriceSmoothingIsActive);
}

/// returns a flag whether the dual price smoothing stabilization is active at the moment
bool ColGenStabilization::solValueSmoothingIsActive()
{
  return (_dualPriceSmoothingIsActive);
}

/// calculates the reduced costs of pure master variables at the SEP point
/// (this is needed for correct calculation of the master dual bound at the SEP point)
void ColGenStabilization::changePureMasterVarsReducedCostUsingSepValues(std::map<VarConstr *, double> & modifiedRedCost)
{
  /// we recalculate the reduced costs of master variables for \pi_sep
  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag == 2)
      {
        Double diff = (*stabIt)->_constrPtr->val() - (*stabIt)->_sepPointVal;
        if (!diff.isZero())
          {
            for (ConstVarConstrPtr2Double::iterator itm = (*stabIt)->_constrPtr->member2coefMap().begin();
                 itm != (*stabIt)->_constrPtr->member2coefMap().end(); ++itm)
              if (itm->first->inCurForm() && (itm->first->isTypeOf(VcId::InstMasterVarMask)))
                {
                  modifiedRedCost[itm->first] -= itm->second * diff;
                  if (printL(6))
                    std::cout << "Var[" << itm->first->name() << "] in const[" << (*stabIt)->_constrPtr->name()
                              << "] of rc_diff[" << diff << "] has coef[" << itm->second
                              << "]  rc= " << modifiedRedCost[itm->first] << std::endl;
                }
          }
      }
}

/// updates the local (to stabilization) information about the stability center
void ColGenStabilization::updateOnLagrBoundChange(Bound const & lagrBound, const int & currentColGenStage,
                                                  const bool & imposeBoundChange)
{
  if (!imposeBoundChange && (_currentColGenStage <= currentColGenStage))
  {
    if (_bestIntermediateBoundRecord >= lagrBound)
      return; /// no improvement in bound
  }
  else
  {
    _currentColGenStage = currentColGenStage;
  }

  _bestIntermediateBoundRecord = lagrBound;

  if (_stabilizationCandConstrList.empty())
    return;

  _incumbentSoValWasUpdated = true;

  std::list<VarConstrStabInfo *>::iterator stabIt;
  for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
  {
    (*stabIt)->_constrPtr->incumbentVal((*stabIt)->valOrSepPointVal());
  }

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _colGenSubProbConfPts.begin();
      spcPt != _colGenSubProbConfPts.end(); spcPt++)
    {
      if (*((*spcPt)->upperBoundPtr()) > 0)
        (*spcPt)->saveBestMastColumnPtrAsIncumbentBest(); /// one of this two is superflous
    }

  /// TO DO: implement saving incumbent mult in a map inside ColGenStabilization
//    for (InstMasterVarPtrSet::const_iterator imvPt = _setOfPureMasterVar.begin();
//        imvPt != _setOfPureMasterVar.end(); imvPt++)
//      if ((*imvPt)->inCurProb())
//        (*imvPt)->saveMultAsIncumbentMult();

  return;
}

/// called each time after the pricing oracle
/// returns true if the mis-price is detected and pricing oracle should be called again
/// 1) updates the parameter alpha in the automatic case
/// 2) detectes if there is a mis-price and updates SEP point in this case
bool ColGenStabilization::updateAfterPricingProblemSolution(int nbAddedNegRedCostCol)
{
  if (_stabilizationCandConstrList.empty())
    return false;

  if ((_param.colGenDualPriceSmoothingAlphaFactor() == 1) && (_nbOfDualPriceSmoothingMisprices == 0))
  {
    calculateAngleAtPricingPoint(true);
    dynamicUpdateSmoothingDualValFactor(nbAddedNegRedCostCol == 0);
  }

  if ((nbAddedNegRedCostCol > 0) || !_dualPriceSmoothingIsActive)
    return false; /// stop misprice sequence

  /// update alpha after a misprice
  if (_param.colGenDualPriceSmoothingAlphaFactor() < 1)
    /// variant for fixed alpha
    _mispriceDualPriceSmoothingAlpha = 1 - (_nbOfDualPriceSmoothingMisprices + 1) * (1 - _baseDualPriceSmoothingAlpha);
  else
    /// variant for dynamic alpha
    _mispriceDualPriceSmoothingAlpha = 1 - (1 - _mispriceDualPriceSmoothingAlpha) * 2;

  _nbOfDualPriceSmoothingMisprices += 1;

  std::list<VarConstrStabInfo *>::iterator stabIt;
  if ((_nbOfDualPriceSmoothingMisprices > _param.colGenDualPriceSmoothingMaxNbOfUpdate())
      || (_mispriceDualPriceSmoothingAlpha <= 0))
  {
    /// return to the Kelley point
    for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
      if ((*stabIt)->_stabilizationParticipationFlag)
        (*stabIt)->_stabilizationParticipationFlag = 1;
    _dualPriceSmoothingIsActive = false;
  }
  else /// update separation point
  {
    for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
      if ((*stabIt)->_stabilizationParticipationFlag)
      {
        (*stabIt)->recomputeSmoothedValue(_mispriceDualPriceSmoothingAlpha);
        if (printL(StabilizationPrintLevel))
          std::cout << "smooth dualSol[" << (*stabIt)->_constrPtr->name() << "] = "
                    << (*stabIt)->_sepPointVal << std::endl;
      }
  }
  _probPtr->updateInDualSol(); /// val values are changed, we need to update _inDualSol

  return true; /// continue misprice sequence
}

void ColGenStabilization::resetOnColGenTermination()
{
  if (_stabilizationCandConstrList.empty())
    return;

  /// return val values to the Kelley point
  for (std::list<VarConstrStabInfo *>::iterator stabIt = _stabilizationCandConstrList.begin();
       stabIt != _stabilizationCandConstrList.end(); ++stabIt)
    if ((*stabIt)->_stabilizationParticipationFlag)
      (*stabIt)->_stabilizationParticipationFlag = 1;

  if (_dualPriceSmoothingIsActive)
    {
      _dualPriceSmoothingIsActive = false;
      _probPtr->updateInDualSol();
    }
}

/// called after each column generation iteration (after SEP point produced no mis-price)
/// 1) updates the set of multi point stabilization constraints
/// 2) desactivates dual price smoothing (temporarily until next iteration)
/// 3) saves the subgradient at the IN point (if the current SEP point becomes the IN one)
/// 4) updates penalty functions
void ColGenStabilization::updateAfterColGenIteration()
{
  if (_stabilizationCandConstrList.empty())
    return;

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::multiPointMode)
    {
      /// order is important!
      addMultiPointStabConstraint();
      checkMultiPointStabConstraints();
    }

  std::list<VarConstrStabInfo *>::iterator stabIt;
  if (_dualPriceSmoothingIsActive)
    {
      /// return val values to the Kelley point
      for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
        if ((*stabIt)->_stabilizationParticipationFlag)
          (*stabIt)->_stabilizationParticipationFlag = 1;
      _dualPriceSmoothingIsActive = false;
      _probPtr->updateInDualSol();
    }

  /// this should be done after column generation iteration,
  /// as smoothed constraints list should be constant over the whole misprice sequence
  if (_incumbentSoValWasUpdated)
    for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
      (*stabIt)->_stabilizationParticipationFlag = 1;

  /// should be done after incumbentDualValIsSet
  if (_incumbentSoValWasUpdated && _param.colGenDualPriceSmoothingBetaFactor() > 0)
      saveNormalizedSubgradientAtIncumbent();

  if (_stabilizationCandArtVarPtrList.empty())
    return;

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
    {
      if ((_stabFunctionCurvature == BapcodInfinity) && (_averDifference > 0))
        {
          /// first setting of curvature, based on the average difference
          /// between the initial dual solution and the current one
          _savedOuterInterval = _savedInnerInterval = _averDifference;
          _stabFunctionCurvature = _averDifference / _param.StabilFuncKappa();
        }
      updatePenaltyFunctionBasedOnCurvature();
    }

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
    {
      if (_innerHalfInterval == 1e12)
        _savedInnerInterval =  -_param.StabilFuncKappa();
      bool updateBounds(false);
      if ((_averDifference > 0) && (_innerHalfInterval == 1e12))
        {
          /// first setting of half intervals, based on the average difference
          /// between the initial dual solution and the current one
          if (_param.colGenStabilizationFunctionType().status() == StabilizationFunctionType::threePiece)
            _maxHalfInterval = _innerHalfInterval = _averDifference * _param.StabilFuncKappa();
          else
            {
              _outerHalfInterval = _averDifference * _param.StabilFuncKappa();
              _maxHalfInterval = _innerHalfInterval = _outerHalfInterval* 0.1;
            }
        }
      else if ((_innerHalfInterval < 1e12) && (_maxInnerArtVarValue == 0))
        {
          _innerHalfInterval /= 2;
        }

      setStabArtVarsCostsAndBounds();
      _probPtr->updateObjCoeffsInForm(_stabilizationCandArtVarPtrList);
      if (updateBounds)
        _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
    }

  return;
}

/// verifies if there are stabilization artificial variables in the current solution of the master
/// if yes, updates the penalty function and asks to repeat the column generation
bool ColGenStabilization::updateOnArtVarsInFinalSolution()
{
  if (!stabVarsInSolution())
    return (false);

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
    {
      _stabFunctionCurvature *= _param.StabilFuncArtVarInSolUpdateFactor();
      updatePenaltyFunctionBasedOnCurvature();
    }
  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
    {
      _innerAngle /= _param.StabilFuncArtVarInSolUpdateFactor();
      _outerAngle /= _param.StabilFuncArtVarInSolUpdateFactor();
      setStabArtVarsCostsAndBounds();
      _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
    }
  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::multiPointMode)
    {
      _innerAngle /= _param.StabilFuncArtVarInSolUpdateFactor();
      setStabArtVarsCostsAndBounds();
      _probPtr->updateBoundsInForm(_stabilizationCandArtVarPtrList);
    }
  return (true);
}

/// saves the stabilization information for the child nodes
/// (or just to restart stabilization "at the same" point after)
StabilizationInfo * ColGenStabilization::recordStabilizationInfo()
{
  StabilizationInfo * stabInfoPtr = new StabilizationInfo();
  if (_baseDualPriceSmoothingAlpha < 0.5)
    stabInfoPtr->stabilizationBySmoothingAutoAlpha = 0.5;
  else
    stabInfoPtr->stabilizationBySmoothingAutoAlpha = _baseDualPriceSmoothingAlpha;
  stabInfoPtr->stabFunctionCurvature = _stabFunctionCurvature;
  stabInfoPtr->innerAngle = _param.StabilFuncInnerAngle();
  stabInfoPtr->outerAngle = _param.StabilFuncOuterAngle();

  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::explicitMode)
    {
      stabInfoPtr->maxHalfInterval = _maxHalfInterval;
      stabInfoPtr->averInnerHalfInterval = _maxHalfInterval * _param.StabilFuncHalfIntervalChildNodeFactor();
      stabInfoPtr->averOuterHalfInterval = stabInfoPtr->averInnerHalfInterval * 10;
    }
  if (_param.colGenProximalStabilizationRule().status() == ColGenProximalStabilizationMode::curvatureMode)
    {
      stabInfoPtr->averInnerHalfInterval = _savedInnerInterval;
      stabInfoPtr->averOuterHalfInterval = _savedOuterInterval;
    }

  stabInfoPtr->stabilityCenter.clear();

  if (!_stabilizationCandConstrList.empty())
    {
      std::list<VarConstrStabInfo *>::iterator stabIt;
      for (stabIt = _stabilizationCandConstrList.begin(); stabIt != _stabilizationCandConstrList.end(); ++stabIt)
        if ((*stabIt)->_stabilizationParticipationFlag
            && ((*stabIt)->_constrPtr->vcIndexStatus() == VcIndexStatus::Active))
        {
          (*stabIt)->_constrPtr->incrParticipation(5);
          if (printL(7))
            std::cout << "ColGenStabilization::recordStabilizationInfo() participation of constr "
                      << (*stabIt)->_constrPtr->name()
                      << " was incremented to " << (*stabIt)->_constrPtr->participation() << std::endl;

          stabInfoPtr->stabilityCenter.push_back(
                  StabCenterPair((*stabIt)->_constrPtr, (*stabIt)->_constrPtr->incumbentVal()));
        }

    }

  return stabInfoPtr;
}

/// called at the set down of the evaluation algorithm
void ColGenStabilization::setDownStab()
{
  deactivate();
}

double ColGenStabilization::curAlphaValue()
{
  if (_dualPriceSmoothingIsActive)
    return _mispriceDualPriceSmoothingAlpha;
  return 0.0;
}
