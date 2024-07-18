/**
 *
 * This file bcModelParameterC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcModelParameterC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastColumnC.hpp"

const int allStatus[] =
{ SolutionStatus::Optimum, SolutionStatus::Infeasible,
  SolutionStatus::Unbounded, SolutionStatus::UnSolved,
  SolutionStatus::PrimalFeasSolFound, SolutionStatus::DualFeasSolFound,
  SolutionStatus::OptimumUnscalInfeas};

const int N = 7;

const std::set<int> _admStatus = std::set<int>(allStatus, allStatus + N);

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                             -- methods of class SolutionMethod
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SolutionMethod::SolutionMethod(const int & stat)
{
  switch (stat)
  {
    case 0:
      _status = none;
      return;
    case 1:
      _status = lpSolver;
      return;
    case 2:
      _status = mipSolver;
      return;
    case 3:
      _status = customSolver;
      return;
    case 4:
      _status = custom2mipSolver;
      return;
    default:
      _status = undefined;
      break;
  }
  return;
}
const SolutionMethod::SMenum & SolutionMethod::status() const
{
  return (_status);
}

int SolutionMethod::getStatusAsInteger() const
{
  switch (_status)
  {
    case none:
      return (0);
    case lpSolver:
      return (1);
    case mipSolver:
      return (2);
    case customSolver:
      return (3);
    case custom2mipSolver:
      return (4);
    default:
      break;
  }
  return (-1);
}

std::ostream& SolutionMethod::print(std::ostream& os) const
{
  switch (_status)
  {
    case none:
      return (os << "none");
    case lpSolver:
      return (os << "lpSolver");
    case mipSolver:
      return (os << "mipSolver");
    case customSolver:
      return (os << "customSolver");
    case custom2mipSolver:
      return (os << "custom2mipSolver");
    default:
      break;
  }
  return (os << "undefined");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//            -- methods of class StabilizationFunctionType
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

StabilizationFunctionType::StabilizationFunctionType(const int & stat)
{
  set(stat);
}

bool StabilizationFunctionType::set(int method)
{
  switch (method)
  {
    case 0:
      _status = none;
      break;
    case 1:
      _status = boxStep;
      break;
    case 2:
      _status = threePiece;
      break;
    case 3:
      _status = fivePiece;
      break;
    case 4:
      _status = bundle;
      break;
    default:
      _status = undefined;
      return false;
  }
  return true;
}

const StabilizationFunctionType::StabilizationFunctionTypeMenum & StabilizationFunctionType::status() const
{
  return (_status);
}

int StabilizationFunctionType::getStatusAsInteger() const
{
  switch (_status)
  {
    case none:
      return (0);
    case boxStep:
      return (1);
    case threePiece:
      return (2);
    case fivePiece:
      return (3);
    case bundle:
      return (4);
    default:
      break;
  }
  return (-1);
}

std::ostream& StabilizationFunctionType::print(std::ostream& os) const
{
  switch (_status)
  {
    case none:
      return (os << "none");
    case boxStep:
      return (os << "box step");
    case threePiece:
      return (os << "3-piece");
    case fivePiece:
      return (os << "5-piece");
    case bundle:
      return (os << "Bundle");
    default:
      break;
  }
  return (os << "undefined");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//            -- methods of class ColGenProximalStabilizationMode
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

ColGenProximalStabilizationMode::ColGenProximalStabilizationMode(const int & stat)
{
  set(stat);
}

bool ColGenProximalStabilizationMode::set(int method)
{
  switch (method)
  {
    case 0:
      _status = curvatureMode;
      break;
    case 1:
      _status = explicitMode;
      break;
    case 2:
      _status = multiPointMode;
      break;
    default:
      _status = undefined;
      return false;
  }
  return true;
}

const ColGenProximalStabilizationMode::ColGenProximalStabilizationModeMenum & ColGenProximalStabilizationMode::status() const
{
  return (_status);
}

int ColGenProximalStabilizationMode::getStatusAsInteger() const
{
  switch (_status)
  {
    case curvatureMode:
      return (0);
    case explicitMode:
      return (1);
    case multiPointMode:
      return (2);
    default:
      break;
  }
  return (-1);
}

std::ostream& ColGenProximalStabilizationMode::print(std::ostream& os) const
{
  switch (_status)
  {
    case curvatureMode:
      return (os << "curvatureMode");
    case explicitMode:
      return (os << "explicitMode");
    case multiPointMode:
      return (os << "multiPointMode");
    default:
      break;
  }
  return (os << "undefined");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                  -- methods of class SmoothingMode
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SmoothingMode::SmoothingMode(const int & stat)
{
  set(stat);
}

bool SmoothingMode::set(int method)
{
  switch (method)
  {
    case 0:
      _status = none;
      break;
    case 1:
      _status = Wentges;
      break;
    case 2:
      _status = Neame;
      break;
    case 3:
      _status = ActiveColDs;
      break;
    default:
      _status = undefined;
      return false;
  }
  return true;
}

const SmoothingMode::SmoothingMenum & SmoothingMode::status() const
{
  return (_status);
}

int SmoothingMode::getStatusAsInteger() const
{
  switch (_status)
  {
    case none:
      return (0);
    case Wentges:
      return (1);
    case Neame:
      return (2);
    case ActiveColDs:
      return (3);
    default:
      break;
  }
  return (-1);
}

std::ostream& SmoothingMode::print(std::ostream& os) const
{
  switch (_status)
  {
    case none:
      return (os << "none");
    case Wentges:
      return (os << "Wentges");
    case Neame:
      return (os << "Neame");
    case ActiveColDs:
      return (os << "ActiveColDs");
    default:
      break;
  }
  return (os << "undefined");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                             -- methods of class SolutionStatus
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SolutionStatus::SolutionStatus()
{
}

SolutionStatus::SolutionStatus(const int & stat)
{
  insert(stat);
  return;
}
SolutionStatus::SolutionStatus(const int & stat1, const int & stat2)
{
  insert(stat1);
  insert(stat2);
  return;
}
SolutionStatus::SolutionStatus(const std::vector<int> & statVect)
{
  insert(statVect.begin(), statVect.end());
  return;
}
std::pair<std::set<int>::iterator, bool> SolutionStatus::insert(
    const int & stat)
{
  if (_admStatus.count(stat))
    return (this->std::set<int>::insert(stat));
  else
    return (this->std::set<int>::insert(Undefined));
}
SolutionStatus::SolutionStatus(const SolutionStatus & right)
{
  clear();
  insert(right.begin(), right.end());
  return;
}
SolutionStatus::SolutionStatus(const int * first, const int * last)
{
  insert(first, last);
  return;
}

SolutionStatus & SolutionStatus::operator=(const SolutionStatus & right)
{
  clear();
  insert(right.begin(), right.end());
  return (*this);
}

SolutionStatus & SolutionStatus::operator=(const int & stat)
{
  clear();
  insert(stat);
  return (*this);
}

bool SolutionStatus::intersects(const SolutionStatus & right) const
{
  for (std::set<int>::const_iterator it = begin(); it != end(); ++it)
    if (right.count(*it))
      return (true);
  return (false);
}

std::string SolutionStatus::stat2string(const int & stat) const
{
  switch (stat)
  {
    case Optimum:
      return ("Optimum");
    case Infeasible:
      return ("Infeasible");
    case Unbounded:
      return ("Unbounded");
    case UnSolved:
      return ("UnSolved");
    case PrimalFeasSolFound:
      return ("PrimalFeasSolFound");
    case DualFeasSolFound:
      return ("DualFeasSolFound");
    case OptimumUnscalInfeas:
      return ("OptimumUnscalInfeas");
    default:
      return ("Undefined");
  }
}

std::ostream& SolutionStatus::print(std::ostream& os) const
{
  if (empty())
    os << "SolutionStatus is empty " << std::endl;
  else
  {
    os << "SolutionStatus includes ";
    for (std::set<int>::const_iterator it = begin(); it != end(); ++it)
      os << "    " << stat2string(*it) << " ,";
    os << std::endl;
  }
  return (os);
}

void SolutionStatus::validate_one(std::string token)
{
  int element = boost::lexical_cast<int>(token);

  /// Only insert element contained into a given interval
  if (element < N && element > -2)
    insert(element);
  else
    throw boost::program_options::invalid_option_value("");
}

MasterInitMode::MasterInitMode(const int & stat)
{
  set(stat);
}

bool MasterInitMode::set(const int & stat)
{
  switch (stat)
  {
    case 0:
      _status = noArtCol;
      break;
    case 1:
      _status = globalArtCol;
      break;
    case 2:
      _status = subProbArtCol;
      break;
    case 3:
      _status = localArtCol;
      break;
    case 4:
      _status = incSolCol;
      break;
    case 5:
      _status = incSolColAndGac;
      break;
    case 6:
      _status = incSolColAndLac;
      break;
    case 7:
      _status = localAndGlobAc;
      break;
    default:
      _status = defaultInit;
      return false;
  }
  return true;
}

std::ostream& MasterInitMode::print(std::ostream& os) const
{
  switch (_status)
  {
    case noArtCol:
      return (os << "noArtCol");
    case globalArtCol:
      return (os << "globalArtCol");
    case subProbArtCol:
      return (os << "subProbArtCol");
    case localArtCol:
      return (os << "localArtCol");
    case incSolCol:
      return (os << "incSolCol");
    case incSolColAndGac:
      return (os << "incSolColAndGac");
    case incSolColAndLac:
      return (os << "incSolColAndLac");
    case localAndGlobAc:
      return (os << "localAndGlobAc");
    case defaultInit:
    default:
      break;
  }
  return (os << "defaultInit");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//                             -- methods of class SelectionStrategy
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SelectionStrategy::SelectionStrategy(const int & stat)
{
  set(stat);
}

SelectionStrategy::SelectionStrategy(const PriorityEnum & stat) :
    _selectedRule(stat)
{
}

bool SelectionStrategy::set(const int & stat)
{
  switch (stat)
  {
    case 0:
      _selectedRule = NotConsideredForSelection;
      break;
    case 1:
      _selectedRule = FirstFound;
      break;
    case 2:
      _selectedRule = HighestPriority;
      break;
    case 3:
      _selectedRule = MostFractional;
      break;
    case 4:
      _selectedRule = LeastFractional;
      break;
    case 5:
      _selectedRule = FracWeightedPriority;
      break;
    case 6:
      _selectedRule = Closest2RoundUp;
      break;
    case 7:
      _selectedRule = Closest2RoundDown;
      break;
    case 8:
      _selectedRule = GuidedSearch;
      break;
    case 9:
      _selectedRule = LeastCost;
      break;
    case 10:
      _selectedRule = LeastReducedCost;
      break;
    case 11:
      _selectedRule = LeastGreedyCost;
      break;
    case 12:
      _selectedRule = LeastSteepestEdgeCost;
      break;
    case 13:
      _selectedRule = LeastPseudoCost;
      break;
    case 15:
      _selectedRule = LeastInfeasibility;
      break;
    case 16:
      _selectedRule = MostViolated;
      break;
    case 17:
      _selectedRule = NotConsideredForIntegralityCheck;
      break;
    default:
      _selectedRule = Undefined;
      return false;
  }
  return true;
}

void SelectionStrategy::initializeIncumbent()
{
  switch (_selectedRule)
  {
    case NotConsideredForSelection:
      _incumbentCriteria = 0;
      break; // unused
    case FirstFound:
      _incumbentCriteria = 0;
      break; // consider only fractional var
    case HighestPriority:
      _incumbentCriteria = 0;
      break; // maxPriority
    case MostFractional:
      _incumbentCriteria = 2;
      break; //  maxFrac
    case LeastFractional:
      _incumbentCriteria = 2;
      break;  // minFrac
    case FracWeightedPriority:
      _incumbentCriteria = BapcodInfinity;
      break; // maxWeightedPriority
    case Closest2RoundUp:
      _incumbentCriteria = 2;
      break;  // minFrac
    case Closest2RoundDown:
      _incumbentCriteria = 2;
      break;  // minFrac
    case GuidedSearch:
      _incumbentCriteria = 2;
      break;  // minFrac
    case LeastCost:
      _incumbentCriteria = BapcodInfinity;
      break; // minCost
    case LeastReducedCost:
      _incumbentCriteria = BapcodInfinity;
      break; // minCost
    case LeastGreedyCost:
      _incumbentCriteria = BapcodInfinity;
      break; // minCost
    case LeastSteepestEdgeCost:
      _incumbentCriteria = BapcodInfinity;
      break; // minCost
    case LeastPseudoCost:
      _incumbentCriteria = BapcodInfinity;
      break; // minCost
    case LeastInfeasibility:
      _incumbentCriteria = 0;
      break; // minInfeasibilities
    case MostViolated:
      _incumbentCriteria = 0;
      break; // maxViolation
    case NotConsideredForIntegralityCheck:
      _incumbentCriteria = 0;
      break; // maxViolation
    default:
      break;
  }
  if (printL(5))
      std::cout << "SelectionStrategy::initializeCriteria(): _incumbentCriteria = "
                << _incumbentCriteria << " selectedRule = " << _selectedRule
                << std::endl;
  return;
}

const int SelectionStrategy::computeCriteriaAndValToCompareToIncumbent(const SolutionVarInfo * varInfo,
                                                                       Double & challengerRoundedValue,
                                                                       const int & printlevel)
{
  int status(-1);
  Double challengerCriteriaValue(0);

  Double floorVal = Dfloor(varInfo->value);
  Double ceilVal = Dceil(varInfo->value);

  switch (_selectedRule)
    {
  case NotConsideredForSelection: // unused
  case FirstFound:
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = -1;
        if (challengerCriteriaValue < _incumbentCriteria)
          {
            if (printL(printlevel))
                std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() varConstr "
                          << varInfo->varPtr->name() << " criteriaValue = "
                          << challengerCriteriaValue << " _incumbentCriteria = "
                          << _incumbentCriteria << " rounded down replace incumbent "
                          << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incumbent
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            status = 0;
            challengerRoundedValue = floorVal;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = -1;
        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
                std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() varConstr "
                          << varInfo->varPtr->name() << " criteriaValue = "
                          << challengerCriteriaValue << " _incumbentCriteria = "
                          << _incumbentCriteria << " rounded up replace incumbent "
                          << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1) /// if round down is not the incumbent
        && (challengerCriteriaValue == _incumbentCriteria))
          {
            status = 0;
            challengerRoundedValue = ceilVal;
          }
      }
    break;

  case HighestPriority: // maxPriority
    challengerCriteriaValue = -1;
    if (varInfo->canRoundDown || varInfo->canRoundUp)
      challengerCriteriaValue = varInfo->priorityLevel;
    if (varInfo->canRoundDown)
      /// rounding down is better, if possible
      challengerRoundedValue = floorVal;
    else
      challengerRoundedValue = ceilVal;

    if (challengerCriteriaValue > _incumbentCriteria) // strictly better
      {
        if (printL(printlevel))
            std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                      << varInfo->varPtr->name() << " criteriaValue = "
                      << challengerCriteriaValue << " _incumbentCriteria = "
                      << _incumbentCriteria << " replace incumbent "
                      << std::endl;
        _incumbentCriteria = challengerCriteriaValue;
        status = 1;
      }
    else if (challengerCriteriaValue == _incumbentCriteria)
      {
        status = 0;
        challengerRoundedValue = floorVal;
      }
    break;

  case MostFractional: //  maxFrac
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = 1 - (varInfo->value - floorVal);

        if (printL(5))
            std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() MostFractional var = "
                      << varInfo->varPtr->name() << " floor = " << floorVal
                      << " challengerCriteriaValue = " << challengerCriteriaValue
                      << std::endl;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
                std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                          << varInfo->varPtr->name() << " criteriaValue = "
                          << challengerCriteriaValue << " _incumbentCriteria = "
                          << _incumbentCriteria << " rounded down replace incumbent "
                          << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incumbent
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            status = 0;
            challengerRoundedValue = floorVal;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = 1 - (ceilVal - varInfo->value);
        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
                std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                          << varInfo->varPtr->name() << " criteriaValue = "
                          << challengerCriteriaValue << " _incumbentCriteria = "
                          << _incumbentCriteria << " ceilVal = " << ceilVal
                          << " rounded up replace incumbent " << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else
            if ((status != 1) /// if round down is not the incumbent
                && (challengerCriteriaValue == _incumbentCriteria))
            {
                status = 0;
                challengerRoundedValue = ceilVal;
            }
      }

    if (printL(5))
        std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() MostFractional varConstr "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue << " status = " << status
                  << " _incumbentCriteria = " << _incumbentCriteria << std::endl;
    break;


  case MostViolated: // maxViolation
  case Closest2RoundDown: // unused
  case LeastSteepestEdgeCost: // unused
  case FracWeightedPriority: // unused
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = (varInfo->value - floorVal) * varInfo->cost;

        if (printL(5))
            std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() FracWeightedPriority var = "
                      << varInfo->varPtr->name() << " floor = " << floorVal
                      << " challengerCriteriaValue = " << challengerCriteriaValue
                      << std::endl;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
              std::cout
                  << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue << " _incumbentCriteria = "
                  << _incumbentCriteria << " rounded down replace incumbent "
                  << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incumbent
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            status = 0;
            challengerRoundedValue = floorVal;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = (ceilVal - varInfo->value) * varInfo->cost;
        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
              std::cout
                  << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() FracWeightedPriority var "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue << " _incumbentCriteria = "
                  << _incumbentCriteria << " ceilVal = " << ceilVal
                  << " rounded up replace incumbent " << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1) /// if round down is not the incumbent
        && (challengerCriteriaValue == _incumbentCriteria))
          {
            status = 0;
            challengerRoundedValue = ceilVal;
          }
      }

    if (printL(5))
      std::cout
          << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() FracWeightedPriority varConstr "
          << varInfo->varPtr->name() << " criteriaValue = "
          << challengerCriteriaValue << " status = " << status
          << " _incumbentCriteria = " << _incumbentCriteria << std::endl;
    break;

  case LeastInfeasibility: // criteriaValue = BapcodInfinity; break; // minInfeasibilities
  {
    if (varInfo->varPtr->isTypeOf(VcId::MastColumnMask) && varInfo->canRoundUp)
    {
      MastColumn *colPtr = static_cast<MastColumn *>(varInfo->varPtr);
      challengerCriteriaValue = 0;
      if (!colPtr->spSol()->resConsumption().empty())
        challengerCriteriaValue = colPtr->spSol()->resConsumption().back()[0] * sqrt(varInfo->value.val());
      if (challengerCriteriaValue > _incumbentCriteria) // strictly better
      {
        if (printL(printlevel))
          std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                    << varInfo->varPtr->name() << " criteriaValue = "
                    << challengerCriteriaValue << " _incumbentCriteria = "
                    << _incumbentCriteria << " ceilVal = " << ceilVal
                    << " rounded up replace incumbent " << std::endl;
        _incumbentCriteria = challengerCriteriaValue;
        status = 1; // round up replace incumbent
        challengerRoundedValue = ceilVal;
      }
      else if (challengerCriteriaValue == _incumbentCriteria)
      {
        status = 0;
        challengerRoundedValue = ceilVal;
      }
    }

    if (printL(5))
      std::cout << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() LeastInfeasibility varConstr "
                << varInfo->varPtr->name() << " criteriaValue = "
                << challengerCriteriaValue << " status = " << status
                << " _incumbentCriteria = " << _incumbentCriteria << std::endl;
    break;
  }
  case LeastFractional: // minFrac
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = varInfo->value - floorVal;

        if (printL(5))
          std::cout
              << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() LeastFractional var = "
              << varInfo->varPtr->name() << " floor = " << floorVal
              << " challengerCriteriaValue = " << challengerCriteriaValue
              << std::endl;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
              std::cout
                  << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue << " _incumbentCriteria = "
                  << _incumbentCriteria << " rounded down replace incumbent "
                  << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incumbent
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            status = 0;
            challengerRoundedValue = floorVal;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = ceilVal - varInfo->value;
        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            if (printL(printlevel))
              std::cout
                  << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue << " _incumbentCriteria = "
                  << _incumbentCriteria << " ceilVal = " << ceilVal
                  << " rounded up replace incumbent " << std::endl;
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1) /// if round down is not the incumbent
        && (challengerCriteriaValue == _incumbentCriteria))
          {
            status = 0;
            challengerRoundedValue = ceilVal;
          }
      }

    if (printL(5))
      std::cout
          << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() LeastFractional varConstr "
          << varInfo->varPtr->name() << " criteriaValue = "
          << challengerCriteriaValue << " status = " << status
          << " _incumbentCriteria = " << _incumbentCriteria << std::endl;
    break;

  case Closest2RoundUp: // minFrac
    if (varInfo->canRoundUp)
      {
        // consider round up if different from round down
        challengerCriteriaValue = ceilVal - varInfo->value;

        if (printL(printlevel))
          std::cout
              << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
              << varInfo->varPtr->name() << " criteriaValue = "
              << challengerCriteriaValue << " _incumbentCriteria = "
              << _incumbentCriteria << " col is suitable" << std::endl;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            challengerRoundedValue = ceilVal;
            status = 0;
          }
      }

    if (varInfo->canRoundDown)
    {
        challengerRoundedValue = floorVal;
        status = 0;
    }
    break;

  case GuidedSearch: // minFrac
    /// consider round down if it is in the sense of defaultVal(), and not dominated by round up
    if ((varInfo->canRoundDown) && (varInfo->value >= varInfo->incumbentValue))
      {
        challengerCriteriaValue = varInfo->value - varInfo->incumbentValue;

        if (printL(printlevel))
          std::cout
              << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
              << varInfo->varPtr->name() << " to round down at val " << floorVal
              << " while incumbent val is " << varInfo->incumbentValue
              << " criteriaValue = " << challengerCriteriaValue
              << " _incumbentCriteria = " << _incumbentCriteria
              << " col is suitable" << std::endl;
        if (challengerCriteriaValue < _incumbentCriteria) /// strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; /// round down replace incument
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            status = 0;
            challengerRoundedValue = floorVal;
          }
      }

    /// now consider round up if different from round down, and not dominated by round down
    if ((varInfo->canRoundUp) && (varInfo->value <= varInfo->incumbentValue))
      {
        challengerCriteriaValue = varInfo->incumbentValue - varInfo->value;

        if (printL(printlevel))
          std::cout
              << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() var "
              << varInfo->varPtr->name() << " to round up at val " << ceilVal
              << " criteriaValue = " << challengerCriteriaValue
              << " _incumbentCriteria = " << _incumbentCriteria
              << " col is suitable" << std::endl;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1)
            && (challengerCriteriaValue == _incumbentCriteria))
          {
            status = 0;
            challengerRoundedValue = ceilVal;
          }
      }
    break;

  case LeastCost: // minCost
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = floorVal - varInfo->value; /// negative val
        challengerCriteriaValue *= varInfo->cost;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incument
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            challengerRoundedValue = floorVal;
            status = 0;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = ceilVal - varInfo->value;
        challengerCriteriaValue *= varInfo->cost;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1)
            && (challengerCriteriaValue == _incumbentCriteria))
          {
            challengerRoundedValue = ceilVal;
            status = 0;
          }
      }
    break;

  case LeastReducedCost: // criteriaValue = varConstrPtr->computeReducedCost(); break; // minCost
    if (varInfo->canRoundDown)
      {
        challengerCriteriaValue = floorVal - varInfo->value; /// negative val
        challengerCriteriaValue *= varInfo->reducedCost;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incumbent
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            challengerRoundedValue = floorVal;
            status = 0;
          }
      }

    if (varInfo->canRoundUp)
      {
        challengerCriteriaValue = ceilVal - varInfo->value;
        challengerCriteriaValue *= varInfo->reducedCost;

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1)
            && (challengerCriteriaValue == _incumbentCriteria))
          {
            challengerRoundedValue = ceilVal;
            status = 0;
          }
      }
    break;

  case LeastGreedyCost: // criteriaValue = varConstrPtr->greedyCost(); break; // minCost
    if (varInfo->canRoundDown)
      {
        if (floorVal >= varInfo->value) ///positive sense
          challengerCriteriaValue = varInfo->varPtr->greedyCost(); /// what is greedyCost???
        else
          /// negative sense
          challengerCriteriaValue = -varInfo->varPtr->greedyCost();

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round down replace incument
            challengerRoundedValue = floorVal;
          }
        else if (challengerCriteriaValue == _incumbentCriteria)
          {
            challengerRoundedValue = floorVal;
            status = 0;
          }
      }

    /// now consider round up if different from round down
    if (varInfo->canRoundUp)
      {
        if (ceilVal >= varInfo->value) ///positive sense
          challengerCriteriaValue = varInfo->varPtr->greedyCost();
        else
          /// negative sense
          challengerCriteriaValue = -varInfo->varPtr->greedyCost();

        if (challengerCriteriaValue < _incumbentCriteria) // strictly better
          {
            _incumbentCriteria = challengerCriteriaValue;
            status = 1; // round up replace incumbent
            challengerRoundedValue = ceilVal;
          }
        else if ((status != 1)
            && (challengerCriteriaValue == _incumbentCriteria))
          {
            challengerRoundedValue = ceilVal;
            status = 0;
          }
      }

    break;

  default:
    break;
    }
  if (printL(printlevel))
    {
      std::cout
        << "SelectionStrategy::computeCriteriaAndValToCompareToImcumbent() varConstr "
        << varInfo->varPtr->name() << " RoundedValue = "
        << challengerRoundedValue << " _selectedRule = ";
      print() << " status = " << status << std::endl;
    }
  return status;
}

const int SelectionStrategy::computeCriteriaToCompareToIncumbent(const SolutionVarInfo * varInfo,
                                                                 const Double & challengerRoundedValue,
                                                                 const int & printlevel)
{
  /// note that the challengerRoundedValue is assumed to be suitable
  /// for residual problem and therefore needs not be tested.

  int status(-1);
  Double challengerCriteriaValue(0);
  Double floorVal = Dfloor(varInfo->value);
  Double ceilVal = Dceil(varInfo->value);

  switch (_selectedRule)
    {
  case NotConsideredForSelection: // unused
  case FirstFound:
    {
      if (challengerRoundedValue > 0)
        {
          challengerCriteriaValue = -1;
          if (challengerCriteriaValue <= _incumbentCriteria)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() var "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded down replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round down replace incument
                    }
                  else /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      if ((_incumbentCriteria == 2) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  case HighestPriority: // maxPriority
    {
      challengerCriteriaValue = -1;
      if (varInfo->canRoundDown || varInfo->canRoundUp)
        challengerCriteriaValue = varInfo->priorityLevel;
      if (challengerCriteriaValue > _incumbentCriteria) // strictly better
        {
          if (printL(printlevel))
            std::cout << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                      << varInfo->varPtr->name() << " criteriaValue = "
                      << challengerCriteriaValue
                      << " _incumbentCriteria = " << _incumbentCriteria
                      << " replace incumbent " << std::endl;
          _incumbentCriteria = challengerCriteriaValue;
          status = 1;
        }
      else if (challengerCriteriaValue == _incumbentCriteria)
        {
          status = 0;
        }
      break;
    }
  case MostFractional: //  maxFrac
    {
      /// consider round down
      if ((floorVal > 0) && (challengerRoundedValue < varInfo->value))
        {
          challengerCriteriaValue = 1 - (varInfo->value - challengerRoundedValue);
          if (challengerCriteriaValue <= _incumbentCriteria)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded down replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round down replace incument
                    }
                  else /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      /// now consider round up if different from round down
      if (challengerRoundedValue > floorVal)
        {
          challengerCriteriaValue = 1 - (challengerRoundedValue - varInfo->value);
          if (challengerCriteriaValue <= _incumbentCriteria) // better or equal to incumbent (aht already includes the round down option)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded up replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round up replace incumbent
                    }
                  else if (status != 1) /// if round down is not the incumbent
                  /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      if ((_incumbentCriteria == 2) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  case MostViolated: // maxViolation
  case Closest2RoundDown: // unused
  case LeastSteepestEdgeCost: // unused
  case FracWeightedPriority: // unused
    {
      /// consider round down
      if ((floorVal > 0) && (challengerRoundedValue < varInfo->value))
        {
          challengerCriteriaValue = (varInfo->value - challengerRoundedValue) * varInfo->cost;
          if (challengerCriteriaValue <= _incumbentCriteria)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded down replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round down replace incument
                    }
                  else /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      /// now consider round up if different from round down
      if (challengerRoundedValue > floorVal)
        {
          challengerCriteriaValue = (challengerRoundedValue - varInfo->value) * varInfo->cost;
          if (challengerCriteriaValue <= _incumbentCriteria) // better or equal to incumbent (aht already includes the round down option)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded up replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round up replace incumbent
                    }
                  else if (status != 1) /// if round down is not the incumbent
                  /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      if ((_incumbentCriteria == BapcodInfinity) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  case LeastInfeasibility:
  {
    challengerCriteriaValue = 0;
    if (varInfo->varPtr->isTypeOf(VcId::MastColumnMask) && varInfo->canRoundUp)
    {
      MastColumn *colPtr = static_cast<MastColumn *>(varInfo->varPtr);
      challengerCriteriaValue = 0;
      if (!colPtr->spSol()->resConsumption().empty())
        challengerCriteriaValue = -colPtr->spSol()->resConsumption().back()[0] * sqrt(varInfo->value.val());
    }
    if (challengerCriteriaValue > _incumbentCriteria) // strictly better
    {
      if (printL(printlevel))
        std::cout << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                  << varInfo->varPtr->name() << " criteriaValue = "
                  << challengerCriteriaValue
                  << " _incumbentCriteria = " << _incumbentCriteria
                  << " replace incumbent " << std::endl;
      _incumbentCriteria = challengerCriteriaValue;
      status = 1;
    }
    else if (challengerCriteriaValue == _incumbentCriteria)
    {
      status = 0;
    }
    break;
  }
  case LeastFractional: // minFrac
    {
      /// consider round down
      if ((floorVal > 0) && (challengerRoundedValue < varInfo->value))
        {
          challengerCriteriaValue = varInfo->value - challengerRoundedValue;
          if (challengerCriteriaValue <= _incumbentCriteria)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded down replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round down replace incument
                    }
                  else /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      /// now consider round up if different from round down
      if (challengerRoundedValue > floorVal)
        {
          challengerCriteriaValue = challengerRoundedValue - varInfo->value;
          if (challengerCriteriaValue <= _incumbentCriteria) // better or equal to incumbent (aht already includes the round down option)
            {
                  if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                    {
                      if (printL(printlevel))
                        std::cout
                            << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                            << varInfo->varPtr->name() << " criteriaValue = "
                            << challengerCriteriaValue
                            << " _incumbentCriteria = " << _incumbentCriteria
                            << " rounded up replace incumbent " << std::endl;
                      _incumbentCriteria = challengerCriteriaValue;
                      status = 1; /// round up replace incumbent
                    }
                  else if (status != 1) /// if round down is not the incumbent
                  /// equal to incumbent
                    {
                      status = 0;
                    }
            }
        }

      if ((_incumbentCriteria == 2) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  case Closest2RoundUp: /// minFrac
    {
      /// consider gap to rounded up val
      challengerCriteriaValue = ceilVal - challengerRoundedValue;
      if (challengerCriteriaValue <= _incumbentCriteria) /// better or equal to incumbent (that already includes the round down option)
        {
              if (printL(printlevel))
                std::cout
                    << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
                    << varInfo->varPtr->name() << " criteriaValue = "
                    << challengerCriteriaValue << " _incumbentCriteria = "
                    << _incumbentCriteria << std::endl;
              if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                {
                  _incumbentCriteria = challengerCriteriaValue;
                  status = 1; // round up replace incumbent
                }
              else // equal to incumbent
                {
                  status = 0;
                }
        }

      if ((_incumbentCriteria == 2) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  case GuidedSearch: /// minFrac
    {
      /// consider gap to default val
      if (challengerRoundedValue < varInfo->incumbentValue)
        challengerCriteriaValue = varInfo->incumbentValue - challengerRoundedValue;
      else
        challengerCriteriaValue = challengerRoundedValue - varInfo->incumbentValue;

      if (challengerCriteriaValue <= _incumbentCriteria)
        {
              if (challengerCriteriaValue < _incumbentCriteria) /// strictly better
                {
                  _incumbentCriteria = challengerCriteriaValue;
                  status = 1; ///  replace incumbent
                }
              else /// equal to incumbent
                {
                  status = 0;
                }
        }

      if ((_incumbentCriteria == 2) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }
      break;
    }

  case LeastCost: /// minCost
    {
      /// consider gap to current val with its sign
      challengerCriteriaValue = challengerRoundedValue - varInfo->value;

      challengerCriteriaValue *= varInfo->cost;
      if (challengerCriteriaValue <= _incumbentCriteria)
        {
              if (challengerCriteriaValue < _incumbentCriteria) /// strictly better
                {
                  _incumbentCriteria = challengerCriteriaValue;
                  status = 1; /// replace incumbent
                }
              else // equal to incumbent
                {
                  status = 0;
                }
        }

      /// no existing incumbent: cannot be compared to incumbent
      if ((_incumbentCriteria == BapcodInfinity) && (status == -1))
        {
          status = 0;
        }

      break;
    }
  case LeastReducedCost: // criteriaValue = varConstrPtr->computeReducedCost(); break; // minCost
    {
      /// consider gap to current val with its sign
      challengerCriteriaValue = challengerRoundedValue - varInfo->value;

      challengerCriteriaValue *= varInfo->reducedCost;
      if (challengerCriteriaValue <= _incumbentCriteria)
        {
          if (challengerCriteriaValue < _incumbentCriteria) /// strictly better
            {
              _incumbentCriteria = challengerCriteriaValue;
              status = 1; /// replace incumbent
            }
          else // equal to incumbent
            {
              status = 0;
            }
        }

      /// no existing incumbent: cannot be compared to incumbent
      if ((_incumbentCriteria == BapcodInfinity) && (status == -1))
        {
          status = 0;
        }

      break;
    }
  case LeastGreedyCost: // criteriaValue = varConstrPtr->greedyCost(); break; // minCost
    {
      /// consider cost per unit of delta val sith it sense
      if (challengerRoundedValue < varInfo->value) /// negative sense
        challengerCriteriaValue = -varInfo->varPtr->greedyCost(); /// what is greedyCost ???
      else
        challengerCriteriaValue = varInfo->varPtr->greedyCost();

      if (challengerCriteriaValue <= _incumbentCriteria)
        {
              if (challengerCriteriaValue < _incumbentCriteria) // strictly better
                {
                  _incumbentCriteria = challengerCriteriaValue;
                  status = 1; // replace incumbent
                }
              else /// equal to incumbent
                {
                  status = 0;
                }
        }

      if ((_incumbentCriteria == BapcodInfinity) && (status == -1)) /// no existing incumbent: cannot be compared to incumbent
        {
          status = 0;
        }

      break;
    }
  default:
  break;
    }
  if (printL(printlevel))
    {
      std::cout
        << "SelectionStrategy::computeCriteriaToCompareToImcumbent() varConstr "
        << varInfo->varPtr->name() << " RoundedValue = "
        << challengerRoundedValue << " _selectedRule = ";
      print() << " status = " << status << std::endl;
    }
  return status;
}

std::ostream& SelectionStrategy::print(std::ostream& os) const
{
  switch (_selectedRule)
  {
    case NotConsideredForSelection:
      return (os << "NotConsideredForSelection");
    case FirstFound:
      return (os << "FirstFound");
    case HighestPriority:
      return (os << "HighestPriority");
    case MostFractional:
      return (os << "MostFractional");
    case LeastFractional:
      return (os << "LeastFractional");
    case FracWeightedPriority:
      return (os << "FracWeightedPriority");
    case Closest2RoundUp:
      return (os << "Closest2RoundUp");
    case Closest2RoundDown:
      return (os << "Closest2RoundDown");
    case GuidedSearch:
      return (os << "GuidedSearch");
    case LeastCost:
      return (os << "LeastCost");
    case LeastReducedCost:
      return (os << "LeastReducedCost");
    case LeastGreedyCost:
      return (os << "LeastGreedyCost");
    case LeastSteepestEdgeCost:
      return (os << "LeastSteepestEdgeCost");
    case LeastPseudoCost:
      return (os << "LeastPseudoCost");
    case LeastInfeasibility:
      return (os << "LeastInfeasibility");
    case MostViolated:
      return (os << "MostViolated");
    case NotConsideredForIntegralityCheck:
      return (os << "NotConsideredForIntegralityCheck");
    default:
      break;
  }
  return (os << "Undefined");
}

