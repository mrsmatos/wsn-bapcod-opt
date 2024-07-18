/**
 *
 * This file bcProbConfigC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcProbConfigC.hpp"

#include <utility>
#include "bcSpVarConstrC.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcModelFormulationC.hpp"
#include "bcModelColGenSpC.hpp"
#include "bcFormC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_data.hpp"
#endif //BCP_RCSP_IS_FOUND

/**
 * Generic code
 */
using namespace std;
using namespace VcIndexStatus;

BranchingEvaluationInfo::BranchingEvaluationInfo(const int & depth_, const Double & dualBound_,
                                                 const Double & lhsFracPart_, const char & direction_,
                                                 BranchingGeneratorHistory * const historyPtr_) :
  depth(depth_), dualBound(dualBound_), lhsFracPart(lhsFracPart_), direction(direction_),
  phaseInfo(), historyPtr(historyPtr_)
{
}

void BranchingEvaluationInfo::addEvalEntry(const int & phaseNumber, const int & childNumber,
                                           const double & lpValue, const double & evalTime)
{
  /// first, we add the evaluation entry
  if (phaseInfo.size() < phaseNumber)
    phaseInfo.resize(phaseNumber, BranchingPhaseEvoluationInfo());
  BranchingPhaseEvoluationInfo & phaseEvalInfo = phaseInfo[phaseNumber - 1];
  if (phaseEvalInfo.size() < childNumber)
    phaseEvalInfo.resize(childNumber);
  phaseEvalInfo.setEvalEntry(childNumber - 1, lpValue - dualBound, evalTime);
  /// second, we update the pseudo-costs
  /// Ruslan : note that for the moment pseudo-costs are valid only for the case of
  /// two child nodes (the code needs to be adapted for the component set branching)
  double unitaryImprovement(0.0);
  if (((childNumber == 1) && (direction == 'U')) || ((childNumber == 2) && (direction == 'L')))
    unitaryImprovement = (lpValue - dualBound) / (1 - lhsFracPart);
  else
    unitaryImprovement = (lpValue - dualBound) / lhsFracPart;
  if (historyPtr->numEvaluations.size() < phaseNumber)
    {
      historyPtr->numEvaluations.resize(phaseNumber);
      historyPtr->pseudoCosts.resize(phaseNumber);
    }
  std::vector<double> & pseudoCosts = historyPtr->pseudoCosts[phaseNumber - 1];
  std::vector<int> & numEvaluations = historyPtr->numEvaluations[phaseNumber - 1];
  if (numEvaluations.size() < childNumber)
    {
      numEvaluations.resize(childNumber, 0);
      pseudoCosts.resize(childNumber, 0.0);
    }
  numEvaluations[childNumber - 1] += 1;
  pseudoCosts[childNumber - 1] = (pseudoCosts[childNumber - 1] * (numEvaluations[childNumber - 1] - 1)
                                  + unitaryImprovement) / numEvaluations[childNumber - 1];

}

void ProbConfig::increasePCSolCount()
{
    _pcSolCount += 1;
}

void ProbConfig::increasePCNodeCount()
{
    _pcNodeCount += 1;
}

ProbConfig::ProbConfig(ProbConfType configType,
       		           Model * modelPtr,
                       std::string genericName,
		               const IndexCell & id,
		               const Bound & primalIncBound,
		               const Bound & dualIncBound,
		               Problem * problemPtr):
        _configType(configType),
        _modelPtr(modelPtr),
        _genericName(std::move(genericName)),
        _ref(problemPtr->ref()),
        _pcSolCount(0),
        _pcNodeCount(0),
        _id(id),
        _primalIncBound(Bound::infPrimalBound(_modelPtr->objectiveSense())),
        _dualIncBound(dualIncBound),
        _cutOffValue(Bound::infPrimalBound(_modelPtr->objectiveSense())),
        _ovfConfPtr(NULL),
        _networkFlowPtr(NULL),
#ifdef BCP_RCSP_IS_FOUND
        _rcspGraphPtr(NULL),
#endif //BCP_RCSP_IS_FOUND
        _defaultGenericVar(NULL),
        _defaultGenericConstr(NULL),
        _dualBoundContrib(0, _modelPtr->objectiveSense()),
        _probPtr(problemPtr),
        _curFormIsInfeasible(false),
        _primalSolPtr(NULL),
        _isPrepared(false)
{
    _probPtr->probConfPtr(this);
    resetCutOffValue(primalIncBound);
    _primalIncBound = primalIncBound;
    return;
}

ProbConfig::~ProbConfig()
{
  _colGenSubProbConfPts.clear();

  if (_probPtr != NULL)
    {
      delete _probPtr;
      _probPtr = NULL;
    }

  for (std::list< Variable *>::const_iterator it = _pcVarPtrList.begin(); it != _pcVarPtrList.end(); ++it)
    {
      VarConstr* varConstr = *it;
      _modelPtr->garbageCollector().erase(varConstr);
      delete varConstr;
    }
  _pcVarPtrList.clear();

  for (std::list< Constraint *>::const_iterator it = _pcConstrPtrList.begin(); it != _pcConstrPtrList.end(); ++it)
    {
      VarConstr * varConstr = *it;
      _modelPtr->garbageCollector().erase(varConstr);
      delete varConstr;
    }
  _pcConstrPtrList.clear();

  for (map<string, GenericVar*>::const_iterator it = _name2GenericVarPtrMap.begin();
       it != _name2GenericVarPtrMap.end(); it++)
    {
      GenericVar* genVar = (*it).second;
      delete genVar;
    }
  _name2GenericVarPtrMap.clear();

  for (map<string, GenericConstr*>::const_iterator it = _name2GenericConstrPtrMap.begin();
       it != _name2GenericConstrPtrMap.end(); it++)
    {
      GenericConstr* genConstr = (*it).second;
      delete genConstr;
    }
  _name2GenericConstrPtrMap.clear();

  std::set<GenericBranchingConstr *, DynamicGenConstrSort>::iterator brGenConstrPtrIt;
  for (brGenConstrPtrIt = _candidateBranchingGenericConstr.begin();
       brGenConstrPtrIt != _candidateBranchingGenericConstr.end(); ++brGenConstrPtrIt)
    delete *brGenConstrPtrIt;
  _candidateBranchingGenericConstr.clear();

  std::set<GenericCutConstr *, DynamicGenConstrSort>::iterator cutGenConstrPtrIt;
  for (cutGenConstrPtrIt = _candidateCutGenericConstr.begin();
       cutGenConstrPtrIt != _candidateCutGenericConstr.end(); ++cutGenConstrPtrIt)
    delete *cutGenConstrPtrIt;
  _candidateCutGenericConstr.clear();

  if (_networkFlowPtr != NULL)
    delete _networkFlowPtr;
  _networkFlowPtr = NULL;

#ifdef BCP_RCSP_IS_FOUND
    if (_rcspGraphPtr != NULL)
        delete _rcspGraphPtr;
    _rcspGraphPtr = NULL;
#endif //BCP_RCSP_IS_FOUND
  return;
}

void ProbConfig::resetCutOffValue(const Bound & incumbentVal) 
{ 
  /// this is needed as resetCutOffValue can be called from a node
  if (incumbentVal >= _primalIncBound)
    return;

  BcObjStatus::MinMaxIntFloat objStatus = _modelPtr->objectiveSense();
  if ((objStatus == BcObjStatus::minInt) || (objStatus == BcObjStatus::minFloat))
  {
    if (objStatus == BcObjStatus::minInt)
      _cutOffValue = floor((double)(incumbentVal - Double::precision)) + param().optimalityGapTolerance() * 2;
    else
      _cutOffValue = incumbentVal - param().optimalityGapTolerance() / 2;
  }
  else
  {
    if (objStatus == BcObjStatus::maxInt)
      _cutOffValue = ceil((double)(incumbentVal + Double::precision)) - param().optimalityGapTolerance() * 2;
    else
      _cutOffValue = incumbentVal + param().optimalityGapTolerance() / 2;
  }
}/// if one finds a dual bound above this custoff value, the solver stops.

bool ProbConfig::updatePrimalIncBound(const Bound & newIncVal)
{
  if (printL(2))
    {
      std::cout << "ProbConfig::updatePrimalIncBound() new primal val = " << newIncVal << std::endl;
    }

  if (newIncVal < _primalIncBound)
    {
      if (printL(1))
        {
          std::cout << "*************************** " << std::endl;

          std::cout << "New Incumbent Solution: " << newIncVal << " < " << _primalIncBound
		    << std::endl;

          std::cout << "*************************** " << std::endl;
        }
      resetCutOffValue(newIncVal);
      _primalIncBound = newIncVal;

      if (printL(2) && (probPtr() != NULL) && (probPtr()->incumbentSol() != NULL))
      {
          std::cout << "ProbConfig::updatePrimalIncBound() new Inc Solution " << _primalIncBound
                    << std::endl;
          probPtr()->incumbentSol()->printVar();
      }

      return (true);
    }

  return (false);
}


bool ProbConfig::updateDualIncBound(const Bound & bd)
{
  if (printL(2))
    {
      std::cout << "ProbConfig::updateDualIncBound() new dual val = " << bd << std::endl;
    }
  if (bd > _dualIncBound)
    {
      /// verification in order not to set dual bound greater than the primal one
      /// (bd > _primalIncBound is possible when using strong branching)
      if (bd > _primalIncBound)
        _dualIncBound = _primalIncBound;
      else
        _dualIncBound = bd;
      return true;
    }
  return false;
}

void ProbConfig::resetDualIncBound()
{
    _dualIncBound = Bound::infDualBound(modelPtr()->objectiveSense());
}

void ProbConfig::resetPrimalIncBound()
{
    _primalIncBound = Bound::infPrimalBound(modelPtr()->objectiveSense());
}

const IndexCell & ProbConfig::id() const
{
  return(_id);
}

const bool ProbConfig::enumeratedStatus() const
{
  if (probPtr() == NULL)
      return false;
  return probPtr()->getEnumeratedStatus();
}

GenericVar * ProbConfig::defaultGenericVarPtr() const
{
  return _defaultGenericVar;
}

GenericConstr * ProbConfig::defaultGenericConstrPtr() const
{
  return _defaultGenericConstr;
}

ColGenSpConf * ProbConfig::colGenSubProbConf(const int & index) const
{
  if ((index < 0) || (index >= static_cast<int>(_colGenSubProbConfPts.size())))
      return NULL;
  return _colGenSubProbConfPts[index];
}


void ProbConfig::createDefaultGenericVarConstr()
{
  _defaultGenericVar = _modelPtr->createGenericVar(this,
						   BcVarConstrType::local2Formulation,
						   "",
						   MultiIndexNames(),
						   'I',
						   0,  /// cost
						   BapcodInfinity, /// ub
						   SelectionStrategy::MostFractional,
						   /// priorities changed by Ruslan
						   11.0,  /// genericBranchingOnAggregateVarPL
						   1.0 /// compBoundSetBranchingPL, this will be adjusted during setupProbConfig
						   );

  _defaultGenericConstr = _modelPtr->createGenericConstr(this,
							 BcVarConstrType::local2Formulation,
							 "",
							 MultiIndexNames(),
							 'E',
							 0,
							 0,
							 true,
							 's');

}

Constraint * ProbConfig::castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddConstraint() should not be called");
  return constrPtr;
}

InstanciatedConstr * ProbConfig::castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddConstraint() should not be called");
  return iconstrPtr;
}

Variable * ProbConfig::castAndAddVariable(Variable * varPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddVariable() should not be called");
  return varPtr;
}

InstanciatedVar * ProbConfig::castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddVariable() should not be called");
  return ivarPtr;
}


const Bound & ProbConfig::dualBoundContrib() const
{
  return(_dualBoundContrib);
}

void ProbConfig::insertInstVar(InstanciatedVar * iVarPtr)
{
  _iVarPts.push_back(iVarPtr);
  
  return;
}

void ProbConfig::insertInstConstr(InstanciatedConstr * iConstrPtr)
{
  _iConstrPts.push_back(iConstrPtr);
  
  return;
}

Problem * ProbConfig::probPtr() const
{
  return(_probPtr);
}

/** 
 * Setup the math programming problem  without generating the formulation matrix yet
 */
void ProbConfig::prepareProbConfig()
{
  if (_isPrepared)
      return;

  _isPrepared = true;

  probPtr()->defineFormulation();

  _pcConstrPtrList.insert(_pcConstrPtrList.end(), _iConstrPts.begin(), _iConstrPts.end());

  _pcVarPtrList.insert(_pcVarPtrList.end(), _iVarPts.begin(), _iVarPts.end());

  probPtr()->addVarSet(_pcVarPtrList, 1, 0);
  probPtr()->addConstrSet(_pcConstrPtrList, 1, 0);
  probPtr()->buildProblem();

  return;
}

/// this function just solves a MIP without applying any decomposition
Solution * ProbConfig::solvePC(bool showOutput)
{
    if (!progStatus().doRun())
        return(NULL);

    if (_probPtr == NULL)
        return(NULL);

    if (_primalSolPtr != NULL)
    {
        delete _primalSolPtr;
        _primalSolPtr = NULL;
    }

    MathProgSolverInterface * MIPinterfacePtr = _probPtr->primalFormulationPtr()->interfacePtr();
    MIPinterfacePtr->setScreenOutput(showOutput);

    Double incumbentValue(_modelPtr->initialPrimalBound());
    if (_modelPtr->objectiveSense() == BcObjStatus::minInt)
        incumbentValue = floor((double)(incumbentValue - Double::precision)) + Double::precision;
    MIPinterfacePtr->setUpperCutOffValue(incumbentValue);

    MIPinterfacePtr->setSolveFromScratch(true);

    Time start;
    /// Return true if a feasible solution to the problem was found, false otherwise
    int maxLevelOfRestriction(0);
    char flag = (probPtr()->solMode().status() == SolutionMethod::lpSolver) ? 'd' : ' ';
    int solverReturnStatus(probPtr()->solveProb(maxLevelOfRestriction, flag, printL(3)));
    double time = start.getElapsedTime_dbl();
    if (printL(3))
        std::cout << "bcTimeMIPSol" << time << std::endl;
    bapcodInit().statistics().incrTimer("bcTimeMIPSol", time);
    bapcodInit().statistics().incrTimer("bcTimeMain", bapcodInit().startTime().getElapsedTime_dbl());

    bapcodInit().statistics().incrValue("bcRecBestInc", _probPtr->primalBound());
    bapcodInit().statistics().incrValue("bcRecBestDb", _probPtr->dualBound());

    if (solverReturnStatus > 0)
        _primalSolPtr = extractCurrentSol();

    MIPinterfacePtr->resetUpperCutOffValue();
    MIPinterfacePtr->resetSolveFromScratch();
    MIPinterfacePtr->setScreenOutput(false);

    return _primalSolPtr;
}

Solution * ProbConfig::getSolution(const VarPtrSet & curSol) 
{
  Solution * curSolPtr = new Solution(this, NULL);

  curSolPtr->includeVarSet(curSol);

  return curSolPtr;
}

Solution * ProbConfig::getSolution(const MasterVarSolution & varPtrList)
{
  Solution * curSolPtr = new Solution(this, NULL);

  for (MasterVarSolution::const_iterator varPt = varPtrList.begin(); varPt !=varPtrList.end(); ++varPt)
    {       
      curSolPtr->includeVar(varPt->first, varPt->second._value, true);
    }

  return curSolPtr;
}

Solution * ProbConfig::getDissagregatedSolution(Solution * curSolPtr) 
{
  return curSolPtr;
}

Solution * ProbConfig::getAggregatedSolution(Solution * curSolPtr) 
{
  return curSolPtr;
}

Solution * ProbConfig::extractCurrentSol()
{
  if (_primalSolPtr != NULL)
    {
      delete _primalSolPtr;
      _primalSolPtr = NULL;
    }

  /// Create Solution from inPrimalLpSol()
  if (!probPtr()->inPrimalLpSol().empty())
    {
      recordCurrentPrimalSol();
      _primalSolPtr = probPtr()->extractIncumbent();

      /// Reset the solution cost (not equal to primalBound()
      _primalSolPtr->resetCost();
    }
 
  return _primalSolPtr;
}

void ProbConfig::recordCurrentPrimalSol()
{
  Solution * intermSolPtr = probPtr()->retrieveCurPrimalLpSol();
  if (intermSolPtr != NULL)
    {
      /// Reset the solution cost (not equal to primalBound()
      intermSolPtr->resetCost();
      probPtr()->recordSolution(intermSolPtr);
    }
  
  return;
}

void ProbConfig::sortTreeOfColClasses()
{
  /// Update dual val cost
  for (ColClassesVector::const_iterator bcitPt = _treeOfColClasses.begin(); bcitPt != _treeOfColClasses.end(); bcitPt++)
    {
      /// Partial reset
      (*bcitPt)->setSigma(-(*bcitPt)->valOrSepPointVal());
      (*bcitPt)->solPtr(NULL);
      (*bcitPt)->associatedPricingSPsolved(false);
      (*bcitPt)->bestReducedCost(Bound(0, _modelPtr->objectiveSense()));

      for (CompSetInstBrConstrPtrSet::const_iterator bcPt =
          (*bcitPt)->setOfPredCSconstrPtr().begin();
          bcPt != (*bcitPt)->setOfPredCSconstrPtr().end(); bcPt++)
        /// Bc constraint dual val
        (*bcitPt)->addToSigma(-(*bcPt)->valOrSepPointVal());
    }

  /**
   * Sort list of classes by increasing depth (= number of precedessors:
   * in case of ties treat first largest dual value class)
   */
  if (printL(5))
    for (ColClassesVector::const_iterator bcitPt = _treeOfColClasses.begin(); bcitPt != _treeOfColClasses.end();
         bcitPt++)
      {
        std::cout << " ColClasses before sorting " << (*bcitPt)->name() << " Lvalue = " << (*bcitPt)->marginLvalue()
                  << " sigma = " << (*bcitPt)->sigma() << " preds : ";
        for (CompSetInstBrConstrPtrSet::iterator it = (*bcitPt)->setOfPredCSconstrPtr().begin();
             it != (*bcitPt)->setOfPredCSconstrPtr().end(); ++it)
          std::cout << " " << (*it)->name();
        std::cout << std::endl;
      }

  std::stable_sort(_treeOfColClasses.begin(), _treeOfColClasses.end(),
                   CompSetInstMastBranchConstr::CSbrGreedyComparator());

  if (printL(5))
    for (ColClassesVector::const_iterator bcitPt = _treeOfColClasses.begin(); bcitPt != _treeOfColClasses.end();
         bcitPt++)
        std::cout << " ColClasses after sorting " << (*bcitPt)->name() << std::endl;

  return;
}

std::ostream& ProbConfig::print(std::ostream& os) const
{
  os << "ProbConfig: " << std::endl;
  os << "   id = " << _id << std::endl;

  if (printL(3))
    {
      os << "  Variables: " << std::endl;
      for (std::list<Variable *>::const_iterator itv = _pcVarPtrList.begin(); 
	   itv != _pcVarPtrList.end(); 
	   itv++)
	{ 
	  os << (*itv)->name() << ", " ;
	}
      
      os << std::endl;
      os << "  Constraints: " << std::endl;
      
      for (std::list<Constraint *>::const_iterator itc = _pcConstrPtrList.begin(); 
	   itc != _pcConstrPtrList.end(); 
	   itc++)
	{ 
	  os << (*itc)->name() << ", " ;
	}
      
      os << std::endl;
    }

  return(os);
}


std::ostream& ProbConfig::printForm(std::ostream& os) const
{
  if  (_probPtr != NULL)
    _probPtr->printForm(os);
  
  return(os);
}

void ProbConfig::nicePrintAllConstraints(std::ostream& os) const
{
  os << "Printing all constrains of ProbConfig " << name() << std::endl;
  os << "Generic constraints : " << std::endl;
  for (std::map< std::string , GenericConstr * >::const_iterator mapIt = _name2GenericConstrPtrMap.begin();
       mapIt != _name2GenericConstrPtrMap.end(); ++mapIt)
    mapIt->second->nicePrintAllConstraints(os);
  os << "Generic cuts : " << std::endl;
  for (std::map< std::string , GenericCutConstr * >::const_iterator mapIt = _name2GenericCutConstrPtrMap.begin();
       mapIt != _name2GenericCutConstrPtrMap.end(); ++mapIt)
    mapIt->second->nicePrintAllConstraints(os);
  os << "Generic branching constraints : " << std::endl;
  for (std::map< std::string , GenericBranchingConstr * >::const_iterator
       mapIt = _name2GenericBranchingConstrPtrMap.begin(); mapIt != _name2GenericBranchingConstrPtrMap.end(); ++mapIt)
    mapIt->second->nicePrintAllConstraints(os);
}

void ProbConfig::insertGenericVar(GenericVar * gvPtr)
{
  _name2GenericVarPtrMap[gvPtr->defaultName()] = gvPtr;
}

void ProbConfig::insertGenericConstr(GenericConstr * gcPtr)
{
  _name2GenericConstrPtrMap[gcPtr->defaultName()] = gcPtr;
}

void ProbConfig::insertGenericCutConstr(GenericCutConstr * gcPtr)
{
 if (printL(5))
    std::cout << "ProbConfig::insertGenericCutConstr: GenericCutConstr = " << gcPtr->defaultName() << std::endl;

  _candidateCutGenericConstr.insert(gcPtr);
  _name2GenericCutConstrPtrMap[gcPtr->defaultName()] = gcPtr;
  return;
}


void ProbConfig::insertGenericBranchingConstr(GenericBranchingConstr * gbcPtr)
{
    if (printL(6))
        std::cout << "ProbConfig::insertGenericBranchingConstr: GenericBranchingConstr = "
                  << gbcPtr->defaultName() << std::endl;

    if ((gbcPtr->priorityLevel() > 0))
    {
        _candidateBranchingGenericConstr.insert(gbcPtr);
        _name2GenericBranchingConstrPtrMap[gbcPtr->defaultName()] = gbcPtr;
    }

    return;
}

GenericVar * ProbConfig::getGenericVar(const std::string & name) const
{
  std::map < std::string, GenericVar * >::const_iterator it = _name2GenericVarPtrMap.find(name);
  if (it != _name2GenericVarPtrMap.end())
      return (it->second);

  return NULL;
}

GenericConstr * ProbConfig::getGenericConstr(const std::string & name) const
{
  std::map < std::string, GenericConstr * >::const_iterator it = _name2GenericConstrPtrMap.find(name);
  
  if (it != _name2GenericConstrPtrMap.end())
      return (it->second);

  return NULL;
}

GenericCutConstr * ProbConfig::getGenericCutConstr(const std::string & name) const
{
  
  if (printL(6)) 
    std::cout << "ProbConfig::getGenericCutConstr search for " << name
              << " in map of size " << _name2GenericCutConstrPtrMap.size() << std::endl;

  std::map < std::string, GenericCutConstr * >::const_iterator it = _name2GenericCutConstrPtrMap.find(name);
  
  if (it != _name2GenericCutConstrPtrMap.end()) 
    return (it->second);

  return NULL;
}

GenericBranchingConstr * ProbConfig::getGenericBranchingConstr(const std::string & name) const
{
  
  if (printL(6)) 
    std::cout << "ProbConfig::getGenericBranchingConstr search for " << name 
	          << " in map of size " << _name2GenericBranchingConstrPtrMap.size() << std::endl;

  std::map < std::string, GenericBranchingConstr * >::const_iterator it = _name2GenericBranchingConstrPtrMap.find(name);
  
  if (it != _name2GenericBranchingConstrPtrMap.end()) 
    return (it->second);

  return NULL;
}

BapcodInit & ProbConfig::bapcodInit() const
{
  return _modelPtr->bapcodInit();
}

BapcodInit* ProbConfig::bapcodInitPtr() const
{
  return _modelPtr->bapcodInitPtr();
}

const ControlParameters& ProbConfig::param() const
{
  return _modelPtr->param();
}

ControlParameters& ProbConfig::param()
{
  return _modelPtr->param();
}

const ProgStatus& ProbConfig::progStatus() const
{
  return _modelPtr->progStatus();
}

ProgStatus& ProbConfig::progStatus()
{
  return _modelPtr->progStatus();
}

void ProbConfig::setMipRequiredStatus(const SolutionStatus & newStatus)
{
    _probPtr->setMIPRequiredStatus(newStatus);
}

#ifdef BCP_RCSP_IS_FOUND
bool ProbConfig::fillRCSPGraph()
{
    delete _rcspGraphPtr;

    _rcspGraphPtr = new bcp_rcsp::GraphData(_id.first());

    if (_networkFlowPtr == nullptr)
    {
        std::cerr << "BaPCod error: cannot build the RCSP graph as a network is not associated with the formulation"
                  << std::endl;
        delete _rcspGraphPtr;
        _rcspGraphPtr = nullptr;
        return false;
    }

    int instVarId = 0;
    for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 's');
         varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); ++varPtrIt)
    {
        auto * instVarPtr = dynamic_cast<InstanciatedVar *>(*varPtrIt);
        if (instVarPtr != nullptr)
        {
            _instVarPts.push_back(instVarPtr);
            _instVarToIdMap[instVarPtr] = instVarId++;
        }
    }

    _rcspGraphPtr->nbElementaritySets = (int)_networkFlowPtr->elementaritySetPts().size();
    _rcspGraphPtr->nbPackingSets = (int)_networkFlowPtr->packingSetPts().size();
    _rcspGraphPtr->nbCoveringSets = (int)_networkFlowPtr->coveringSetPts().size();

    int nbElementaritySets = (int)_networkFlowPtr->elementaritySetPts().size();

    std::vector<const NetworkSet *> elementaritySets(nbElementaritySets, nullptr);
    for (std::vector<NetworkSet *>::const_iterator elemSetPtrIt = _networkFlowPtr->elementaritySetPts().begin();
         elemSetPtrIt != _networkFlowPtr->elementaritySetPts().end(); ++elemSetPtrIt)
    {
        const NetworkSet * elemSetPtr = *elemSetPtrIt;
        int elemSetId = elemSetPtr->id();
        if ((elemSetId < 0) || (elemSetId >= nbElementaritySets))
        {
            std::cerr << "BaPCod error: error: IDs of elementarity sets should be between 0 and size-1" << std::endl;
            delete _rcspGraphPtr;
            _rcspGraphPtr = nullptr;
            return false;
        }
        elementaritySets[elemSetId] = elemSetPtr;
    }
    for (int elemSetAlgIt = 0; elemSetAlgIt < nbElementaritySets; ++elemSetAlgIt)
    {
        if (elementaritySets[elemSetAlgIt] == nullptr)
        {
            std::cerr << "BaPCod error : elementary set with ID = " << elemSetAlgIt << " is not defined" << std::endl;
            delete _rcspGraphPtr;
            _rcspGraphPtr = nullptr;
            return false;
        }
    }

    _rcspGraphPtr->packSetNeighbourhoodForRank1CutSeparation.resize(_rcspGraphPtr->nbPackingSets);
    for (int packSetId = 0; packSetId < _rcspGraphPtr->nbPackingSets; ++packSetId)
    {
        const NetworkSet * packSetPtr = _networkFlowPtr->packingSetPts()[packSetId];
        if (!packSetPtr->cutSeparationNeighbourhood().empty())
        {
            for (const auto packSetPtrInNeighbourhood : packSetPtr->cutSeparationNeighbourhood())
                _rcspGraphPtr->packSetNeighbourhoodForRank1CutSeparation[packSetId].push_back(
                        packSetPtrInNeighbourhood->id());
        }
    }

    if (!_networkFlowPtr->sourceList().empty())
        _rcspGraphPtr->sourceVertexId
            = _networkFlowPtr->netVertexPtr(_networkFlowPtr->sourceList().front())->id();

    if (!_networkFlowPtr->sourceList().empty())
        _rcspGraphPtr->sinkVertexId
            = _networkFlowPtr->netVertexPtr(_networkFlowPtr->sinkList().front())->id();

    for (auto * netResourcePtr : _networkFlowPtr->sideResourceConstrList())
    {
        _rcspGraphPtr->resources.emplace_back(netResourcePtr->id());
        bcp_rcsp::ResourceData & resource = _rcspGraphPtr->resources.back();
        if (netResourcePtr->main())
        {
          if (!netResourcePtr->disposable())
          {
            std::cerr << "BaPCod error : a main resource should always be disposable" << std::endl;
            delete _rcspGraphPtr;
            _rcspGraphPtr = nullptr;
            return false;
          }
          resource.type = bcp_rcsp::ResourceData::Main;
          resource.step = netResourcePtr->step();
        }
        else if (netResourcePtr->base())
        {
          resource.type = bcp_rcsp::ResourceData::Base;
          resource.step = netResourcePtr->step();
        }
        else if (netResourcePtr->dependent())
        {
          resource.type = bcp_rcsp::ResourceData::Dependent;
          resource.baseResId = netResourcePtr->baseResourceId();
          resource.initAccumResConsFunction = netResourcePtr->initAccumResConsFunction();
          if (netResourcePtr->associatedVar() != nullptr)
          {
              resource.associatedVarIds.push_back(_instVarToIdMap[netResourcePtr->associatedVar()]);
              _resIdToVarIdMap[netResourcePtr->id()] = netResourcePtr->associatedVar();
          }
        }
        else if (netResourcePtr->selection())
        {
            resource.type = bcp_rcsp::ResourceData::Selection;
            int nbValues = netResourcePtr->nbValuesForSelection();
            resource.associatedVarIds.resize(nbValues, -1);
            for (auto & pair : netResourcePtr->associatedVarsForSelection())
                if (pair.first >= 0 && pair.first < nbValues)
                    resource.associatedVarIds[pair.first] = _instVarToIdMap[pair.second];
        }
        else if (netResourcePtr->disposable())
        {
          resource.type = bcp_rcsp::ResourceData::SecondaryDisposable;
        }
        else
        {
          resource.type = bcp_rcsp::ResourceData::SecondaryNonDisposable;
        }
        if (!netResourcePtr->selection() && (netResourcePtr->associatedVar() != nullptr))
        {
            resource.associatedVarIds.push_back(_instVarToIdMap[netResourcePtr->associatedVar()]);
            _resIdToVarIdMap[netResourcePtr->id()] = netResourcePtr->associatedVar();
        }
    }

    for (auto resId : _networkFlowPtr->nonDisposableSpecialResourceIds())
        _rcspGraphPtr->nonDisposableBinaryResIds.push_back(resId);

    for (const auto & triple : _networkFlowPtr->permanentRyanFosterConstraints())
        _rcspGraphPtr->permanentRyanFosterConstraints.push_back(triple);

    _rcspGraphPtr->elemSetsDistanceMatrix = _networkFlowPtr->elemSetsDistanceMatrix();

    int maxVertexId = -1;
    for (lemon::ListDigraph::NodeIt lemonVertex(_networkFlowPtr->digraph()); lemonVertex != lemon::INVALID;
         ++lemonVertex)
    {
        NetworkVertex * netVertPtr = _networkFlowPtr->netVertexPtr(lemonVertex);
        if (netVertPtr->isFake())
            continue;

        _rcspGraphPtr->vertices.emplace_back(netVertPtr->id());
        bcp_rcsp::VertexData & vertex = _rcspGraphPtr->vertices.back();
        if (maxVertexId < vertex.id)
            maxVertexId = vertex.id;

        std::set<int> elemSetIds;
        for (auto * elemSetPtr : netVertPtr->elementaritySetPts())
            elemSetIds.insert(elemSetPtr->id());
        vertex.elemSetIds = std::vector<int>(elemSetIds.begin(), elemSetIds.end());

        std::set<int> enumElemSetIds;
        for (auto * elemSetPtr : netVertPtr->enumerationOnlyElementaritySetPts())
            enumElemSetIds.insert(elemSetPtr->id());
        vertex.enumElemSetIds = std::vector<int>(enumElemSetIds.begin(), enumElemSetIds.end());

        std::set<int> packSetIds;
        for (auto * packSetPtr : netVertPtr->packingSetPts())
            packSetIds.insert(packSetPtr->id());
        vertex.packSetIds = std::vector<int>(packSetIds.begin(), packSetIds.end());

        std::set<int> covSetIds;
        for (auto * covSetPtr : netVertPtr->coveringSetPts())
            covSetIds.insert(covSetPtr->id());
        vertex.covSetIds = std::vector<int>(covSetIds.begin(), covSetIds.end());

        vertex.name = netVertPtr->name();

        for (auto * resPtr : _networkFlowPtr->sideResourceConstrList())
        {
            int resId = resPtr->id();
            vertex.accumResConsumptionLB[resId] = resPtr->vertexConsumptionLB(lemonVertex);
            vertex.accumResConsumptionUB[resId] = resPtr->vertexConsumptionUB(lemonVertex);
            if (resPtr->selection())
                vertex.incompatibleValuesForSelection[resId] = resPtr->incompatibleValuesForSelection(lemonVertex);
        }

        for (auto & pair : netVertPtr->specResConsumptionBounds())
            vertex.accumBinResConsBounds[pair.first] = pair.second;

        for (auto * elemSetPtr : netVertPtr->inMemoryOfElemSets())
            vertex.ngNeighbourhood.push_back(elemSetPtr->id());
    }

    int maxArcId = -1;
    for (lemon::ListDigraph::ArcIt lemonArc(_networkFlowPtr->digraph()); lemonArc != lemon::INVALID; ++lemonArc)
    {
        NetworkArc * netArcPtr = _networkFlowPtr->netArcPtr(lemonArc);
        if (netArcPtr->isFake())
            continue;
        _rcspGraphPtr->arcs.emplace_back(netArcPtr->id());
        bcp_rcsp::ArcData & arc = _rcspGraphPtr->arcs.back();
        if (maxArcId < arc.id)
            maxArcId = arc.id;
        arc.tailVertexId = netArcPtr->tailVertexPtr()->id();
        arc.headVertexId = netArcPtr->headVertexPtr()->id();

        std::set<int> elemSetIds;
        for (auto * elemSetPtr : netArcPtr->elementaritySetPts())
            elemSetIds.insert(elemSetPtr->id());
        arc.elemSetIds = std::vector<int>(elemSetIds.begin(), elemSetIds.end());

        std::set<int> enumElemSetIds;
        for (auto * elemSetPtr : netArcPtr->enumerationOnlyElementaritySetPts())
            enumElemSetIds.insert(elemSetPtr->id());
        arc.enumElemSetIds = std::vector<int>(enumElemSetIds.begin(), enumElemSetIds.end());

        std::set<int> packSetIds;
        for (auto * packSetPtr : netArcPtr->packingSetPts())
            packSetIds.insert(packSetPtr->id());
        arc.packSetIds = std::vector<int>(packSetIds.begin(), packSetIds.end());

        std::set<int> covSetIds;
        for (auto * covSetPtr : netArcPtr->coveringSetPts())
            covSetIds.insert(covSetPtr->id());
        arc.covSetIds = std::vector<int>(covSetIds.begin(), covSetIds.end());

        arc.pureCost = netArcPtr->cost();

        for (auto * resPtr : _networkFlowPtr->sideResourceConstrList())
        {
            int resId = resPtr->id();
            arc.accumResConsumptionLB[resId] = resPtr->arcConsumptionLB(lemonArc);
            arc.accumResConsumptionUB[resId] = resPtr->arcConsumptionUB(lemonArc);
            arc.resConsumption[resId] = resPtr->arcConsumption(lemonArc);
            if (resPtr->dependent() || resPtr->base())
              arc.resConsumptionPWLfunction[resId] = resPtr->arcConsumptionPWLfunction(lemonArc);
            if (resPtr->dependent())
            {
              arc.depResWaitingTimeCoeff[resId] = resPtr->arcDepResWaitingTimeCoeff(lemonArc);
              arc.accumResConsPWLlowerBoundFunction[resId] = resPtr->arcConsumptionPWLlowerBoundFunction(lemonArc);
              arc.accumResConsPWLupperBoundFunction[resId] = resPtr->arcConsumptionPWLupperBoundFunction(lemonArc);
            }
            if (resPtr->selection())
                arc.incompatibleValuesForSelection[resId] = resPtr->incompatibleValuesForSelection(lemonArc);
        }

        for (auto & pair : netArcPtr->specResConsumption())
            arc.binaryResConsumption[pair.first] = pair.second;

        for (auto * elemSetPtr : netArcPtr->inMemoryOfElemSets())
            arc.ngNeighbourhood.push_back(elemSetPtr->id());

        auto mappingIt = netArcPtr->varToCoeffMaps().begin();
        for (auto & pair : *mappingIt)
            arc.varIdToCostAndCoeff[_instVarToIdMap[pair.first]] = std::make_pair(pair.first->costrhs(),pair.second);

        for (++mappingIt; mappingIt != netArcPtr->varToCoeffMaps().end(); ++mappingIt)
        {
            arc.alternativeMappings.emplace_back();
            for (auto & pair : *mappingIt)
                arc.alternativeMappings.back()[_instVarToIdMap[pair.first]]
                    = std::make_pair(pair.first->costrhs(),pair.second);
        }

        arc.name = netArcPtr->name();
    }

    if (!_rcspGraphPtr->generatePreprocessedInfo())
    {
        delete _rcspGraphPtr;
        _rcspGraphPtr = nullptr;
        return false;
    }
    return true;
}
#endif //BCP_RCSP_IS_FOUND
