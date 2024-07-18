/**
 *
 * This file bcModelFormulationC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <utility>

#include "bcModelingLanguageC.hpp"
#include "bcModelNetworkFlow.hpp"
#include "bcProbConfigC.hpp"
#include "bcModelC.hpp"
#include "bcFormC.hpp"

BcFormulation::BcFormulation(ProbConfig * pConfPointer) :
    _probConfPtr(pConfPointer)
{
}

BcFormulation::BcFormulation(const BcFormulation & that) :
    _probConfPtr(that._probConfPtr)
{
}

BcFormulation::BcFormulation(BcModel & model, const std::string & name) :
    _probConfPtr(((Model *)model)->createProbConf(name, MultiIndex()))
{
}


BcFormulation::~BcFormulation()
{
}

bool BcFormulation::isDefined() const
{
    return _probConfPtr != NULL;
}

bool BcFormulation::isMaster() const
{
    return isDefined() && _probConfPtr->isTypeOf(PcId::MasterMask);
}

bool BcFormulation::isColGenSp() const
{
    return isDefined() && _probConfPtr->isTypeOf(PcId::ColGenSpConfMask);
}

BcFormulation BcFormulation::master() const
{
  if (!isColGenSp())
      return BcFormulation(NULL);
  return BcFormulation(_probConfPtr->mastConfPtr());
}

const std::list<BcFormulation> & BcFormulation::colGenSubProblemList() const
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }
    return _probConfPtr->colGenSubProbFormList();
}

ProbConfig * BcFormulation::probConfPtr() const
{
  return _probConfPtr;
}

const MultiIndex & BcFormulation::id() const
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }
    return _probConfPtr->id();
}

const std::string & BcFormulation::name() const
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }
    return _probConfPtr->probPtr()->name();
}

void BcFormulation
     ::getRyanAndFosterBranchingConstrsList(std::list<BcRyanAndFosterBranchConstr> & ryanAndFosterBranchConstrList) const
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    _probConfPtr->probPtr()->getActiveRyanAndFosterBranchingConstraints(ryanAndFosterBranchConstrList);
}

void BcFormulation::getActiveSoftConflictCutsList(std::list<BcSoftConflictsCut> & softConflictsCutsList) const
{
    MasterConf * mastConfPtr = NULL;
    if (isMaster())
        mastConfPtr = static_cast<MasterConf *>(_probConfPtr);
    else
        mastConfPtr = _probConfPtr->mastConfPtr();

    GenericCutConstr * _genSoftConflictCutConstrPtr = mastConfPtr->getGenericCutConstr("TLCC");
    if (_genSoftConflictCutConstrPtr == NULL)
        return;

    const IndexCell2InstancConstrPtrMap & constrPtrMap = _genSoftConflictCutConstrPtr->indexCell2InstancConstrPtrMap();
    for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
        if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
            && it->second->isTypeOf(VcId::SoftConflictsCutConstrMask))
        {
            SoftConflictsCut * cutPtr = static_cast<SoftConflictsCut *>(it->second);
            softConflictsCutsList.push_back(BcSoftConflictsCut(cutPtr));
        }
}

void BcFormulation::getColumnsInPrimalLpSolution(std::vector<std::pair<double, BcSolution>> & columnsInSol) const
{
    columnsInSol.clear();
    MasterConf * mastConfPtr = nullptr;
    if (isMaster())
        mastConfPtr = static_cast<MasterConf *>(_probConfPtr);
    else
        mastConfPtr = _probConfPtr->mastConfPtr();

    if (mastConfPtr == nullptr || mastConfPtr->probPtr() == nullptr)
        return;

    for (auto * varPtr : mastConfPtr->probPtr()->inPrimalLpSol())
        if (varPtr->isTypeOf(VcId::MastColumnMask))
        {
            auto * colPtr = static_cast<MastColumn *>(varPtr);
            columnsInSol.emplace_back(colPtr->val(), BcSolution(colPtr->spSol()));
        }
}

void BcFormulation::getCustomNonLinearActiveCutsList(std::list<BcCustomNonLinearCut> & cutsList) const
{
    MasterConf * mastConfPtr = NULL;
    if (isMaster())
        mastConfPtr = static_cast<MasterConf *>(_probConfPtr);
    else
        mastConfPtr = _probConfPtr->mastConfPtr();

    std::set<GenericCutConstr *, DynamicGenConstrSort>::const_iterator genCutPtrIt;
    for (genCutPtrIt = mastConfPtr->candidateCutGenericConstr().begin();
         genCutPtrIt != mastConfPtr->candidateCutGenericConstr().end(); ++genCutPtrIt)
    {
        const GenericCustomNonLinearCutConstr * genCustNonLinCutConstrPtr
                = dynamic_cast<GenericCustomNonLinearCutConstr *>(*genCutPtrIt);
        if (genCustNonLinCutConstrPtr == NULL)
            continue;

        const IndexCell2InstancConstrPtrMap & constrPtrMap = (*genCutPtrIt)->indexCell2InstancConstrPtrMap();
        for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin();
             it != constrPtrMap.end(); ++it)
            if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
                && it->second->isTypeOf(VcId::CustomNonLinearCutConstrMask))
            {
                CustomNonLinearCut * cutPtr = static_cast<CustomNonLinearCut *>(it->second);
                cutsList.push_back(BcCustomNonLinearCut(cutPtr));
            }
    }
}

void BcFormulation::resetConstrRHS()
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    _probConfPtr->probPtr()->hardResetConstraintRHS('s');
}

void BcFormulation::resetObjective()
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    _probConfPtr->probPtr()->hardResetObjective('s');
}

void BcFormulation::update()
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    _probConfPtr->prepareProbConfig();

    _probConfPtr->probPtr()->updateProblem();

    _probConfPtr->resetDualIncBound();
    _probConfPtr->resetPrimalIncBound();
}

void BcFormulation::priorityLevel(const double & prLevel)
{
  if (!isDefined())
  {
    if (printL(6))
	  std::cout << "BaPCod info :  Model BcFormulation == NULL" << std::endl;
    return;
  }

  if (isColGenSp())
    static_cast<ColGenSpConf *>(_probConfPtr)->priorityLevel(prLevel);
}

double BcFormulation::zeroReducedCostThreshold() const
{
    if (isColGenSp())
        return - static_cast<ColGenSpConf *>(_probConfPtr)->fixedDualCost();

    return 0.0;
}

bool BcFormulation::rollbackPointSavedStatus() const
{
    if (isColGenSp())
        return static_cast<ColGenSpConf *>(_probConfPtr)->getRollbackPointSavedStatus();

    return false;
}

bool BcFormulation::enumeratedStatus() const
{
  if (_probConfPtr == NULL)
    return false;
  return _probConfPtr->enumeratedStatus();
}

const BcNetwork BcFormulation::network()
{
  if (!isDefined() || (_probConfPtr->networkFlowPtr() == NULL))
    {
      std::cerr << "BapCod error in BcFormulation::network(): network is not defined" << std::endl;
      exit(1);
    }
  return BcNetwork(_probConfPtr->networkFlowPtr());
}

#ifdef BCP_RCSP_IS_FOUND
const bcp_rcsp::GraphData * BcFormulation::RCSPGraph()
{
    if (!isDefined())
    {
        std::cerr << "BapCod error in BcFormulation::RCSPGraph(): network is not defined" << std::endl;
        exit(1);
    }
    return _probConfPtr->rcspGraphPtr();
}
#endif

void BcFormulation::setFixedCost(const double & value)
{
    if (isColGenSp())
        static_cast<ColGenSpConf *>(_probConfPtr)->fixedCost(value);
}

const BcFormulation & BcFormulation::attach(BcSolverOracleFunctor * oraclePtr)
{
    if (isColGenSp())
        static_cast<ColGenSpConf *>(_probConfPtr)->probPtr()->solverOracleFunctorPtr(oraclePtr);
    return *this;
}

const BcFormulation & BcFormulation::attach(BcFracSolBasedHeuristicFunctor * functorPtr)
{
    if (isMaster())
        static_cast<MasterConf *>(_probConfPtr)->probPtr()->fracSolBasedHeuristicFunctorPtr(functorPtr);
    return *this;
}

const BcFormulation & BcFormulation::attach(BcMasterHeuristicFunctor * functorPtr)
{
    if (isMaster())
        static_cast<MasterConf *>(_probConfPtr)->probPtr()->masterHeuristicFunctorPtr(functorPtr);
    return *this;
}

const BcFormulation & BcFormulation::attach(BcEnumSolBasedHeuristicFunctor * functorPtr)
{
    if (isMaster())
        static_cast<MasterConf *>(_probConfPtr)->probPtr()->enumSolBasedHeuristicFunctorPtr(functorPtr);
    return *this;
}

const BcFormulation & BcFormulation::attach(BcDivingFixingFunctor * functorPtr)
{
    if (isMaster())
        static_cast<MasterConf *>(_probConfPtr)->probPtr()->divingFixingFunctorPtr(functorPtr);
    return *this;
}

const BcFormulation & BcFormulation::operator<=(const double & ubOnNbOfUsedSol)
{
  if (!isDefined())
  {
    if (printL(6))
      std::cout << "BaPCod info :  Model BcFormulation == NULL" << std::endl;
  }
  else
    _probConfPtr->upperBoundPtr(new Double(ubOnNbOfUsedSol));
  return *this;
}

const BcFormulation & BcFormulation::operator>=(const double & lbOnNbOfUsedSol)
{
    if (!isDefined())
    {
        if (printL(6))
            std::cout << "BaPCod info :  Model BcFormulation == NULL" << std::endl;
    }
    else
        _probConfPtr->lowerBoundPtr(new Double(lbOnNbOfUsedSol));
    return *this;
}

const BcFormulation & BcFormulation::operator==(const double & nbOfUsedSol)
{
    if (!isDefined())
    {
        if (printL(6))
            std::cout << "BaPCod info :  Model BcFormulation == NULL" << std::endl;
    }
    else
    {
        _probConfPtr->lowerBoundPtr(new Double(nbOfUsedSol));
        _probConfPtr->upperBoundPtr(new Double(nbOfUsedSol));
    }
    return *this;
}

BcSolution BcFormulation::solve(bool printForm, bool printOutput)
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    _probConfPtr->prepareProbConfig();
    if(printForm)
    {
        _probConfPtr->printForm();
    }
    return BcSolution(_probConfPtr->solvePC(printOutput));
}

SolutionStatus BcFormulation::getStatus()
{
  if (!isDefined() || !_probConfPtr->isPrepared())
    return SolutionStatus::UnSolved;
  
  return _probConfPtr->probPtr()->formulationPtr()->status();
}

std::ostream& BcFormulation::print(std::ostream& os) const
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }

    return  _probConfPtr->printForm(os);
}

std::ostream& operator<<(std::ostream& os, const BcFormulation & that)
{
  if (!that.isDefined())
    return os << "FormulationPtr::operator<<: undefined FormulationPtr" << std::endl;

  if (printL(5))
    std::cout << "operator<<(BcFormulation) " << std::endl;

  that._probConfPtr->prepareProbConfig();

  return that._probConfPtr->print(os);
}

bool BcFormulation::currentNodeIsRoot()
{
    if (!isDefined())
    {
        std::cerr << "BaPCod error : Model BcFormulation == NULL" << std::endl;
        exit(1);
    }
    if ((_probConfPtr->probPtr() == NULL) && (_probConfPtr->probPtr()->curNodePtr() == NULL))
    {
        std::cerr << "BaPCod error : cannot get the current node pointer in BcFormulation::currentNodeIsRoot"
                  << std::endl;
        exit(1);
    }
    return _probConfPtr->probPtr()->curNodePtr()->isRoot();
}

bool BcFormulation::debugSolutionIsValidAtThisNode()
{
    ProbConfig * probConfPtr = _probConfPtr;
    if (!isMaster())
    {
        if (!isColGenSp())
            return false;
        else
            probConfPtr = master().probConfPtr();
    }

    MasterConf * mastConfPtr = static_cast<MasterConf *>(probConfPtr);

    Solution * debugSolPtr = mastConfPtr->getDebugSolution();
    const Node * curNodePtr = probConfPtr->probPtr()->curNodePtr();

    if (curNodePtr != nullptr)
    {
        if ((debugSolPtr == NULL) || !curNodePtr->debugSolutionAtThisNode())
            return false;

        /// if the current best solution value is equal to the debug solution value, then the debug solution
        /// may be cut off, so we should not verify it
        if (curNodePtr->nodeIncIpPrimalBound() < debugSolPtr->cost() + Double::precision)
            return false;
    }

    return true;
}

void BcFormulation::initializeWithColumns(BcSolution & sol, const bool activeColumns)
{
    if (!isMaster())
        return;

    _probConfPtr->prepareProbConfig();

    Solution * curSolPtr = sol._solutionPtr;
    Solution * masterSolPtr = new Solution(_probConfPtr);

    while (curSolPtr != NULL)
    {
        if (!curSolPtr->solVarValMap().empty())
        {
            ProbConfig * pcPtr = curSolPtr->probConfPtr();
            pcPtr->recordSubproblemSolution(curSolPtr, false, 2, masterSolPtr);
        }
        curSolPtr = curSolPtr->_nextSolPtr;
    }

    if (activeColumns)
        _probConfPtr->recordInitialActiveSetOfColumns(masterSolPtr);
    else
        _probConfPtr->recordInitialInactiveSetOfColumns(masterSolPtr);
}

void BcFormulation::initializeWithSolution(BcSolution & sol)
{
    if (!isMaster())
        return;

    _probConfPtr->prepareProbConfig();

    Solution * curSolPtr = sol._solutionPtr;
    Solution * probConfSolPtr = new Solution(_probConfPtr);

    while (curSolPtr != NULL)
    {
        if (!curSolPtr->solVarValMap().empty())
        {
            ProbConfig * pcPtr = curSolPtr->probConfPtr();
            /// changed by Ruslan to make Split variant work
            pcPtr->recordSubproblemSolution(curSolPtr, false, 2, probConfSolPtr);
        }
        curSolPtr = curSolPtr->_nextSolPtr;
    }

    _probConfPtr->updatePrimalIncSolution(probConfSolPtr);
    delete probConfSolPtr;
}

const BcFormulation & BcFormulation::operator+=(BcSolution & sol)
{
    initializeWithSolution(sol);
    return *this;
}

BcFormulationArray::BcFormulationArray(BcModel & modelPointer, std::string name) :
    _modelPtr(modelPointer), _name(std::move(name)), _curForm(NULL)
{
}

BcFormIndex BcFormulationArray::operator[](const int & index)
{
  return BcFormIndex(this, MultiIndex(index));
}

BcFormulation & BcFormulationArray::getElement(const MultiIndex & multiIndexId)
{
  if (_curForm.isDefined())
    {
      if (_curForm.id() == multiIndexId)
	    return _curForm;

      _curForm._probConfPtr = NULL;
    }
  return _curForm;
}

BcFormulation & BcFormulationArray::createElement(const MultiIndex & multiIndexId)
{
  if (_curForm._probConfPtr != NULL)
    {
      if (_curForm.id() == multiIndexId)
	    return _curForm;
    }

  _curForm = BcFormulation(_modelPtr->createProbConf(_name, multiIndexId));

  return _curForm;
}

BcFormulationArray::~BcFormulationArray()
{
}


BcFormIndex::BcFormIndex(BcFormulationArray * modFormPtr, const MultiIndex & indexArray) :
    std::pair<BcFormulationArray *, MultiIndex>(modFormPtr, indexArray)
{
}

BcFormIndex::BcFormIndex(BcFormulationArray * modFormPtr) :
    std::pair<BcFormulationArray *, MultiIndex>(modFormPtr, MultiIndex(0))
{
}

BcFormIndex::~BcFormIndex()
{
}

BcFormIndex & BcFormIndex::operator[](const int & index)
{
  second.operator+=(index);
  return *this;
}

bool BcMasterHeuristicFunctor::operator() (IN_ BcFormulation spPtr,
                                           IN_ std::vector<std::pair<BcSolution, double> > & colsInFixedSolution,
                                           OUT_ BcSolution & primalSol)
{
  return false;
}

bool BcFracSolBasedHeuristicFunctor::operator() (IN_ const BcFracSolBasedHeuristicFunctorInput & in,
                                                 OUT_ BcSolution & primalSol)
{
  return false;
}

void BcDivingFixingFunctor::colCutGenTerminated(IN_ BcFormulation spPtr,
                                                IN_ std::vector<std::pair<BcSolution, double> > & colsInFixedSolution,
                                                IN_ std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                                OUT_ BcSolution & primalSol)
{
}

void BcDivingFixingFunctor::operator() (const BcDivingFixingFunctorInput & in, BcDivingFixingFunctorOutput & out)
{
}

bool BcEnumSolBasedHeuristicFunctor::operator()(IN_ BcFormulation spPtr,
                                                IN_ double incPrimalIpBound,
                                                IN_ const std::vector<BcSolution> & enumSolution,
                                                OUT_ std::vector<std::pair<int, double> > & solution,
                                                OUT_ BcSolution & addSolution)
{
    return false;
}

bool BcSolverOracleFunctor::setupNode(BcFormulation spPtr, const BcSolverOracleInfo * infoPtr)
{
  return false;
}

BcSolverOracleInfo * BcSolverOracleFunctor::recordSolverOracleInfo(const BcFormulation spPtr)
{
  return NULL;
}

int BcSolverOracleFunctor::getMessageIdToCutGeneration() const
{
  return PricingSolverCutsMessage::noMessage;
}

int BcSolverOracleFunctor::getNumberOfEnumeratedSolutions() const
{
  return -1; /// by default, the subproblem is not enumerated
}

bool BcSolverOracleFunctor::getDebugSolution(BcFormulation spPtr, BcSolution & primalSol)
{
    return true;
}

bool BcSolverOracleFunctor::setDebugSolution(const std::vector<std::vector<int> > & ids, bool vertexBased)
{
    return false;
}

void BcSolverOracleFunctor::reducedCostFixingAndEnumeration(IN_ BcFormulation spPtr,
                                                            IN_ const int & enumerationMode,
                                                            IN_ const double & threshold)
{
}

bool BcSolverOracleFunctor::getEnumeratedStatus() const
{
    return false;
}


void BcSolverOracleFunctor::checkEnumeratedSolutions(IN_ BcFormulation spPtr,
                                                     IN_ const std::vector<Solution *> & solPts,
                                                     OUT_ std::vector<bool> & solIsEnumerated)
{

}

void BcSolverOracleFunctor::getEnumeratedSolutions(IN_ BcFormulation spPtr,
                                                   IN_ const int & maxNumberOfSolutions,
                                                   OUT_ BcSolution & enumeratedSol,
                                                   OUT_ std::vector<double> & reducedCosts)
{
}


bool BcSolverOracleFunctor::isProperSolution(IN_ BcFormulation spPtr,
                                             IN_ BcSolution & solution)
{
  return true;
}

bool BcSolverOracleFunctor::solSatisfiesCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                                            IN_ const BcSolution & solution)
{
  return true;
}

bool BcSolverOracleFunctor
     ::improveCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                  IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                  IN_ const bool & masterConverged)
{
  return false;
}

bool BcSolverOracleFunctor
     ::drawPrimalSolutionToDotFile(IN_ BcFormulation spPtr,
                                   IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                   IN_ const std::string & filename)
{
    return false;
}

void BcSolverOracleFunctor::columnGenerationTerminated(IN_ BcFormulation spPtr, bool afterRedCostFixing,
                                                       int nodeOrder, int nodeDepth, int cutSepRound, double dualBound,
                                                       double elapsedTime, bool masterConverged)
{
    return;
}


bool BcSolverOracleFunctor::lightenCurrentSpRelaxation(IN_ BcFormulation spPtr, const int & masterConverged,
                                                       const int & callMode)
{
  return false;
}

bool BcSolverOracleFunctor::operator() (IN_ BcFormulation spPtr,
                                        IN_ int colGenPhase,
                                        OUT_ double & objVal,
                                        OUT_ double & dualBound,
                                        OUT_ BcSolution & primalSol)
{
    if (printL(-1))
        std::cout << "BaPCod WARNING : BcSolverOracleFunctor::operator() SHOULD NOT BE CALLED" << std::endl;
    return false;
}

/// for backward compatibility
bool BcSolverOracleFunctor::operator()(IN_ BcFormulation spPtr,
				                       IN_ OUT_ double & objVal,
                                       IN_ OUT_ double & primalBound,
				                       IN_ OUT_ double & dualBound,
                                       OUT_ BcSolution & primalSol,
                                       OUT_ BcDualSolution & dualSol,
                                       IN_ OUT_ int & phaseOfStageApproach)
{
    bool returnValue = operator()(spPtr, phaseOfStageApproach, objVal, dualBound, primalSol);
    primalBound = objVal;
    return returnValue;
}
