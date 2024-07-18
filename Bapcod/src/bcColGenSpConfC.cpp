/**
 *
 * This file bcColGenSpConfC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcColGenSpConfC.hpp"
#include "bcFormC.hpp"
#include "bcMasterConfC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcModelColGenSpC.hpp"
#include "bcModelC.hpp"

/**
 * Generic code
 */
using namespace std;
using namespace VcIndexStatus;

ColGenSpConf::ColGenSpConf(const std::string genericName,
			               const IndexCell & id,
			               MasterConf * mastConfPtr,
			               const Double & fixedCost,
			               const bool & implicitlyFixCardinality,
			               Double * upperBoundPtr,
			               Double * lowerBoundPtr,
			               const Double & defaultDualVal4UbConstr,
			               const Double & defaultDualVal4LbConstr,
			               Problem * problemPtr):
  ProbConfig(ProbConfig::colGenSp,
	         mastConfPtr->modelPtr(),
	         genericName,
	         id,
    	     Bound::infPrimalBound(mastConfPtr->modelPtr()->objectiveSense()),
	         Bound::infDualBound(mastConfPtr->modelPtr()->objectiveSense()),
	         problemPtr),
  _objStatus(mastConfPtr->modelPtr()->objectiveSense()),
  _mastConfPtr(mastConfPtr),
  _misColPtr(NULL),
  _spConfHasClassInducedSpVarBounds(false),
  _implicitlyFixCardinality(implicitlyFixCardinality),
/// @todo replace the following two by simple double
  _upperBoundPtr(upperBoundPtr),
  _lowerBoundPtr(lowerBoundPtr),
  _defaultDualVal4UbConstr(defaultDualVal4UbConstr),
  _defaultDualVal4LbConstr(defaultDualVal4LbConstr),
  _lowerBoundMastConstrPtr(NULL),
  _upperBoundMastConstrPtr(NULL),
  _fixedCost(fixedCost),
  _spRootReducedCost(0),
  _toSplit(param().SplitColIntoDissagregateSpVarInHeadIn),
  _rollbackPointSaved(false),
  _fixedDualCost(0),
  _target(BapcodInfinity, _objStatus),
  _mult(1),
  _priorityLevel(0.5),
  _curBestMastColumnPtr(NULL),
  _incBestMastColumnPtr(NULL)
{
  if (mastConfPtr != NULL)
    mastConfPtr->insertColGenSpConf(this);

  return;
}



ColGenSpConf::~ColGenSpConf()
{
  if (_upperBoundPtr != NULL)
    delete _upperBoundPtr;

  if (_lowerBoundPtr != NULL)
    delete _lowerBoundPtr;

  if (_primalSolPtr != NULL)
    delete _primalSolPtr;

  _iVarPts.clear();

  return;
}

Constraint * ColGenSpConf::castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ColGenSpConf::castAndAddConstraint() should not be called");
  return constrPtr;
}

InstanciatedConstr * ColGenSpConf::castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ColGenSpConf::castAndAddConstraint() should not be called");
  return iconstrPtr;
}

Variable * ColGenSpConf::castAndAddVariable(Variable * varPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ColGenSpConf::castAndAddVariable() should not be called");
  return varPtr;
}

InstanciatedVar * ColGenSpConf::castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately)
{
  SubProbVariable * spVarPtr = NULL;
  if (ivarPtr->isTypeOf(VcId::SubProbVariableMask))
    {
      spVarPtr = (SubProbVariable *)(ivarPtr);
    }

  if( printL(7))
    std::cout << " sp InstanciatedVar " << ivarPtr;

  if (spVarPtr == NULL) /// cast the variable
    {
      spVarPtr = new SubProbVariable(ivarPtr, _mastConfPtr);
      if (ivarPtr) delete (ivarPtr);
      else
	{
	  std::cerr << " ERROR: Item *it is NULL\n";
	  abort();
	}
    }

  if (printL(6))
    std::cout << "cast instantiated var of sp " << name()
	      << " into SubProbVariable : " << spVarPtr->name()
	      << std::endl;

  if( printL(7))
    std::cout << " cast variables to subprob var " << spVarPtr;

  if (insertImmediately)
    {
      probPtr()->addVar(spVarPtr);
    }
  else
    {
      if( printL(6))
	     std::cout << " cast variables to subprob var : add in _pcVarPtrList " << spVarPtr->name() << std::endl;
      _pcVarPtrList.push_back(spVarPtr);
    }

  return spVarPtr;
}

void ColGenSpConf::prepareProbConfig()
{
  if (!progStatus().doRun()) return;

  if ((param().ovfSolMode().status() ==  SolutionMethod::none)
      && (param().colGenSubProbSolMode().status() ==  SolutionMethod::none))
    return;

  if (_isPrepared) return;
    _isPrepared = true;

  if (printL(5))
    std::cout << " ColGenSpConf::prepareProbConfig() " << name() << " _name2GenericVarPtrMap.size() "
	          << _name2GenericVarPtrMap.size() << std::endl;

  bool branchOnNumberOfColumns(false);
  double maxCompBoundSetBranchingPL(-1);
  for (std::map < std::string, GenericVar * >::const_iterator it = _name2GenericVarPtrMap.begin();
      it != _name2GenericVarPtrMap.end(); ++it)
    {
      if (printL(5))
	    std::cout << " setup col gen BranchingConstr for " << it->first << std::endl;

      it->second->setupGenericBranchingConstr();

      if (it->second->compBoundSetGenBranchConstrPtr() != NULL)
        {
          if (maxCompBoundSetBranchingPL < it->second->compBoundSetGenBranchConstrPtr()->priorityLevel())
            maxCompBoundSetBranchingPL = it->second->compBoundSetGenBranchConstrPtr()->priorityLevel();
          if (!it->second->indexCell2InstancVarPtrMap().empty())
            branchOnNumberOfColumns = true;
        }
    }

  /// if component set branching is not used for variables of this subproblem,
  /// we remove component set branching for the default generic var
  /// otherwise, we give it the maximum priority, as branching on the columns number is
  /// done only for component set branching of the default generic var (to avoid repetition)
  if (defaultGenericVarPtr()->compBoundSetGenBranchConstrPtr() != NULL)
    {
      if (branchOnNumberOfColumns)
        defaultGenericVarPtr()->compBoundSetGenBranchConstrPtr()->priorityLevel(maxCompBoundSetBranchingPL + 1);
      else
        _mastConfPtr->candidateBranchingGenericConstr().erase(defaultGenericVarPtr()->compBoundSetGenBranchConstrPtr());
    }

  probPtr()->defineFormulation();

  _pcConstrPtrList.insert(_pcConstrPtrList.end(), _iConstrPts.begin(), _iConstrPts.end());

  std::list<InstanciatedVar *> tempListOfVar2Upcast(_iVarPts.begin(), _iVarPts.end());

  /// Cast variables to subprob var
  for (std::list<InstanciatedVar *>::const_iterator it = tempListOfVar2Upcast.begin(); it != tempListOfVar2Upcast.end();
       ++it)
    {
      castAndAddVariable(*it);
    }

  probPtr()->addVarSet(_pcVarPtrList, 1, 0);
  probPtr()->addConstrSet(_pcConstrPtrList, 1, 0);
  probPtr()->buildProblem();
  
  if (_networkFlowPtr != NULL)
  {
      _networkFlowPtr->generateArcInfo();
#ifdef BCP_RCSP_IS_FOUND
      fillRCSPGraph();
#endif //BCP_RCSP_IS_FOUND
  }

  /// should be done after _networkFlowPtr->generateArcInfo();
  if (!probPtr()->prepareSolverOracleFunctor())
      progStatus().setStat(ProgStatus::terminate);

  return;
}


Double * ColGenSpConf::upperBoundPtr() const
{
  return(_upperBoundPtr);
}

Double * ColGenSpConf::lowerBoundPtr() const
{
  return(_lowerBoundPtr);
}

void ColGenSpConf::upperBoundPtr(Double * ubPtr)
{
  if (_upperBoundPtr != NULL)
    delete _upperBoundPtr;
  _upperBoundPtr = ubPtr;
}

void ColGenSpConf::lowerBoundPtr(Double * lbPtr)
{
  if (_lowerBoundPtr != NULL)
    delete _lowerBoundPtr;
  _lowerBoundPtr = lbPtr;
}

const Double & ColGenSpConf::defaultDualVal4UbConstr() const
{
  return _defaultDualVal4UbConstr;
}

const Double & ColGenSpConf::defaultDualVal4LbConstr() const
{
  return _defaultDualVal4LbConstr;
}


void  ColGenSpConf::lowerBoundMastConstrPtr(InstMastConvexityConstr * mccPtr)
{
  _lowerBoundMastConstrPtr = mccPtr;
  return;
}

void  ColGenSpConf::upperBoundMastConstrPtr(InstMastConvexityConstr * mccPtr)
{
  _upperBoundMastConstrPtr = mccPtr;
  return;
}
InstMastConvexityConstr * ColGenSpConf::lowerBoundMastConstrPtr()  const
{

  return _lowerBoundMastConstrPtr;
}

InstMastConvexityConstr * ColGenSpConf::upperBoundMastConstrPtr()  const
{
  return _upperBoundMastConstrPtr;
}


MissingColumn * ColGenSpConf::misColPtr() const
{
  return _misColPtr;
}

void ColGenSpConf::misColPtr(MissingColumn * misCptr)
{
  _misColPtr = misCptr;
  return;
}

const bool & ColGenSpConf::spConfHasClassInducedSpVarBounds() const
{
  return _spConfHasClassInducedSpVarBounds;
}

void ColGenSpConf::spConfHasClassInducedSpVarBounds(const bool & f)
{
  _spConfHasClassInducedSpVarBounds = f;
}

MasterConf * ColGenSpConf::mastConfPtr() const
{
  return _mastConfPtr;
}

void ColGenSpConf::fixedCost(const Double & b)
{
  _fixedCost = b;
}

void ColGenSpConf::fixedDualCost(const Double & b)
{
  _fixedDualCost = b;
}

const Double & ColGenSpConf::spRootReducedCost() const
{
  return _spRootReducedCost;
}

const Double & ColGenSpConf::mult() const
{
  return _mult;
}

/// Reset min Cost to min Sp value and target according du dual master value
void ColGenSpConf::updateTarget(const bool currentlyPerformingPhaseI)
{
  const int printlevel = 5;

  /// Reset target
 _fixedDualCost = Bound(0, _objStatus);

 switch (_mastConfPtr->probPtr()->solMode().status())
   {
   case SolutionMethod::lpSolver:
   case SolutionMethod::mipSolver:
   case SolutionMethod::customSolver:
   case SolutionMethod::custom2mipSolver:
     {
       if (_lowerBoundMastConstrPtr != NULL)
	 {
	   _fixedDualCost += _lowerBoundMastConstrPtr->valOrSepPointVal();

	   if (printL(printlevel))
	     std:: cout << "ColGenSpConf::updateTarget()  after lowerBoundMastConstr = "
			<<  _fixedDualCost
			<< std::endl;
	 }

       if (_upperBoundMastConstrPtr != NULL)
	 {
	   _fixedDualCost += _upperBoundMastConstrPtr->valOrSepPointVal();

	   if (printL(printlevel))
	     std:: cout << "ColGenSpConf::updateTarget()  after upperBoundMastConstr = "
			<<  _fixedDualCost
			<< std::endl;
	 }

       break;
     }
   case SolutionMethod::none:
     break;
   case SolutionMethod::undefined:
     bapcodInit().check(true,
			"ColGenSpConf::updateTarget(: ERROR undefined solution method");
     break;
   }

 /// check if node can be pruned based on primal bound a ancestor Node dual bound

 _target = Bound(-_fixedDualCost, _objStatus);

 Double dualBoundAdjustmentForNonAccountedConstraints(0);
 /// add fixedCost;
 if (currentlyPerformingPhaseI)
   {
     dualBoundAdjustmentForNonAccountedConstraints -= fixedCost();
     if (printL(5))
       std:: cout << "ColGenSpConf::updateTarget()  dualBoundAdjustmentForNonAccountedConstraints = "
		  <<  dualBoundAdjustmentForNonAccountedConstraints
		  << std::endl;
   }

 Double contribMaxCompSetConstr(0);

 /// Check for other master constraints
 for ( ConstrPtrSet::const_iterator cPt = _mastConfPtr->probPtr()->inDualSol().begin();
       cPt != _mastConfPtr->probPtr()->inDualSol().end(); ++cPt)
   {
     if ((*cPt)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
       {
           CompSetInstMastBranchConstr * csIbcPtr = static_cast<CompSetInstMastBranchConstr *>(*cPt);
           if (csIbcPtr->ColGenSpConfPtr() != this)
               continue;
           if (contribMaxCompSetConstr < - csIbcPtr->valOrSepPointVal())
               contribMaxCompSetConstr =  - csIbcPtr->valOrSepPointVal();
           if (printL(5))
               std:: cout << "ColGenSpConf::updateTarget() after constr " << csIbcPtr->name()
                          << " contribMaxCompSetConstr = " <<  contribMaxCompSetConstr << std::endl;
           continue;
	}

     if (!(*cPt)->isTypeOf(VcId::InstMasterConstrMask))
 	    continue;

      if ((*cPt)->type() == 'S')
	    /// constraint is a Subproblem constraints generated in the master
	    /// by a row-and-column generation method
	    /// should be ignored
	    continue;

     /// Check non linear master constraints
     if ((*cPt)->isTypeOf(VcId::Base4NonLinearConstraintMask))
       {
	    dualBoundAdjustmentForNonAccountedConstraints = BapcodInfinity;
	    break;
       }
    }

  dualBoundAdjustmentForNonAccountedConstraints += contribMaxCompSetConstr;

  _target += dualBoundAdjustmentForNonAccountedConstraints;

  return;
}

/// Set var cost and bounds, target, minCost,
bool ColGenSpConf::updateConf(bool inPurePhaseOne)
{
  /**
      set problem configuration before each iteration of the column generation procedure and,
      within a column generation iteration, for each new pricing subproblem defined by branching constraints.
      Modification involve setting the objective function and new bounds on subproblem varaibles 
      (that call for new constraint checking).
   */

  Time start;

  probPtr()->resetSolution();
  probPtr()->clearRecordedSol();

  /// Test if problem infeasible
  if (probPtr()->probInfeasibleFlag())
    {
      if (printL(3))
        std::cout << "ColGenSpConf::updateConf():  infeasible colGenSP " << name() << std::endl;
      return(true);
    }

  if ((_upperBoundMastConstrPtr != NULL) && (_upperBoundMastConstrPtr->curRhs() <= 0))
    {
      if (printL(3))
        std::cout << "ColGenSpConf::updateConf(): no more solution to generate  from  colGenSP" << name() << std::endl;
      /// No need to attempt generating columns
      return(false);
    }

  /// Compute _minCost, _maxCost
  if (probPtr()->updateProbForColGen(inPurePhaseOne))
    return(true);

  probPtr()->setDualBound(probPtr()->minCost());
  probPtr()->setPrimalLpBound(target());


  if (printL(3))
    {
      std::cout << "ColGenSpConf::updateConf(): colGenSP name = " << name()
                << ", bestPossibleReducedCost==dual_bound = " << probPtr()->dualBound()
                << ", minCost= " << probPtr()->minCost() << ", target = " << target()
                << ", CutOffValue==primal_Bound==target = " << probPtr()->primalBound()
                << ", (if dual_bound < primal_Bound, solve SP: otherwise no hope to get a neg red cost column)"
                << std::endl;
    }

  bapcodInit().statistics().incrTimer("bcTimeSpUpdateProb", start.getElapsedTime_dbl());

  return(false);
}

void ColGenSpConf::getSolPC(bool inPurePhaseOne)
{
  _primalSolPtr = probPtr()->retrieveCurPrimalLpSol();

  if (printL(5))
      std::cout << " ColGenSpConf::getSolPC() found sol after solvePC primalSol = " << *_primalSolPtr << std::endl;

  /// prepare constraints generated in the oracle for insertion
  int insertionLevel = 1;
  for (ConstrPtrSet::const_iterator sPtr = probPtr()->inDualSol().begin(); sPtr != probPtr()->inDualSol().end();
       sPtr++)
  {
      (*sPtr)->problemPtr(_probPtr);
      checkConstraint4Insertion(*sPtr, insertionLevel);
  }

  _curBestMastColumnPtr = recordSubproblemSolution(_primalSolPtr, inPurePhaseOne, 1);

  /// added by Ruslan : verification of the correctness of objVal returned by the user
  /// (I think it is very important to have it)
  if ((_curBestMastColumnPtr != NULL) && (param().colGenDualPriceSmoothingAlphaFactor() == 0)
      && treeOfColClasses().empty())
    {
      long long int scaleFactor = param().SafeDualBoundScaleFactor();
      Double realObjVal = fixedDualCost() + fixedCost();
      if (scaleFactor > 0)
          realObjVal += probPtr()->objVal()._val / scaleFactor;
      else
          realObjVal += probPtr()->objVal();
      if  (_curBestMastColumnPtr->computeReducedCost() != realObjVal)
        {
            std::cerr << "BaPCod WARNING : objVal = " << probPtr()->objVal()
                      << " does not correspond to the cost of the best column = "
                      << _curBestMastColumnPtr->reducedCost() - fixedDualCost() - fixedCost() << std::endl;
            if (printL(-1))
              std::cout << "BaPCod WARNING : objVal = " << probPtr()->objVal()
                        << " does not correspond to the cost of the best column = "
                        << _curBestMastColumnPtr->reducedCost() - fixedDualCost() - fixedCost() << " ("
                        << realObjVal - _curBestMastColumnPtr->reducedCost()
                        << "), subproblem id = " << id().first()
                        << (enumeratedStatus() ? " (enum.)" : " (not enum.)") << std::endl
                        << "objVal = " << probPtr()->objVal() << ", fixedDualCost = "
                        << fixedDualCost() << ", fixedCost = " << fixedCost()
                        << ", reduced cost = " << _curBestMastColumnPtr->reducedCost() << std::endl;

          /* left for debugging */
          _curBestMastColumnPtr->spSol()->printOrderedSolution();
          double dp(_curBestMastColumnPtr->cost());
          for (ConstVarConstrPtr2Double::iterator itm = _curBestMastColumnPtr->member2coefMap().begin();
               itm != _curBestMastColumnPtr->member2coefMap().end(); ++itm)
            if (itm->first->inCurForm() && !itm->first->val().isZero())
              {
                dp += itm->first->val() * itm->first->membCoef(_curBestMastColumnPtr);
                std::cout << "Var[" << _curBestMastColumnPtr->name() << "] in const["
                          << itm->first->name() << "] of val[" << itm->first->val()
                          << "] has coef[" << itm->first->membCoef(_curBestMastColumnPtr) << "]  rc= "
                          << dp << std::endl;

              }
          dp = 0;
          for (VarPtr2DoubleMap::const_iterator itm = _curBestMastColumnPtr->spSol()->solVarValMap().begin();
               itm != _curBestMastColumnPtr->spSol()->solVarValMap().end(); ++itm)
            if (!itm->first->curCost().isZero())
              {
                dp += itm->first->curCost() * itm->second;
                std::cout << "Sp var[" << itm->first->name()
                          << "] with cost[" << itm->first->costrhs()
                          << "] with rc[" << itm->first->curCost()
                          << "] has coef[" << itm->second << "]  rc= " << dp << std::endl;
              }
          /* */

          if (printL(2))
            {
              const std::vector<int> & orderedIds = _curBestMastColumnPtr->spSol()->orderedIds();
              if (!orderedIds.empty() && (_networkFlowPtr != NULL))
                {
                  std::vector<int>::const_iterator arcIdIt = orderedIds.begin();
                  int vertId = _networkFlowPtr->netArcPtr(*arcIdIt)->tailVertexPtr()->id();
                  std::cout << "Column ordered vertex ids : " << vertId;
                  while (arcIdIt != orderedIds.end())
                    {
                      vertId = _networkFlowPtr->netArcPtr(*arcIdIt)->headVertexPtr()->id();
                      std::cout << " -> " << vertId;
                      ++arcIdIt;
                    }
                  std::cout << std::endl;
                }
              std::cout << *_curBestMastColumnPtr;
            }
        }
    }

    /**
     * We return all solutions found except the best solution recorded above
     */
    if (printL(5))
        std::cout << " ColGenSpConf::getSolPC(): MultiColGeneration : number of col to record = "
                  << probPtr()->recordedSolList().size() << std::endl;

    /// Ruslan : changed inserting level to 3 (unsuitable)
    /// real insertion with appropriate level will be done in insertColumnsInMaster()
    insertionLevel = param().Search4NegRedCostColInInactivePool ? 2 : 3;
    for (std::list <Solution * >::const_iterator spSolPt = probPtr()->recordedSolList().begin();
         spSolPt != probPtr()->recordedSolList().end(); spSolPt++)
    {
        recordSubproblemSolution(*spSolPt, inPurePhaseOne, insertionLevel);
    }
}

/// updates conf, solves the pricing suproblem and record its solutions, return the best solution found
Solution * ColGenSpConf::solvePC(int & maxLevelOfSubProbRestriction, bool inPurePhaseOne)
{
  if (curFormIsInfeasible()) return NULL;

  bapcodInit().check(mastConfPtr() == NULL, "ColGenSpConf::solvePC() masterConf should be defined");
  if (printL(2))
    std::cout << "ColGenSpConf::solvePC(): mastConfPtr()->treeOfColClasses().size() = "
              << mastConfPtr()->treeOfColClasses().size() << std::endl;

  if (spConfHasClassInducedSpVarBounds())
  {
    if (printL(5))
      std::cout << "recallMemorisedBounds() is called" << std::endl;

    for (VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 's');
        spVarPt != probPtr()->probVarSet().end(Active, 's'); spVarPt++)
      (*spVarPt)->recallMemorisedBounds();


    //Added by Romain: Comment this part if you don't need artificial var or dynamic var.
    for (VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'd');
        spVarPt != probPtr()->probVarSet().end(Active, 'd'); spVarPt++)
      (*spVarPt)->recallMemorisedBounds();

    for (VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'a');
        spVarPt != probPtr()->probVarSet().end(Active, 'a'); spVarPt++)
      (*spVarPt)->recallMemorisedBounds();

    if (updateConf(inPurePhaseOne)) // problem infeasible
    {
      if (printL(1))
        std::cout << "ColGenSpConf::solvePC() : updated subproblem is infeasible  " << std::endl;

      return(NULL);
    }
  }

  Time start;
  /**
   * probPtr()->solveProb returns true if a feasible solution to the problem was found,
   * false otherwise (i.e., if SP is infeasible)
   */
  int solverReturnStatus = probPtr()->solveProb(maxLevelOfSubProbRestriction, ' ', printL(3));

  if (printL(3))
    std::cout << "ColGenSpConf solveProb ReturnStatus = "
	      << solverReturnStatus
	      << std::endl;

  bapcodInit().statistics().incrCounter("bcCountSpSol");
  bapcodInit().statistics().incrTimer("bcTimeSpMPsol",  start.getElapsedTime_dbl());
  Bound rootDualBound(probPtr()->dualBound(), _objStatus);
  if (printL(5))
    std::cout << "  ColGenSpConf::solvePC() rootDualBound = " << rootDualBound << std::endl;

  if (printL(2) && (solverReturnStatus <= 0))
    std::cout << "ColGenSpConf::solve: no solution found for SP" << std::endl;

  std::list<Solution *> garbageCollection;
  /// @todo use smart pointeurs instead of garbageCollector

  /// Retrieve implicit solution \delta, compute mult and dualBoundContrib
  computeSpDualBoundContrib();
  Bound bestRedCost(spRootReducedCost(), _objStatus);

  /**
   * Retrieve explicit solution
   * Create Solution from inPrimalLpSol()
   */
  if (solverReturnStatus > 0)
  {
    getSolPC(inPurePhaseOne);
  }
  else if (solverReturnStatus < 0)
  {
    /// if reduced cost fixing is used, then the fact that a subproblem did not find any solutions
    /// does not mean that it is infeasible
    if ((maxLevelOfSubProbRestriction == 0) && (param().ReducedCostFixingThreshold() == 0))
      probPtr()->setProbStatus(SolutionStatus::Infeasible);

    if (printL(1))
    {
      std::cout << "  ColGenSpConf::solvePC() bestRedCost found at pricing tree root = " << bestRedCost
                << " Infeasible " << std::endl;
    }

    return(NULL);
  }
  else //(SpSolverReturnStatus == 0)
  {
    probPtr()->setProbStatus(SolutionStatus::UnSolved);
    if (printL(1))
    {
      std::cout << "  ColGenSpConf::solvePC() bestRedCost found at pricing tree root = " << bestRedCost
		        << " no sol found " << std::endl;
    }

    return(NULL);
  }

  Solution * bestSolPtr = _primalSolPtr->clone();

  if (!(treeOfColClasses().empty()))
  {
    sortTreeOfColClasses();
  }

  if (spRootReducedCost().negative(param().BapCodReducedCostTolerance))
  {
    if (printL(3))
      std::cout << "  ColGenSpConf::solvePC() found neg red cost SP sol at root: redCost =  " <<  spRootReducedCost()
                << ", solPtrRoot = " << *_primalSolPtr << std::endl;

    if (!(treeOfColClasses().empty()))
    {
      correctDualboundContribAndBestSol(rootDualBound, bestRedCost, bestSolPtr);
    }

    /// Reset normal curLb and Ub on SP var
    if (spConfHasClassInducedSpVarBounds())
    {
      if (printL(5))
        std::cout << "recallMemorisedBounds() is called" << std::endl;

      for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 's');
          spVarPt != probPtr()->probVarSet().end(Active, 's'); spVarPt++)
      {
        (*spVarPt)->recallMemorisedBounds();
      }

      //Added by Romain: Comment this part if you don't want to use dynamic or artificial var.
      for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'd');
          spVarPt != probPtr()->probVarSet().end(Active, 'd'); spVarPt++)
      {
        (*spVarPt)->recallMemorisedBounds();
      }

      for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'a');
          spVarPt != probPtr()->probVarSet().end(Active, 'a'); spVarPt++)
      {
        (*spVarPt)->recallMemorisedBounds();
      }

    }

    /// Negative red cost column found at the root of pricing SP tree
    if (printL(1))
    {
      std::cout << "  ColGenSpConf::solvePC() bestRedCost found = " << bestRedCost << std::endl;
    }

    if (_primalSolPtr != NULL)
        delete _primalSolPtr;
    _primalSolPtr = bestSolPtr;  /// in this case  = bestSolPtr;

    return (_primalSolPtr); /// in this case  = bestSolPtr;
  }
  else
  {
    if (printL(5))
    {
      std::cout << "  ColGenSpConf::solvePC() found NO neg red cost SP  sol at root: redCost =  "
                <<  spRootReducedCost() << ", solPtrRoot = " << *_primalSolPtr << std::endl;
    }
  }


  if (!(treeOfColClasses().empty()))
  {

    Solution * solPtrRoot = _primalSolPtr->clone();

    Double localBestPossibleRedCost(0), updatedRootRedCost(0);

    for (ColClassesVector::const_iterator bcPt = treeOfColClasses().begin();
         bcPt != treeOfColClasses().end(); bcPt++)
    {
      if (printL(5))
        std::cout << "ColGenSpConf::solvePC(): consider col class " << (*bcPt)->name() << std::endl;

      if ((*bcPt)->associatedPricingSPsolved())
      {
        if (printL(5))
          std::cout << "ColGenSpConf::solvePC(): associatedPricingSPsolved" << std::endl;

        /// A pred node sol solves this SP
        continue;
      }

      /// Check if predecessor class solution does not solve that class problem
      if ((*bcPt)->dirPredCSconstrPtr() == NULL)
      {
        /// Predecessor is root node class Q
        localBestPossibleRedCost = spRootReducedCost() + (*bcPt)->val();

        if (printL(5))
          std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name() << " has direct pred = the root node SP"
                    << std::endl;

        if ((*bcPt)->CBsatisfiedBySol(solPtrRoot))
        {
          Bound dualBd = rootDualBound; // + (*bcPt)->val();  // FVCHECK
          (*bcPt)->recSol(solPtrRoot, dualBd, Bound(localBestPossibleRedCost, _objStatus));
          if (printL(5))
            std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                      << " direct pred solution does satisfy coump bounds" << std::endl;
        }
        else
	  {
	    if (printL(5))
	      std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                    << " direct pred solution does not satisfy coump bounds" << std::endl;
	  }
      }
      /// Predecessor is a regular node class Q_k
      else
      {
        if (printL(5))
          std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name() << " has direct pred "
                    << (*bcPt)->dirPredCSconstrPtr()->name() << std::endl;

        CompSetInstMastBranchConstr * predCSconstrPtr((*bcPt)->dirPredCSconstrPtr());
        do
        {
          if (predCSconstrPtr->associatedPricingSPsolved())
          {
            localBestPossibleRedCost = predCSconstrPtr->reducedCost() + (*bcPt)->val();

            if ((*bcPt)->CBsatisfiedBySol(predCSconstrPtr->solPtr()))
            {
              /**
               * Record reduced cost and SP sol
               * sol of relaxed problem solves the more constrained problem
               */

	          Bound dualBd = predCSconstrPtr->pricingSpZetaVal();
              (*bcPt)->recSol(predCSconstrPtr->solPtr(), dualBd, Bound(localBestPossibleRedCost, _objStatus));

              if (printL(5))
                std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                          << " has direct pred sol that does  satisfy coump bounds" << std::endl;
            }
            else if (printL(5))
              std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                        << " has direct pred sol that does not satisfy coump bounds" << std::endl;
            break;
          }
          else
          {
            if (printL(5))
              std::cout << " predecessor problem has not been solved even though we used breadth first search"
                        << std::endl;
            bapcodInit().require(predCSconstrPtr != NULL,
                                 "ColGenSpConf::solvePC(): should have treated above the case"
                                 "where precedessor is the root of the class tree");
            predCSconstrPtr = predCSconstrPtr->dirPredCSconstrPtr();

            /// Predecessor of predecessor is the root
            if (predCSconstrPtr == NULL)
            {
              localBestPossibleRedCost = spRootReducedCost() + (*bcPt)->val();

              if (printL(5))
                std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                          << " has direct pred = the root node SP" << std::endl;

              if ((*bcPt)->CBsatisfiedBySol(solPtrRoot))
              {
		        Bound dualBd = rootDualBound; // + (*bcPt)->val(); // FVCHECK
                (*bcPt)->recSol(solPtrRoot, Bound(rootDualBound + (*bcPt)->val(), _objStatus),
                                Bound(localBestPossibleRedCost, _objStatus));
                if (printL(5))
                  std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                            << " pred solution does  satisfy coump bounds" << std::endl;
              }
              else if (printL(5))
                std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name()
                          << " pred solution does not satisfy coump bounds" << std::endl;
              break;
            }
          }
        } while (true);

      }

      if ((*bcPt)->associatedPricingSPsolved())
      {
        /// Check if predessor class value gives rise to find a neg red cost col
        if (printL(5))
          std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name() << " is solved " << std::endl;

        if ((*bcPt)->reducedCost().negative(param().BapCodReducedCostTolerance))
        {
          if (printL(5))
            std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name() << " has neg reduced cost "
                      << (*bcPt)->reducedCost() << std::endl;
          break;
        }

        /// SP^k solved go to next one
        continue;
      }
      /**
       * Test if pricing SP has a chance to yield a negative reduced cost col
       * error in dualbound computation ?
       */
      if (printL(5))
        std::cout << "ColGenSpConf::solvePC(): class bestPossibleRedCost = " << localBestPossibleRedCost <<  std::endl;

      if (!param().PriceAllSubproblemsForBestDualBound())
      /// if true, one prices columns over subproblem class even if the reduced cost is expected to be positive;
      /// this allow to get tigher dual bound; this is needed in the case a node has integer solution --
      /// i.e. cannot be separated by branching  -- and ther is an optimality gap
      {

        if (!(localBestPossibleRedCost.negative(param().BapCodReducedCostTolerance)))
        {
          if (printL(5))
            std::cout << "ColGenSpConf::solvePC(): class cannot yield neg red cost column " <<  std::endl;
          continue;
        }
      }

      bool SPinfeasible(false);

            /**
             * Reset memorized SP var bounds only for var concerned by class:
             * Reset normal curLb and Ub on SP var
             */
            for (ComponentSequence::const_iterator cbPt = (*bcPt)->compBoundSet().begin();
                 cbPt != (*bcPt)->compBoundSet().end(); cbPt++)
                (cbPt)->varPtr()->recallMemorisedBounds();

            /// Record tightest spVarBd
            VarPtr2DoubleMap spVarLbMap;
            VarPtr2DoubleMap spVarUbMap;
            spVarLbMap.clear();  spVarUbMap.clear();

            for (ComponentSequence::const_iterator cbPt = (*bcPt)->compBoundSet().begin();
                 cbPt != (*bcPt)->compBoundSet().end(); cbPt++)
            {
                /// Set lower bound
                if ((cbPt)->sign() == 'G')
                {
                    if (!spVarLbMap.count((cbPt)->varPtr())) spVarLbMap[(cbPt)->varPtr()] =  (cbPt)->val();
                    else if ((cbPt)->val() > spVarLbMap[(cbPt)->varPtr()])
                        /// Record highest LB
                        spVarLbMap[(cbPt)->varPtr()] =  (cbPt)->val();
                }
                /// Set upper bound because ((cbPt)->sign() == 'L')
                else
                {
                    if (!spVarUbMap.count((cbPt)->varPtr()))
                        spVarUbMap[(cbPt)->varPtr()] =  (cbPt)->val();
                    else if ((cbPt)->val() < spVarUbMap[(cbPt)->varPtr()])
                        // Record lowest UB
                        spVarUbMap[(cbPt)->varPtr()] =  (cbPt)->val();
                }
            }

            /**
             * Set lower bound ignoring local lower bound
             * that is valid for class not restricted
             * by specific upper bound
             */
            for (VarPtr2DoubleMap::const_iterator spLbPt = spVarLbMap.begin();
                 spLbPt != spVarLbMap.end(); spLbPt++)
            {
                spLbPt->first->curLb(spLbPt->second);
                if (printL(5))
                    std::cout << " set var " << spLbPt->first->name() << " bounds in [" << spLbPt->first->curLb()
                              << " , " << spLbPt->first->curUb() << "]" << std::endl;
            }

            /** Set upper bound, ignoring local upper bound
             * that is valid for class not restricted
             * by specific upper bound
             */
            for (VarPtr2DoubleMap::const_iterator spUbPt = spVarUbMap.begin();
                 spUbPt != spVarUbMap.end(); spUbPt++)
            {
                spUbPt->first->curUb(spUbPt->second);
                if (printL(5))
                    std::cout << " set var " << spUbPt->first->name() << " bounds in [" << spUbPt->first->curLb()
                              << " , " << spUbPt->first->curUb() << "]" << std::endl;
            }

            /// Need to reset subproblem ?
            probPtr()->clearRecordedSol();
            if (SPinfeasible)
            {
                bapcodInit().check(1, "ColGenSpConf::solve: SP infeasible");
            }
            else
            {
                bapcodInit().check(updateConf(inPurePhaseOne), "ColGenSpConf::solve: SP should not infeasible");
                /// Solve class pricing SP with new var bounds
                Time start4;
                /// Return true if a feasible solution to the problem was found, false otherwise
                solverReturnStatus = probPtr()->solveProb(maxLevelOfSubProbRestriction, ' ', printL(3));
                bapcodInit().statistics().incrCounter("bcCountSpSol");
                bapcodInit().statistics().incrTimer("bcTimeSpMPsol",  start4.getElapsedTime_dbl());

                if (printL(2) && (solverReturnStatus <= 0))
                    std::cout << "ColGenSpConf::solve: no solution found for SP" << std::endl;
            }

            /// Retrieve explicit solution x
            updatedRootRedCost = spRootReducedCost() - rootDualBound + probPtr()->dualBound();
            localBestPossibleRedCost = updatedRootRedCost - (*bcPt)->sigma();

            if (printL(5))
                std::cout << "constraint " << (*bcPt)->name() << std::endl
                << "  spRootReducedCost() = " << spRootReducedCost() << std::endl
                << "  rootDualBound = " << rootDualBound << std::endl
                << "  probPtr()->dualBound() = " << probPtr()->dualBound() << std::endl
                << "  updatedRootRedCost = " << updatedRootRedCost << std::endl
                << "  (*bcPt)->sigma() = " <<   (*bcPt)->sigma() << std::endl
                << "  localBestPossibleRedCost = " << localBestPossibleRedCost << std::endl
			    << "solverReturnStatus = " << solverReturnStatus << std::endl;

            /**
             * Create Solution from inPrimalLpSol()
             * && !probPtr()->inPrimalLpSol().empty() empty sol is acceptable
             */
            if (solverReturnStatus > 0)
            {
                getSolPC(inPurePhaseOne);
                bapcodInit().require((*bcPt)->CBsatisfiedBySol(_primalSolPtr),
                                     "ColGenSpConf::solvePC() pricing SP solution should satisfy CS bound" );
            }

            /// Record  pricing SP solution
            (*bcPt)->recSol(_primalSolPtr, Bound(probPtr()->dualBound(), _objStatus),
                            Bound(localBestPossibleRedCost, _objStatus));

            /**
             * Reset SP var bounds, do it after col record
             * otherwise col might be unsuitable
             */
            for (ComponentSequence::const_iterator cbPt = (*bcPt)->compBoundSet().begin();
                 cbPt != (*bcPt)->compBoundSet().end(); cbPt++)
            {
                if (printL(6))
                    std::cout << " var " << (cbPt)->varPtr()->name() << " bounds in [" << (cbPt)->varPtr()->curLb()
                              << " , " << (cbPt)->varPtr()->curUb() << "] Before reset" << std::endl;

                (cbPt)->varPtr()->recallMemorisedBounds();

                if (printL(6))
                    std::cout << " var " << (cbPt)->varPtr()->name() << " bounds in [" << (cbPt)->varPtr()->curLb()
                              << " , " << (cbPt)->varPtr()->curUb() << "] After reset" << std::endl;
            }

            /// Break if a neg red col has been identified
            if (localBestPossibleRedCost.negative(param().BapCodReducedCostTolerance))
            {
                if (printL(5))
                    std::cout << "ColGenSpConf::solvePC(): col class " << (*bcPt)->name() << " has neg reduced cost "
                              << (*bcPt)->reducedCost() << std::endl;

                /// Goto LABEL endOfListingSp
                break;
            }

            // Attempt to see if _primalSolPtr solves successor pricing SP
            for (ColClassesVector::const_iterator nextbcPt = bcPt + 1; nextbcPt != treeOfColClasses().end(); nextbcPt++)
            {
                /// Variable (*bcPt) is a predecessor of (*nextbcPt)
                if ((*nextbcPt)->setOfPredCSconstrPtr().count(*bcPt))
                {

                    if ((*nextbcPt)->CBsatisfiedBySol(_primalSolPtr))
                    {
                        localBestPossibleRedCost = updatedRootRedCost - (*nextbcPt)->sigma();

                        /// Reduced cost computed compare to root reduced cost
                        (*nextbcPt)->recSol(_primalSolPtr, Bound(probPtr()->dualBound(), _objStatus),
                                            Bound(localBestPossibleRedCost, _objStatus));

                        if (localBestPossibleRedCost.negative(param().BapCodReducedCostTolerance))
                        {
                            if (printL(5))
                                std::cout << "ColGenSpConf::solvePC(): col class " << (*nextbcPt)->name()
                                          << " has neg reduced cost " << (*nextbcPt)->reducedCost() << std::endl;
                            break;
                        }
                    }

                }
            }

            /// Break if a neg red col has been identified
            if (localBestPossibleRedCost.negative(param().BapCodReducedCostTolerance))
            {
                if (printL(5))
                    std::cout << "ColGenSpConf::solvePC(): col class "  << (*bcPt)->name() << " has neg reduced cost "
                              << (*bcPt)->reducedCost() << std::endl;
                /// Goto LABEL endOfListingSp
                break;
            }
            else
            {
                if (printL(5))
                    std::cout << "ColGenSpConf::solvePC(): class found no solution column " <<  std::endl;
            }

        }
        /// Label endOfListingSp

        correctDualboundContribAndBestSol(rootDualBound, bestRedCost, bestSolPtr);

    } // (!(treeOfColClasses().empty())

    if (printL(7))
    {
        std::cout <<" AFTER probPtr()->solveProb()" << std::endl;
        if (probPtr()->formulationPtr() != NULL)
            probPtr()->formulationPtr()->printForm();
    }

    if (spConfHasClassInducedSpVarBounds())
    {
        if (printL(5))
            std::cout << "recallMemorisedBounds() is called" << std::endl;

        for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 's');
            spVarPt != probPtr()->probVarSet().end(Active, 's'); spVarPt++)
        {
            (*spVarPt)->recallMemorisedBounds();
        }

        //Added by Romain: Comment this part if you don't want to use dynamic or artificial var.
        for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'd');
            spVarPt != probPtr()->probVarSet().end(Active, 'd'); spVarPt++)
        {
            (*spVarPt)->recallMemorisedBounds();
        }

        for(VarIndexManager::iterator spVarPt = probPtr()->probVarSet().begin(Active, 'a');
            spVarPt != probPtr()->probVarSet().end(Active, 'a'); spVarPt++)
        {
            (*spVarPt)->recallMemorisedBounds();
        }
    }

    if (printL(1))
    {
        std::cout << "  ColGenSpConf::solvePC() bestRedCost found = " << bestRedCost << std::endl;
    }

        if (_primalSolPtr != NULL) delete _primalSolPtr;
        _primalSolPtr = bestSolPtr;  /// in this case  = bestSolPtr;

    probPtr()->resetSolution();

    return _primalSolPtr;
}

bool ColGenSpConf::cannotGenerateAnyMoreCol()
{

    if (_upperBoundMastConstrPtr != NULL)
    {
        if (_upperBoundMastConstrPtr->curRhs() <= 0)
        {
            if (printL(3))
                std::cout << "ColGenSpConf::cannotGenerateAnyMoreCol(): _upperBoundMastConstrPtr->curRhs() ="
                          << _upperBoundMastConstrPtr->curRhs() << std::endl;

            /// No need to attempt generating columns
            return true;
        }
    }
    return false;

}

bool ColGenSpConf::cardinalityIsFixed()
{
    if (_implicitlyFixCardinality) return true;

    if (_lowerBoundMastConstrPtr == NULL)
        return false;

    if (_upperBoundMastConstrPtr == NULL)
        return false;

    return (_lowerBoundMastConstrPtr->rhs() == _upperBoundMastConstrPtr->rhs());
}

bool ColGenSpConf::needNotGenerateAnyMoreCol()
{

    if (_lowerBoundMastConstrPtr == NULL)
        return true;

    /// No suproblem solution required in the master
    if (_lowerBoundMastConstrPtr->curRhs() <= 0)
        return true;

    return false;

}

bool ColGenSpConf::performSpRelaxationLigntening(const bool & masterConverged, const int & callMode)
{
  if (!probPtr()->solverOracleFunctorDefined())
    return false;
    
  updateTarget(false);
  updateConf(false);
  return probPtr()->lightenCurrentSpRelaxation(masterConverged, callMode);
}

/// returns the number of generated columns
int ColGenSpConf::genNewColumns(bool currentlyPerformingPhaseI, int & maxLevelOfSubProbRestriction,
                                int & allAddedColumns, int & addedNegRedCostColumns)
{
  _curBestMastColumnPtr = NULL;     // stabilization: make sure that no previous column is used as the best (Artur)

  const int returnFlagNeedNotGenerateAnyMoreCol = 0;
  const int returnFlagSpIsInfeasible = -1;
  const int returnFlagCannotGenerateAnyMoreCol = -2;

  bapcodInit().statistics().incrCounter("bcCountCgSpSolverCall");

  /// Important to reset _dualBoundContrib to zero in case
  /// there are not computed due to (minCost() >= target())
  _dualBoundContrib = Bound(0, _objStatus);

  if (cannotGenerateAnyMoreCol())
    return returnFlagCannotGenerateAnyMoreCol;

  /// Compute target
  updateTarget(currentlyPerformingPhaseI);

  /// Reset var bounds, var cost, sp minCost
  if (updateConf(currentlyPerformingPhaseI))
    {
      if (printL(3))
        std::cout << "ColGenSpConf::genNewCol(): SP infeasible" << std::endl;

      computeSpDualBoundContrib();

      if (needNotGenerateAnyMoreCol())
        return returnFlagNeedNotGenerateAnyMoreCol;

      return returnFlagSpIsInfeasible;
    }

  /// switch off the reduced cost estimation when stabilization is applied
  if ((param().colGenDualPriceSmoothingAlphaFactor() == 0)
      && (param().colGenStabilizationFunctionType().getStatusAsInteger() == 0)
      && !bapcodInit().param().SplitColIntoDissagregateSpVar()
      && (param().mastInitMode().status() != MasterInitMode::noArtCol)
      && !spConfHasClassInducedSpVarBounds())    /// added by Ruslan: TO DO verify if it is ok
    {
      /// check whether need to call sp sol
      switch (_mastConfPtr->probPtr()->solMode().status())
        {
        case SolutionMethod::lpSolver:
        case SolutionMethod::mipSolver:
        case SolutionMethod::customSolver:
        case SolutionMethod::custom2mipSolver:
        {
          Bound dualPrimalBoundDiff(probPtr()->dualBound() - probPtr()->primalBound(), _objStatus);
          if (!dualPrimalBoundDiff.negative(param().BapCodReducedCostTolerance()))
            {
              /// problem solved to optimality
              if (printL(3))
                std::cout << "ColGenSpConf::genNewCol(): a priori dualBound =" << probPtr()->dualBound()
                          << " >= a priori primalBound  =" << probPtr()->primalBound()
                          << " equal to target = " << target() << std::endl;

              computeSpDualBoundContrib();
              return returnFlagNeedNotGenerateAnyMoreCol;
            }
          break;
        }
        case SolutionMethod::undefined:
        {
          bapcodInit().check(true, "ColGenSpConf::genNewCol: ERROR undefined solution method");
          break;
        }
        case SolutionMethod::none:
        default:
        {
          break;
        }
        }
    }

  /// Solve sub-problem: call oracle
  Time start1;
  Solution * solPtr = solvePC(maxLevelOfSubProbRestriction, currentlyPerformingPhaseI);
  bapcodInit().statistics().incrTimer("bcTimeSpSol",  start1.getElapsedTime_dbl());

  clearSolutions();
  
  /// Empty solution
  if (solPtr == NULL)
    {
      if (printL(2))
        std::cout << "ColGenSpConf::genNewCol(): NEW col is the empty column" << std::endl;

      return (probPtr()->probStatus().intersects(SolutionStatus::Infeasible) ? -1 : 0);
    }
    
  insertConstraintsInMaster();

  insertColumnsInMaster(allAddedColumns, addedNegRedCostColumns);

  return 0;
}

void ColGenSpConf::clearSolutions()
{
  if (_primalSolPtr != NULL)
    delete _primalSolPtr;
  _primalSolPtr = NULL;
}

/// Record other SP solutions as extra columns
/// Check whether master columns already exist in pool
MastColumn * ColGenSpConf::checkColumn4Insertion(MastColumn * colPtr, bool inPurePhaseOne, const int & insertionLevel)
{
  if (printL(5))
    std::cout << "ColGenSpConf::checkColumn4Insertion(): Test New column " << colPtr->name()
              << ", insertionLevel = " << insertionLevel << std::endl;

  if (insertionLevel > 0)
    {
      VcIndexStatus::VcStatus colStatus =  VcIndexStatus::Undefined;
      Variable * tempColPtr = _mastConfPtr->probPtr()->probVarSet().findPtr(colPtr);
      /// if UseColumnsPool is false, if the same column exists already, we use it only if it is active
      /// if column is not active, its membership might be wrong
      if ((tempColPtr != NULL) && (param().UseColumnsPool() || (tempColPtr->vcIndexStatus() == VcIndexStatus::Active)))
        {
          MastColumn * cPtr =  static_cast< MastColumn * >(tempColPtr);

          colStatus = tempColPtr->vcIndexStatus();

          if (printL(2))
            std::cout << "ColGenSpConf::checkColumn4Insertion(): New column " << colPtr->name()
                      << "  already exists as " << cPtr->name() << " colsize " << cPtr->spSol()->size()
                      << "  colStatus = " << colStatus << std::endl;
          /**
           * Reset its multiplicity to that of generated column needed in cutting planes approaches
           */
          cPtr->mult(colPtr->mult());
          delete colPtr;
          colPtr = cPtr;

          switch(_mastConfPtr->probPtr()->solMode().status())
            {
            case SolutionMethod::lpSolver:
            case SolutionMethod::mipSolver:
            case SolutionMethod::customSolver:
            case SolutionMethod::custom2mipSolver:
            {
              if ((insertionLevel == 1) && (colStatus == VcIndexStatus::Active))
                {
                  if (!param().SplitColIntoDissagregateSpVar())
                  {
                    Double redCost = colPtr->computeReducedCost();
                    if (printL(0) && redCost.negative(10 * param().BapCodReducedCostTolerance()))
                    {
                        std::cout << "BaPCod warning : Existing Column in master should have non-negative reduced cost"
                                  << redCost << std::endl;
                    }
                  }
                }
              break;
            }
            case SolutionMethod::none:
              break;
            case SolutionMethod::undefined:
              bapcodInit().check(true, "ColGenSpConf::checkColumn4Insertion(): ERROR undefined solution method");
              break;
            }
        }
      else
	    {
	      bapcodInit().statistics().incrCounter("bcCountCol");


          if (printL(2))
            std::cout << "ColGenSpConf::checkColumn4Insertion(): New column " << colPtr->name()
		              << "  does not already exists " << std::endl;

          colPtr->setAggregateVariable(colPtr); /// moved here from the constructor of MastColumn by Ruslan
          colPtr->recordInMembershipOfSubProbVar();
        }
    }

  /// this will add the column to the problem and set the membership (if the column is not already in the problem),
  /// insertion to the formulation (at level 1) will be done using _tempColPtrList4Insertion
  if (colPtr->vcIndexStatus() != VcIndexStatus::Active)
    {
      if ((insertionLevel <= 2) && param().Search4NegRedCostColInInactivePool())
        _mastConfPtr->probPtr()->addVar(colPtr, 2);
      else
        /// we never make a column inactive if parameter Search4NegRedCostColInInactivePool is false,
        /// as inactive columns are never cleaned up
        _mastConfPtr->probPtr()->addVar(colPtr, 3);
      colPtr->resetCost(inPurePhaseOne); /// we need to reset cost due to pure phase one
    }

  if (printL(3))
    colPtr->print();

  if (param().GenerateProperColumns())
    bapcodInit().require(colPtr->suitableForResidualProb(),
                         "ColGenSpConf::checkColumn4Insertion(): generated Mast Column should be suitable "
                         "if it is a SP solution, CHECK that the oracle enforces bounds on subproblem variables");

    switch(_mastConfPtr->probPtr()->solMode().status()) {
       case SolutionMethod::lpSolver:
       case SolutionMethod::mipSolver:
       case SolutionMethod::customSolver:
       case SolutionMethod::custom2mipSolver:
       {
           if ( ((insertionLevel == 1) || param().InsertAllGeneratedColumnsInFormRatherThanInPool())
               && (colPtr->vcIndexStatus() != VcIndexStatus::Active))
           {
             colPtr->incrParticipation(17);
             _tempColPtrList4Insertion.push_back(colPtr);
             if (printL(2))
                   std::cout << "ColGenSpConf::checkColumn4Insertion(): insertionLevel = " << insertionLevel
                             << ", column registered for direct inclusion " << colPtr->name() << std::endl;
           }
           break;
       }
       case SolutionMethod::none:
       {
           break;
       }
       case SolutionMethod::undefined:
       {
           bapcodInit().check(true, "ColGenSpConf::checkColumn4Insertion(): ERROR undefined solution method");
           break;
       }
   }
  return(colPtr);
}


InstanciatedConstr * ColGenSpConf::checkConstraint4Insertion(Constraint * constrPtr, const int & insertionLevel)
{
  if (constrPtr->isTypeOf(VcId::InstanciatedConstrMask))
    {
      return checkConstraint4Insertion((InstanciatedConstr*) constrPtr, insertionLevel);
    }
  else if (printL(3))
    {
      std::cout << "ColGenSpConf::checkConstraint4Insertion(Constraint *) UNDEFINED: ";
      constrPtr->print();
    }
  return NULL;

}

InstanciatedConstr * ColGenSpConf::checkConstraint4Insertion(InstanciatedConstr * iconstrPtr,
                                                             const int & insertionLevel)
{
  if (printL(3))
    std::cout << "ColGenSpConf::checkConstraint4Insertion(InstanciatedConstr *) check constraint  "
	          << iconstrPtr->name() << " insertionLevel = " << insertionLevel
	          << " is MasterMask? " << iconstrPtr->probConfPtr()->isTypeOf(PcId::MasterMask) << std::endl;

  if ((insertionLevel > 0) && true)
    {
      iconstrPtr = _mastConfPtr->castAndAddConstraint(iconstrPtr);

      if (printL(3))
	    std::cout << "ColGenSpConf::checkConstraint4Insertion(InstanciatedConstr *) check constraint  != NULL "
		          << (iconstrPtr != NULL) << std::endl;

      if (iconstrPtr != NULL)
	    {
	      if (printL(3))
	        iconstrPtr->print();

	      _tempMastConstrPtrList4Insertion.push_back(iconstrPtr);
	      if (printL(3))
	        {
              std::cout << "ColGenSpConf::checkConstraint4Insertion(InstanciatedConstr *) push constraint  "
                        << iconstrPtr->name();
		      iconstrPtr->print() << std::endl;

	        }
	    }

      return(iconstrPtr);
    }

  return NULL;

}

MastColumn * ColGenSpConf::recordSubproblemSolution(Solution * spSolPtr, bool inPurePhaseOne,
                                                    const int & insertionLevel, Solution * masterSolPtr,
                                                    bool changeEnumeratedFlag)
{
#ifdef BC_MORE_TIMERS
    Time start;
#endif

    /**
     * Define column from sp solution
     * Call application specific column constructor
     * Constraint membership is set in addVar if col does not already exists
     */

    if (printL(3))
        std::cout << " RecordSubproblemSolution for ColGenSpConf "  << name() << " with insertionLevel "
                  << insertionLevel << std::endl;

    MastColumn * colPtr(NULL);

    if (spSolPtr == NULL) return NULL;


    if (param().SplitColIntoDissagregateSpVar && _toSplit)
    {
        for (VarPtr2DoubleMap::const_iterator it = spSolPtr->solVarValMap().begin();
             it != spSolPtr->solVarValMap().end(); ++it)
        {
            Solution * newSolPtr = new Solution(this);
            newSolPtr->cost(it->first->costrhs());
            newSolPtr->includeVar(it->first, 1, false);

            colPtr = new MastColumn(_mastConfPtr, this, newSolPtr, it->first->name());
            if (probPtr()->curNodePtr() != NULL)
                colPtr->treatOrderId(probPtr()->curNodePtr()->treatOrder());

            if (printL(3))
                std::cout << "NEWLY GENERATED Dissagr Mast Column: " << colPtr->name() << std::endl;
            /// The column pointer can be changed within checkColumn4Insertion() if the column already exists
            colPtr = checkColumn4Insertion(colPtr, inPurePhaseOne, insertionLevel);
            /// Needs to be done after  checkColumn4Insertion, as colPtr can change in  checkColumn4Insertion
            if (masterSolPtr != NULL)
                masterSolPtr->includeVar(colPtr, spSolPtr->multiplicity(), true);
        }
    }
    else
    {
        spSolPtr->resetCost();
        colPtr = new MastColumn(_mastConfPtr, this, spSolPtr);
        if (changeEnumeratedFlag)
            colPtr->spSol()->enumeratedFlag(!spSolPtr->enumeratedFlag());
        if (probPtr()->curNodePtr() != NULL)
            colPtr->treatOrderId(probPtr()->curNodePtr()->treatOrder());
        if (printL(3))
            std::cout << "NEWLY GENERATED Mast Column: " << colPtr->name() << ", insertionLevel = "
                      << insertionLevel << std::endl;

        /// The column pointer can be changed within checkColumn4Insertion() if the column already exists
        colPtr = checkColumn4Insertion(colPtr, inPurePhaseOne, insertionLevel);
        /// Needs to be done after  checkColumn4Insertion, as colPtr can change in  checkColumn4Insertion
        if (masterSolPtr != NULL)
            masterSolPtr->includeVar(colPtr, spSolPtr->multiplicity(), true);

        return colPtr;
    }

#ifdef BC_MORE_TIMERS
    bapcodInit().statistics().incrTimer("bcTimeRecordColumns",  start.getElapsedTime_dbl());
#endif

    return NULL;
}

int ColGenSpConf::insertAllColumnsInMaster()
{
  int counter = 0;
  for (std::list<MastColumn *>::iterator colPt = _tempColPtrList4Insertion.begin();
       colPt != _tempColPtrList4Insertion.end(); ++colPt)
    {
      counter += _mastConfPtr->probPtr()->addVar((*colPt), 1, 1);
    }
  if (printL(0))
    std::cout << "Added " << counter << " columns in the formulation " << std::endl;
  clearColPtrList4Insertion();
  return counter;
}

void ColGenSpConf::clearColPtrList4Insertion()
{
  for (std::list<MastColumn *>::iterator colPt = _tempColPtrList4Insertion.begin();
       colPt != _tempColPtrList4Insertion.end(); ++colPt)
    (*colPt)->decrParticipation(11);
  _tempColPtrList4Insertion.clear();
  if (printL(5))
      std::cout << "BaPCod info :  _tempColPtrList4Insertion has been cleared " << std::endl;
}

void ColGenSpConf::insertColumnsInMaster(int & allAddedColumns, int & addedNegRedCostColumns)
{
    int counter(0);

    for (std::list<MastColumn *>::iterator colPt = _tempColPtrList4Insertion.begin();
         colPt != _tempColPtrList4Insertion.end(); ++colPt)
    {
        switch(_mastConfPtr->probPtr()->solMode().status())
        {
            case SolutionMethod::lpSolver :
            case SolutionMethod::mipSolver :
            case SolutionMethod::customSolver :
            case SolutionMethod::custom2mipSolver :
            {
                (*colPt)->computeReducedCost();
                bool redCostIsNegative = (param().SafeDualBoundScaleFactor() > 0 ? (*colPt)->reducedCost()._val < 0
                                          : (*colPt)->reducedCost().negative(param().BapCodReducedCostTolerance()));
                if (printL(5))
                    std::cout << "Column Reduced Cost = " << (*colPt)->reducedCost() << std::endl;
                if (redCostIsNegative)
                {
                    int numAddedVars = _mastConfPtr->probPtr()->addVar((*colPt), 1,
                                                                       (param().SplitColIntoDissagregateSpVar ? 2 : 1));
                    allAddedColumns += numAddedVars;
                    addedNegRedCostColumns += numAddedVars;

                    if (printL(5))
                        std::cout << "counter = " << counter << std::endl;
                    if (printL(2))
                    {
                        std::cout << "ColGenSpConf::insertColumnsInMaster(): NEW col " << (*colPt)->name()
                                  << " has NEG RED COST (" << (*colPt)->reducedCost()
                                  << ") and therefore inserted in the formulation" << std::endl;
                        (*colPt)->spSol()->shortPrint();
                        (*colPt)->spSol()->printOrderedSolution();
                        std::cout << std::endl;
                    }
                }
                else if (param().InsertNewNonNegColumnsDirectlyInFormRatherThanInPool())
                {
                    allAddedColumns += _mastConfPtr->probPtr()->addVar((*colPt), 1,
                                                                       (param().SplitColIntoDissagregateSpVar ? 2 : 1));

                    if (printL(2))
                    {
                        std::cout << "ColGenSpConf::insertColumnsInMaster(): NEW col " << (*colPt)->name()
                                  << " has non neg red cost (" << (*colPt)->reducedCost()
                                  << ") but is inserted directly in the formulation" << std::endl;
                        (*colPt)->spSol()->shortPrint();
                        (*colPt)->spSol()->printOrderedSolution();
                        std::cout << std::endl;
                    }
                }
                else if (param().Search4NegRedCostColInInactivePool())
                {
                    if (printL(2))
                        std::cout << "ColGenSpConf::insertColumnsInMaster(): NEW col " << (*colPt)->name()
                                  << " is made inactive" << std::endl;
                    _mastConfPtr->probPtr()->addVar((*colPt), 2, 0);
                }
                else
                {
                    if (printL(2))
                        std::cout << "NEW col " << (*colPt)->name()
                                  << " is NOT inserted (vcIndex = " << (*colPt)->vcIndexStatus()
                                  << "), redCost = " << (*colPt)->reducedCost() << std::endl ;
                }
                break;
            }
            case SolutionMethod::none :
                break;
            case SolutionMethod::undefined :
                bapcodInit().check(true,
                                   "ColGenSpConf::insertColumnsInMaster(): ERROR undefined solution method");
                break;
        }
    }

    clearColPtrList4Insertion();
}


int ColGenSpConf::insertConstraintsInMaster()
{
  if (printL(2))
	      std::cout << "ColGenSpConf::insertConstraintsInMaster(): _tempMastConstrPtrList4Insertion.size() =  "
			        << _tempMastConstrPtrList4Insertion.size() << std::endl;
  int counter(0);

  for (std::list<InstanciatedConstr *>::iterator constrPt = _tempMastConstrPtrList4Insertion.begin();
       constrPt != _tempMastConstrPtrList4Insertion.end(); ++constrPt)
    {
      switch(_mastConfPtr->probPtr()->solMode().status())
        {
        case SolutionMethod::lpSolver :
        case SolutionMethod::mipSolver :
        case SolutionMethod::customSolver :
        case SolutionMethod::custom2mipSolver :
          {
            /// Test of no negative reduced cost is not a valid stopping criteria
            if (printL(2))
	          std::cout << "ColGenSpConf::insertConstraintsInMaster(): has inserted NEW master constraint "
			            << (*constrPt)->name() << std::endl;
            _mastConfPtr->probPtr()->addConstr((*constrPt), 1,  2);
            counter++;
            break;
          }
        case SolutionMethod::none :
          break;
        case SolutionMethod::undefined :
	      bapcodInit().check(true,
                             "ColGenSpConf::insertConstraintsInMaster(): ERROR undefined solution method");
	      break;
        }
    }

  _tempMastConstrPtrList4Insertion.clear();

  return(counter);
}

/**
 * Retrieve implicit solution \delta, compute mult and dualBoundContrib
 * Value dualBoundContrib includes a correction of the master objective value to eliminate the input
 * of the convexity constraints
 */
void ColGenSpConf::computeSpDualBoundContrib()
{
  long long int scaleFactor = param().SafeDualBoundScaleFactor();

  _dualBoundContrib = Bound(0, _objStatus);
  _spRootReducedCost = 0;

  Double pricingSpSolutionValue = 0;
  if (probPtr()->probInfeasibleFlag())
    {
      pricingSpSolutionValue = BapcodInfinity;
      if (printL(5))
        std::cout << " subProb Infeasible: individual spDualBdContrib = " << pricingSpSolutionValue << std::endl;
    }
  else
    {
      if (scaleFactor > 0)
          pricingSpSolutionValue = fixedCost() + probPtr()->dualBound()._val / scaleFactor;
      else
          pricingSpSolutionValue =  zero(fixedCost() + probPtr()->dualBound(),
                                         param().BapCodReducedCostTolerance());

      if (printL(5))
        std::cout << " subProb Feasible: individual spDualBdContrib = " << pricingSpSolutionValue
                  << " upper convexity constr contrib = "
                  << (_upperBoundMastConstrPtr != NULL ? _upperBoundMastConstrPtr->valOrSepPointVal()
                                                       : Double::staticZero)
                  << " lower convexity constr contrib = "
                  << (_lowerBoundMastConstrPtr != NULL ? _lowerBoundMastConstrPtr->valOrSepPointVal()
                                                       : Double::staticZero)
		          << " fixedCost() = " << fixedCost() << " fixedDualCost() = " << fixedDualCost() << std::endl;
    }

  _spRootReducedCost = pricingSpSolutionValue + fixedDualCost();

  if (printL(5))
      std::cout << "fixedCost = "  << fixedCost() << std::endl
                << "  probPtr()->dualBound() = " << probPtr()->dualBound() << std::endl
                << "  fixedCost + probPtr()->dualBound() = " << fixedCost() + probPtr()->dualBound() << std::endl
                << "  zero(fixedCost() + probPtr()->dualBound(), param().BapCodReducedCostTolerance) = "
                << zero(fixedCost() + probPtr()->dualBound(), param().BapCodReducedCostTolerance())
                << std::endl
                << "  pricingSpSolutionValue " << pricingSpSolutionValue
                << " <0? " << (pricingSpSolutionValue.negative()) << std::endl
                << "  spRootReducedCost = " << _spRootReducedCost
                << " <0? " << (_spRootReducedCost.negative()) << std::endl;

    switch (_mastConfPtr->probPtr()->solMode().status())
    {
        case SolutionMethod::lpSolver:
        case SolutionMethod::mipSolver:
        case SolutionMethod::customSolver:
        case SolutionMethod::custom2mipSolver:
        {
            double dualValUB = 0;
            double dualValLB = 0;
            if (_upperBoundMastConstrPtr != NULL)
            {
                if (scaleFactor > 0)
                    /// floor or ceil??
                    dualValUB = floor(_upperBoundMastConstrPtr->valOrSepPointVal()._val * scaleFactor) / scaleFactor;
                else
                    dualValUB = _upperBoundMastConstrPtr->valOrSepPointVal();
            }
            if (_lowerBoundMastConstrPtr != NULL)
            {
                if (scaleFactor > 0)
                    dualValLB = floor(_lowerBoundMastConstrPtr->valOrSepPointVal()._val * scaleFactor) / scaleFactor;
                else
                    dualValLB = _lowerBoundMastConstrPtr->valOrSepPointVal();
            }
            bool pricingSpSolutionValueIsNegative = (scaleFactor > 0 ? pricingSpSolutionValue._val < 0
                                                                     : pricingSpSolutionValue.negative());

            if (pricingSpSolutionValueIsNegative) {
                if (_upperBoundMastConstrPtr != NULL) {
                    _mult = _upperBoundMastConstrPtr->curRhs();
                    _dualBoundContrib += (pricingSpSolutionValue + dualValUB) * _mult;
                } else {
                    _mult = BapcodInfinity;
                    /// Multiplying by pricingSpSolutionValue gives the correct sign
                    _dualBoundContrib += pricingSpSolutionValue * _mult;
                }
                /// also count correction in dual objective value
                if (_lowerBoundMastConstrPtr != NULL) {
                    _dualBoundContrib += dualValLB * _lowerBoundMastConstrPtr->curRhs();
                }
            }
            else {
                if (_lowerBoundMastConstrPtr != NULL) {
                    _mult = _lowerBoundMastConstrPtr->curRhs();
                    _dualBoundContrib += (pricingSpSolutionValue + dualValLB) * _mult;
                } else {
                    _mult = 0; // - BapcodInfinity;
                    //_dualBoundContrib += _mult * pricingSpSolutionValue;
                }
                /// also count correction in dual objective value
                if (_upperBoundMastConstrPtr != NULL) {
                    _dualBoundContrib += dualValUB * _upperBoundMastConstrPtr->curRhs();
                }
            }
            break;
        }
        case SolutionMethod::none:
            break;
        case SolutionMethod::undefined:
            bapcodInit().check(true,
                               "ColGenSpConf::updateTarget(: ERROR undefined solution method");
            break;
    }


  if (printL(5))
    {
      if (_upperBoundMastConstrPtr != NULL)
        std::cout << "ub constr name " << _upperBoundMastConstrPtr->name()
                  << " val " << _upperBoundMastConstrPtr->valOrSepPointVal()
                  << " rhs " << _upperBoundMastConstrPtr->curRhs() << std::endl;

      if (_lowerBoundMastConstrPtr != NULL)
        std::cout << "constraint name " << _lowerBoundMastConstrPtr->name()
                  << " val " << _lowerBoundMastConstrPtr->valOrSepPointVal()
                  << " rhs " << _lowerBoundMastConstrPtr->curRhs() << std::endl;

      std::cout << "  pricingSpSolutionValue = " << pricingSpSolutionValue << " mult " << _mult
                << " SP_dualBoundContrib " << _dualBoundContrib << std::endl;
    }

  return;
}

void ColGenSpConf::correctDualboundContribAndBestSol(const Bound & rootDualBound, Bound & bestRedCost,
                                                     Solution * & bestSolPtr)
{
  /// Compute corrected dualbound contrib and select best reduced cost column
  Double dbObjContribOfBSBR(0);
  Double dbCorrectionsOfDeltaPricingSpZetaVal(0);
  Double replacement(0);
  Double verifDbCorrections(0);

  /**
   * Remove dualVal contrib of branching constraints and
   * update dualContrib with solved SP from the treeOfColClasses
   */
  for (ColClassesVector::const_iterator bcitPt = treeOfColClasses().begin();
      bcitPt != treeOfColClasses().end(); bcitPt++)
  {
    dbObjContribOfBSBR += (*bcitPt)->val() * (*bcitPt)->curRhs();

    if (printL(5))
      std::cout << "constraint " << (*bcitPt)->name() << std::endl
                << "  val = " << (*bcitPt)->val() << std::endl
                << "  sigma = " << (*bcitPt)->sigma() << std::endl
                << "  reducedCost = " << (*bcitPt)->reducedCost() << std::endl
                << "  zeta = " << (*bcitPt)->pricingSpZetaVal() << std::endl
                << "  rhs = " <<  (*bcitPt)->curRhs() << std::endl
                << "  marginLvalue = " <<(*bcitPt)->margLvalue4DualBd() << std::endl
                << "current dbObjContribOfBSBR = " << dbObjContribOfBSBR << std::endl;

    if (!((*bcitPt)->associatedPricingSPsolved()))
      continue;


    replacement = (*bcitPt)->pricingSpZetaVal() - rootDualBound;
    dbCorrectionsOfDeltaPricingSpZetaVal += replacement * (*bcitPt)->margLvalue4DualBd();
    if (printL(5))
      std::cout << "  replacement * marginalL = " << replacement
                << " * " << (*bcitPt)->margLvalue4DualBd()  << std::endl
                << "current correction = " << dbCorrectionsOfDeltaPricingSpZetaVal << std::endl;

    /// And select best reduced cost column
    if ((*bcitPt)->reducedCost() < bestRedCost)
    {
      bestRedCost = (*bcitPt)->reducedCost();
      if (bestSolPtr != NULL)
      {
        delete bestSolPtr;
        bestSolPtr = NULL;
      }
      bestSolPtr = (*bcitPt)->solPtr();
      bestSolPtr = bestSolPtr->clone();
    }
  }

  _dualBoundContrib += dbObjContribOfBSBR + dbCorrectionsOfDeltaPricingSpZetaVal;

  if (printL(5))
  {
    std::cout << " bestRedCost " << bestRedCost << ", new _dualBoundContrib = " << _dualBoundContrib << std::endl;

    if (bestSolPtr != NULL)
      std::cout << " bestSol " << *bestSolPtr;
  }

  return;
}

const Double & ColGenSpConf::fixedCost()  const
{
  return _fixedCost;
}

const Bound & ColGenSpConf::target()  const
{
  return _target;
}

const Double & ColGenSpConf::fixedDualCost()  const
{
  return _fixedDualCost;
}


void ColGenSpConf::target(const Bound & t)
{
  _target = t;

  return;
}

std::ostream& ColGenSpConf::print(std::ostream& os) const
{
  os << "ColGenSpConf "  << std::endl;

  ProbConfig::print(os);

  os << "  _fixedCost = " <<  _fixedCost << std::endl;
  os << "  _fixedDualCost = " <<  _fixedDualCost << std::endl;
  os << "  _target = " <<  _target << std::endl;

  if (_lowerBoundMastConstrPtr != NULL)
    os << "  cur LB = " << _lowerBoundMastConstrPtr->curRhs() << std::endl;

  if (_upperBoundMastConstrPtr != NULL)
    os << "  cur UB = " << _upperBoundMastConstrPtr->curRhs() << std::endl;

  return(os);
}

