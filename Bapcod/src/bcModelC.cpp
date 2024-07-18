/**
 *
 * This file bcModelC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcModelC.hpp"
#include "bcDoubleC.hpp"
#include "bcModelObjectiveC.hpp"
#include "bcModelPointerC.hpp"
#include "bcNodeC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcVarConstrC.hpp"

/**
 * Generic code
 */
using namespace std;

Model::Model(BapcodInit * bapcodInitPtr,
	         const std::string & modelName,
	         const BcObjStatus::MinMaxIntFloat & objectiveSense):
  _bapcodInitPtr(bapcodInitPtr),
  _solutionFoundCallbackPtr(NULL),
  _bcModelPtr(new BcModel(this)),
  _modelMasterCnt(0),
  _modelOvfCnt(0),
  _modelColGenSpCnt(0),
  _modelConstrCnt(0),
  _modelVarCnt(0),
  _modelMastColCnt(0),
  _modelProblemCnt(0),
  _modelBrConstrGenCnt(0),
  _modelGenVarConstrCnt(0),
  _sideProbConfPts(),
  _masterConfPtr(NULL),
  _curColGenSpConfPtr(NULL),
  _mapOfColGenSpConf(),
  _objectiveSense(objectiveSense),
  _modelIsSetup(false),
  _garbageCollector(),
  _initialPrimalBound(Bound::infPrimalBound(objectiveSense)),
  _initialDualBound(Bound::infDualBound(objectiveSense)),
  _modelName(modelName),
  _artVarCost(0.0)
{
    _bapcodInitPtr->statistics().incrValue("bcRecBestInc", Bound::infPrimalBound(objectiveSense));

    _bapcodInitPtr->modelPtr(this);
}

Model::~Model()
{
    for (auto * sideProbConfPtr : _sideProbConfPts)
        delete sideProbConfPtr;
    _sideProbConfPts.clear();

    if (_masterConfPtr != NULL)
    {
        if (_masterConfPtr->ovfConfPtr() != NULL)
            delete  _masterConfPtr->ovfConfPtr();
        delete  _masterConfPtr;
    }
    _masterConfPtr = NULL;

    _bcModelPtr->clearModel();   // by Artur/Romain
    delete _bcModelPtr;

    if (_solutionFoundCallbackPtr != NULL)
        delete _solutionFoundCallbackPtr;

    if (_bapcodInitPtr != NULL)
        _bapcodInitPtr->modelPtr(NULL);

    if (printL(2))
        std::cout << "Destructed the model" << std::endl;
    return;
}

void Model::increaseModelConstrCnt()
{
  _modelConstrCnt += 1;
}

void Model::increaseModelVarCnt()
{
    _modelVarCnt += 1;
}

void Model::increaseModelMastColCnt()
{
    _modelMastColCnt += 1;
}

void Model::increaseModelBrConstrGenCnt()
{
    _modelBrConstrGenCnt += 1;
}

void Model::increaseModelProblemCnt()
{
    _modelProblemCnt += 1;
}

void Model::increaseModelGenVarConstrCnt()
{
    _modelGenVarConstrCnt += 1;
}

const BcObjStatus::MinMaxIntFloat & Model::objectiveSense() const
{
    return _objectiveSense;
}

void Model::objectiveSense(const BcObjStatus::MinMaxIntFloat & newObjectiveSense)
{
    _objectiveSense = newObjectiveSense;
}

void Model::setup()
{
  if (_modelIsSetup)
      return;

  _modelIsSetup = true;
  prepareModel();

  return;
}

Solution * Model::enumerateAllColumns(int & nbEnumColumns)
{
    setup();

    int saveParamValue = param().RCSPmaxNumOfLabelsInHeurEnumeration();
    param().RCSPmaxNumOfLabelsInHeurEnumeration(param().RCSPmaxNumOfLabelsInEnumeration());

    Solution * solPtr = _masterConfPtr->enumerateAllColumns(nbEnumColumns);

    param().RCSPmaxNumOfLabelsInHeurEnumeration(saveParamValue);

    return solPtr;
}

Solution * Model::solve()
{

 if (printL(1))
   std::cout << "NEXT PROBLEM " << _modelName << std::endl;

  setup();

  if (!progStatus().doRun())
  {
     std::cerr << "BaPCod error : cannot build the model" << std::endl;
     return NULL;
  }

  Solution * solutionPtr(NULL);

  if ((ovfConfPtr() != NULL) && (param().ovfSolMode().status() !=  SolutionMethod::none))
  {
    solutionPtr = ovfConfPtr()->solvePC();

    _masterConfPtr->updatePrimalIncBound(ovfConfPtr()->primalIncBound());
    _masterConfPtr->updateDualIncBound(ovfConfPtr()->dualIncBound());
  }

  if ((masterConfPtr() != NULL) && (param().masterSolMode().status() !=  SolutionMethod::none))
  {
    solutionPtr = _masterConfPtr->solvePC();
  }

  if (!gapSmallerThanTol(_masterConfPtr->dualIncBound(), _masterConfPtr->primalIncBound(), param()))
  {
    bapcodInit().statistics().incrCounter("bcFailToSolveModel");
  }

  return(solutionPtr);
}


void Model::mastConfPtr(MasterConf * MCptr)
{
    _masterConfPtr = MCptr;
}

void Model::prepareModel()
{
  if (_artVarCost == 0.0)
  {
      if (printL(-1))
          std::cout << "BaPCod WARNING : artificial variable cost is not set, setting it to 1e+6" << std::endl;
      std::cerr << "BaPCod WARNING : artificial variable cost is not set, setting it to 1e+6" << std::endl;
      _artVarCost = 1e+6;
  }
  /// Need to prepare master and subproblems first
  if (_masterConfPtr != NULL)
    {
      _masterConfPtr->prepareProbConfig();

      /// Ovf conf is indeed generated from master and ColGenSpConf
      if (ovfConfPtr() != NULL)
	    ovfConfPtr()->prepareProbConfig();
    }
}

std::ostream& Model::print(std::ostream& os) const
{
  os << "Model: " << std::endl;
  os << "  modelName: " << _modelName << std::endl;
  if (ovfConfPtr() != NULL)
    ovfConfPtr()->print(os);

  if (masterConfPtr() != NULL)
    masterConfPtr()->print(os);

  return(os);
}

std::ostream& Model::printSol(std::ostream& os) const
{
  os << "Model: " << _modelName << std::endl;
  os << "MASTER SOL" << std::endl;
  if (masterConfPtr() != NULL)
    masterConfPtr()->print(os);

  os << "OVF SOL" << std::endl;

  if (ovfConfPtr() != NULL) ovfConfPtr()->print(os);

  return(os);
}


BapcodInit & Model::bapcodInit() const
{
  return *_bapcodInitPtr;
}

BapcodInit* Model::bapcodInitPtr() const
{
  return _bapcodInitPtr;
}

const ProgStatus& Model::progStatus() const
{
  return _bapcodInitPtr->progStatus();
}

PlainVarConstrPointerSet& Model::garbageCollector()
{
  return _garbageCollector;
}

ProgStatus& Model::progStatus()
{
  return _bapcodInitPtr->progStatus();
}

const BcModel & Model::bcModel() const
{
  return (*_bcModelPtr);
}

BcModel & Model::bcModel()
{
  return (*_bcModelPtr);
}


void Model::setArtCostValue(const Double & costAprioriEvaluation)
{
  _artVarCost = costAprioriEvaluation;
}

Variable * Model::recordSubproblemSolution(ColGenSpConf * cgSpConfPtr, Solution * solPtr,
					                       const int & insertionLevel) const
{
  return cgSpConfPtr->recordSubproblemSolution(solPtr, false, insertionLevel);
}


OvfConf * Model::createOvfConf(const std::string & genericName, const IndexCell & id)
{
  std::string aname(genericName);

  if (id == MultiIndex())
    {
      MultiIndex newid(_modelOvfCnt++);
      newid.appendRef2name(aname);
    }
  else
    {
      id.appendRef2name(aname);
    }

  if (ovfConfPtr() == NULL)
    {
        OvfConf * ptr = new OvfConf(this,
                                    new MipProblem(modelProblemCnt(),
                                                   param().MasterMipSolverRightHandSideZeroTol(),
                                                   param().MasterMipSolverReducedCostTolerance(),
                                                   _objectiveSense,
                                                   param().ovfSolMode,
                                                   aname,
                                                   param().RequiredSolStatForOvf,
                                                   param().PreprocessorOnForOvf,
                                                   param().ProbingOnForOvf,
                                                   param().RequiredSolStatForOvf,
                                                   param().PreprocessorOnForOvf,
                                                   param().ProbingOnForOvf,
                                                   param().AutomaticCuttingPlanesOnForOvf,
                                                   param().SolverSelectForOvf));

      if (!ovfConfPtr(ptr))
          delete ptr;
    }

  return ovfConfPtr();
}


MasterConf * Model::createMasterConf(const std::string & genericName, const IndexCell & id)
{
    if (_masterConfPtr == NULL)
    {
        std::string aname(genericName);
        if (id == MultiIndex())
        {
            MultiIndex newid(_modelMasterCnt++);
            newid.appendRef2name(aname);
        }
        else
        {
            id.appendRef2name(aname);
        }

        Problem * probPtr = NULL;

        if (param().masterSolMode().status() == SolutionMethod::mipSolver)
        {

            probPtr = new MipProblem(modelProblemCnt(),
                                     param().MasterMipSolverRightHandSideZeroTol(),
                                     param().MasterMipSolverReducedCostTolerance(),
                                     _objectiveSense, param().masterSolMode(),
                                     aname,
                                     param().RequiredSolStatForMast(),
                                     false,
                                     false,
                                     param().RequiredSolStatForMast(),
                                     param().PreprocessorOnForMast(),
                                     param().ProbingOnForMast(),
                                     param().AutomaticCuttingPlanesOnForMast(),
                                     param().SolverSelectForMast());
        }
        else
        {
            probPtr = new Problem(modelProblemCnt(),
                                  param().MasterMipSolverRightHandSideZeroTol(),
                                  param().MasterMipSolverReducedCostTolerance(),
                                  _objectiveSense,
                                  param().masterSolMode(),
                                  aname,
                                  param().RequiredSolStatForMast(),
                                  param().PreprocessorOnForMast(),
                                  false,
                                  param().SolverSelectForMast());
        }
        increaseModelProblemCnt();

        _masterConfPtr = new MasterConf(this, probPtr, _initialPrimalBound, _initialDualBound);

        if (param().ovfSolMode().status() !=  SolutionMethod::none)
            createOvfConf("ovf", IndexCell());

        _masterConfPtr->createDefaultGenericVarConstr();
    }
    return _masterConfPtr;
}


ProbConfig * Model::createProbConf(const std::string & genericName, const IndexCell & id)
{

  std::string aname(genericName);

  if (id == MultiIndex())
    {
      MultiIndex newid(_modelOvfCnt++);
      newid.appendRef2name(aname);
    }
  else
    {
      id.appendRef2name(aname);
    }

  Problem * problemPtr = NULL;

  problemPtr = new MipProblem(modelProblemCnt(),
                              param().MasterMipSolverRightHandSideZeroTol(),
                              param().MasterMipSolverReducedCostTolerance(),
                              _objectiveSense,
			                  param().masterSolMode(), //SolutionMethod::mipSolver,
			                  aname,
                              param().RequiredSolStatForOvf(),
                              param().PreprocessorOnForOvf(),
                              param().ProbingOnForOvf(),
                              param().RequiredSolStatForOvf(),
                              param().PreprocessorOnForOvf(),
                              param().ProbingOnForOvf(),
                              param().AutomaticCuttingPlanesOnForOvf(),
                              param().SolverSelectForOvf());
  increaseModelProblemCnt();


  ProbConfig * pcPtr = new ProbConfig(ProbConfig::sideProb,
				                      this,
				                      genericName,
				                      id,
				                      Bound::infPrimalBound(_objectiveSense),
				                      Bound::infDualBound(_objectiveSense),
				                      problemPtr);

  _sideProbConfPts.push_back(pcPtr);

  return pcPtr;
}

ColGenSpConf * Model::createColGenSubproblem(const std::string & spname,
						                     const MultiIndex & id,
						                     const bool & implicitlyFixCardinality,
						                     const Double & upperBound,
						                     const Double & lowerBound,
						                     const Double & fixedCost,
						                     const Double & defaultDualVal)
{
  std::string aspname(spname);

  if (id == MultiIndex())
    {
      MultiIndex newid(_modelColGenSpCnt++);
      newid.appendRef2name(aspname);
    }

  Double initDualVal(_artVarCost);
  if (defaultDualVal > 1)
      initDualVal = defaultDualVal;

  Double * maxNbOfSpSolUsedPtr = new Double(upperBound);
  Double * minNbOfSpSolUsedPtr = new Double(lowerBound);

    _curColGenSpConfPtr = new ColGenSpConf(spname,
                                           id,
                                           _masterConfPtr,
                                           fixedCost,
                                           implicitlyFixCardinality,
                                           maxNbOfSpSolUsedPtr,
                                           minNbOfSpSolUsedPtr,
                                           initDualVal,
                                           - initDualVal,
                                           new MipProblem(modelProblemCnt(),
                                                          param().ColGenSpMipSolverRightHAndSideZeroTol(),
                                                          param().ColGenSpMipSolverReducedCostTolerance(),
                                                          _objectiveSense,
                                                          param().colGenSubProbSolMode,
                                                          aspname, 0,
                                                          param().PreprocessorOnForColGenSp,
                                                          param().ProbingOnForColGenSp,
                                                          param().RequiredSolStatForColGenSp,
                                                          param().PreprocessorOnForColGenSp,
                                                          param().ProbingOnForColGenSp,
                                                          param().AutomaticCuttingPlanesOnForColGenSp,
                                                          param().SolverSelectForColGenSp));

    increaseModelProblemCnt();
    _curColGenSpConfPtr->createDefaultGenericVarConstr();

    if (printL(6))
        std::cout << "createColGenSubproblem(): confPtr = " << aspname << std::endl;

  _mapOfColGenSpConf[spname][id] = _curColGenSpConfPtr;

  return _curColGenSpConfPtr;
}

GenericVar * Model::createGenericVar(ProbConfig * probConfigPtr,
					                 const BcVarConstrType::BcVcType & vctype,
					                 const std::string & name,
					                 const MultiIndexNames & multiIndexNames,
					                 const char & type,
					                 const Double & cost,
					                 const Double & ub,
					                 const SelectionStrategy & priorityRule,
					                 const Double & genericBranchingOnAggregateVarPL,
					                 const Double & compBoundSetBranchingPL,
					                 const char & flag, /// 's' for static, 'd' for dynamic
					                 const char & sense,
                                     int firstIndexMax, int secondIndexMax, int thirdIndexMax)
{
  GenericVar * gvPtr = new GenericVar(this, vctype, probConfigPtr, name, multiIndexNames, type, cost, ub,
                                      priorityRule, genericBranchingOnAggregateVarPL, compBoundSetBranchingPL, flag,
                                      firstIndexMax, secondIndexMax, thirdIndexMax);

  if (printL(5))
    std::cout << " Model::createGenericVar() = " << gvPtr << std::endl;

  if (probConfigPtr != NULL)
    {
      if (printL(5))
       	std::cout << " inserted "  << std::endl;
      probConfigPtr->insertGenericVar(gvPtr);
    }

  gvPtr->defaultLb((sense == 'P'? 0: -BapcodInfinity));
  gvPtr->defaultGlobalLb((sense == 'P'? 0: -BapcodInfinity));

  gvPtr->defaultSense(sense);

  return gvPtr;
}

GenericConstr * Model::createGenericConstr(ProbConfig * probConfigPtr,
                                           const BcVarConstrType::BcVcType & vctype,
                                           const std::string & name,
                                           const MultiIndexNames & multiIndexNames,
                                           const char & sense,
                                           const Double & rhs,
                                           const Double & val,
                                           const bool & toBeUsedInPreprocessing,
                                           const char & flag,
                                           const char & type,
                                           const char & kind,
                                           const SelectionStrategy & priority,
                                           const Double & priorityLevel)
/// type  'E' for explicit (placed in the MIP formulation),
/// 'I' for implicit (only used for preprocessing),
/// 'S' for constraints defining a subsystem in column generation for extended formulation approach
{
  if (printL(6))
    std::cout << " Model::createGenericConstr() : GenConstr =  " << name << std::endl;

  GenericConstr * gcPtr =  new GenericConstr(this, vctype, probConfigPtr, name, multiIndexNames, sense, rhs,
                                             priority, priorityLevel, toBeUsedInPreprocessing, flag);

  if (probConfigPtr != NULL)
    probConfigPtr->insertGenericConstr(gcPtr);

  Double defVal(val);
  /// if stabilization is used, we do not change initialization of default incumbent values for constraints
  if ((val == 0) && (param().colGenStabilizationFunctionType().getStatusAsInteger() == 0)
      && (param().colGenDualPriceSmoothingAlphaFactor() == 0))
    {
      if (sense == 'L')
        defVal = _artVarCost;
      else
        defVal = -_artVarCost;
    }
  gcPtr->defaultVal(defVal);

  gcPtr->defaultType(type);
  gcPtr->defaultKind(kind);
  return gcPtr;
}

GenericCutConstr * Model::createGenericCut(ProbConfig * probConfigPtr,
						                   const std::string & name,
                                           const char & type,
						                   const SelectionStrategy & separationPriorityRule,
						                   const Double & nonRootPriorityLevel,
						                   const Double & rootPriorityLevel,
						                   const char & sense,
						                   const Double & rhs,
						                   const bool & toBeUsedInPreprocessing)
/// type  'C' for core (required for the IP formulation,
/// 'F' for facultative (only helpfull for the LP formulation),
/// 'S' for constraints defining a subsystem in column generation for extended formulation approach
{
  GenericCutConstr * gcPtr =  new GenericCutConstr(this, _masterConfPtr, name, type, separationPriorityRule,
						                           nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing);

  if (probConfigPtr != NULL)
    gcPtr->probConfPtr(probConfigPtr);

  gcPtr->defaultName(name);
  gcPtr->defaultSense(sense);
  gcPtr->defaultCostRhs(rhs);
  gcPtr->defaultFlag('d');
  if (type == 'S')
    gcPtr->defaultVal(0);
  else
    gcPtr->defaultVal((sense == 'L' ? _artVarCost : - _artVarCost));
  return gcPtr;
}

GenericCustomNonLinearCutConstr *
Model::createGenericCustomNonLinearCutConstr(BcCustomNonLinearCutArrayFunctor * functorPtr,
                                             ProbConfig * probConfigPtr,
                                             const std::string & name,
                                             const char & type,
					                         const SelectionStrategy & separationPriorityRule,
                                             const Double & nonRootPriorityLevel,
						                     const Double & rootPriorityLevel,
                                             const char & sense,
					                         const Double & rhs)
{
  GenericCustomNonLinearCutConstr * gcPtr
    =  new GenericCustomNonLinearCutConstr(this, _masterConfPtr, name, type, separationPriorityRule,
                                           nonRootPriorityLevel, rootPriorityLevel, functorPtr);
  if (probConfigPtr != NULL)
    gcPtr->probConfPtr(probConfigPtr);

  gcPtr->defaultName(name);
  gcPtr->defaultSense(sense);
  gcPtr->defaultCostRhs(rhs);
  gcPtr->defaultFlag('d');
  gcPtr->defaultVal((sense == 'L' ? _artVarCost : - _artVarCost));
  return gcPtr;
}

GenericSoftConflictsCutConstr *
Model::createGenericSoftConflictsCutConstr(BcSoftConflictsCutArrayFunctor * functorPtr,
                                           ProbConfig * probConfigPtr,
                                           const std::string & name,
                                           const char & type,
                                           const SelectionStrategy & separationPriorityRule,
                                           const Double & nonRootPriorityLevel,
                                           const Double & rootPriorityLevel,
                                           const char & sense,
                                           const Double & rhs)
{
    GenericSoftConflictsCutConstr * gcPtr
            =  new GenericSoftConflictsCutConstr(this, _masterConfPtr, name, type, separationPriorityRule,
                                                 nonRootPriorityLevel, rootPriorityLevel, functorPtr);
    if (probConfigPtr != NULL)
        gcPtr->probConfPtr(probConfigPtr);

    gcPtr->defaultName(name);
    gcPtr->defaultSense(sense);
    gcPtr->defaultCostRhs(rhs);
    gcPtr->defaultFlag('d');
    gcPtr->defaultVal((sense == 'L' ? _artVarCost : - _artVarCost));
    return gcPtr;
}

 GenericBranchingConstr * Model::createGenericBranching(ProbConfig * probConfigPtr,
                                                       const std::string & name,
							                           const char & type,
							                           const SelectionStrategy & separationPriorityRule,
                                                       const Double & nonRootPriorityLevel,
						                               const Double & rootPriorityLevel,
							                           const char & sense,
                                                       const Double & rhs,
							                           const bool & toBeUsedInPreprocessing)
/// type  'C' for core (required for the IP formulation,
/// 'F' for facultative (only helpfull for the LP formulation),
/// 'S' for constraints defining a subsystem in column generation for extended formulation approach
{
  GenericBranchingConstr * gcPtr =  new GenericBranchingConstr(this, _masterConfPtr, name,
                                                               separationPriorityRule, nonRootPriorityLevel,
                                                               rootPriorityLevel, toBeUsedInPreprocessing);

  if (probConfigPtr != NULL)
    {
       gcPtr->probConfPtr(probConfigPtr);
    }

  gcPtr->defaultName(name);
  gcPtr->defaultSense(sense);
  gcPtr->defaultCostRhs(rhs);
  gcPtr->defaultFlag('d');
  if (type == 'S')
    gcPtr->defaultVal(0);
  else
    gcPtr->defaultVal((sense == 'L' ? _artVarCost : - _artVarCost));

  return gcPtr;
}

InstanciatedConstr * Model::createConstraint(ProbConfig * probConfPtr,
                                             GenericConstr * genConstrPtr,
                                             const MultiIndex & id,
                                             const Double & rhs,
                                             const char & sense,
                                             const Double & val,
                                             const std::string & name,
                                             const bool & toBeUsedInPreprocessing,
                                             const bool & considerAsEqualityInPreprocessing)
{


  if (printL(6))
      std::cout << " Model::addConstraint: adding  constraint Name = " <<  name
			   << " GenConstrName = " << genConstrPtr->defaultName() << std::endl;

  InstanciatedConstr * iConstrPtr =  genConstrPtr->getConstrPtr(id);

  if (iConstrPtr == NULL)
    {
      std::string cname(name);

      id.appendRef2name(cname, genConstrPtr->multiIndexNames());
      if (probConfPtr != NULL) probConfPtr->id().appendRef2name(cname);

      iConstrPtr =  genConstrPtr->newInstanciation(IndexCell(id), probConfPtr, cname, rhs, sense,
                                                   genConstrPtr->defaultType(), genConstrPtr->defaultKind(),
                                                   genConstrPtr->defaultFlag(), val, BapcodInfinity, - BapcodInfinity,
                                                   'U', 1.0, true, toBeUsedInPreprocessing,
                                                   considerAsEqualityInPreprocessing);

      if (probConfPtr != NULL)
		  probConfPtr->insertPureConstr(iConstrPtr);

      genConstrPtr->addConstrPtr2MultiIndexMap(id, iConstrPtr);

    }

  return iConstrPtr;

}

InstanciatedConstr * Model::createConstraint(ProbConfig * probConfPtr,
                                             GenericConstr * genConstrPtr,
                                             const MultiIndex & id,
                                             const Double & rhs,
                                             const char & sense,
                                             const Double & val)
{
  return createConstraint(probConfPtr, genConstrPtr, id, rhs, sense, val, genConstrPtr->defaultName(),
                          genConstrPtr->toBeUsedInPreprocessing(), genConstrPtr->considerAsEqualityInPreprocessing());
}

InstanciatedConstr * Model::createConstraint(ProbConfig * probConfPtr,
                                             GenericConstr * genConstrPtr,
                                             const MultiIndex & id,
                                             const Double & rhs,
                                             const char & sense)
{
  return createConstraint(probConfPtr, genConstrPtr, id, rhs, sense, genConstrPtr->defaultVal());
}

InstanciatedConstr * Model::createConstraint(ProbConfig * probConfPtr,
                                             GenericConstr * genConstrPtr,
                                             const MultiIndex & id,
                                             const Double & rhs)
{
  return createConstraint(probConfPtr, genConstrPtr, id, rhs, genConstrPtr->defaultSense());
}


InstanciatedConstr * Model::createConstraint(ProbConfig * probConfPtr,
                                             GenericConstr * genConstrPtr,
                                             const MultiIndex & id)
{
  return createConstraint(probConfPtr, genConstrPtr, id, genConstrPtr->defaultCostRhs());
}


InstanciatedVar * Model::createVariable(ProbConfig * probConfPtr,
					                    GenericVar * genVarPtr,
					                    const MultiIndex & id,
					                    const Double & cost,
					                    const char & type,
					                    const std::string & name,
					                    const Double & ub,
					                    const Double & lb,
					                    const Double & priority,
					                    const char & directive,
					                    const Double & globalUb,
					                    const Double & globalLb,
					                    const char & kind,
					                    const char & sense,
					                    const char & flag)
{
  InstanciatedVar * iVarPtr =  genVarPtr->getVarPtr(id);

  if (iVarPtr == NULL)
    {
      if (printL(6))
          std::cout << " Model::addVariable: adding VariableName = " <<  name << " id " << id
			        << " GenVarName = " << genVarPtr->defaultName() << " lb = " << lb
			        << " ub = " << ub << " globalLb = " << globalLb << " globalUb = " << globalUb << std::endl;
      std::string vname("");

	  vname = name;
	  id.appendRef2name(vname, genVarPtr->multiIndexNames());
	  if (probConfPtr != NULL)
          probConfPtr->id().appendRef2name(vname);

      iVarPtr =  genVarPtr->newInstanciation(IndexCell(id), probConfPtr, vname, cost, sense, type, kind, ub, lb,
					                         flag, directive, priority, 0, globalUb, globalLb,true);

      if (probConfPtr != NULL)
          probConfPtr->insertPureVar(iVarPtr);

      genVarPtr->addVarPtr2MultiIndexMap(id, iVarPtr);

      if (flag == 'd')
          probConfPtr->castAndAddVariable(iVarPtr, true);

    }

  return iVarPtr;

}

InstanciatedVar * Model::createVariable(ProbConfig * probConfPtr,
                                        GenericVar * genVarPtr,
                                        const MultiIndex & id,
                                        const Double & cost,
                                        const char & type,
                                        const Double & ub,
                                        const Double & lb,
                                        const Double & priority,
                                        const char & directive,
                                        const Double & globalUb,
                                        const Double & globalLb,
                                        const char & kind,
                                        const char & sense,
                                        const char & flag)
{
  return createVariable(probConfPtr, genVarPtr, id, cost, type, genVarPtr->defaultName(), ub, lb, priority, directive,
                        globalUb, globalLb, kind, sense, flag);

}

InstanciatedVar * Model::createVariable(ProbConfig * probConfPtr,
                                        GenericVar * genVarPtr,
                                        const MultiIndex & id,
                                        const Double & cost,
                                        const char & type)
{
  return createVariable(probConfPtr, genVarPtr, id, cost, type, genVarPtr->defaultName(), genVarPtr->defaultUb(),
                        genVarPtr->defaultLb(),1, genVarPtr->defaultDirective(), genVarPtr->defaultGlobalUb(),
                        genVarPtr->defaultGlobalLb(), genVarPtr->defaultKind(), genVarPtr->defaultSense(),
                        genVarPtr->defaultFlag());
}

InstanciatedVar * Model::createVariable(ProbConfig * probConfPtr,
                                        GenericVar * genVarPtr,
                                        const MultiIndex & id,
                                        const Double & cost)
{
  return createVariable(probConfPtr, genVarPtr, id, cost, genVarPtr->defaultType());
}


InstanciatedVar * Model::createVariable(ProbConfig * probConfPtr,
					                    GenericVar * genVarPtr,
					                    const MultiIndex & id)
{
  return createVariable(probConfPtr, genVarPtr, id, genVarPtr->defaultCostRhs());
}

bool Model::addCoefficient(Constraint * constrPtr, Variable * varPtr, const Double & coef)
{

  if (constrPtr == NULL)
      return true;
  if (varPtr == NULL)
      return true;

  if (printL(6))
      std::cout << "Model::addCoefficient: constraintName = " <<  constrPtr->name()
			    <<  " variableName = " << varPtr->name() << " coef = " << coef << std::endl;

  constrPtr->includeMember(varPtr, coef, true);

  return false;
}

MasterConf * Model::master() const
{
//  if (_masterConfPtr == NULL)
//    {
//      std::cout <<"BaPCod info : Model::master(): _masterConfPtr == NULL" << std::endl;
//    }
  return  _masterConfPtr;
}

OvfConf * Model::ovf() const
{
  if (_masterConfPtr == NULL)
    {
      std::cerr << "ERROR Model::ovf(): _masterConfPtr == NULL" << std::endl;
      exit(1);
    }

  return  _masterConfPtr->ovfConfPtr();
}

ColGenSpConf * Model::getColGenSubproblem(const std::string & name, const MultiIndex & id)
{
  if (_mapOfColGenSpConf.empty())
    {
      return NULL;
    }

  if (_curColGenSpConfPtr != NULL)
    {
      if ((_curColGenSpConfPtr->genericName() == name) && (_curColGenSpConfPtr->id() == id))
	return _curColGenSpConfPtr;
    }

  MapColGenSpConfPtrByNameAndMultiIndex::const_iterator itName = _mapOfColGenSpConf.find(name);
  if (itName == _mapOfColGenSpConf.end())
    {
      return NULL;
    }

  MapColGenSpConfPtrByMultiIndex::const_iterator itId = itName->second.find(id);

  if (itId == itName->second.end())
    {
      if (printL(6))
          std::cout << "  Model::getColGenSubproblem(" << id <<  ") could not be found in table" << std::endl;
      return NULL;
    }
  if (printL(6))
      std::cout << "  Model::getColGenSubproblem(" << id <<  ") was found in table" << std::endl;

  _curColGenSpConfPtr = itId->second;

  return _curColGenSpConfPtr;
}

void Model::getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions)
{
    for (auto & namePair : _mapOfColGenSpConf)
        for (auto & indexPair : namePair.second)
            indexPair.second->probPtr()->getEnumeratedSolutions(enumSolutions);
}

 bool Model::checkIfSolutionIsFeasibleUsingCallback(Solution *solPtr) const
 {
     if (_solutionFoundCallbackPtr == NULL)
         return true;
     Solution * dissagrSolPtr = NULL;
     if (_solutionFoundCallbackPtr->disagregateSolutionBeforeProcessing())
         dissagrSolPtr = _masterConfPtr->getDissagregatedSolution(solPtr);
     else
         dissagrSolPtr = solPtr->clone();

     bool solIsFeasible = (*_solutionFoundCallbackPtr)(BcSolution(dissagrSolPtr));

     if (!solIsFeasible && printL(0))
         solPtr->printDetailedSolution(std::cout);

     dissagrSolPtr->deleteSolutionsChain();
     delete dissagrSolPtr;

     return solIsFeasible;
 }



