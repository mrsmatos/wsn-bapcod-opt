/**
 *
 * This file bcGenVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcGlobalException.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcFormC.hpp"
#include "bcModelMasterC.hpp"
#include "bcModelCutConstrC.hpp"

/**
 * Generic code
 */
using namespace std;

/**
 * Sorting procedure to select among cuts of branching constraints: 
 * by genericVarConstr _priorityLevel than by constraint PriorityRule 
 * and priority of violation or fract val.
 */
bool CutSeparationPriorityComp::operator()(InstanciatedConstr * a, InstanciatedConstr * b) const
{
  if (a->genVarConstrPtr() != b->genVarConstrPtr())
  {
    a->bapcodInit().require(a->genVarConstrPtr() != NULL,
                            "CutSeparationPriorityComp::operator() a->genVarConstr should be != NULL");

    a->bapcodInit().require(b->genVarConstrPtr() != NULL,
                            "CutSeparationPriorityComp::operator() b->genVarConstr should be != NULL");

    if (a->genVarConstrPtr()->priorityLevel() > b->genVarConstrPtr()->priorityLevel())
      return true;

    if (a->genVarConstrPtr()->priorityLevel() < b->genVarConstrPtr()->priorityLevel())
      return false;
  }

  /// Constraints are from the same generic class or different class but with the same priority level
  if (a->genVarConstrPtr()->priorityRule() == b->genVarConstrPtr()->priorityRule())
  {
    switch (a->genVarConstrPtr()->priorityRule())
    {
      case SelectionStrategy::FirstFound:
      {
        return (a->ref() < b->ref());
        break;
      }
      case SelectionStrategy::HighestPriority:
      {

        return (a->priority() > b->priority());
        break;
      }
      case SelectionStrategy::MostViolated:
      {
        return (a->violation() > b->violation());
        break;
      }
      default:
        break;
    }
  }

  /// Default case
  return VcRefST(a, b);
}

/***************************************************************************
 ***************   TODO: Methods for GenericVarConstr   ********************
 ***************************************************************************/

GenericVarConstr::GenericVarConstr(Model * modelPtr,
                                   const BcVarConstrType::BcVcType & gvcType,
                                   ProbConfig * probConfPtr,
                                   const std::string & name,
                                   const MultiIndexNames & multiIndexNames,
                                   const SelectionStrategy & priorityRule,
                                   const Double & priorityLevel,
                                   const bool & toBeUsedInPreprocessing) :
    _gvcType(gvcType), _probConfigPtr(probConfPtr), _ref(modelPtr->modelGenVarConstrCnt()), _model(modelPtr),
    _priorityRule(priorityRule), _priorityLevel(priorityLevel), _toBeUsedInPreprocessing(toBeUsedInPreprocessing),
    _considerAsEqualityInPreprocessing(false), _defaultName(name),
    _defaultType('C'), _defaultKind('E'), _defaultFlag('s'), _defaultSense('P'),
    _defaultDirective('U'), _defaultCostRhs(0), _defaultUb(99999),
    _defaultLb(0), _defaultVal(0), _multiIndexNames(multiIndexNames), _dimension(-1)
{
  modelPtr->increaseModelGenVarConstrCnt();
  switch (gvcType)
  {
    case BcVarConstrType::globalSumOfLocals:
      break;
    case BcVarConstrType::globalCloneOfLocals:
      break;
    case BcVarConstrType::local2Formulation:
    default:
    {
      if (probConfPtr == NULL)
      {
        std::cerr << " GenericVarConstr without probConfPtr specification should be global " << std::endl;
        exit(1);
      }
    }
  };

  return;
}

GenericVarConstr::~GenericVarConstr()
{
}

void GenericVarConstr::priorityLevel(const Double & pl)
{
  if (pl <= 0)
    _priorityRule = SelectionStrategy::NotConsideredForSelection;

  _priorityLevel = pl;
}

const Double & GenericVarConstr::genericCost(const InstanciatedVar * const ivarconstrPtr) const
{
  return (ivarconstrPtr->Variable::costrhs());
}

const Double & GenericVarConstr::genericRhs(const InstanciatedConstr * const ivarconstrPtr) const
{
  return (ivarconstrPtr->Constraint::rhs());
}

BcModel & GenericVarConstr::bcModel()
{
  if (modelPtr() == NULL)
  {
    throw GlobalException("GenericVarConstr::modelPtr: _modelPtr is null");
  }
  return (modelPtr()->bcModel());
}

const BcModel & GenericVarConstr::bcModel() const
{
  if (modelPtr() == NULL)
  {
    throw GlobalException("GenericVarConstr::modelPtr: _modelPtr is null");
  }
  return (modelPtr()->bcModel());
}

std::ostream& GenericVarConstr::print(std::ostream& os) const
{
  return (os << "GenericVarConstr" << std::endl);
}

BapcodInit & GenericVarConstr::bapcodInit() const
{
  return _model->bapcodInit();
}

BapcodInit * GenericVarConstr::bapcodInitPtr() const
{
  return _model->bapcodInitPtr();
}

const ControlParameters & GenericVarConstr::param() const
{
  return _model->param();
}

ControlParameters & GenericVarConstr::param()
{
  return _model->param();
}

/***************************************************************************
 *********************   TODO: Methods for GenericVar   ********************
 ***************************************************************************/

GenericVar::GenericVar(Model * modelPtr,
                       const BcVarConstrType::BcVcType & gvctype,
		               ProbConfig * probConfPtr,
                       const std::string & name,
                       const MultiIndexNames & multiIndexNames,
                       const char & type, const Double & cost,
                       const Double & ub,
                       const SelectionStrategy & branchingPriorityRule,
                       const Double & genericBranchingOnAggregateVarPL,
                       const Double & compBoundSetBranchingPL,
                       const char & flag,
                       int firstIndexMax,int secondIndexMax,int thirdIndexMax) :
    GenericVarConstr(modelPtr,
                     gvctype,
                     probConfPtr,
                     name,
                     multiIndexNames,
                     branchingPriorityRule,
                     genericBranchingOnAggregateVarPL,
                     false),
    _genericBranchingOnAggregateVarPL(genericBranchingOnAggregateVarPL),
    _compBoundSetBranchingPL(compBoundSetBranchingPL),
    _ryanFosterBranchingPL(-1),
    _brCmpPtr(NULL),
    _compBoundSetGenBranchConstrPtr(NULL),
    _target4MostFractionalBranchingPriority(0.5),
    _isDirectAccessForInstanciatedVars(false)
{
  if (modelPtr == NULL)
  {
    std::cout << "GenericVar::GenericVar(): model * must be defined"
              << std::endl;
  }
  defaultType(type);
  defaultKind('E');
  defaultCostRhs(cost);
  defaultUb(ub);
  defaultLb(0);
  defaultGlobalUb(BapcodInfinity);
  defaultGlobalLb(0);
  defaultVal(0);
  defaultFlag(flag);
  defaultDirective('U');

  //if at least one max index is given, it means that we want to use direct access for instanciated vars
  if (firstIndexMax != -1)
  {
    _isDirectAccessForInstanciatedVars = true;
    if(secondIndexMax != -1)
    {
      if(thirdIndexMax != -1)
      {
        _varPtrVec3D = vector<vector<vector<InstanciatedVar*> > > (firstIndexMax, 
                           vector<vector<InstanciatedVar*> > (secondIndexMax,
                                    vector<InstanciatedVar*> (thirdIndexMax)));        
      }      
      else
      {
        _varPtrVec2D = vector<vector<InstanciatedVar*> > (firstIndexMax,
                            vector<InstanciatedVar*> (secondIndexMax));  
      }
    }
    else
    {
      _varPtrVec1D = vector<InstanciatedVar*> (firstIndexMax);
    }
  }

  return;
}

void GenericVar::target4MostFractionalBranchingPriority(const Double & target)
{
  _target4MostFractionalBranchingPriority = target;
}

void GenericVar::setupGenericBranchingConstr()
{
  if (printL(5))
    std::cout
      << " GenericVar::setupGenericBranchingConstr() : GenericVar " << defaultName() 
      << " branchingPriorityRule = "
      << priorityRule() << std::endl 
      << " genericBranchingOnAggregateVar priority level "
      << _genericBranchingOnAggregateVarPL
      << " compBoundSetBranching priority level " 
      << _compBoundSetBranchingPL
      << " ryanFosterBranching priority level "
      << _ryanFosterBranchingPL
      << std::endl;

  if ((defaultType() != 'B') && (defaultType() != 'I'))
    {
       priorityRule(SelectionStrategy::NotConsideredForSelection);
    }

  /// Setup the generic branching constraint on this variable
  else if (priorityRule() != SelectionStrategy::NotConsideredForSelection)
    {
      if (printL(5))
	std::cout
          << " var is for Branching -> GenericVar branchingPriorityRule = "
          << priorityRule() << std::endl;

      GenVarGenBranchConstr * locVarBranchingConstrPtr(NULL);

      /// *Be carrefull*, Order is important
      if (_genericBranchingOnAggregateVarPL.positive())
	{
          locVarBranchingConstrPtr = new GenVarGenBranchConstr(modelPtr(), modelPtr()->masterConfPtr(),
              this, prioritySelectionRule(), _genericBranchingOnAggregateVarPL);

	  if (printL(3))
	    std::cout << " genericBranchingOnAggregateVar "
		      << locVarBranchingConstrPtr << std::endl;

	  _locVarBranchingConstrPtrList.push_back(locVarBranchingConstrPtr);
	}

      if (_compBoundSetBranchingPL.positive())
	{
          _compBoundSetGenBranchConstrPtr = new CompBoundSetGenBranchConstr(modelPtr(), this, prioritySelectionRule(),
                                                                            _compBoundSetBranchingPL);
	  _locVarBranchingConstrPtrList.push_back(_compBoundSetGenBranchConstrPtr);
	}

      if (_ryanFosterBranchingPL.positive())
        {
          locVarBranchingConstrPtr = new RyanAndFosterGenBranchConstr(modelPtr(), this, prioritySelectionRule(),
                                                                      _ryanFosterBranchingPL);

          if (printL(3))
            std::cout << " ryanFosterBranching " << locVarBranchingConstrPtr
                      << std::endl;

          _locVarBranchingConstrPtrList.push_back(locVarBranchingConstrPtr);
        }

      switch (priorityRule())
	{
      case SelectionStrategy::FirstFound:
        _brCmpPtr = new BrVarPriorityCalcAndComp_FirstFound;
        break;
      case SelectionStrategy::HighestPriority:
        _brCmpPtr = new BrVarPriorityCalcAndComp_HighestPriority;
        break;
      case SelectionStrategy::FracWeightedPriority:
        _brCmpPtr = new BrVarPriorityCalcAndComp_FracWeightedPriority;
        break;
      case SelectionStrategy::MostFractional:
        _brCmpPtr = new BrVarPriorityCalcAndComp_MostFractional;
        break;
      case SelectionStrategy::LeastFractional:
        _brCmpPtr = new BrVarPriorityCalcAndComp_LeastFractional;
        break;
      case SelectionStrategy::Closest2RoundUp:
        _brCmpPtr = new BrVarPriorityCalcAndComp_Closest2RoundUp;
        break;
      case SelectionStrategy::Closest2RoundDown:
        _brCmpPtr = new BrVarPriorityCalcAndComp_Closest2RoundDown;
        break;
      case SelectionStrategy::NotConsideredForSelection:
        break;
      default:
        bapcodInit().check(
            true,
            "GenericVar::GenericVar(): error cannot separate fract sol on var that is not for branching");
    }

  }

  return;
}

GenericVar::~GenericVar()
{
  if (_brCmpPtr != NULL)
    delete _brCmpPtr;
}

void GenericVar::addVarPtr2MultiIndexMap(const MultiIndex & id, InstanciatedVar * ivarPtr)
{
  if(_isDirectAccessForInstanciatedVars)
  {
    if (id.endPosition == 1)
    {
      _varPtrVec1D[id.first()] = ivarPtr;
    } else if (id.endPosition == 2)
    {
      _varPtrVec2D[id.first()][id.second()] = ivarPtr;
    } else if (id.endPosition == 3)
    {
      _varPtrVec3D[id.first()][id.second()][id.third()] = ivarPtr;
    }   
  }
  else
  {
    _multiIndex2VarPtrMap[id] = ivarPtr;
  }   
}

InstanciatedVar * GenericVar::getVarPtr(const MultiIndex & id) const
{
  if(_isDirectAccessForInstanciatedVars)
    {
      if (id.endPosition == 1)
        {
          return _varPtrVec1D[id.first()];
        }
      else if (id.endPosition == 2)
        {
          return _varPtrVec2D[id.first()][id.second()];
        }
      else if (id.endPosition == 3)
        {
          return _varPtrVec3D[id.first()][id.second()][id.third()];
        }
      return NULL;
    }
    
  MapInstVarPtrByMultiIndex::const_iterator it = _multiIndex2VarPtrMap.find(id);
  if (it == _multiIndex2VarPtrMap.end())
    {
      return NULL;
    }
  return (it->second);
}

bool GenericVar::consecutiveVarWhenBrOnCBS(InstanciatedVar * var1Ptr, InstanciatedVar * var2Ptr) const
{
  return true;
}

bool GenericVar::genericCount(const InstanciatedConstr * const iconstrPtr, const InstanciatedVar * const ivarPtr) const
{
  if (iconstrPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    return (false);

  return (iconstrPtr->genVarConstrPtr()->genericCount(iconstrPtr, ivarPtr));
}

const LpCoef GenericVar::genericCoef(const InstanciatedConstr * const iconstrPtr,
                                     const InstanciatedVar * const ivarPtr) const
{
  if (iconstrPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
    return LpCoef::ZeroCoef;

  return (iconstrPtr->genVarConstrPtr()->genericCoef(iconstrPtr, ivarPtr));
}

void GenericVar::buildMembership(InstanciatedVar * ivPtr)
{
  if (printL(6))
    std::cout << "GenericVar::buildMembership was called" << std::endl;
  /// it is enough to set the flag to false in GenericConstr
  // ivPtr->presetMembership(false); do not change the flag; leave it to what it is in InstanciatedVar
}

std::ostream& GenericVar::print(std::ostream& os) const
{
  return (os << "GenericVar" << std::endl);
}

InstanciatedVar * GenericVar::newInstanciation(const IndexCell & id,
					                           ProbConfig * probConfigPtr,
					                           const std::string & name,
					                           const Double & costrhs,
					                           const char & sense,
					                           const char & type,
					                           const char & kind,
					                           const Double & upperBound,
					                           const Double & lowerBound,
					                           const char & flag,
					                           const char & directive,
					                           const Double & priority,
					                           const Double & val,
					                           const Double & globalUb,
					                           const Double & globalLb,
					                           const bool & presetMembership)
{
  InstanciatedVar * instVarPtr( NULL);

  if (printL(6))
    std::cout << " GenericVar::newInstanciation(): name = " << name << std::endl;

  if (bapcodInit().testLevel() >= 2)
    {
      instVarPtr = checkIfInstanciationAlreadyExist(id); /// iVar does not already exists

      if (instVarPtr != NULL)
	{
	  throw GlobalException(
				"GenericVar::newInstanciation(): error instanciation should not already exists",
				true);
	}
    }

  if (probConfigPtr != NULL)
    {

      switch (probConfigPtr->configType())
	{
	case ProbConfig::ovf:
	  {
	    {
	      throw GlobalException(
				    "GenericVar::newInstanciation(): error instanciation : InstOvfVar not defined yet",
				    true);
	    }
	    break;
	  }
	case ProbConfig::master:
	  {
	    MasterConf * masterConfPtr = static_cast<MasterConf *>(probConfigPtr);
	    instVarPtr = new InstMasterVar(id, 
					   this, 
					   masterConfPtr, 
					   name, 
					   costrhs,
					   sense, 
					   type, 
					   kind, 
					   upperBound,
					   lowerBound, 
					   flag, 
					   directive, 
					   priority,
					   val, 
					   globalUb, 
					   globalLb,
					   presetMembership);
	    break;
	  }
	case ProbConfig::colGenSp:
	  {
	    ColGenSpConf * cgSpConfPtr = static_cast<ColGenSpConf *>(probConfigPtr);
	    instVarPtr = new SubProbVariable(cgSpConfPtr->mastConfPtr(), 
					     id, 
					     this,
					     probConfigPtr, 
					     name, 
					     costrhs, 
					     sense,
					     type, 
					     kind, 
					     upperBound, 
					     lowerBound,
					     flag, 
					     directive, 
					     priority, 
					     val,
					     globalUb, 
					     globalLb, 
					     presetMembership);
	    break;
	  }
	default:
	  {
	    instVarPtr = new InstanciatedVar(id, 
					     this, 
					     probConfigPtr, 
					     name, 
					     costrhs,
					     sense, 
					     type, 
					     kind, 
					     upperBound,
					     lowerBound, 
					     flag, 
					     directive, 
					     priority,
					     val, 
					     globalUb, 
					     globalLb,
					     presetMembership);
	  }
	}
    }
  else
    {
      instVarPtr = new InstanciatedVar(id, 
				       this, 
				       probConfigPtr, 
				       name, 
				       costrhs,
				       sense, 
				       type, 
				       kind, 
				       upperBound, 
				       lowerBound,
				       flag, 
				       directive, 
				       priority, 
				       val, 
				       globalUb,
				       globalLb, 
				       presetMembership);
    }
  if (printL(6))
    std::cout << "GenericVar::createNewInstanciation() created " << instVarPtr->name() << std::endl;

  return instVarPtr;
}

InstanciatedVar * GenericVar::checkIfInstanciationAlreadyExist(const IndexCell & id)
{
  IndexCell2InstancVarPtrMap::const_iterator it;

  it = _indexCell2InstancVarPtrMap.find(id);
  if (it != _indexCell2InstancVarPtrMap.end())
    {
      if (printL(6))
	    std::cout << "checkIfInstanciationAlreadyExist exists  " << it->second->name() << ", id = " << id << std::endl;

      return it->second;
    }

  return NULL;

}

void GenericVar::recordInstanciation(InstanciatedVar * iPtr)
{
  _indexCell2InstancVarPtrMap[iPtr->id()] = iPtr;
}

void GenericVar::deleteInstanciation(InstanciatedVar * iPtr)
{
  _indexCell2InstancVarPtrMap.erase(iPtr->id());

  return;
}

/***************************************************************************
 *********************   TODO: Methods for GenericConstr   *****************
 ***************************************************************************/

GenericConstr::GenericConstr(Model * modelPtr,
                             const BcVarConstrType::BcVcType & gvcType,
                             ProbConfig * probConfPtr, const std::string & name,
                             const MultiIndexNames & multiIndexNames,
                             const char & sense, const Double & rhs,
                             const SelectionStrategy & priority,
                             const Double & priorityLevel,
                             const bool & toBeUsedInPreprocessing,
                             const char & flag) :
    GenericVarConstr(modelPtr, gvcType, probConfPtr, name, multiIndexNames,
                     priority, //SelectionStrategy::NotConsideredForSelection,
                     priorityLevel, // -1: negative priority level  => not for cut/br separation,
                     toBeUsedInPreprocessing),
    _addConstrFunctorPtr(NULL)
{
  if (modelPtr == NULL)
  {
    std::cout << "GenericConstr::GenericConstr(): model * must be defined"
              << std::endl;
  }
  defaultSense(sense);
  defaultCostRhs(rhs);
  defaultFlag(flag);

  _multiIndex2ConstrPtrMap.max_load_factor(0.1);

  return;
}

GenericConstr::~GenericConstr()
{
}

void GenericConstr::addConstrPtr2MultiIndexMap(const MultiIndex & id, InstanciatedConstr * iconstrPtr)
{
  _multiIndex2ConstrPtrMap[id] = iconstrPtr;
}

InstanciatedConstr * GenericConstr::getConstrPtr(const MultiIndex & id) const
{
  MapInstConstrPtrByMultiIndex::const_iterator it = _multiIndex2ConstrPtrMap.find(id);
  if (it == _multiIndex2ConstrPtrMap.end())
  {
    return NULL;
  }
  return (it->second);
}

void GenericConstr::addConstrFunctorPtr(BcAddConstrFunctor * addConstrRoutinePtr)
{
    _addConstrFunctorPtr = addConstrRoutinePtr;
}

InstanciatedConstr * GenericConstr::newInstanciation(const IndexCell & id,
                                                     ProbConfig * probConfigPtr,
                                                     const std::string& name,
                                                     const Double& rhs,
                                                     const char& sense,
                                                     const char& type,
                                                     const char& kind,
                                                     const char& flag,
                                                     const Double& val,
                                                     const Double& upperBound,
                                                     const Double& lowerBound,
                                                     const char & directive,
                                                     const Double & priority,
                                                     const bool & presetMembership,
                                                     const bool & toBeUsedInPreprocessing,
                                                     const bool & considerAsEqualityInPreprocessing)
{
  if (printL(6))
    std::cout << " GenericConstr::newInstanciation(): name = " << name << std::endl;

  InstanciatedConstr * instConstrPtr( NULL);

  if (bapcodInit().testLevel() >= 2)
  {
    instConstrPtr = checkIfInstanciationAlreadyExist(id);

    if (instConstrPtr != NULL)
    {
      throw GlobalException(
          "GenericConstr::newInstanciation(): error instanciation should not already exists",
          true);
    }
  }

  if (probConfigPtr != NULL)
  {
    switch (probConfigPtr->configType())
    {
      case ProbConfig::master:
      {
        MasterConf * masterConfPtr = static_cast<MasterConf *>(probConfigPtr);
        if (masterConfPtr != NULL)
        {

          instConstrPtr = new InstMasterConstr(id, this, probConfigPtr, name,
                                               rhs, sense, type, kind, flag,
                                               val, upperBound, lowerBound,
                                               directive, priority,
                                               presetMembership,
                                               toBeUsedInPreprocessing,
                                               considerAsEqualityInPreprocessing);
        }
        break;
      }
      case ProbConfig::colGenSp:
      case ProbConfig::ovf:
      default:
      {
        instConstrPtr = new InstanciatedConstr(id, this, probConfigPtr, name,
                                               rhs, sense, type, kind, flag,
                                               val, upperBound, lowerBound,
                                               directive, priority,
                                               presetMembership,
                                               toBeUsedInPreprocessing,
                                               considerAsEqualityInPreprocessing);
      }
    }
  }
  else
  {

    instConstrPtr = new InstanciatedConstr(id, this, probConfigPtr, name, rhs,
                                           sense, type, kind, flag, val,
                                           upperBound, lowerBound, directive,
                                           priority, presetMembership,
                                           toBeUsedInPreprocessing,
                                           considerAsEqualityInPreprocessing);
  }

  if (printL(6))
    std::cout << "GenericConstr::createNewInstanciation() created "
              << instConstrPtr->name() << std::endl;

  return instConstrPtr;
}

InstanciatedConstr * GenericConstr::checkIfInstanciationAlreadyExist(const IndexCell & id)
{
  IndexCell2InstancConstrPtrMap::const_iterator it = _indexCell2InstancConstrPtrMap.find(id);
  if (it != _indexCell2InstancConstrPtrMap.end())
  {
    if (printL(6))
      std::cout << "checkIfInstanciationAlreadyExist exists  " << it->second->name() << std::endl;

    return it->second;
  }

  return NULL;
}

void GenericConstr::recordInstanciation(InstanciatedConstr * iPtr)
{
  _indexCell2InstancConstrPtrMap[iPtr->id()] = iPtr;

  return;
}

void GenericConstr::deleteInstanciation(InstanciatedConstr * iPtr)
{
  _multiIndex2ConstrPtrMap.erase(iPtr->id());
  _indexCell2InstancConstrPtrMap.erase(iPtr->id());

  return;
}

bool GenericConstr::genericCount(const InstanciatedConstr * const iconstrPtr,
                                 const InstanciatedVar * const ivarPtr) const
{

  return (genericCoef(iconstrPtr, ivarPtr)).first;
}

void GenericConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  if (printL(6))
    std::cout << "GenericConstr::buildMembership has been called" << std::endl;
}

std::ostream& GenericConstr::print(std::ostream& os) const
{
  return (os << "GenericConstr" << std::endl);
}

void GenericConstr::nicePrintAllConstraints(std::ostream& os) const
{
  os << "Printing all constraints of GenericConstr " << defaultName() << std::endl;
  for (IndexCell2InstancConstrPtrMap::const_iterator mapIt = _indexCell2InstancConstrPtrMap.begin();
       mapIt != _indexCell2InstancConstrPtrMap.end(); ++mapIt)
  {
    mapIt->second->nicePrint(os);
  }
}

/***************************************************************************
 ************   TODO: Methods for Base4NonLinearGenericConstr   ************
 ***************************************************************************/

bool Base4NonLinearGenericConstr::genericMastColumnCount(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (_nlGenericConstrPtr != NULL)
    return(_nlGenericConstrPtr->genericCount(icPtr,colPtr->spSol()->solVarValMap()));

  return(false);
}

const LpCoef Base4NonLinearGenericConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (_nlGenericConstrPtr != NULL)
    return(_nlGenericConstrPtr->genericCoef(icPtr,colPtr->spSol()->solVarValMap()));

  return LpCoef::ZeroCoef;
}


/***************************************************************************
 *****************   TODO: Methods for NonLinearGenericConstr   ************
 ***************************************************************************/

const LpCoef NonLinearGenericConstr::genericCoef(const InstanciatedConstr * const iconstrPtr,
                                                 const InstanciatedVar * const ivarPtr) const
{
  bapcodInit().check(true,
	                "NonLinearGenericConstr::genericCoef() should not be called", ProgStatus::run);
  
  return LpCoef::ZeroCoef;
}

std::ostream& NonLinearGenericConstr::print(std::ostream& os) const
{
  return(os << "NonLinearGenericConstr" << std::endl);
}

InstanciatedConstr * NonLinearGenericConstr::newInstanciation(const IndexCell & id, 
							      ProbConfig* probConfigPtr, 
							      const std::string& name, 
							      const Double& rhs, 
							      const char& sense, 
							      const char& type, 
							      const char& kind, 
							      const char& flag, 
							      const Double& val, 
							      const Double& upperBound, 
							      const Double& lowerBound, 
							      const char & directive, 
							      const Double & priority, 
							      const bool & presetMembership, 
							      const bool & toBeUsedInPreprocessing)
{
  InstanciatedConstr * instConstrPtr( NULL);

  if (bapcodInit().testLevel() >= 2)
    {
      instConstrPtr = checkIfInstanciationAlreadyExist(id); /// iConstr does not already exists

      bapcodInit().require(instConstrPtr == NULL,
                           "GenericConstr::newInstanciation(): error instanciation should not already exsit");
    }


  MasterConf * masterConfPtr = dynamic_cast<MasterConf * >(probConfigPtr);
  if (masterConfPtr != NULL)
    {

      instConstrPtr = new NonLinearInstMastConstr(id,
						  this, 
						  probConfigPtr,              
						  name,                
						  rhs,              
						  sense,                
						  type,                 
						  kind,                 
						  flag,                 
						  val,                  
						  upperBound, 
						  lowerBound,  
						  directive,            
						  priority); 
    }
  else
    {
      instConstrPtr = new NonLinearInstConstr(id,
					      this, 
					      probConfigPtr,              
					      name,                
					      rhs,              
					      sense,                
					      type,                 
					      kind,                 
					      flag,                 
					      val,                  
					      upperBound, 
					      lowerBound,  
					      directive,            
					      priority); 
    }

  if (printL(5))
    std::cout << "GenericConstr::createNewInstanciation() create " << instConstrPtr->name() << std::endl;

  return instConstrPtr;
}

/***************************************************************************
 *****************   TODO: Methods for DynamicGenericConstr   **************
 ***************************************************************************/

DynamicGenericConstr::DynamicGenericConstr(Model * modelPtr,
                                           ProbConfig * probConfPtr,
                                           const std::string & name,
                                           const char & type,
                                           const SelectionStrategy & priorityRule,
                                           const Double & nonRootPriorityLevel,
                                           const Double & rootPriorityLevel,
                                           const bool & toBeUsedInPreprocessing) :
  GenericConstr(modelPtr, BcVarConstrType::local2Formulation,
                probConfPtr, name, MultiIndexNames(), 'G', 0, priorityRule,
                nonRootPriorityLevel, toBeUsedInPreprocessing, 'd'),
  _oracleDefined(false), _type(type), _constrPrototypes(), _rootPriorityLevel(rootPriorityLevel)
{
}

DynamicGenericConstr::~DynamicGenericConstr()
{
  std::list<InstanciatedConstr *>::iterator instConstrPtrIt;
  for (instConstrPtrIt = _constrPrototypes.begin();
       instConstrPtrIt != _constrPrototypes.end(); ++instConstrPtrIt)
  {
    (*instConstrPtrIt)->decrParticipation(10);
    if (printL(-1) && ((*instConstrPtrIt)->participation() != 0))
    {
      std::cout << "BaPCod warning : prototype constraint participation is not zero at destruction" << std::endl;
    }
    delete *instConstrPtrIt;
  }
  _constrPrototypes.clear();

  return;
}

void DynamicGenericConstr::rootPriorityLevel(const double & value)
{
    _rootPriorityLevel = Double(value);
}

bool DynamicGenConstrSort::operator()(DynamicGenericConstr * a, DynamicGenericConstr * b) const
{
  if (a->priorityLevel() > b->priorityLevel())
    return true;

  if (a->priorityLevel() < b->priorityLevel())
    return false;

  return (a->ref() < b->ref());
}

bool DynamicGenConstrSortAtRoot::operator()(DynamicGenericConstr * a, DynamicGenericConstr * b) const
{
  if (a->rootPriorityLevel() > b->rootPriorityLevel())
    return true;

  if (a->rootPriorityLevel() < b->rootPriorityLevel())
    return false;

  return (a->ref() < b->ref());
}

void DynamicGenericConstr::nicePrintAllConstraints(std::ostream& os) const
{
  if (_constrPrototypes.empty())
  {
    os << "Separation of DynamicGenericConstr " << defaultName() << " is based on an oracle" << std::endl;
  }
  else
  {
    os << "Printing all prototype constraints of DynamicGenericConstr " << defaultName() << std::endl;
    for (std::list<InstanciatedConstr *>::const_iterator listIt = _constrPrototypes.begin();
         listIt != _constrPrototypes.end(); ++listIt)
    {
      (*listIt)->nicePrint(os);
    }
  }
}

/***************************************************************************
 *****************   TODO: Methods for GenericCutConstr   ******************
 ***************************************************************************/

GenericCutConstr::GenericCutConstr(Model * modelPtr, ProbConfig * probConfPtr,
                                   const std::string & name, const char & type,
                                   const SelectionStrategy & priorityRule,
                                   const Double & nonRootPriorityLevel,
                                   const Double & rootPriorityLevel,
                                   const bool & toBeUsedInPreprocessing) :
    DynamicGenericConstr(modelPtr, probConfPtr, name, type, priorityRule,
                         nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing),
    _cutSeparationFunctorPtr(NULL)
{
  if (probConfPtr != NULL)
    probConfPtr->insertGenericCutConstr(this);

  return;
}

GenericCutConstr::~GenericCutConstr()
{
  return;
}

void GenericCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  if (printL(6))
    std::cout << "GenericCutConstr::buildMembership has been called"
              << std::endl;
}

void GenericCutConstr::cutSeparationFunctorPtr(BcCutSeparationFunctor * cutSepfunctPtr)
{
    _cutSeparationFunctorPtr = cutSepfunctPtr;
    _oracleDefined = true;
}

bool GenericCutConstr::prepareSeparation()
{
  return true;
}

void GenericCutConstr::cutSeparationRoutine(const VarPtrSet & curSol,
                                            std::multiset<InstanciatedConstr *,
                                                          CutSeparationPriorityComp> & generatedCutConstrSet)
{
  Double maxViolation = param().BapCodCutViolationTolerance();
  int generatedCutCounter(0);

  if (!_oracleDefined)
  {

    for (std::list<InstanciatedConstr *>::const_iterator constrIt = _constrPrototypes.begin();
         constrIt != _constrPrototypes.end(); constrIt++)
    {
      if (printL(6))
        std::cout << "constrPrototypes = " << (*constrIt)->name() << std::endl;

#ifdef BC_MORE_TIMERS
      Time start;
#endif
      Double lhs((*constrIt)->computeLhs(curSol));
#ifdef BC_MORE_TIMERS
      bapcodInit().statistics().incrTimer("bcTimeSepComputLHS", start.getElapsedTime_dbl());
#endif

      switch (priorityRule())
      {
        case SelectionStrategy::FirstFound:
        {
          (*constrIt)->computeViolation(lhs);
          if ((*constrIt)->violation() > param().BapCodCutViolationTolerance())
          {
            if (printL(6))
              std::cout
                  << "GenericCutConstr::cutSeparationRoutine() FirstFound cut/br constraint "
                  << (*constrIt)->name() << "  violation = "
                  << (*constrIt)->violation() << std::endl;

            generatedCutConstrSet.insert(*constrIt);
            generatedCutCounter++;
          }
          break;
        }
        case SelectionStrategy::MostFractional:
        {
          (*constrIt)->computeViolation(lhs);

          if ((*constrIt)->violation() < maxViolation)
            continue;

          if (printL(6))
            std::cout
                << "GenericCutConstr::cutSeparationRoutine() new  MostFractional cut/br constraint "
                << (*constrIt)->name() << "  violation = "
                << (*constrIt)->violation() << std::endl;

          generatedCutConstrSet.insert(*constrIt);
          maxViolation = (*constrIt)->violation();

          break;
        }
        default:
          std::cout
              << "GenericCutConstr::cutSeparationRoutine(): priorityRule not well defined"
              << std::endl;
      }
      if (printL(6))
        std::cout << "CutConstraint " << (*constrIt)->name() << " violation "
                  << (*constrIt)->violation() << std::endl;
    }
    return;
  }

  if (probConfPtr() == NULL)
    return;

  Solution * tempSolution = probConfPtr()->getSolution(curSol);
  Solution * primalSolPtr = probConfPtr()->getAggregatedSolution(tempSolution);
    
  if (printL(5))
    std::cout << "GenericCutConstr::cutSeparationRoutine: primalSol "
              << primalSolPtr << std::endl;

  std::list<BcConstr> cutList;

  BcSolution modelSol(primalSolPtr);

  generatedCutCounter = (*_cutSeparationFunctorPtr)(BcFormulation(modelPtr()->masterConfPtr()),
                                                    modelSol, maxViolation.val(),
                                                    cutList);
  modelPtr()->masterConfPtr()->addVariablesToForm();
  
  if (printL(5))
    std::cout << "GenericCutConstr::cutSeparationRoutine: generated CutConstraint "
              << generatedCutCounter << std::endl;

  if (generatedCutCounter >= 1)
    {
      for (std::list<BcConstr>::iterator cutIt = cutList.begin();
           cutIt != cutList.end(); ++cutIt)
        {
          if (printL(5))
            std::cout << "CutConstraint " << (*cutIt) << std::endl;
          generatedCutConstrSet.insert((InstanciatedConstr *) (*cutIt));
        }
    }

  delete tempSolution;
  primalSolPtr->deleteSolutionsChain();
  delete primalSolPtr;
  return;
}

void GenericCutConstr::cutSeparationBasedOnFixedSol(const VarPtr2DoubleMap & oldPartialSol,
                                                    const VarPtr2DoubleMap & fixedSol,
                                                    std::multiset < InstanciatedConstr * ,
                                                                    CutSeparationPriorityComp > & generatedCutConstrSet)
{
   return;
}

/***************************************************************************
 ******   TODO: Methods for GenericCustomNonLinearCutConstr   **************
 ***************************************************************************/

GenericCustomNonLinearCutConstr
::GenericCustomNonLinearCutConstr(Model * modelPtr,
                                  ProbConfig * probConfPtr,
                                  const std::string & name,
                                  const char & type,
                                  const SelectionStrategy & priorityRule,
                                  const Double & nonRootPriorityLevel,
                                  const Double & rootPriorityLevel,
                                  BcCustomNonLinearCutArrayFunctor * userFunctorPtr):
  GenericCutConstr(modelPtr, probConfPtr, name, type, SelectionStrategy::MostFractional,
                   nonRootPriorityLevel, rootPriorityLevel, false),
  Base4NonLinearGenericConstr(NULL), _userFunctorPtr(userFunctorPtr)
{
}

const LpCoef GenericCustomNonLinearCutConstr
             ::genericMastColumnCoef(InstanciatedConstr * icPtr,
                                     MastColumn * colPtr) const
{
  if (!icPtr->isTypeOf(VcId::CustomNonLinearCutConstrMask))
    return LpCoef(false, 0.0);
    
  return getMastColumnCoeff(static_cast<CustomNonLinearCut *>(icPtr), colPtr);
}

void GenericCustomNonLinearCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
    iconstrPtr->presetMembership(true);
    return;
}

void GenericCustomNonLinearCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> & generatedCutConstrSet)
{
  if (probConfPtr() == NULL)
    return;
  
  Solution * solPtr = probConfPtr()->getSolution(curSol);
  Solution * primalSolPtr = probConfPtr()->getAggregatedSolution(solPtr);
    
  if (printL(5))
    std::cout << "GenericCustomNonLinearCutConstr::cutSeparationRoutine: primalSol "
              << primalSolPtr << std::endl;

  std::list<BcCustomNonLinearCut> cutList;

  BcSolution projectedSol(primalSolPtr);

  std::list<std::pair<double, BcSolution> > columnsInSol;
  for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
      {
        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
        columnsInSol.push_back(std::make_pair(colPtr->val(), BcSolution(colPtr->spSol())));
      }

  int generatedCutCounter = _userFunctorPtr->cutSeparationRoutine(BcFormulation(modelPtr()->masterConfPtr()),
                                                                  projectedSol, columnsInSol,
                                                                  param().BapCodCutViolationTolerance(), cutList);
  delete solPtr;
  primalSolPtr->deleteSolutionsChain();
  delete primalSolPtr;
    
  if (printL(5))
    std::cout << "GenericCustomNonLinearCutConstr::cutSeparationRoutine: generated CutConstraint "
              << generatedCutCounter << std::endl;

  if (generatedCutCounter >= 1)
    {
      for (std::list<BcCustomNonLinearCut>::iterator cutIt = cutList.begin();
           cutIt != cutList.end(); ++cutIt)
        {
          if (printL(5))
            std::cout << "CutConstraint " << ((InstanciatedConstr *) *cutIt)  << std::endl;
          generatedCutConstrSet.insert((InstanciatedConstr *) *cutIt);
        }
    }
}

const LpCoef GenericCustomNonLinearCutConstr::getMastColumnCoeff(CustomNonLinearCut * cutPtr,
                                                                 MastColumn * colPtr) const
{
  double coeff = _userFunctorPtr->getCoefficient(BcCustomNonLinearCut(cutPtr),
                                                 BcSolution(colPtr->spSol()));
  if (coeff == 0.0)
    return LpCoef::ZeroCoef;
  return LpCoef(true, coeff);
}
    
GenericCustomNonLinearCutConstr::~GenericCustomNonLinearCutConstr()
{
}

/***************************************************************************
 ********   TODO: Methods for  GenericSoftConflictsCutConstr  **************
 ***************************************************************************/

GenericSoftConflictsCutConstr
::GenericSoftConflictsCutConstr(Model * modelPtr,
                                ProbConfig * probConfPtr,
                                const std::string & name,
                                const char & type,
                                const SelectionStrategy & priorityRule,
                                const Double & nonRootPriorityLevel,
                                const Double & rootPriorityLevel,
                                BcSoftConflictsCutArrayFunctor * softConflictsCutSepFunctorPtr):
        GenericCutConstr(modelPtr, probConfPtr, name, type, SelectionStrategy::MostFractional,
                         nonRootPriorityLevel, rootPriorityLevel, false),
        Base4NonLinearGenericConstr(NULL), _softConflictsCutSepFunctorPtr(softConflictsCutSepFunctorPtr),
        _genIndicVarPts()
{
}

GenericSoftConflictsCutConstr::~GenericSoftConflictsCutConstr()
{
}

bool GenericSoftConflictsCutConstr::prepareSeparation()
{
  for (std::vector< ColGenSpConf * >::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
       cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
  {
    _genIndicVarPts[*cgSpConfPtrIt] = modelPtr()->createGenericVar(*cgSpConfPtrIt, BcVarConstrType::local2Formulation,
                                                                   defaultName() + "V",
                                                                   MultiIndexNames('i', 'j'),
                                                                   'B', 0, 1);
  }
  return true;
}

const LpCoef GenericSoftConflictsCutConstr::genericMastColumnCoef(InstanciatedConstr * icPtr,
                                                                  MastColumn * colPtr) const
{
  if (!icPtr->isTypeOf(VcId::SoftConflictsCutConstrMask))
    return LpCoef(false, 0.0);

  return getMastColumnCoeff(static_cast<SoftConflictsCut *>(icPtr), colPtr);
}

void GenericSoftConflictsCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}


void GenericSoftConflictsCutConstr::updateSubprobemsWithIndicatorVarAndConstr(const std::list<BcConstr> & cutList)
{
  if (param().colGenSubProbSolMode().status() == SolutionMethod::customSolver)
    return;

  /// subproblems are solved by MIP, thus we add conflict indicator variables and linking constraints in subproblems
  for (std::list<BcConstr>::const_iterator constrIt = cutList.begin(); constrIt != cutList.end(); ++constrIt)
  {
    InstanciatedConstr * icPtr = (InstanciatedConstr *)(*constrIt);
    if (icPtr->isTypeOf(VcId::SoftConflictsCutConstrMask))
    {
      SoftConflictsCut * cutPtr = static_cast<SoftConflictsCut *>(icPtr);
      std::map<ColGenSpConf *, std::list<SubProbVariable *> > newIndicVarPts;
      std::map<ColGenSpConf *, std::list<Constraint *> > newIndicConstrPts;

      for (std::vector< ColGenSpConf * >::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
           cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
      {
        newIndicVarPts.insert(std::make_pair(*cgSpConfPtrIt, std::list<SubProbVariable *>()));
        newIndicConstrPts.insert(std::make_pair(*cgSpConfPtrIt, std::list<Constraint *>()));
      }

      /// now for each conflict we add the constraint linking the variables in conflict and the indicator constraint
      int conflictOrd = 0;
      std::vector<std::pair<SubProbVariable *, SubProbVariable *> >::const_iterator pairIt;
      for (pairIt = cutPtr->conflicts().begin(); pairIt != cutPtr->conflicts().end(); ++pairIt, ++conflictOrd)
      {
        ColGenSpConf * cgSpConfPtr = pairIt->first->cgSpConfPtr();
        GenericVar * indicGenVarPtr = _genIndicVarPts[cgSpConfPtr];

        if (cutPtr->cutType() == 0)
        {
          /// we first verify whether this subproblem already has the indicator variable
          MultiIndex conflictId(pairIt->first->id().first(), pairIt->second->id().first());
          InstanciatedVar *ivPtr = indicGenVarPtr->checkIfInstanciationAlreadyExist(conflictId);
          if (ivPtr == NULL)
          {
            std::list<Constraint *> &newIndicConstrList = newIndicConstrPts[cgSpConfPtr];
            std::list<SubProbVariable *> &newIndicVarList = newIndicVarPts[cgSpConfPtr];
            ivPtr = modelPtr()->createVariable(cgSpConfPtr, indicGenVarPtr, conflictId, 0, 'C');
            ivPtr = cgSpConfPtr->castAndAddVariable(ivPtr, false);
            SubProbVariable *spVarPtr = static_cast<SubProbVariable *>(ivPtr);
            newIndicVarList.push_back(spVarPtr);

            InstanciatedConstr * indicConstrPtr = modelPtr()->createConstraint(cgSpConfPtr,
                                                                              cgSpConfPtr->defaultGenericConstrPtr(),
                                                                              conflictId, 1, 'L', 0,
                                                                              defaultName() + "C");

            indicConstrPtr->includeMember(ivPtr, -1, false);
            indicConstrPtr->includeMember(pairIt->first, 1, false);
            indicConstrPtr->includeMember(pairIt->second, 1, false);
            newIndicConstrList.push_back(indicConstrPtr);
          }
        }
        else /// cutType == 1
        {
          MultiIndex indicVarId(cutPtr->id().first());
          InstanciatedVar * ivPtr = indicGenVarPtr->checkIfInstanciationAlreadyExist(indicVarId);
          if (ivPtr == NULL)
          {
            std::list<SubProbVariable *> & newIndicVarList = newIndicVarPts[cgSpConfPtr];
            ivPtr = modelPtr()->createVariable(cgSpConfPtr, indicGenVarPtr, indicVarId, 0, 'B');
            ivPtr = cgSpConfPtr->castAndAddVariable(ivPtr, false);
            SubProbVariable *spVarPtr = static_cast<SubProbVariable *>(ivPtr);
            newIndicVarList.push_back(spVarPtr);
          }

          std::list<Constraint *> &newIndicConstrList = newIndicConstrPts[cgSpConfPtr];
          MultiIndex indicConstrId(cutPtr->id().first(), pairIt->first->id().first(), pairIt->second->id().first());
          InstanciatedConstr * indicConstrPtr = modelPtr()->createConstraint(cgSpConfPtr,
                                                                             cgSpConfPtr->defaultGenericConstrPtr(),
                                                                             indicConstrId, 1, 'L', 0,
                                                                             defaultName() + "C");
          indicConstrPtr->includeMember(ivPtr, -1, false);
          indicConstrPtr->includeMember(pairIt->first, 1, false);
          indicConstrPtr->includeMember(pairIt->second, 1, false);
          newIndicConstrList.push_back(indicConstrPtr);
        }
      }

      /// pushing the indicator variables to the subproblem formulations
      for (std::map<ColGenSpConf *, std::list<SubProbVariable *> >::iterator mapIt = newIndicVarPts.begin();
           mapIt != newIndicVarPts.end(); ++ mapIt)
      {
        mapIt->first->probPtr()->addVarSet(mapIt->second, 1, 2);
      }

      /// pushing the indicator constraints to the subproblem formulations
      for (std::map<ColGenSpConf *, std::list<Constraint *> >::iterator mapIt = newIndicConstrPts.begin();
           mapIt != newIndicConstrPts.end(); ++ mapIt)
      {
        mapIt->first->probPtr()->addConstrSet(mapIt->second, 1, 2);
      }
    }
  }
}

void GenericSoftConflictsCutConstr
      ::cutSeparationBasedOnFixedSol(const VarPtr2DoubleMap & oldPartialSol, const VarPtr2DoubleMap & fixedSol,
                                     std::multiset < InstanciatedConstr * ,
                                                     CutSeparationPriorityComp > & generatedCutConstrSet)
{
  if (probConfPtr() == NULL)
    return;

  std::list<BcConstr> cutList;

  std::list<std::pair<double, BcSolution> > columnsInOldFixedSol;
  for (VarPtr2DoubleMap::const_iterator mapIt = oldPartialSol.begin(); mapIt != oldPartialSol.end(); ++mapIt)
    if (mapIt->first->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * colPtr = static_cast<MastColumn *>(mapIt->first);
      columnsInOldFixedSol.push_back(std::make_pair(mapIt->second, BcSolution(colPtr->spSol())));
    }

  std::list<std::pair<double, BcSolution> > columnsInNewFixedSol;
  for (VarPtr2DoubleMap::const_iterator mapIt = fixedSol.begin(); mapIt != fixedSol.end(); ++mapIt)
    if (mapIt->first->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * colPtr = static_cast<MastColumn *>(mapIt->first);
      columnsInNewFixedSol.push_back(std::make_pair(mapIt->second, BcSolution(colPtr->spSol())));
    }

  int generatedCutCounter
          = _softConflictsCutSepFunctorPtr->cutSeparationBasedOnFixedSol(BcFormulation(modelPtr()->masterConfPtr()),
                                                                         columnsInOldFixedSol, columnsInNewFixedSol,
                                                                         cutList);

  updateSubprobemsWithIndicatorVarAndConstr(cutList);

  if (printL(5))
    std::cout << "GenericSoftConflictsCutConstr::cutSeparationBasedOnFixedSol: generated CutConstraint "
              << generatedCutCounter << std::endl;

  if (generatedCutCounter >= 1)
  {
    for (std::list<BcConstr>::iterator cutIt = cutList.begin(); cutIt != cutList.end(); ++cutIt)
    {
      if (printL(5))
        std::cout << "CutConstraint " << ((InstanciatedConstr *) *cutIt)  << std::endl;
      generatedCutConstrSet.insert((InstanciatedConstr *) *cutIt);
    }
  }

}

void GenericSoftConflictsCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> & generatedCutConstrSet)
{
  if (probConfPtr() == NULL)
    return;

  std::list<BcConstr> cutList;

  std::list<std::pair<double, BcSolution> > columnsInSol;
  for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
      columnsInSol.push_back(std::make_pair(colPtr->val(), BcSolution(colPtr->spSol())));
    }

  std::list<std::pair<double, BcSolution> > columnsInFixedSol;
  const VarPtr2DoubleMap & partialSolution = probConfPtr()->probPtr()->partialSolution();
  for (VarPtr2DoubleMap::const_iterator mapIt = partialSolution.begin(); mapIt != partialSolution.end(); ++mapIt)
  {
    if (mapIt->first->isTypeOf(VcId::MastColumnMask))
    {
      MastColumn * colPtr = static_cast<MastColumn *>(mapIt->first);
      columnsInFixedSol.push_back(std::pair<double, BcSolution>(mapIt->second, BcSolution(colPtr->spSol())));
    }
  }


  int generatedCutCounter = _softConflictsCutSepFunctorPtr->cutSeparationRoutine(BcFormulation(modelPtr()->masterConfPtr()),
                                                                                 columnsInFixedSol, columnsInSol,
                                                                                 cutList);

  updateSubprobemsWithIndicatorVarAndConstr(cutList);


  if (printL(5))
    std::cout << "GenericSoftConflictsCutConstr::cutSeparationRoutine: generated CutConstraint "
              << generatedCutCounter << std::endl;

  if (generatedCutCounter >= 1)
  {
    for (std::list<BcConstr>::iterator cutIt = cutList.begin(); cutIt != cutList.end(); ++cutIt)
    {
      if (printL(5))
        std::cout << "CutConstraint " << ((InstanciatedConstr *) *cutIt)  << std::endl;
      generatedCutConstrSet.insert((InstanciatedConstr *) *cutIt);
    }
  }
}

const LpCoef GenericSoftConflictsCutConstr::getMastColumnCoeff(SoftConflictsCut * cutPtr,
                                                               MastColumn * colPtr) const
{

  if (cutPtr->cutType() == 0)
  {
    /// the coefficient is the number of conflicts in the column
    int numConflicts = 0;
    std::vector<std::pair<SubProbVariable *, SubProbVariable *> >::const_iterator pairIt;
    for (pairIt = cutPtr->conflicts().begin(); pairIt != cutPtr->conflicts().end(); ++pairIt)
    {
      SubProbVariable * const firstSpVarPtr = pairIt->first;
      SubProbVariable * const secondSpVarPtr = pairIt->second;
      if (colPtr->cgSpConfPtr() == firstSpVarPtr->cgSpConfPtr())
      {
        const VarPtr2DoubleMap & colValMap = colPtr->spSol()->solVarValMap();
        if ((colValMap.find(firstSpVarPtr) != colValMap.end()) && (colValMap.find(secondSpVarPtr) != colValMap.end()))
          numConflicts += 1;
      }
    }
    if (numConflicts > 0)
       return LpCoef(true, (double) numConflicts);
  }
  else /// cutType = 1
  {
    /// the coefficient is 1 if at least one conflict is detected in the column
    std::vector<std::pair<SubProbVariable *, SubProbVariable *> >::const_iterator pairIt;
    for (pairIt = cutPtr->conflicts().begin(); pairIt != cutPtr->conflicts().end(); ++pairIt)
    {
      SubProbVariable *const firstSpVarPtr = pairIt->first;
      SubProbVariable *const secondSpVarPtr = pairIt->second;
      if (colPtr->cgSpConfPtr() == firstSpVarPtr->cgSpConfPtr())
      {
        const VarPtr2DoubleMap &colValMap = colPtr->spSol()->solVarValMap();
        if ((colValMap.find(firstSpVarPtr) != colValMap.end()) && (colValMap.find(secondSpVarPtr) != colValMap.end()))
          return LpCoef(true, 1.0);
      }
    }
  }

  return LpCoef::ZeroCoef;
}
