/**
 *
 * This file bcGenBranchingConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcDoubleC.hpp"
#include "bcGenBranchingConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcModelBranchingConstrC.hpp"

using namespace std;

bool BranchingSeparationPriorityComp::operator()(BranchingConstrGenerator * bcgA, BranchingConstrGenerator * bcgB) const
{
  if (bcgA->genericBrConstrPtr() != bcgB->genericBrConstrPtr())
    {
      if (bcgA->genericBrConstrPtr() == NULL)
        {
          throw GlobalException("BranchingSeparationPriorityComp::operator() a->genVarConstr should be != NULL", true);
        }
      if (bcgB->genericBrConstrPtr() == NULL)
        {
          throw GlobalException("BranchingSeparationPriorityComp::operator() b->genVarConstr should be != NULL", true);
        }

      if (bcgA->genericBrConstrPtr()->priorityLevel() > bcgB->genericBrConstrPtr()->priorityLevel())
        return true;

      if (bcgA->genericBrConstrPtr()->priorityLevel() < bcgB->genericBrConstrPtr()->priorityLevel())
        return false;
    }

  /**
   * Constraints are from the same generic class
   * or different class but with the same priority level
   */
  if (bcgA->genericBrConstrPtr()->priorityRule() == bcgB->genericBrConstrPtr()->priorityRule())
    {
      switch (bcgA->genericBrConstrPtr()->priorityRule())
        {
        case SelectionStrategy::FirstFound:
        {
          return (bcgA->ref() < bcgB->ref());
          break;
        }
        case SelectionStrategy::HighestPriority:
        {
          if (bcgA->priority() > bcgB->priority())
            return true;
          if (bcgA->priority() < bcgB->priority())
            return false;
          break;
        }
        case SelectionStrategy::MostFractional:
        {
          if (bcgA->fracPart() > bcgB->fracPart())
            return true;
          if (bcgA->fracPart() < bcgB->fracPart())
            return false;
          break;
        }
        case SelectionStrategy::LeastFractional:
        {
          if (bcgA->fracPart() < bcgB->fracPart())
            return true;
          if (bcgA->fracPart() > bcgB->fracPart())
            return false;
          break;
        }
        case SelectionStrategy::Closest2RoundUp:
        {
          if (bcgA->uFracPart() < bcgB->uFracPart())
            return true;
          if (bcgA->uFracPart() > bcgB->uFracPart())
            return false;
          break;
        }
        case SelectionStrategy::Closest2RoundDown:
        {
          if (bcgA->lFracPart() < bcgB->lFracPart())
            return true;
          if (bcgA->lFracPart() > bcgB->lFracPart())
            return false;
          break;
        }
        case SelectionStrategy::FracWeightedPriority:
        {
          if ((bcgA->fracPart() * bcgA->priority()) > (bcgB->fracPart() * bcgB->priority()))
            return true;
          if ((bcgA->fracPart() * bcgA->priority()) < (bcgB->fracPart() * bcgB->priority()))
            return false;
          break;
        }
        case SelectionStrategy::Undefined:
        case SelectionStrategy::NotConsideredForSelection:
        case SelectionStrategy::GuidedSearch:
        case SelectionStrategy::LeastCost:
        case SelectionStrategy::LeastReducedCost:
        case SelectionStrategy::LeastGreedyCost:
        case SelectionStrategy::LeastSteepestEdgeCost:
        case SelectionStrategy::LeastPseudoCost:
        case SelectionStrategy::LeastInfeasibility:
        case SelectionStrategy::MostViolated:
        case SelectionStrategy::NotConsideredForIntegralityCheck:
        default:
          break;
        } // switch
    }

  /// Default
  if (bcgA->defaultComparisonFactor() != bcgB->defaultComparisonFactor())
    return (bcgA->defaultComparisonFactor() < bcgB->defaultComparisonFactor());

  return (bcgA->ref() < bcgB->ref());
}

/**
 * Methods of class BranchingConstrGenerator
 */
BranchingConstrGenerator::BranchingConstrGenerator(GenericBranchingConstr * gbcPtr,
                                                   const char & priorityDir,
                                                   const Double & candidateLhs,
                                                   InstanciatedConstr * prototypeInstConstrPtr,
                                                   const std::string & description) :
  _genericBrConstrPtr(gbcPtr), _ref(gbcPtr->modelPtr()->modelBrConstrGenCnt()),
  _description(description), _direction(priorityDir),
  _prototypeInstConstrPtr(prototypeInstConstrPtr),
  _candidateLhs(candidateLhs), _childNbCounter(0)
{
    gbcPtr->modelPtr()->increaseModelBrConstrGenCnt();
  if (_prototypeInstConstrPtr != NULL)
    {
      /// we do not need artificial varaibles associated to prototype, we delete them
      if (_prototypeInstConstrPtr->posLocalArtVarPtr() != NULL)
        {
          delete _prototypeInstConstrPtr->posLocalArtVarPtr();
          _prototypeInstConstrPtr->posLocalArtVarPtr(NULL);
        }
      if (_prototypeInstConstrPtr->negLocalArtVarPtr() != NULL)
        {
          delete _prototypeInstConstrPtr->negLocalArtVarPtr();
          _prototypeInstConstrPtr->negLocalArtVarPtr(NULL);
        }
      _prototypeInstConstrPtr->deleteStabInfoPtr();
      
      _prototypeInstConstrPtr->incrParticipation(6);
    }
}

BranchingConstrGenerator::BranchingConstrGenerator(const BranchingConstrGenerator & that) :
  _genericBrConstrPtr(that._genericBrConstrPtr), _ref(0), _description(that._description),
  _direction(that._direction), _prototypeInstConstrPtr(that._prototypeInstConstrPtr),
  _candidateLhs(that._candidateLhs), _childNbCounter(0)
{
  _ref = that._genericBrConstrPtr->modelPtr()->modelBrConstrGenCnt();
  that._genericBrConstrPtr->modelPtr()->increaseModelBrConstrGenCnt();
  if (_prototypeInstConstrPtr != NULL)
    _prototypeInstConstrPtr->incrParticipation(7);
}

BranchingConstrGenerator::~BranchingConstrGenerator()
{
  if (_prototypeInstConstrPtr != NULL)
    {
      _prototypeInstConstrPtr->decrParticipation(5);
      if (_prototypeInstConstrPtr->participation() == 0)
        delete _prototypeInstConstrPtr;
    }
}

BapcodInit & BranchingConstrGenerator::bapcodInit() const
{ return _genericBrConstrPtr->bapcodInit();}

std::ostream & BranchingConstrGenerator::print(std::ostream& os) const
{
  os << "BranchingConstrGenerator" << std::endl;
  os << "   direction = " << _direction << std::endl;
  os << "   candidateLhs " << _candidateLhs << std::endl;
  os << "   childNbCounter " << _childNbCounter << std::endl;
  if (_prototypeInstConstrPtr != NULL)
    os << "   constr = " << _prototypeInstConstrPtr->name() << std::endl;

  return (os);
}

void BranchingConstrGenerator::nicePrint(std::ostream& os) const
{
  os << _description << " (lhs=" << _candidateLhs << ")";
}

void BranchingConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curMastSol)
{
  _candidateLhs = 0;

  if ((_prototypeInstConstrPtr == NULL)
      || !_prototypeInstConstrPtr->isTypeOf(VcId::InstMasterConstrMask))
  {
    if (printL(5))
      std::cout << "BranchingConstrGenerator::computeLhs = 0 as (_prototypeInstConstrPtr == NULL) "
                << " or _prototypeInstConstrPtr is not an instantiated master constraint" << std::endl;
    _candidateLhs = 0;
    return;
  }

  InstMasterConstr * instMastConstrPtr = static_cast<InstMasterConstr *>(_prototypeInstConstrPtr);
  ///@todo: should setConstr, otherwise membCoef are not defined instead
  for (SolutionVarInfoPtrList::const_iterator infoIt = curMastSol.begin();
      infoIt != curMastSol.end(); infoIt++)
  {    
      bool varIsMastColumn = (*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask);

      if (varIsMastColumn)
	    {
	      MastColumn * colPtr = static_cast<MastColumn *> ((*infoIt)->varPtr);
          
          for (VarPtr2DoubleMap::const_iterator mapIt = colPtr->spSol()->solVarValMap().begin();
               mapIt != colPtr->spSol()->solVarValMap().end(); ++mapIt)
            {
              SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(mapIt->first);
              MapSubProbVariablePtr2Double::const_iterator spVarIt
                = instMastConstrPtr->subProbVarMember2coefMap().find(spVarPtr);
              if (spVarIt != instMastConstrPtr->subProbVarMember2coefMap().end())
                _candidateLhs += mapIt->second * spVarIt->second * (*infoIt)->value;
            }
        }
      else /// not master column
        {
          ConstVarConstrPtr2Double::const_iterator mapIt = instMastConstrPtr->member2coefMap().find((*infoIt)->varPtr);
          if (mapIt != instMastConstrPtr->member2coefMap().end())
            _candidateLhs += mapIt->second * (*infoIt)->value;
        }
	}

  if (printL(5))
    std::cout << "BranchingConstrGenerator::computeLhs = " << _candidateLhs << std::endl;
  return;
}

bool BranchingConstrGenerator::nextNodeBrConstr(Node * parentNodePtr,
                                                std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList,
                                                const ConstrPtrSet & existingMasterBranchingConstr)
{
  /// Node branching constraint defined last is treated first (LIFO)
  nextBranchingConstrPtrList.clear();
  bool success(true);

  int ancestorNodeRef(-1);
  if (parentNodePtr != NULL)
    ancestorNodeRef = parentNodePtr->ref();

  if (printL(5))
    std::cout << "BranchingConstrGenerator::nextNodeBrConstr ancestorNodeRef = "
              << ancestorNodeRef << std::endl;

  switch (_direction)
  {
    case 'U':
    {
      if (_childNbCounter == 0)
        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G', nextBranchingConstrPtrList);
      else if (_childNbCounter == 1)
        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L', nextBranchingConstrPtrList);
      else
        success = false;

      return (success);
    }
    default:
    {
      if (_childNbCounter == 0)
        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dfloor(_candidateLhs), 'L', nextBranchingConstrPtrList);
      else if (_childNbCounter == 1)
        instanciateBrConstr(ancestorNodeRef, ++_childNbCounter, Dceil(_candidateLhs), 'G', nextBranchingConstrPtrList);
      else
        success = false;

      return (success);
    }
  }

  /// Not path to this instruction
  return (success);
}

void BranchingConstrGenerator::instanciateBrConstr(const int & parentNodeNb, const int & childNb, const Double & rhs,
                                                   const char & sense,
                                                   std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
  if (printL(5))
    std::cout << "BranchingConstrGenerator::instanciateBrConstr() BranchingConstrBaseType  = " << std::endl;

  /**
   * Upcasting to type InstMasterBranchingConstr of an InstanciatedConstr or an InstMasterConstr
   * emanting from genericBrConstrPtr 
   */
  std::string name("BC");
  if (_prototypeInstConstrPtr != NULL)
  {

    if (printL(5))
      std::cout << "BranchingConstrGenerator::instanciateBrConstr() BranchingConstrBaseType  = "
                << _prototypeInstConstrPtr->name() << std::endl;
    name = name + _prototypeInstConstrPtr->name();
  }

  /// changed by Ruslan, we always call the full constructor of the branching constraint,
  /// calling partial construtor as below in case (imcPtr != NULL) causes problems
  /// multi-index is started by -1, -1, -1 in order to avoid having the same index as the prototype branching
  /// constraints belonging to the same generic branching constraint
  /// TO DO : we need to have two generic constraints, one for prototype br. constraints, another one for real
  ///         branching constraints
  BranchingConstrBaseType * newConstrPtr
    = new BasicConstrInstMastBranchingConstr(MultiIndex(-1, -1, -1, parentNodeNb, childNb), _genericBrConstrPtr,
                                             _genericBrConstrPtr->modelPtr()->masterConfPtr(),
                                             _prototypeInstConstrPtr, _description,
                                             name + "p" + parentNodeNb + "c" + childNb,
                                             rhs, sense, ' ', 'E', 'd');

  if (printL(5))
    newConstrPtr->print(cout);

  nextBranchingConstrPtrList.push_back(newConstrPtr);

  return;

}

GenAggrSubProbVarBranchingConstr
::GenAggrSubProbVarBranchingConstr(Model * modelPtr,
			                       ProbConfig * probConfPtr,
			                       const std::string & name,
                                   const std::string & genVarName,
                                   const double & targetFraction,
                                   const int & numIgnoredIndices,
                                   const SelectionStrategy & priorityRule,
                                   const Double & nonRootPriorityLevel,
                                   const Double & rootPriorityLevel,
			                       const bool & toBeUsedInPreprocessing) :
  GenericBranchingConstr(modelPtr, probConfPtr, name, priorityRule,
                         nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing),
  _genVarName(genVarName), _targetFraction(targetFraction), _numIgnoredIndices(numIgnoredIndices),
  _numGeneratedBrConstrs(0), _branchingVarPriorityFunctorPtr(NULL)
{
}

GenAggrSubProbVarBranchingConstr::~GenAggrSubProbVarBranchingConstr()
{
}

void GenAggrSubProbVarBranchingConstr
     ::branchingVarPriorityFunctorPtr(BcVarBranchingPriorityFunctor * branchVarPriorityFuncPtr)
{
    _branchingVarPriorityFunctorPtr = branchVarPriorityFuncPtr;
}

void GenAggrSubProbVarBranchingConstr
     ::branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                         const MasterColSolution & curListOfMasterCol,
                                         const int & maxNumOfCandidates,
		  			                     BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  if (probConfPtr() == NULL)
    return;

  std::map<MultiIndex, double> index2valueMap;
  std::pair<std::map<MultiIndex, double>::iterator, bool> mapInsertResult;

  for (MasterVarSolution::const_iterator solIt = curListOfMastAndSubprobVar.begin();
       solIt != curListOfMastAndSubprobVar.end(); ++solIt)
    if (solIt->first->isTypeOf(VcId::SubProbVariableMask)
        && (solIt->first->genVarPtr()->defaultName() == _genVarName))
      {
        MultiIndex varId(static_cast<const InstanciatedVar *>(solIt->first)->id(), _numIgnoredIndices);
        double varValue = solIt->second._value;
        mapInsertResult = index2valueMap.insert(std::make_pair(varId, varValue));
        if (!mapInsertResult.second)
          (mapInsertResult.first)->second += varValue;
      }

  std::set<std::pair<double, MultiIndex > > sortedIdsByPriority;

  std::map<MultiIndex, double>::iterator ivMapIt;
  for (ivMapIt = index2valueMap.begin(); ivMapIt != index2valueMap.end(); ++ivMapIt)
    {
      double intPart;
      double fracPart = std::modf(ivMapIt->second, &intPart);
      if ((fracPart < param().BapCodIntegralityTolerance()) || (fracPart > 1 - param().BapCodIntegralityTolerance()))
        continue;
      
      double priority = std::abs(fracPart - _targetFraction);
      if (_branchingVarPriorityFunctorPtr != NULL)
        priority = - (*_branchingVarPriorityFunctorPtr)(ivMapIt->first, intPart, fracPart);
      sortedIdsByPriority.insert(std::make_pair(priority, ivMapIt->first));
    }
    
  std::set<std::pair<double, MultiIndex > >::iterator setIt;
  int candCounter = 0;
  for (setIt = sortedIdsByPriority.begin(); setIt != sortedIdsByPriority.end(); ++setIt)
    {
      if (++candCounter > maxNumOfCandidates)
        break;
      
      MultiIndex constrId(_numGeneratedBrConstrs++);
      InstMasterConstr * prototypeConstrPtr = new InstMasterConstr(constrId, this, probConfPtr(), defaultName(),
                                                                   0, defaultSense());
      
      for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
           cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
        {
          GenericVar * genVarPtr = (*cgSpConfPtrIt)->getGenericVar(_genVarName);
          if (genVarPtr == NULL)
            continue;
          
          if (_numIgnoredIndices == 0)
            {
              InstanciatedVar * iVarPtr = genVarPtr->getVarPtr(setIt->second);
              if (iVarPtr == NULL)
                continue;

              prototypeConstrPtr->includeMember(iVarPtr, 1.0, true);
            }
          else
            {
              for (IndexCell2InstancVarPtrMap::const_iterator mapIt = genVarPtr->indexCell2InstancVarPtrMap().begin();
                   mapIt != genVarPtr->indexCell2InstancVarPtrMap().end(); ++mapIt)
                {
                  MultiIndex varId(mapIt->first, _numIgnoredIndices);
                  if (varId == setIt->second)
                    prototypeConstrPtr->includeMember(mapIt->second, 1.0, true);
                }
            }
        }
      
      std::string description(defaultName());
      BranchingConstrGenerator * gentorPtr = new BranchingConstrGenerator(this, 'U', 0.0, prototypeConstrPtr,
                                                                          setIt->second.appendRef2name(description));
      generatedBrConstrGeneratorSet.insert(gentorPtr);
    }
    
  return;
}

/**
 * Methods of class GenericBranchingConstr
 */
GenericBranchingConstr::GenericBranchingConstr(Model * modelPtr,
                                               ProbConfig * probConfPtr,
                                               const std::string & name,
                                               const SelectionStrategy & prioritySelectionRule,
                                               const Double & nonRootPriorityLevel,
                                               const Double & rootPriorityLevel,
                                               const bool & toBeUsedInPreprocessing) :
  DynamicGenericConstr(modelPtr, probConfPtr, name, 'C', prioritySelectionRule,
                       nonRootPriorityLevel, rootPriorityLevel, toBeUsedInPreprocessing),
  _branchingHistory(), _branchingSeparationFunctorPtr(NULL)
{
  if (probConfPtr != NULL)
    probConfPtr->insertGenericBranchingConstr(this);

  return;
}

GenericBranchingConstr::~GenericBranchingConstr()
{
  std::vector<BranchingEvaluationInfo *>::iterator evalInfoPtrIt;
  for (BranchingGeneratorHistoryMap::iterator mapIt = _branchingHistory.begin();
       mapIt != _branchingHistory.end(); ++mapIt)
    {
      for (evalInfoPtrIt = mapIt->second.evaluationsInfo.begin();
           evalInfoPtrIt != mapIt->second.evaluationsInfo.end(); ++evalInfoPtrIt)
        delete *evalInfoPtrIt;
      delete mapIt->first;
    }
  _branchingHistory.clear();

  return;
}

void GenericBranchingConstr::branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                               const int & maxNumOfCandidates,
                                                               BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  return;
}

double GenericBranchingConstr::computeLhs(const MasterVarSolution & curListOfMastVarInSol,
                                          const InstanciatedConstr * instConstrPtr)
{
  double lhs = 0.0;

  if ((instConstrPtr == NULL) || !instConstrPtr->isTypeOf(VcId::InstMasterConstrMask))
    return lhs;

  const InstMasterConstr * instMastConstrPtr = static_cast<const InstMasterConstr *>(instConstrPtr);
  for (MasterVarSolution::const_iterator pairIt = curListOfMastVarInSol.begin();
       pairIt != curListOfMastVarInSol.end(); pairIt++)
  {
    bool isSubprobVar = pairIt->first->isTypeOf(VcId::SubProbVariableMask);
    if (isSubprobVar)
    {
      SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(pairIt->first);
      MapSubProbVariablePtr2Double::const_iterator mapIt;
      mapIt = instMastConstrPtr->subProbVarMember2coefMap().find(spVarPtr);
      if (mapIt != instMastConstrPtr->subProbVarMember2coefMap().end())
        lhs += mapIt->second * pairIt->second._value;
    }
    else /// pure master variable
    {
      ConstVarConstrPtr2Double::const_iterator mapIt;
      mapIt = instMastConstrPtr->member2coefMap().find(pairIt->first);
      if (mapIt != instMastConstrPtr->member2coefMap().end())
        lhs += mapIt->second * pairIt->second._value;
    }
  }

  return lhs;
}

void GenericBranchingConstr::branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                               const MasterColSolution & curListOfMasterCol,
                                                               const int & maxNumOfCandidates,
                                                               BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  if (probConfPtr() == NULL)
    return;

  Solution *solPtr = probConfPtr()->getSolution(curListOfMastAndSubprobVar);
  Solution *primalSolPtr = probConfPtr()->getAggregatedSolution(solPtr);

  if (printL(5))
    std::cout << "GenericBranchingConstr::branchingSeparationFindCandidates: primalSol " << primalSolPtr << std::endl;

  double curFractionalLhs(0);
  std::list<std::pair<BcConstr, std::string> > returnBrConstrList;

  bool branchConstrHasBeenGenerated = false;
  if (_branchingSeparationFunctorPtr == NULL) /// prototype constraints are pre-specified by the user
  {
    std::vector<std::pair<double, std::pair<IndexCell, InstanciatedConstr *> > > prototypeConstrPtsVector;
    for (std::list<InstanciatedConstr *>::iterator listIt = _constrPrototypes.begin();
         listIt != _constrPrototypes.end(); ++listIt)
    {
      /// we calculate the left-hand side of each candidate
      double intPart, lhs = computeLhs(curListOfMastAndSubprobVar, *listIt);
      double lhsFracPart = std::modf(lhs, &intPart);
      if ((lhsFracPart > param().BapCodIntegralityTolerance())
          && (lhsFracPart < 1 - param().BapCodIntegralityTolerance()))
      {
        prototypeConstrPtsVector.push_back(std::make_pair(std::abs(lhsFracPart - 0.5),
                                                 std::make_pair((*listIt)->id(), *listIt)));
      }
    }
    std::stable_sort(prototypeConstrPtsVector.begin(), prototypeConstrPtsVector.end());
    int nbCandidates = 0;
    std::vector<std::pair<double, std::pair<IndexCell, InstanciatedConstr *> > >::iterator pairIt;
    for (pairIt = prototypeConstrPtsVector.begin();
         (pairIt != prototypeConstrPtsVector.end()) && (nbCandidates < maxNumOfCandidates); ++pairIt, ++nbCandidates)
    {
       std::string description(defaultName());
       if (_constrPrototypes.size() > 1)
         description = pairIt->second.first.appendRefWithBrackets2name(description);
       returnBrConstrList.push_back(std::make_pair(BcConstr(pairIt->second.second), description));
    }
    branchConstrHasBeenGenerated = !returnBrConstrList.empty();
  }
  else
  {
      BcSolution projectedSol(primalSolPtr);
      std::list<std::pair<double, BcSolution> > columnsInSol;

      for (auto & pair : curListOfMasterCol)
          columnsInSol.push_back(std::make_pair(pair.second._value, BcSolution(pair.first->spSol())));

      branchConstrHasBeenGenerated = (*_branchingSeparationFunctorPtr)(BcFormulation(probConfPtr()), projectedSol,
                                                                       columnsInSol, maxNumOfCandidates,
                                                                       returnBrConstrList);
  }

  delete solPtr;
  primalSolPtr->deleteSolutionsChain();
  delete primalSolPtr;

  if (printL(5))
    std::cout << "GenericBranchingConstr::branchingSeparationFindCandidates: branchConstrhasBeenGenerated "
              << branchConstrHasBeenGenerated << std::endl;

  if (branchConstrHasBeenGenerated)
  {
    for (std::list<std::pair<BcConstr, std::string> >::iterator bcit = returnBrConstrList.begin();
        bcit != returnBrConstrList.end(); ++bcit)
      {
        BranchingConstrGenerator * gentorPtr = new BranchingConstrGenerator(this, 'U', curFractionalLhs,
                                                                            (InstanciatedConstr *) (bcit->first),
                                                                            bcit->second);
        generatedBrConstrGeneratorSet.insert(gentorPtr);
      }

    if ((int) generatedBrConstrGeneratorSet.size() > maxNumOfCandidates)
    {
      if (printL(5))
      {
        std::cout
            << "branchingSeparationFindCandidates : remove last  candidate of list of size = "
            << generatedBrConstrGeneratorSet.size() << std::endl;
      }
      delete *(--(generatedBrConstrGeneratorSet.end()));
      generatedBrConstrGeneratorSet.erase(--(generatedBrConstrGeneratorSet.end()));
    }

  }

  return;
}


bool GenericBranchingConstr::prepareSeparation()
{
  /// we copy all constraints specified by the user to the set of prototype constrains
  for (IndexCell2InstancConstrPtrMap::iterator mapIt = _indexCell2InstancConstrPtrMap.begin();
       mapIt != _indexCell2InstancConstrPtrMap.end(); ++mapIt)
  {
    insertConstrPrototypes(mapIt->second);
    /// this is needed in order to keep every prototype constraint "alive" during the whole execution,
    /// as otherwise it may be deleted when a branching generator is deleted;
    mapIt->second->incrParticipation(16);
  }
  /// we now delete all prototype constraints from the generic constr maps
  _indexCell2InstancConstrPtrMap.clear();
  _multiIndex2ConstrPtrMap.clear();

  return true;
}

