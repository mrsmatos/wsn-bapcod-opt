/**
 *
 * This file bcMasterConfC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMasterConfC.hpp"
#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcModelC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcGenMastVarConstrC.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4GenChildrenInBranching.hpp"

using namespace std;
using namespace VcIndexStatus;

MasterConf::MasterConf(Model * modelPtr,
		               Problem * problemPtr,
                       const Bound & initialPrimalBound,
		               const Bound & initialDualBound) :
                    ProbConfig(ProbConfig::master, modelPtr,  "master",
                               IndexCell(), initialPrimalBound, initialDualBound, problemPtr),
                    _enumSolutionPtr(NULL),
                    _debugSolutionPtr(NULL),
                    _masterCommons4ProblemSetup(*this),
                    _masterCommons4EvalAlg(*this),
                    _masterCommons4GenChildNodes(*this),
                    _masterCommons4PrimalHeuristic(*this),
                    _nonStabStaticArtVarPtrList(),
                    _initialSetOfActiveColumns(NULL),
                    _initialSetOfInactiveColumns(NULL),
                    _bestPrimalBound(initialPrimalBound),
                    _bestDualBound(initialDualBound),
                    _bapTreeNodes(),
                    _subTreeSizeByDepth()
{
}

MasterConf::~MasterConf()
{
    if (_enumSolutionPtr != NULL)
    {
        _enumSolutionPtr->deleteSolutionsChain();
        delete _enumSolutionPtr;
    }

    if (_debugSolutionPtr != NULL)
    {
        _debugSolutionPtr->deleteSolutionsChain();
        delete _debugSolutionPtr;
    }

    for (size_t i = 0; i < _bapTreeNodes.size(); ++i)
        delete _bapTreeNodes[i];
    _bapTreeNodes.clear();

    /// we delete primal solution
    if (_primalSolPtr != NULL)
    {
        Solution * currentPtr = _primalSolPtr;

        while(currentPtr != NULL)
        {
            Solution * tempPtr = currentPtr;
            currentPtr = currentPtr->nextSolPtr();
            delete tempPtr;
        }

        _primalSolPtr = NULL;
    }

  /// we first delete from memory all dynamic variables and constraints
  _probPtr->deactivateAndRemoveAllVarsAndConstrsFromMemory();

  while (!_nonStabStaticArtVarPtrList.empty())
    {
      delete _nonStabStaticArtVarPtrList.back();
      _nonStabStaticArtVarPtrList.pop_back();
    }

  for (std::list< Constraint *>::iterator it = _pcConstrPtrList.begin(); it != _pcConstrPtrList.end(); )
    {
      Constraint* constrPtr = *it;
      if (constrPtr->isTypeOf(VcId::InstMastConvexityConstrMask))
        {
          _modelPtr->garbageCollector().erase(constrPtr);
          it = _pcConstrPtrList.erase(it);
          delete constrPtr;
        }
      else
      {
        it++;
      }
    }

  while (!_convexityGenConstrPtrList.empty())
    {
      ConvexityGenConstr * convexityGenConstr = _convexityGenConstrPtrList.back();

      _convexityGenConstrPtrList.pop_back();
      delete  convexityGenConstr;
    }


  while (!_colGenSubProbConfPts.empty())
    {
      ColGenSpConf * colGenSp = _colGenSubProbConfPts.back();
      _colGenSubProbConfPts.pop_back();
      delete colGenSp;
    }
}

std::ostream& MasterConf::print(std::ostream& os) const
{
  os << "MasterConf: " << std::endl;
  ProbConfig::print(os);

  if (printL(6))
    {
      os << "  pure Master Variables: " << std::endl;

      for (InstMasterVarPtrSet::const_iterator it = _setOfPureMasterVar.begin();
           it != _setOfPureMasterVar.end(); ++it)
        os << "  var name = " << (*it)->name() << std::endl;

      os << "  col gen subProblems: " << std::endl;
      for (std::vector<ColGenSpConf *>::const_iterator cgIt = colGenSubProbConfPts().begin();
           cgIt != colGenSubProbConfPts().end(); ++cgIt)
        (*cgIt)->print(os);
    }

  if (probPtr() != NULL)
    probPtr()->print(os);

  return (os);
}

Solution * MasterConf::getAggregatedSolution(Solution * masterSolPtr)
{
  if (masterSolPtr == NULL)
    return NULL;

  const int printlevel = 6;

  Solution * pureMasterSolPtr = new Solution(this, NULL);

  pureMasterSolPtr->cost(masterSolPtr->cost());

  Solution * previousSolptr = pureMasterSolPtr;

  MasterColSolution listOfFractMastCol;

  const VarPtr2DoubleMap & mastColMap = masterSolPtr->solVarValMap();
  VarPtr2DoubleMap mastSolMap;

  for (VarPtr2DoubleMap::const_iterator it = mastColMap.begin(); it != mastColMap.end(); ++it)
    {
      if (it->first->isTypeOf(VcId::MastColumnMask))
        {
          ValueRecord rec(it->second);
          listOfFractMastCol.push_back(it->first, rec);
        }
      else
        {
          pureMasterSolPtr->includeVar(it->first, it->second, false);
          pureMasterSolPtr->cost( pureMasterSolPtr->cost() + it->second * it->first->curCost());
        }
    }

  resetColGenSpListOfFractMastCol(listOfFractMastCol);

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = colGenSubProbConfPts().begin();
       spcPt != colGenSubProbConfPts().end(); ++spcPt)
    {
      VarPtr2DoubleMap curAggregateSol;

      /**
       * _listOfFractMastCol  includes integer variables up to value 1.999; columns are sorted in ILO
       */
      for (MasterColSolution::const_iterator colPt = (*spcPt)->listOfFractMastColInColGenSp().begin();
           colPt != (*spcPt)->listOfFractMastColInColGenSp().end(); colPt++)
        {
          colPt->first->fillMapOfIntSpVar(curAggregateSol, colPt->second._value);
        } /// ColPt

      Solution * projectedSPSolPtr = new Solution(*spcPt, previousSolptr);
      previousSolptr = projectedSPSolPtr;
      Double cumCost(0);

      for (VarPtr2DoubleMap::const_iterator spit = curAggregateSol.begin();
           spit != curAggregateSol.end(); spit++)
        {
          if (printL(printlevel))
            std::cout << "curAggregateSol[" << spit->first->name() << "] = " << spit->second << std::endl;

          projectedSPSolPtr->includeVar(spit->first, spit->second,false);
        }
    }

  return pureMasterSolPtr;
}


void MasterConf::resetColGenSpListOfFractMastCol(const MasterColSolution & curListOfFractMastCol)
{
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = colGenSubProbConfPts().begin();
       spcPt != colGenSubProbConfPts().end(); ++spcPt)
    {
      (*spcPt)->listOfFractMastColInColGenSp().clear();
    }

  for (MasterColSolution::const_iterator colPt = curListOfFractMastCol.begin();
       colPt != curListOfFractMastCol.end(); colPt++)
    {
	  colPt->first->cgSpConfPtr()->listOfFractMastColInColGenSp().push_back(colPt->first, colPt->second);
    }
}

Solution * MasterConf::getDissagregatedSolution(Solution * masterSolPtr)
{
  if (masterSolPtr == NULL)
    return NULL;

  const int printlevel = 5;

  Solution * pureMasterSolPtr = new Solution(this, NULL);
  pureMasterSolPtr->cost(masterSolPtr->cost());

  Solution * previousSolptr = pureMasterSolPtr;
  int floorPart;

  /// records fractional columns to aggregate them into sp solutions
  MasterColSolution listOfFractMastCol;

  const VarPtr2DoubleMap & mastColMap = masterSolPtr->solVarValMap();

  for (VarPtr2DoubleMap::const_iterator it = mastColMap.begin(); it != mastColMap.end(); ++it)
    {
      if (it->first->isTypeOf(VcId::MastColumnMask))
        {
          ValueRecord rec(it->second, param().BapCodIntegralityTolerance());
          if (rec._isFractional)
            {
              listOfFractMastCol.push_back(it->first, rec);
            }
          floorPart = rec._value - rec._lfracValue;
          if (floorPart > param().BapCodIntegralityTolerance())
            {
              Solution * projectedSPSolPtr = new Solution(it->first->cgSpConfPtr(), previousSolptr);
              previousSolptr = projectedSPSolPtr;

              projectedSPSolPtr->cost(it->first->curCost());
              for (VarPtr2DoubleMap::const_iterator spit = it->first->spSol()->solVarValMap().begin();
                   spit != it->first->spSol()->solVarValMap().end(); spit++)
                {
                  if (printL(printlevel))
                    std::cout << "curDisaggregateSol[" << spit->first->name()
                              << "] = " << spit->second << std::endl;

                  projectedSPSolPtr->includeVar(spit->first, spit->second, false);
                }
              /// save ordered solution if exists
#ifdef BCP_RCSP_IS_FOUND
              projectedSPSolPtr->setRcspSolPtr(it->first->spSol()->copyRcspSolPtr());
#endif
              const std::vector<int> & orderedIds = it->first->spSol()->orderedIds();
              if (!orderedIds.empty())
                {
                  for (std::vector<int>::const_iterator arcIdIt = orderedIds.begin();
                       arcIdIt != orderedIds.end(); ++arcIdIt)
                    projectedSPSolPtr->addToOrderedIds(*arcIdIt);
                  const std::vector<std::vector<double> > & resConsumption = it->first->spSol()->resConsumption();
                  for (std::vector<std::vector<double> >::const_iterator resConsIt = resConsumption.begin();
                       resConsIt != resConsumption.end(); ++resConsIt)
                    projectedSPSolPtr->addToResConsumption(*resConsIt);
                }

              projectedSPSolPtr->multiplicity(floorPart);
            }
        }
      else /// pure master var
        {
          if (printL(printlevel))
            std::cout << "getDissagregatedSolution(): pure mast var " << it->first->name()
                      << " val = " << it->second << std::endl;
          pureMasterSolPtr->includeVar(it->first, it->second, false);
        }
    }

  resetColGenSpListOfFractMastCol(listOfFractMastCol);

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = colGenSubProbConfPts().begin();
       spcPt != colGenSubProbConfPts().end(); ++spcPt)
    {
      std::vector< std::pair < MastColumn *, ValueRecord > > listOfColInLexicographicOrder;

      for (MasterColSolution::const_iterator it = (*spcPt)->listOfFractMastColInColGenSp().begin();
           it != (*spcPt)->listOfFractMastColInColGenSp().end(); ++it)
          listOfColInLexicographicOrder.push_back(*it);

      stable_sort(listOfColInLexicographicOrder.begin(), listOfColInLexicographicOrder.end(),
                  LexicographicMastColValSorting());

      std::map<int, VarPtr2DoubleMap> curDisaggregateSol;
      if (printL(printlevel))
        std::cout << "getDissagregatedSolution(): considering subproblem " << (*spcPt)->name() << std::endl;

      Double cumVal(0);
      Double cumCost(0);
      int tindex(0);
      MastColumn * mastColPtr(NULL);
      const Solution * latestSolPtr(NULL);

      std::vector< std::pair < MastColumn *, ValueRecord > >::const_iterator colPt;
      for (colPt = listOfColInLexicographicOrder.begin(); colPt != listOfColInLexicographicOrder.end(); colPt++)
      /**
       * Try to take convex combination of frac col to get an integer solution
       */
        {
          mastColPtr = colPt->first;
          if (printL(printlevel))
            {
              std::cout << "col[ " << mastColPtr->name() << " ] has fractional part "
                        << colPt->second._lfracValue << " cumVal = " << cumVal << std::endl;

              mastColPtr->printColVector();
              mastColPtr->spSol()->printOrderedSolution();
            }

          Double tmpColVal(colPt->second._lfracValue);
          Double lambdaR(0);
          while (tmpColVal > 0)
            {
              /// Slice is filled
              if (cumVal == (tindex + 1))
                {
                  Solution * projectedSPSolPtr = new Solution(*spcPt, previousSolptr);
                  previousSolptr = projectedSPSolPtr;
                  projectedSPSolPtr->cost(cumCost);
                  cumCost = 0;
                  for (VarPtr2DoubleMap::const_iterator spit = curDisaggregateSol[tindex].begin();
                       spit != curDisaggregateSol[tindex].end(); spit++)
                    {
                      if (printL(printlevel))
                        std::cout << "curDisaggregateSol[" << tindex << "][" << spit->first->name()
                                  << "] = " << spit->second << std::endl;

                      projectedSPSolPtr->includeVar(spit->first, spit->second, false);
                    }

                  /// save ordered solution if exists
#ifdef BCP_RCSP_IS_FOUND
                 projectedSPSolPtr->setRcspSolPtr(latestSolPtr->copyRcspSolPtr());
#endif
                 const std::vector<int> & orderedIds = latestSolPtr->orderedIds();
                 if (!orderedIds.empty())
                   {
                     for (std::vector<int>::const_iterator arcIdIt = orderedIds.begin();
                          arcIdIt != orderedIds.end(); ++arcIdIt)
                      projectedSPSolPtr->addToOrderedIds(*arcIdIt);
                      const std::vector<std::vector<double> > & resConsumption = mastColPtr->spSol()->resConsumption();
                      for (std::vector<std::vector<double> >::const_iterator resConsIt = resConsumption.begin();
                        resConsIt != resConsumption.end(); ++resConsIt)
                      projectedSPSolPtr->addToResConsumption(*resConsIt);
                   }

                  /// Goto next slice
                  tindex++;
                }


              if (cumVal + lambdaR > tindex + 1)
                lambdaR = (tindex + 1 - cumVal);
              else
                lambdaR = tmpColVal;

              mastColPtr->fillAggregateSol(curDisaggregateSol[tindex], lambdaR);
              latestSolPtr = mastColPtr->spSol();
              tmpColVal -= lambdaR;
              cumVal += lambdaR;
              cumCost += lambdaR * mastColPtr->curCost();

              if (printL(printlevel))
                std::cout << "getDissagregatedSolution(): tindex = " << tindex
                          << " cumVal = " << cumVal << std::endl;
            }
        } /// ColPt

      if (cumVal > param().BapCodIntegralityTolerance()) /// check last slice (the others have been check on creation)
        {
          Solution * projectedSPSolPtr = new Solution(*spcPt, previousSolptr);
          previousSolptr = projectedSPSolPtr;
          projectedSPSolPtr->cost(cumCost);
          cumCost = 0;

          for (VarPtr2DoubleMap::const_iterator spit = curDisaggregateSol[tindex].begin();
               spit != curDisaggregateSol[tindex].end(); spit++)
            {
              if (printL(printlevel))
                std::cout << "curDisaggregateSol[" << tindex << "][" << spit->first->name()
                          << "] = " << spit->second << std::endl;

              projectedSPSolPtr->includeVar(spit->first, spit->second, false);
            }

          /// save ordered solution if exists
#ifdef BCP_RCSP_IS_FOUND
          projectedSPSolPtr->setRcspSolPtr(latestSolPtr->copyRcspSolPtr());
#endif
          const std::vector<int> & orderedIds = latestSolPtr->orderedIds();
          if (!orderedIds.empty())
            {
              for (std::vector<int>::const_iterator arcIdIt = orderedIds.begin();
                   arcIdIt != orderedIds.end(); ++arcIdIt)
              projectedSPSolPtr->addToOrderedIds(*arcIdIt);
              const std::vector<std::vector<double> > & resConsumption = mastColPtr->spSol()->resConsumption();
              for (std::vector<std::vector<double> >::const_iterator resConsIt = resConsumption.begin();
                resConsIt != resConsumption.end(); ++resConsIt)
              projectedSPSolPtr->addToResConsumption(*resConsIt);
            }
          
        }
    }

  return pureMasterSolPtr;
}

void MasterConf::recordInitialActiveSetOfColumns(Solution * solPtr)
{
  if (solPtr != NULL)
    _initialSetOfActiveColumns = solPtr;
}

void MasterConf::recordInitialInactiveSetOfColumns(Solution * solPtr)
{
  if (solPtr != NULL)
    _initialSetOfInactiveColumns = solPtr;
}

/// function replaced by Ruslan
void MasterConf::updatePrimalIncSolution(const Bound & primalBound, Solution * solPtr)
{
  if (solPtr == nullptr)
    return;

  solPtr->resetCost();
  if (updatePrimalIncBound(primalBound))
    {
      if (_primalSolPtr != nullptr)
        delete _primalSolPtr;
      if (!_modelPtr->checkIfSolutionIsFeasibleUsingCallback(solPtr))
        {
          std::cerr << "Error: new incumbent solution is infeasible as determined by "
                       "the user specified callback" << std::endl;
          exit(1);
        }
      _primalSolPtr = solPtr->clone();
    }
}

void MasterConf::updateCurValidDualBound(Node * nodePtr)
{
  ///  New node dual bound
  if (nodePtr->updateSubtreeDualBound(nodePtr->nodeIncIpDualBound()))
    {
      ///  Root node
      if (nodePtr->depth() == 0)
        {
          updateDualIncBound(nodePtr->subtreeDualBound());
          return;
        }
    }

  bool dualBoundUpdated(false);
  Bound worstLowerBound(0, modelPtr()->objectiveSense());

  for (Node * bbNodePtr = nodePtr; bbNodePtr->father() != NULL; bbNodePtr = bbNodePtr->father())
    {
      /// temporarily for testing
      if (bbNodePtr->subtreeDualBound() < bbNodePtr->father()->subtreeDualBound())
      {
        if (printL(-1))
          std::cout << "BaPCod WARNING : dual bound " << bbNodePtr->subtreeDualBound()
                    << " of node N_" << bbNodePtr->branchAndPriceOrder()
                    << " is worse than dual bound " << bbNodePtr->father()->subtreeDualBound()
                    << " of its father N_" << bbNodePtr->father()->branchAndPriceOrder() << std::endl;
        std::cerr << "BaPCod WARNING : dual bound of a node is worse that dual bound of its father" << std::endl;
      }

      /// Compute min lower bound amongst child
      worstLowerBound = bbNodePtr->subtreeDualBound();

      for (std::list<Node *>::const_iterator sNpt = bbNodePtr->father()->sons().begin();
           sNpt != bbNodePtr->father()->sons().end(); sNpt++)
        {
          if ((*sNpt)->subtreeDualBound() < worstLowerBound)
              worstLowerBound = (*sNpt)->subtreeDualBound();
          if (printL(5))
            std::cout << "MasterConf::updateCurValidDualBound() brother node  has dual bound "
                      << (*sNpt)->subtreeDualBound() << " while best dual bound is "
                      << worstLowerBound << std::endl;
        }

      /// temporarily for testing
      if (worstLowerBound < bbNodePtr->father()->subtreeDualBound())
      {
        if (printL(-1))
          std::cout << "BaPCod WARNING : min dual bound " << worstLowerBound
                    << " of children of N_" << bbNodePtr->father()->branchAndPriceOrder()
                    << " is worse than its dual bound " << bbNodePtr->father()->subtreeDualBound() << std::endl;
        std::cerr << "BaPCod WARNING : minimum dual bound of children is worse that dual of the node" << std::endl;
      }

      if (printL(2))
        {
          if (worstLowerBound > bbNodePtr->father()->subtreeDualBound())
            std::cout << "MasterConf::updateCurValidDualBound(): Subtree Dual Bound Update at Node "
                      << bbNodePtr->father()->ref() << std::endl
                      << "   It had a dual bound  = "
                      << bbNodePtr->father()->subtreeDualBound() << std::endl
                      << "   It is updated to new bound = " << worstLowerBound
                      << std::endl;
        }

      dualBoundUpdated = bbNodePtr->father()->updateSubtreeDualBound(worstLowerBound);

      ///  Root node not reached
      if (bbNodePtr->father()->depth() != 0)
        dualBoundUpdated = false;
    }

  if (dualBoundUpdated)
    updateDualIncBound(worstLowerBound);

  return;
}

bool MasterConf::updateDualIncBound(const Bound & newDualBound)
{
  if (ProbConfig::updateDualIncBound(newDualBound))
    {
      if (dualIncBound() > bestDualBound())
        bestDualBound(dualIncBound());

      bapcodInit().statistics().incrValue("bcRecBestDb", dualIncBound());
      if (printL(1))
        {
          std::cout << "******************************* " << std::endl;

          std::cout << "New Dual Bound for the problem: " << dualIncBound()
              << " while Primal Bd is " << primalIncBound() << std::endl;

          std::cout << "******************************* " << std::endl;

          printTime(bapcodInit().startTime().getElapsedTime(), std::cout);
        }
      return true;
    }
  return false;
}

std::ostream& MasterConf::printInfoBeforeSolvingNode(const Node * nodePtr, int numPrimOpenNodes, int numSecOpenNodes,
                                                     int estimTreeSize, std::ostream& os) const
{
  os << "************************************************************************************************" << std::endl;

  os << "**** " ;
  nodePtr->nicePrint(os);
  os << ", global bounds : [ " << dualIncBound() << " , " << primalIncBound() << " ], ";
  printTime(bapcodInit().startTime().getElapsedTime(), os);

  os << "**** " << numPrimOpenNodes;
  if (numSecOpenNodes > 0)
    os << " (+" << numSecOpenNodes << ")";
  os << " open nodes, ";
  if (estimTreeSize > 0)
    os << "ETS : " << estimTreeSize << ", ";
  probPtr()->printDynamicVarConstrStats(os); //, true);
  os << std::endl;

  os << "************************************************************************************************" << std::endl;

  return (os);
}

bool MasterConf::updatePrimalIncBound(const Bound & newIncVal)
{
  bool status(ProbConfig::updatePrimalIncBound(newIncVal));

  if (status)
    {
      if (printL(-1))
      {
        std::cout << "New model incumbent solution " << newIncVal << ", ";
        printTime(bapcodInit().startTime().getElapsedTime(), std::cout);
      }
      bapcodInit().statistics().incrValue("bcRecBestInc", _primalIncBound);
    }
  return status;
}
