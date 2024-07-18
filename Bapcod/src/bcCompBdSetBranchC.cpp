/**
 *
 * This file bcCompBdSetBranchC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
//#include "bcBendersCutC.hpp"

#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcGenBranchingConstrC.hpp"

using namespace std;

/// Methods of class CompBdSetBranchConstrGenerator

CompBdSetBranchConstrGenerator
::CompBdSetBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
                                 const ComponentSequence & compBoundSet,
                                 const char & priorityDir) :
  BranchingConstrGenerator(gbcPtr, priorityDir), _compBoundSet(compBoundSet),
  _currentCompBoundSet(compBoundSet)
{
  /// we generate the description string now
  std::stringstream sstream;
  sstream << (_compBoundSet._cgSpConfPtr != NULL ? _compBoundSet._cgSpConfPtr->name() : "undefined subProb.");
  if (!_compBoundSet.empty())
    {
      sstream << " with ";
      for (ComponentSequence::const_iterator it = _compBoundSet.begin(); it != _compBoundSet.end(); ++it)
        {
          if (it != _compBoundSet.begin())
            sstream << ", ";
          sstream << it->_varPtr->name();
          switch (it->_sign) {
            case 'G':
              sstream << " >= ";
              break;
            case 'L':
              sstream << " <= ";
              break;
            }
          sstream << (*it)._val;
        }
    }
  _description = sstream.str();
}

CompBdSetBranchConstrGenerator
::CompBdSetBranchConstrGenerator(const CompBdSetBranchConstrGenerator & that):
  BranchingConstrGenerator(that), _compBoundSet(that._compBoundSet),
  _currentCompBoundSet(that._compBoundSet)
{
}

CompBdSetBranchConstrGenerator::~CompBdSetBranchConstrGenerator()
{
}

const Double CompBdSetBranchConstrGenerator::priority() const
{
  if (_compBoundSet.empty()) /// branching on thenumber of comumns associated to the subproblem
  {
    return 100; /// arbritrary high priority
    /// @todo define the max priority for the subproblme and use it here.
  }
  return (_compBoundSet.lastCP().varPtr()->priority());
}

std::ostream & CompBdSetBranchConstrGenerator::print(std::ostream& os) const
{
  BranchingConstrGenerator::print(os);
  os << "CompBdSetBranchConstrGenerator" << std::endl;
  for (ComponentSequence::const_iterator it = _compBoundSet.begin();
       it != _compBoundSet.end(); ++it)
    os << *it;

  return (os);
}

void CompBdSetBranchConstrGenerator::nicePrint(std::ostream& os) const
{
  os << _description << " (" << _compBoundSet.fracWeight() << ")";
}

void CompBdSetBranchConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curSol)
{
  for (ComponentSequence::iterator it = _compBoundSet.begin(); it != _compBoundSet.end(); ++it)
  {
    it->cardinality(0);
    it->complementCard(0);
  }
  Double curFracWeight(0);

  if (printL(6))
    std::cout << "CompBdSetBranchConstrGenerator::computeLhs" << std::endl;

  for (SolutionVarInfoPtrList::const_iterator infoIt = curSol.begin(); infoIt != curSol.end(); infoIt++)
  {
    if (!(*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
      continue;

    MastColumn * mastColumnPtr = static_cast<MastColumn *> ((*infoIt)->varPtr);

    if (!(_compBoundSet.cgSpConfPtr() == mastColumnPtr->cgSpConfPtr()))
      /// Column's SP is not concerned by componentSet on which we branch
      continue;

    curFracWeight += (*infoIt)->value;

    for (ComponentSequence::iterator it = _compBoundSet.begin(); it != _compBoundSet.end(); ++it)
    {
      /*
       * Check whether column belongs to nested subclass
       * defined by all previous component bound
       * else goto next column
       */
      if (printL(6))
      {
        std::cout << " col[" << mastColumnPtr->name() << "] of val = " << (*infoIt)->value << std::endl;

        std::cout << " spvar " << it->varPtr()->name() << " bound = " << it->val() << "   MCcoef = "
                  << mastColumnPtr->spVarVal(it->varPtr()) << std::endl;
      }

      if (it->satisfiedBy(mastColumnPtr->spVarVal(it->varPtr())))
      {
        it->cardinality(it->cardinality() + (*infoIt)->value);
        if (printL(6))
          std::cout << " added to cur card of " << it->varPtr()->name()
          << " card = " << it->cardinality() << std::endl;
      } else
      {
        it->complementCard(it->complementCard() + (*infoIt)->value);
        if (printL(6))
          std::cout << " added to cur complcard of " << it->varPtr()->name()
                    << " complcard = " << it->complementCard() << std::endl;

        // Goto next Column
        break;
      }
    }
  }

  _compBoundSet.fracWeight(curFracWeight);
  _compBoundSet.roundFracWeight();
  _currentCompBoundSet = _compBoundSet;
  _candidateLhs = curFracWeight;

  return;
}

/** @internal
 * Recursive method with major side effect on "global" 
 * variable _compBoundSet (IN OUT controlling variable)
 */
bool CompBdSetBranchConstrGenerator::nextNodeBrConstr(Node * parentNodePtr,
                                                      std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList,
                                                      const ConstrPtrSet & existingMasterBranchingConstr)
{
  int printLevel = 5;

  /// Node branching constraint defined last is treated first (LIFO)
  nextBranchingConstrPtrList.clear();

  if (printL(printLevel))
    std::cout
          << "CompBdSetBranchConstrGenerator::nextNodeBrConstr() _currentCompBoundSet = " << _currentCompBoundSet;

  int ancestorNodeRef(-1);
  if (parentNodePtr != NULL)
    ancestorNodeRef = parentNodePtr->ref();

  /// Branching on cardinality constraint
  if (_currentCompBoundSet.empty())
  {
    if (_currentCompBoundSet.allBranchesGenerated())
      return (false);

    if (printL(1))
      std::cout
            << "CompBdSetBranchConstrGenerator::nextNodeBrConstr() branching on cardinality constraint,"
               " empty component bound sequence " << std::endl;
  }

  if (_childNbCounter == 0)
  {
    if ((_direction == 'U') && (_currentCompBoundSet.activeSense() == 'L'))
      _currentCompBoundSet.complement();

    /// Clone
    instanciateBrConstr(ancestorNodeRef, ++_childNbCounter,
                        _currentCompBoundSet, nextBranchingConstrPtrList);
    if (printL(printLevel))
      {
        std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr(): constraint ";
        nextBranchingConstrPtrList.back()->shortPrint(std::cout);
        std::cout << " created (_childNbCounter == 0)" << std::endl;
      }
  } else
  {
    /**
     * Check whether existing br constraint is active:
     * see proposition 4 of paper gbrs
     */
    if (_childNbCounter > 1)
    {

      for (ConstrPtrSet::const_iterator constrPt = existingMasterBranchingConstr.begin();
           constrPt != existingMasterBranchingConstr.end(); constrPt++)
      {
        if ((*constrPt)->vcIndexStatus() == VcIndexStatus::Active)
        {
          if ((*constrPt)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
          {
            CompSetInstMastBranchConstr * csIbcPtr = static_cast<CompSetInstMastBranchConstr *> ((*constrPt));
            if (csIbcPtr->enforces(_currentCompBoundSet))
              return (false);
          }
        }
      }
    }

    /// When (_activeSense == 'L'), second node is rounded up (complemented twice)
    _currentCompBoundSet.complement();

    /// Clone
    instanciateBrConstr(ancestorNodeRef, ++_childNbCounter,
                        _currentCompBoundSet, nextBranchingConstrPtrList);
    if (printL(printLevel))
      {
        std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr(): constraint ";
        nextBranchingConstrPtrList.back()->shortPrint(std::cout);
        std::cout << " created (_childNbCounter > 0)" << std::endl;
      }

    if (!_currentCompBoundSet.empty())
    {
      if (printL(printLevel))
        std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr() _currentCompBoundSet.pop since size = "
                  << _currentCompBoundSet.size() << std::endl;
      _currentCompBoundSet.pop_back();
      _currentCompBoundSet.roundFracWeight();

      if (printL(printLevel))
        std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr() _currentCompBoundSet.empty = "
                  << _currentCompBoundSet.empty() << " compBoundSet._cgSpConfPtr->cardinalityIsFixed() = "
                  << _currentCompBoundSet._cgSpConfPtr->cardinalityIsFixed() << std::endl;

      /// Ruslan : this caused bug in GlassCuttingWithDefects/tests/nonRegressionTests/testBugBranching2
      /// therefore "return false;" is changed to  "_currentCompBoundSet.allBranchesGenerated(true);"
      if ((_currentCompBoundSet.empty()) && (_currentCompBoundSet._cgSpConfPtr->cardinalityIsFixed()))
        _currentCompBoundSet.allBranchesGenerated(true);
      //return false;

    } else
    {
      /// As the two nodes have already been generated
      _currentCompBoundSet.allBranchesGenerated(true);
    }
  }

  if (!(*(nextBranchingConstrPtrList.rbegin()))->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
  {
    return (true);
  }

  CompSetInstMastBranchConstr * curCsBcPtr
    = static_cast<CompSetInstMastBranchConstr *> (*(nextBranchingConstrPtrList.rbegin()));

  /** 
   * Check whether branching constraint is feasible; 
   * else node should be pruned by infeasibility 
   * and it might as well not be generated.
   */
  if (curCsBcPtr->incompatibleCsBd())
  {
    if (printL(printLevel))
      std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr "
                   "DO NOT RECORD THIS NODE as  BC  makes the problem infeasible " << std::endl;

    nextBranchingConstrPtrList.clear();

    /// Recursive call
    return nextNodeBrConstr(parentNodePtr, nextBranchingConstrPtrList, existingMasterBranchingConstr);
  }

  /**
   * Else, check whether ancestor node 
   * has a direct sucessor that defines 
   * the same constraint (see Observation 6 of gbrs paper)
   */

  if (printL(printLevel))
    {
      std::cout << "Trying to add branching constraint ";
      curCsBcPtr->shortPrint();
      std::cout << std::endl;
    }

  if (parentNodePtr == NULL)
    /// Current node is dummy
    return (true);

  if (parentNodePtr->isRoot())
    /// Current node is the root
    return (true);

  Node * curSonPtr = parentNodePtr;
  Node * bbNodePtr = parentNodePtr->father();

  /// I.e., while the current node is not the root of the B-a-B tree
  while (bbNodePtr != NULL)
  {
    if (bbNodePtr->sons().empty())
    {
      std::cout
              << "CompBdSetBranchConstrGenerator::nextNodeBrConstr(): localNodeBrConstrList() empty"
              << std::endl;
      exit(1);
    }

    for (std::list<Node *>::const_iterator sonNodePt =
         bbNodePtr->sons().begin(); sonNodePt != bbNodePtr->sons().end();
         sonNodePt++)
    {
      if (*sonNodePt == curSonPtr)
        continue;
      for (std::list<BranchingConstrBaseType *>::const_iterator it =
           (*sonNodePt)->localNodeBrConstrList().begin();
           it != (*sonNodePt)->localNodeBrConstrList().end(); ++it)
      {
        /// Branching constraint is of type other than CompSetInstMastBranchConstr
        if (!(*it)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
        {
          continue;
        }
        CompSetInstMastBranchConstr * cousinCsBcPtr =
                static_cast<CompSetInstMastBranchConstr *> (*it);

        /// Branching constraint exists in the direct sucessor of an ancestor
        if (*curCsBcPtr == *cousinCsBcPtr)
        {
            if (printL(printLevel))
              {
                std::cout << "CompBdSetBranchConstrGenerator::nextNodeBrConstr DO NOT RECORD THIS BC ";
                curCsBcPtr->shortPrint();
                std::cout << " which is identical to ";
                cousinCsBcPtr->shortPrint();
                std::cout << " of node " << (*sonNodePt)->branchAndPriceOrder() << std::endl;
              }

          nextBranchingConstrPtrList.clear();

          /// Recursive call
          return (nextNodeBrConstr(parentNodePtr, nextBranchingConstrPtrList, existingMasterBranchingConstr));
        }
      }
    }

    /// Goto next node up in the tree
    bbNodePtr = bbNodePtr->father();
  }

  return (true);
}

void CompBdSetBranchConstrGenerator
     ::instanciateBrConstr(const int & parentNodeNb, const int & childNb, const ComponentSequence & compBdSeq,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
  std::string name("BC");

  BranchingConstrBaseType * mBCptr(NULL);

  char kindMBC('E'); /// Implicit branching constraint seems not faster

  mBCptr = new CompSetInstMastBranchConstr(compBdSeq,IndexCell(), _genericBrConstrPtr,
                                           _genericBrConstrPtr->modelPtr()->masterConfPtr(),
                                           name + "cs", compBdSeq.classCardinality(), kindMBC);
  if (printL(6))
    {
      std::cout << " new CompSetInstMastBranchConstr " << std::endl;
      mBCptr->print(std::cout);
    }

  nextBranchingConstrPtrList.push_back(mBCptr);

  return;
}

CompBoundSetGenBranchConstr::CompBoundSetGenBranchConstr(Model * modelPtr,
                                                         GenericVar * genVarPtr,
                                                         const SelectionStrategy & priorityRule,
                                                         const Double & priorityLevel) :
    GenVarGenBranchConstr(modelPtr, modelPtr->masterConfPtr(), genVarPtr, priorityRule, priorityLevel),
    Base4NonLinearGenericConstr(NULL)
{
  toBeUsedInPreprocessing(false);
  
  return;
}

CompBoundSetGenBranchConstr::~CompBoundSetGenBranchConstr()
{
  return;
}

void CompBoundSetGenBranchConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
}

void CompBoundSetGenBranchConstr::updateGeneratedBrConstrGeneratorSet(const ComponentSequence & curCPSet,
                                                                      const ComponentBound & lastCP,
                                                                      const int & curAdditionalNbOfCpBd,
                                                                      const bool & pushAdditionalCP,
                                                                      const Double & cardinality,
                                                                      CompBoundSetGenBranchConstr * curCompBoundSetGenBranchConstr,
                                                                      BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  ComponentSequence curCandidate(curCPSet);
  curCandidate.additionalNbOfCpBd(curAdditionalNbOfCpBd);

  if (pushAdditionalCP)
  {
    curCandidate.push_back(lastCP);
    curCandidate.roundFracWeight();
  }
    else if (!curCandidate.empty())
    {
        curCandidate.rbegin()->cardinality(cardinality);
    }

  if (printL(6))
  {
    std::cout
            << " CompBoundSetGenBranchConstr::updateGeneratedBrConstrGeneratorSet(): NEW candidate = "
            << curCandidate << std::endl;
  }

    CompBdSetBranchConstrGenerator * compBdSetBranchConstrGenPtr;
    if (curCandidate.empty())
    {
        compBdSetBranchConstrGenPtr = new CompBdSetBranchConstrGenerator(curCompBoundSetGenBranchConstr, curCandidate,
                                                                         'U');
    }
    else
    {
        compBdSetBranchConstrGenPtr = new CompBdSetBranchConstrGenerator(curCompBoundSetGenBranchConstr, curCandidate,
                                                                         curCandidate.rbegin()->varPtr()->directive());
    }
    generatedBrConstrGeneratorSet.insert(compBdSetBranchConstrGenPtr);

  int candListMaxSize = 1;
  if (param().StrongBranchingPhaseOne().active())
    candListMaxSize = param().StrongBranchingPhaseOne().maxNumOfCandidates();

  if ((int) generatedBrConstrGeneratorSet.size() > candListMaxSize)
  {
    if (printL(6))
    {
      std::cout << " CompBoundSetGenBranchConstr::updateGeneratedBrConstrGeneratorSet(): "
                << "remove last  candidate of list of size = " << generatedBrConstrGeneratorSet.size() << std::endl;
    }
    generatedBrConstrGeneratorSet.erase(--(generatedBrConstrGeneratorSet.end()));
  }

  return;

}

void CompBoundSetGenBranchConstr
     ::branchingSeparationFindCandidates(const MasterColSolution & listOfFractMastCol,
                                         const int & maxNumOfCandidates,
                                         BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  /// Columns in curSol are assumed to be sorted in ILO
  MasterConf * mastConfPtr = modelPtr()->masterConfPtr();
  bapcodInit().check(
                     mastConfPtr == NULL,
                     "CompBoundSetGenBranchConstr::branchingSeparationRoutine(); masterConf should be defined");

  if (printL(6))
  {
    std::cout
            << "CompBoundSetGenBranchConstr::separationRoutine: nb of fract mast var "
            << listOfFractMastCol.size() << std::endl;
  }

  /// Implement separation independently for each subproblem
  for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPt = mastConfPtr->colGenSubProbConfPts().begin();
       cgSpConfPt != mastConfPtr->colGenSubProbConfPts().end(); ++cgSpConfPt)
  {
    if (genVarPtr()->probConfPtr() != *cgSpConfPt)
      /// Not a SP concerned by the gen var on which we branch
      continue;

    ComponentSequence curClassCompBoundSet(*cgSpConfPt);
    const MasterColSolution & curListOfFractMastCol = (*cgSpConfPt)->listOfFractMastColInColGenSp();

    /// Show the probconf associated to current GenVar
    if (printL(6))
    {
      std::cout << "genVar used in CompBoundSet is included in probConfig name = "
              << (*cgSpConfPt)->name() << std::endl;
    }

    VarPtrSet candVarForCompBoundBranching;

    for (IndexCell2InstancVarPtrMap::const_iterator ivPt = genVarPtr()
         ->indexCell2InstancVarPtrMap().begin();
         ivPt != genVarPtr()->indexCell2InstancVarPtrMap().end(); ivPt++)
    {
      if (ivPt->second->isCandForBranching()) //((*ivPt)->type() == 'B') || ((*ivPt)->type() == 'I'))
      {
        candVarForCompBoundBranching.insert(ivPt->second);
        if (printL(6))
        {
          std::cout
                  << "CompBoundSetGenBranchConstr::separationRoutine: cand sp var "
                  << ivPt->second->name() << std::endl;
        }
      }
    }

    //int minAdditionalNbOfCpBd(candVarForCompBoundBranching.size() + 1);
    Double totFracColVal(-1); // default value if not yet computed: it shall be computed in first call to exploreFrac or separateFra

    if (printL(5))
      {
          std::cout << "Current treeOfColClasses of " << (*cgSpConfPt)->name() << " : "
              << std::endl;
          for (ColClassesVector::const_iterator it =
              (*cgSpConfPt)->treeOfColClasses().begin();
              it != (*cgSpConfPt)->treeOfColClasses().end(); ++it)
            {
              std::cout << "Class " << (*it)->name();
              (*it)->shortPrint();
              std::cout << std::endl;
              for (ComponentSequence::const_iterator csIt = (*it)->compBoundSet().begin();
                  csIt != (*it)->compBoundSet().end(); ++csIt)
                std::cout << "   " << *csIt;
            }
      }

    /**
     * Select branching set that yields smallest number
     * of additional cpBd fixing and hence b-a-p nodes
     */
    if (_genVarPtr == _genVarPtr->probConfPtr()->defaultGenericVarPtr())
        separateCardinality(curListOfFractMastCol, curClassCompBoundSet, generatedBrConstrGeneratorSet);
    else
        exploreFrac((*cgSpConfPt)->treeOfColClasses(),
                    curListOfFractMastCol,
                    candVarForCompBoundBranching,
                    curClassCompBoundSet, totFracColVal, 0,
                    generatedBrConstrGeneratorSet);

  }
  
  /// if highest priority rule is used, we leave only highest priority candidates
  /// coming from this generic branching constraint
  if (priorityRule() == SelectionStrategy::HighestPriority)
    {
      Double highestPriority(0.0);
      for (BranchGeneratorsSet::iterator brConstrGenPtrIt = generatedBrConstrGeneratorSet.begin();
           brConstrGenPtrIt != generatedBrConstrGeneratorSet.end(); ++brConstrGenPtrIt)
        if (((*brConstrGenPtrIt)->genericBrConstrPtr() == this) && ((*brConstrGenPtrIt)->priority() > highestPriority))
          {
            highestPriority = (*brConstrGenPtrIt)->priority();
          }
      for (BranchGeneratorsSet::iterator brConstrGenPtrIt = generatedBrConstrGeneratorSet.begin();
           brConstrGenPtrIt != generatedBrConstrGeneratorSet.end();)
        {
          if (((*brConstrGenPtrIt)->genericBrConstrPtr() == this)
              && ((*brConstrGenPtrIt)->priority() < highestPriority))
            {
              delete *brConstrGenPtrIt;
              generatedBrConstrGeneratorSet.erase(brConstrGenPtrIt++);
            }
          else
            brConstrGenPtrIt++;
        }
          
    }
  
  return; //-- not path to this instruction
} //-- CompBoundSetGenBranchConstr::separationRoutine

/**
 * Takes columns of curListOfMasterCol
 * and returns them in a sorted list
 * ILOsortedListOfMastCol in the lexicographic
 * order dictated by the class partition of
 * current branching constaints treeOfColClasses,
 * whose component bound already treated a
 * memorised in curClassCompBoundSet
 */
void CompBoundSetGenBranchConstr::ILOsortMastColumn(const ColClassesVector & treeOfColClasses,
                                                    const MasterColSolution & curListOfMasterCol,
                                                    ComponentSequence & curClassCompBoundSet,
                                                    MasterColSolution & ILOsortedListOfMastCol)
{
  if (curListOfMasterCol.empty())
    return;

  if (treeOfColClasses.empty())
  {
    if (printL(6))
    {
      for (MasterColSolution::const_iterator it =
           curListOfMasterCol.begin();
           it != curListOfMasterCol.end();
           ++it)
      {
        std::cout << " remaining MastColumn col = "
                << it->first->name()
                << std::endl;
      }
    }

    /// Sort in lexicographic order allows to catch hiden integer solution

    std::vector< std::pair < MastColumn *, ValueRecord > > listOfColInLexicographicOrder;

    for (MasterColSolution::const_iterator it =
         curListOfMasterCol.begin(); it != curListOfMasterCol.end(); ++it)
      listOfColInLexicographicOrder.push_back(*it);

    stable_sort(listOfColInLexicographicOrder.begin(),
                listOfColInLexicographicOrder.end(),
                LexicographicMastColValSorting());

    if (printL(6))
      for (std::vector< std::pair < MastColumn *, ValueRecord > >::const_iterator it =
           listOfColInLexicographicOrder.begin();
           it != listOfColInLexicographicOrder.end(); ++it)
        std::cout << " Before: partial sorted list of MastColumn = "
              << it->first->name() << std::endl;

    ///@todo check that columns are approprietely sorted by subproblem first.

    ILOsortedListOfMastCol.insert(ILOsortedListOfMastCol.end(),
                                  listOfColInLexicographicOrder.begin(),
                                  listOfColInLexicographicOrder.end());
    //setOfColInLexicographicOrder.end());

    if (printL(6))
      for (MasterColSolution::const_iterator it = ILOsortedListOfMastCol
           .begin(); it != ILOsortedListOfMastCol.end(); ++it)
        std::cout << " ILOsortedListOfMastCol = " << it->first->name() << std::endl;
    return;
  }

  /// Select any br constraint
  CompSetInstMastBranchConstr * bcPtr = treeOfColClasses.front();

  /// Move up next cp bound
  ComponentSequence::const_iterator bcCsIt;
  bcCsIt = bcPtr->compBoundSet().begin();

  /// Pass previously observed cp bound
  for (ComponentSequence::const_iterator bsIt = curClassCompBoundSet.begin();
       bsIt != curClassCompBoundSet.end(); ++bsIt)
  {
    if (bcCsIt == bcPtr->compBoundSet().end())
    {
      std::cout
              << "error ILOsortMastColumn: bcCsIt == bcPtr->compBoundSet()->end()"
              << std::endl;
      exit(1);
    }

    ++bcCsIt;
  }
  if (bcCsIt == bcPtr->compBoundSet().end())
  {
    std::cout
            << "error ILOsortMastColumn: bcCsIt == bcPtr->compBoundSet()->end()"
            << std::endl;
    exit(1);
  }

  ComponentBound curCpBd(*bcCsIt);

  /// Take round-up bound
  if (curCpBd.sign() == 'L')
    curCpBd.complement();

  if (printL(6))
    std::cout << "ILOsortMastColumn curCpBd = " << curCpBd << std::endl;

  /// Select next component bound et reset spVarInvolvedInFracMastCol if need be
  MasterColSolution colSatisfying;
  MasterColSolution colNotSatisfying;

  for (MasterColSolution::const_iterator colPt =
       curListOfMasterCol.begin();
       colPt != curListOfMasterCol.end();
       colPt++)
  {
    if (printL(6))
      std::cout << "ILOsortMastColumn test col " << colPt->first->name()
      << std::endl;

    /// Col satisfy cp bound
    //Commented by Romain: dynamic_cast is not necessary: colPt is already a MastColumn*
    //    if (curCpBd.satisfiedBy(
    //        (dynamic_cast<const MastColumn * const >(*colPt))->spVarVal(
    //            curCpBd.varPtr())))
    if (curCpBd.satisfiedBy(colPt->first->spVarVal(curCpBd.varPtr())))
    {
      if (printL(6))
        std::cout << "ILOsortMastColumn satisfied " << std::endl;

      colSatisfying.push_back(colPt->first, colPt->second);
    } else
    {
      if (printL(6))
        std::cout << "ILOsortMastColumn not satisfied " << std::endl;

      colNotSatisfying.push_back(colPt->first, colPt->second);
    }
  }

  /// Filter branching constraint
  ColClassesVector remainingSatBrConstr;
  ColClassesVector remainingNotSatBrConstr;
  for (ColClassesVector::const_iterator it =
       treeOfColClasses.begin(); it != treeOfColClasses.end(); ++it)
  {
    /// beaware: bcCsIt2 intially overlapped bcCsIt
    ComponentSequence::const_iterator bcCsIt2, bcCsItNext;
    bcCsIt2 = (*it)->compBoundSet().begin();

    /// Pass previously observed cp bound
    for (ComponentSequence::const_iterator bsIt = curClassCompBoundSet.begin();
         bsIt != curClassCompBoundSet.end(); bsIt++)
    {
      if (bcCsIt2 == (*it)->compBoundSet().end())
      {
        std::cout << "error bcCsIt2 == (*it)->compBoundSet()->end()"
                << std::endl;
        exit(1);
      }

      bcCsIt2++;
    }
    if (bcCsIt2 != (*it)->compBoundSet().end())
    {
      bcCsItNext = bcCsIt2;
      bcCsItNext++;

      /// If their are still others behind)
      if (bcCsItNext != (*it)->compBoundSet().end())
      {
        /// Check next one
        if (*bcCsIt2 == curCpBd)
          remainingSatBrConstr.push_back(*it);
        else
          remainingNotSatBrConstr.push_back(*it);
      }
    }
  }

  if (!colSatisfying.empty())
  {
    ComponentSequence newClassCompBoundSet(curClassCompBoundSet);
    newClassCompBoundSet.reset();

    /// Add cp bound
    newClassCompBoundSet.push_back(curCpBd);

    /**
     * remaing columns in class shall be sorted in lexicographic
     * order in last step of recursive calls
     */
    ILOsortMastColumn(remainingSatBrConstr,
                      colSatisfying,
                      newClassCompBoundSet,
                      ILOsortedListOfMastCol);
  }

  if (!colNotSatisfying.empty())
  {
    curCpBd.complement();

    /// Add cp bound
    curClassCompBoundSet.push_back(curCpBd);

    /**
     * Remaing columns in class shall be sorted in lexicographic
     * order in last step of recursive calls
     */
    ILOsortMastColumn(remainingNotSatBrConstr,
                      colNotSatisfying,
                      curClassCompBoundSet,
                      ILOsortedListOfMastCol);
  }

  return;
}

/** 
 * 
 * 
 * @param treeOfColClasses C 
 * @param curListOfFractMastCol F
 * @param candVarForCompBoundBranching I 
 * @param curClassCompBoundSet S, p is implicitly given by the size of S
 * @param totFracColVal Records the overall value of the fractional columns in the current master solution
 * @param additionalNbOfCpBd  Incremented  in each recursive call to separate
 * @param recordedBestBranchingSet Record of the selected component set for branching along with charateristics
 */
void CompBoundSetGenBranchConstr::exploreFrac(const ColClassesVector & treeOfColClasses,
                                              const MasterColSolution & curListOfFractMastCol,
                                              VarPtrSet & candVarForCompBoundBranching,
                                              ComponentSequence & curClassCompBoundSet,
                                              Double & totFracColVal,
                                              int additionalNbOfCpBd,
                                              BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  int printLevel(5);

  if (curListOfFractMastCol.empty())
  {
    if (printL(printLevel))
      std::cout << "CompBoundSetGenBranchConstr::exploreFrac() received no frac col" << std::endl;

    return;
  }

  if (additionalNbOfCpBd > int(genVarPtr()->indexCell2InstancVarPtrMap().size()))
  {
    std::cout << "CompBoundSetGenBranchConstr::exploreFrac() is cycling in recursive calls" << std::endl;
    exit(1);
  }

  int nbOfBranchToTest(2);
  bool fullSearchForBestBranchingSet(true); /// changed by Ruslan to true
  ComponentBound curCpBd;

  /**
   * Step 1: if one reach a leaf node of 
   * the nested tree of column classes, call separate
   */
  if (treeOfColClasses.empty())
  {
      separateFrac(curListOfFractMastCol, candVarForCompBoundBranching,curClassCompBoundSet,totFracColVal,
                   additionalNbOfCpBd,generatedBrConstrGeneratorSet);

      return;
  }

  /** 
   * Step 2.a: Check whether the current class 
   * $Q(curClassCompBoundSet)$ has fractional value
   */
  Double totalFracWeight(0);
  for (MasterColSolution::const_iterator colPt = curListOfFractMastCol.begin(); colPt != curListOfFractMastCol.end();
       colPt++)
  {
    totalFracWeight += colPt->second._lfracValue;
    if (printL(printLevel))
      std::cout << "exploreFrac(): col set includes in ILO order "
                << colPt->first->name() << " of temp val " << colPt->second._lfracValue
                << " totalFracWeight = " << totalFracWeight << std::endl;
  }

  curClassCompBoundSet.fracWeight(totalFracWeight);

  /**
   * Step a bis: Check whether the current column 
   * class $Q(curClassCompBoundSet)$ has fractional value
   * To try set it to true
   */
  bool testFractionalityOfCurClass = _genVarPtr == _genVarPtr->probConfPtr()->defaultGenericVarPtr();


  /**
   * If false do not branch on current convexity 
   * constraint curClassCompBoundSet.empty());
   */
  if (testFractionalityOfCurClass)
  {
    Double totalCardinalityOfFractCol(0);
    Double fracPart(Dfrac(totalFracWeight));
    if (fracPart > param().BapCodCutViolationTolerance)
    {
      if (printL(printLevel))
        std::cout << "exploreFrac: curClassCompBoundSet yield fract weight = " << fracPart << std::endl;

      updateGeneratedBrConstrGeneratorSet(curClassCompBoundSet, curClassCompBoundSet.lastCP(), additionalNbOfCpBd,
                                          false, totalFracWeight, this,generatedBrConstrGeneratorSet);
      /**
       * Either do not try to split further a class that is fractional:
       * or try only one branch down: nbOfBranchToTest--;
       * or let's consider other original variables for branching
       */
      return;
    }
  }

  /**
   * Step 2.b: Select the next component in ILO order
   * Select any br constraint
   */

  CompSetInstMastBranchConstr * bcPtr = treeOfColClasses.front();

  /// Move up next cp bound
  ComponentSequence::const_iterator bcCsIt;
  bcCsIt = bcPtr->compBoundSet().begin();


  /// Pass previously observed cp bound
  for (ComponentSequence::const_iterator bsIt = curClassCompBoundSet.begin();
       bsIt != curClassCompBoundSet.end(); bsIt++)
  {
    if (bcCsIt == bcPtr->compBoundSet().end())
    {
      std::cout << "error bcCsIt == bcPtr->compBoundSet()->end()" << std::endl;
      exit(1);
    }

    bcCsIt++;
  }

  if (bcCsIt == bcPtr->compBoundSet().end())
  {
    std::cout << "error bcCsIt == bcPtr->compBoundSet()->end()" << std::endl;
    exit(1);
  }

  curCpBd = *bcCsIt;

  /// Take round-up bound
  if (curCpBd.sign() == 'L')
    curCpBd.complement();

  if (printL(printLevel))
    std::cout << "exploreFrac curCpBd = " << curCpBd << std::endl;

  /**
   * Step 2.c and d: Check whether the current class 
   * $Q(curClassCompBoundSet, curCpBd)$ has fractional value, 
   * while preparing the partition
   */
  Double weightOfSatisfyingCol(0);
  MasterColSolution colSatisfying;
  MasterColSolution colNotSatisfying;

  SortMastColPerDecreasingSpVal comparator(curCpBd.varPtr()); //(Issam)
  std::multiset<pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal> ILOsortedSatisfyingCol(comparator);
  int nbOfSatisfyingFracCol(0);
  int nbOfNotSatisfyingFracCol(0);

  for (MasterColSolution::const_iterator colPt = curListOfFractMastCol.begin(); 
       colPt != curListOfFractMastCol.end(); 
       colPt++)
    {

      if (printL(printLevel))
	std::cout << "exploreFrac test col " << colPt->first->name()
		  << " at tmp val = " << colPt->second._lfracValue << std::endl;

      /// Col statisfy cp bound
      if (curCpBd.satisfiedBy(colPt->first->spVarVal(curCpBd.varPtr())))
	{
	  colSatisfying.push_back(colPt->first, colPt->second);
	  colPt->first->fillSpIndicatorColMultiset(ILOsortedSatisfyingCol,
						   curCpBd.varPtr(),
						   colPt->second);
	  if (printL(printLevel))
	    std::cout << "exploreFrac satisfied by col for which sp val = "
		      << colPt->first->spVarVal(curCpBd.varPtr()) << std::endl;

	  weightOfSatisfyingCol += colPt->second._lfracValue;
	  if (colPt->second._isFractional)
	    nbOfSatisfyingFracCol++;
	} 
      else
	{
	  if (printL(printLevel))
	    std::cout << "exploreFrac not satisfied by col for which sp val = "
		      << colPt->first->spVarVal(curCpBd.varPtr()) << std::endl;

	  colNotSatisfying.push_back(colPt->first, colPt->second);
	  if (colPt->second._isFractional)
	    nbOfNotSatisfyingFracCol++;
	}
    }

  Double colClassWeightFracPart(Dfrac(weightOfSatisfyingCol));
  int nbOfCol(curListOfFractMastCol.size());
  int nbOfSatisfyingCol(colSatisfying.size());
  int nbOfNotSatisfyingCol(colNotSatisfying.size());
  bool curColClassHasFractionalValue(colClassWeightFracPart > param().BapCodCutViolationTolerance);

  if (printL(printLevel))
    std::cout << " exploreFrac: nbOfCol = " << nbOfCol
	      << " nbOfSatisfyingCol = " << nbOfSatisfyingCol
	      << " nbOfSatisfyingFracCol = " << nbOfSatisfyingFracCol
	      << " nbOfNotSatisfyingCol = " << nbOfNotSatisfyingCol
	      << " ILOsortedSatisfyingCol.size() = "
	      << ILOsortedSatisfyingCol.size() << " nbOfNotSatisfyingFracCol = "
	      << nbOfNotSatisfyingFracCol << " weightOfSatisfyingCol = "
	      << weightOfSatisfyingCol << " colClassWeightFracPart = "
	      << colClassWeightFracPart << " curColClassHasFractionalValue = "
	      << curColClassHasFractionalValue << std::endl;

  /**
   * In integer case  (curCpBd.varPtr()->type() == 'I'), 
   * column class can be integer while column class 
   * is not but fractional set must be fractional
   * in both the original var and column class weight
   */
  if (curColClassHasFractionalValue)
  {
    /// Check whether the original var is fractional
    int tindex(0);
    Double cumVal(0);
    std::vector<Double> itemDissagregVect;
    itemDissagregVect.push_back(0);
    itemDissagregVect[tindex] = 0;
    std::multiset< std::pair < MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal>::const_iterator mastColPt;
    for (mastColPt = ILOsortedSatisfyingCol.begin(); mastColPt != ILOsortedSatisfyingCol.end(); mastColPt++)
    {
      Double itemVal(mastColPt->first->spVarVal(curCpBd.varPtr()));
      if (printL(printLevel))
      {
        std::cout << "exploreFrac: col[ " << mastColPt->first->name()
                << " ] has fractional part " << mastColPt->second._lfracValue
                << " and item " << curCpBd.varPtr()->name() << " with val "
                << itemVal << std::endl;
      }

      Double tmpColVal(mastColPt->second._lfracValue);
      Double lambdaR(0);
      while (tmpColVal > 0)
      {
        if (cumVal == (tindex + 1))
        {
          tindex++;
          itemDissagregVect.push_back(0);
          itemDissagregVect[tindex] = 0;
        }
        lambdaR = tmpColVal;
        if (cumVal + lambdaR > tindex + 1)
          lambdaR = (tindex + 1 - cumVal);

        itemDissagregVect[tindex] += itemVal * lambdaR;
        tmpColVal -= lambdaR;
        cumVal += lambdaR;
        if (printL(printLevel))
          std::cout << "exploreFrac: itemDissagregVect[" << tindex << "] "
                << itemDissagregVect[tindex] << " cumVal = " << cumVal
                << std::endl;
      }
    }
    /// Cumval denote the cumulative frac weight of columns containting spVar
    int lastIndex(tindex);
    Double componentVal(itemDissagregVect[lastIndex]);

    /**
     * Only check the last index value that
     * determines the boundary of the current column class;
     * indeed, branching on intermediate index would
     * correspond to partitionnng the current class;
     * this can only be done in separateFrac, after
     * the existing hierarchy of column classes
     * has been explored for fractinality
     */
    Double spVarFracPart(Dfrac(componentVal));

    /// Must be fractional inthe original var
    bool orignalVarValIsFractional(
                                   spVarFracPart > param().BapCodCutViolationTolerance);
    componentVal.Cceil();
    bool candBdIsEqualToLastClassCpBd(componentVal == curCpBd.val());
    if (printL(printLevel))
      std::cout << "exploreFrac: last index = " << lastIndex << ", "
            << " test componentVal = " << componentVal
            << " curCpBd.val() = " << curCpBd.val() << " spVarFracPart = "
      << spVarFracPart << " orignalVarValIsFractional = "
            << orignalVarValIsFractional
            << " candBdIsEqualToLastClassCpBd = "
            << candBdIsEqualToLastClassCpBd << std::endl;

    /// Must be fractional in both the original var and column class weight
    if (orignalVarValIsFractional && candBdIsEqualToLastClassCpBd
        && curColClassHasFractionalValue)
    {
      if (printL(printLevel))
        std::cout << "exploreFrac: curClassCompBoundSet yield fact weight = "
              << spVarFracPart << std::endl;

        curCpBd.cardinality(weightOfSatisfyingCol);
      updateGeneratedBrConstrGeneratorSet(curClassCompBoundSet, curCpBd,
                                          additionalNbOfCpBd, true, 0.0, this,
                                          generatedBrConstrGeneratorSet);
      /**
       * Either do not try to split further a class that is fractional:
       * or try only one branch down: nbOfBranchToTest--;
       * or let's consider other original variables for branching
       */
      if (!fullSearchForBestBranchingSet)
        return;
    } else
    {
      /// If curColClassHasFractionalValue explore further the current class
      if (printL(printLevel))
        std::cout << "class is not fractional; " << std::endl;
    }
  }

  /// Else prepare partition
  int nbFracColRequired(((curCpBd.varPtr()->type() == 'I') ? 1 : 2));
  bool testFirstBranch(nbOfSatisfyingFracCol >= nbFracColRequired);
  bool testSecondBranch(nbOfNotSatisfyingFracCol >= nbFracColRequired);
  if (testFirstBranch && testSecondBranch && (nbOfBranchToTest <= 1))
  {

    /// Go for larger fractional set
     if (nbOfSatisfyingCol > nbOfNotSatisfyingCol) testSecondBranch = false;
     else testFirstBranch = false;
  }

  /// Step 2.e: test first side of the partition
  if (testFirstBranch && (nbOfBranchToTest >= 1))
  {

    if (printL(printLevel))
      std::cout << "exploreFrac IN FIRST newCpBd = " << curCpBd << std::endl;

    /// Filter branching constraint
    ColClassesVector remainingBrConstr;
    for (ColClassesVector::const_iterator it =
         treeOfColClasses.begin(); it != treeOfColClasses.end(); ++it)
    {
      /// Attention, bcCsIt2 cover bcCsIt initialy
      ComponentSequence::const_iterator bcCsIt2;
      bcCsIt2 = (*it)->compBoundSet().begin();

      /// Pass previously observed cp bound
      for (ComponentSequence::const_iterator bsIt =
           curClassCompBoundSet.begin(); bsIt != curClassCompBoundSet.end();
           bsIt++)
      {
        if (bcCsIt2 == (*it)->compBoundSet().end())
        {
          std::cout << "error bcCsIt2 == (*it)->compBoundSet()->end()"
                  << std::endl;
          exit(1);
        }

        bcCsIt2++;
      }
      if (bcCsIt2 != (*it)->compBoundSet().end())
      {
        /// Check next one
        if (*bcCsIt2 == curCpBd)
        {
          if (++bcCsIt2 != (*it)->compBoundSet().end())
          {
            /// Their are stil others behind
            remainingBrConstr.push_back(*it);
          }
        }
      }
    }

    /// Add cp bound
    curClassCompBoundSet.push_back(curCpBd);
    curClassCompBoundSet.roundFracWeight();

    /// Reset set of SpVar for Branching
    if (curCpBd.varPtr()->ub() == 1)
    {
      /// Remove it
      candVarForCompBoundBranching.erase(curCpBd.varPtr());

      exploreFrac(remainingBrConstr,
                  colSatisfying,
                  candVarForCompBoundBranching,
                  curClassCompBoundSet,
                  totFracColVal, additionalNbOfCpBd,
                  generatedBrConstrGeneratorSet);

      /// Reinsert it
      candVarForCompBoundBranching.insert(curCpBd.varPtr());
    } else
    {
      exploreFrac(remainingBrConstr,
                  colSatisfying,
                  candVarForCompBoundBranching,
                  curClassCompBoundSet,
                  totFracColVal, additionalNbOfCpBd,
                  //recordedBestBranchingSet,
                  generatedBrConstrGeneratorSet);
    }

    /// Remove cp bound
    curClassCompBoundSet.pop_back();
    curClassCompBoundSet.roundFracWeight();
  }

  /// Step 2.F: Explore further complement of current class
  if (testSecondBranch && (nbOfBranchToTest >= 1))
  {
    curCpBd.complement();
    if (printL(printLevel))
      std::cout << "exploreFrac IN SECOND newCpBd = " << curCpBd << std::endl;

    /// Filter branching constraint
    ColClassesVector remainingBrConstr;
    for (ColClassesVector::const_iterator it =
         treeOfColClasses.begin(); it != treeOfColClasses.end(); ++it)
    {
      /// Attention, bcCsIt2 cover bcCsIt initialy
      ComponentSequence::const_iterator bcCsIt2;
      bcCsIt2 = (*it)->compBoundSet().begin();

      /// Pass previously observed cp bound
      for (ComponentSequence::const_iterator bsIt =
           curClassCompBoundSet.begin(); bsIt != curClassCompBoundSet.end();
           bsIt++)
      {
        if (bcCsIt2 == (*it)->compBoundSet().end())
        {
          std::cout << "error bcCsIt2 == (*it)->compBoundSet()->end()"
                  << std::endl;
          exit(1);
        }

        bcCsIt2++;
      }
      if (bcCsIt2 != (*it)->compBoundSet().end())
      {
        /// Check next one
        if (*bcCsIt2 == curCpBd)
        {
          if (++bcCsIt2 != (*it)->compBoundSet().end())
          {
            /// Their are still others behind
            remainingBrConstr.push_back(*it);
          }
        }
      }
    }

    curClassCompBoundSet.push_back(curCpBd);
    curClassCompBoundSet.roundFracWeight();

    /// Reset set of SpVar for Branching
    if (curCpBd.varPtr()->ub() == 1)
    {
      // Remove it
      candVarForCompBoundBranching.erase(curCpBd.varPtr());
      exploreFrac(remainingBrConstr,
                  colNotSatisfying,
                  candVarForCompBoundBranching,
                  curClassCompBoundSet,
                  totFracColVal,
                  additionalNbOfCpBd,
                  //recordedBestBranchingSet,
                  generatedBrConstrGeneratorSet);

      // Reinsert it
      candVarForCompBoundBranching.insert(curCpBd.varPtr());
    } else
    {

      exploreFrac(remainingBrConstr,
                  colNotSatisfying,
                  candVarForCompBoundBranching,
                  curClassCompBoundSet,
                  totFracColVal,
                  additionalNbOfCpBd,
                  //recordedBestBranchingSet,
                  generatedBrConstrGeneratorSet);
    }

    curClassCompBoundSet.pop_back();
    curClassCompBoundSet.roundFracWeight();
  }

  return;
}

void CompBoundSetGenBranchConstr::separateCardinality(const MasterColSolution & curListOfFractMastCol,
                                                      ComponentSequence & curClassCompBoundSet,
                                                      BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
    int localPrintLevel(5);
    if (curListOfFractMastCol.empty())
        return;
    Double totalFracWeight(0);
    for (MasterColSolution::const_iterator colPt = curListOfFractMastCol
         .begin(); colPt != curListOfFractMastCol.end(); colPt++)
    {
        totalFracWeight += colPt->second._lfracValue;
        if (printL(localPrintLevel))
            std::cout << "separateCardinality: col set includes in ILO order "
            << colPt->first->name() << " of temp val " << colPt->second._lfracValue
            << " totalFracWeight = " << totalFracWeight << std::endl;
        
    }
    Double fracPart(Dfrac(totalFracWeight));
    if (fracPart > param().BapCodCutViolationTolerance)
    {
        if (printL(localPrintLevel))
            std::cout << "separateCardinality: curClassCompBoundSet yields fract weight = "
            << fracPart << std::endl;
        
        updateGeneratedBrConstrGeneratorSet(curClassCompBoundSet,
                                            curClassCompBoundSet.lastCP(),
                                            0, false, totalFracWeight,
                                            this, generatedBrConstrGeneratorSet);
        
    }
}

void CompBoundSetGenBranchConstr::separateFrac(const MasterColSolution & curListOfFractMastCol,
                                               //const Solution & fractMastColInCurSol,
                                               VarPtrSet & candVarForCompBoundBranching,
                                               ComponentSequence & curClassCompBoundSet,
                                               Double & totFracColVal,
                                               int additionalNbOfCpBd,
                                               //RecordedBestCompSeq & recordedBestBranchingSet,
                                               BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{

  int localPrintLevel(5);
  if (printL(localPrintLevel))
    std::cout << "separateFrac(): curListOfFractMastCol "
          << curListOfFractMastCol.size() << std::endl;

  int nbOfBranchToTest(2);
  bool fullSearchForBestBranchingSet(true);
  ComponentBound curCpBd;

  /// Step A: Check whether the current column set has fractional columns
  if (curListOfFractMastCol.empty())
    return;

  /// Columns in fractMastColInCurSol are assumed to be sorted in ILO
  Double totalFracWeight(0);
  for (MasterColSolution::const_iterator colPt = curListOfFractMastCol.begin(); colPt != curListOfFractMastCol.end();
       colPt++)
  {
    totalFracWeight += colPt->second._lfracValue;
    if (printL(localPrintLevel))
      std::cout << "separateFrac: col set includes in ILO order "
            << colPt->first->name() << " of temp val " << colPt->second._lfracValue
            << " totalFracWeight = " << totalFracWeight << std::endl;

  }

  /**
   * Default value of totFracColVal is -1 so as 
   * to be intitialized on the first call trough this test
   */
  if (totFracColVal < 0)
    totFracColVal = totalFracWeight;

  curClassCompBoundSet.fracWeight(totalFracWeight);

  /**
   * Step a bis: Check whether the current column 
   * class $Q(curClassCompBoundSet)$ has fractional value
   */

  bool testFractionalityOfCurClass(true);

  /**
   * If false do not branch on current convexity 
   * constraint curClassCompBoundSet.empty());
   */
  if (testFractionalityOfCurClass)
  {
    Double totalCardinalityOfFractCol(0);
    Double fracPart(Dfrac(totalFracWeight));
    if (fracPart > param().BapCodCutViolationTolerance)
    {
      if (printL(localPrintLevel))
        std::cout << "separateFrac: curClassCompBoundSet yields fract weight = "
              << fracPart << std::endl;

      updateGeneratedBrConstrGeneratorSet(curClassCompBoundSet,
                                          curClassCompBoundSet.lastCP(),
                                          additionalNbOfCpBd,
                                          false, totalFracWeight,
                                          this,
                                          generatedBrConstrGeneratorSet);

      /**
       * Either do not try to split further a class that is fractional:
       * or try only one branch down: nbOfBranchToTest--;
       * or let's consider other original variables for branching
       */
      return;
    }
  }

  /// Further split the list of fractional columns using a component bound
  if (candVarForCompBoundBranching.empty())
    /// no more splitting is possible
    return;

  /**
   * Step b: Compute aggregate  values $(v_i)_{i \in I}$
   * enumerate on all spVar
   */
  std::map<Variable *,
          std::multiset<pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal>,
          Variable::PrioritySort> fractColInvolvingSpVar;

  /**
   * For each SpVar, we compute the list of column including this SpVar;
   * columns in list are sorted in lexicographic order because curListOfFractMastCol is.
   * columns in curListOfFractMastCol are assumed to be sorted in ILO
   * Columns need to be sorted in lexicographic order 
   * of this component value which is done by SortMastColPerDecreasingSpVal
   */
  for (MasterColSolution::const_iterator colPt = curListOfFractMastCol
       .begin(); colPt != curListOfFractMastCol.end(); colPt++)
  {
    colPt->first->fillSpIndicatorMap(fractColInvolvingSpVar,
                                     colPt->second,
                                     candVarForCompBoundBranching);
  }

  /**
   * Step c: Detect fractional $v_i$ if any 
   * and STEP d: Identify discriminating components
   */
  int tindex(0);
  Double cumVal(0);
  Variable * bestVar4FurtherSplit(NULL);
  Double bestSepIndex(0);
  bool keepOnlyBestVar4FurtherSplit = (additionalNbOfCpBd >= 0);
  ComponentBound bestCpBd4FurtherSplit;
  int bestNbOfSatisfyingFracCol(0);
  int bestNbOfNotSatisfyingFracCol(0);
  std::list<FurtherSplitCandidate> furtherSplitCandidates;
  bool foundCPforSepInThisCall2separate(false);

  /// Last component in Component Bound Set
  Variable * lastCptr(NULL);

  /// iterate on SP variables on which we could branch
  for (std::map<Variable *,
       std::multiset< std::pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal>,
       Variable::PrioritySort>::const_iterator spVarPt = fractColInvolvingSpVar
       .begin(); spVarPt != fractColInvolvingSpVar.end(); spVarPt++)
  {
    /**
     * Test if Variable is candidtate for being chosen
     * next to varPtr when branching on component bound set
     */
    if ((priorityRule() == SelectionStrategy::HighestPriority))
      {
        std::vector<InstanciatedVar *> varPts;
        curClassCompBoundSet.allCompLbVarPts(varPts);
        bool consecutiveToOneVar = false;
        for (std::vector<InstanciatedVar *>::iterator varPtrIt = varPts.begin(); varPtrIt != varPts.end(); ++ varPtrIt)
          if (spVarPt->first->consecutive2varWhenBrOnCBS(*varPtrIt))
            {
              consecutiveToOneVar = true;
              break;
            }
        if (!consecutiveToOneVar)
          continue;
      }
    else
      {
        if (!spVarPt->first->consecutive2varWhenBrOnCBS(curClassCompBoundSet.lastCompLbVarPtr()))
          continue;
      }

    if (printL(localPrintLevel))
      std::cout << "separateFrac: fractColInvolvingSpVar["
            << spVarPt->first->name() << "] " << std::endl;

    /**
     * We apply the procedure of Table 1 retricted to
     * the current SP var to identify a fractional
     * var value on which to branch.
     */
    std::vector<Double> itemDissagregVect;
    cumVal = 0;
    tindex = 0;
    itemDissagregVect.push_back(0);
    itemDissagregVect[tindex] = 0;
    if (spVarPt->second.empty())
    {
      std::cerr << "separateFrac: set of columns associated to sp var should not be empty"
                << std::endl;
      exit(1);
    }

    /**
     * spVarPt->second contains the list of
     * fract columns involving this SP var.
     * They need to be sorted in lexicographic
     * order of this component value
     */
    Double minItemVal(99999999);
    bool minItemValIsInitialized(false);
    bool thereExistsAtLeast2DifferentItemVal(false);

    for (std::multiset<pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal>::const_iterator mastColPt =
         spVarPt->second.begin(); mastColPt != spVarPt->second.end();
         mastColPt++)
    {
      Double itemVal(mastColPt->first->spVarVal(spVarPt->first));
      if (!minItemValIsInitialized)
      {
        minItemValIsInitialized = true;
        minItemVal = itemVal;
      } else if (minItemVal > itemVal)
      {
        minItemVal = itemVal;
        thereExistsAtLeast2DifferentItemVal = true;
      }
      if (printL(localPrintLevel))
      {
        std::cout << "col[ " << mastColPt->first->name()
                << " ] has fractional part " << mastColPt->second._lfracValue
                << " and item " << spVarPt->first->name() << " with val "
                << itemVal << " minItemVal " << minItemVal
                << " thereExistsAtLeast2DifferentItemVal "
                << thereExistsAtLeast2DifferentItemVal << std::endl;

        //(*mastColPt)->printColVector();
      }
      Double tmpColVal(mastColPt->second._lfracValue);
      Double lambdaR(0);
      while (tmpColVal > 0)
      {
        if (cumVal == (tindex + 1))
        {
          tindex++;
          itemDissagregVect.push_back(0);
          itemDissagregVect[tindex] = 0;
        }

        lambdaR = tmpColVal;
        if (cumVal + lambdaR > tindex + 1)
          lambdaR = (tindex + 1 - cumVal);

        itemDissagregVect[tindex] += itemVal * lambdaR;
        tmpColVal -= lambdaR;
        cumVal += lambdaR;
        if (printL(localPrintLevel))
          std::cout << "separateFrac: itemDissagregVect[" << tindex << "] "
                << itemDissagregVect[tindex] << " cumVal = " << cumVal
                << std::endl;
      } // while (tmpColVal > 0)
    } // for mastColPt

    // fv 19/06/2013: induce failing to find branching constraints
    // if (!thereExistsAtLeast2DifferentItemVal) continue;
    /// var is not a candidate to discriminate between columns

    /**
     * cumVal denote the cumulative frac
     * weight of columns containting spVar
     * take value equal to median value,
     * i.e. itemDissagregVect[Dfloor(cumval/2))]
     * go get median value (= 1 in binary case, = median in integer case)
     */
    int lastIndex(tindex);

    int index = (spVarPt->first->type() == 'B' ? lastIndex : 0);

    if (printL(localPrintLevel))
      std::cout << "separateFrac: lastIndex = " << lastIndex
            << " starting index = " << index << std::endl;

    for (; index <= lastIndex; ++index)
    {
      Double componentVal(itemDissagregVect[index]);
      Double spVarFracPart(Dfrac(componentVal));

      bool orignalVarValIsFractional(spVarFracPart > param().BapCodCutViolationTolerance);

      componentVal.Cceil();
      ComponentBound candCpBd(spVarPt->first, componentVal, 'G');
      int nbOfSatisfyingFracCol(0);
      int nbOfNotSatisfyingFracCol(0);
      Double totalCardinalityOfFractCol(0);

      /// Must also have fractional column class weight
      if (printL(localPrintLevel))
        std::cout << "separateFrac:  check CpBd  " << candCpBd;

      for (std::multiset<  std::pair< MastColumn *, ValueRecord >, SortMastColPerDecreasingSpVal>::const_iterator
           mastColPt = spVarPt->second.begin(); mastColPt != spVarPt->second.end(); ++mastColPt)
      {
        if (printL(localPrintLevel))
          std::cout << "separateFrac:   check col  " << mastColPt->first->name()
          << std::endl;
        if (candCpBd.satisfiedBy(mastColPt->first->spVarVal(candCpBd.varPtr())))
        {
          totalCardinalityOfFractCol += mastColPt->second._lfracValue;
          candCpBd.cardinality(candCpBd.cardinality() + mastColPt->second._value);
          if (mastColPt->second._isFractional)
          {
            ++nbOfSatisfyingFracCol;
            if (printL(localPrintLevel))
              std::cout << "separateFrac:  compare  SP var "
                    << spVarPt->first->name()
              << " add to nbOfSatisfyingFracCol" << std::endl;
          }
        } else
        {
            candCpBd.complementCard(candCpBd.complementCard() + mastColPt->second._value);
            if (mastColPt->second._isFractional)
              {
                ++nbOfNotSatisfyingFracCol;
                if (printL(localPrintLevel))
                  std::cout << "separateFrac:  compare  SP var "
                  << spVarPt->first->name()
                  << " add to nbOfNotSatisfyingFracCol" << std::endl;
              }
        }
      } // for mastColPt
      Double totalCardinalityFractPart(Dfrac(totalCardinalityOfFractCol));
      bool curColClassHasFractionalValue(
                                         totalCardinalityFractPart > param().BapCodCutViolationTolerance);

      /**
       * For int sp var we must make sure that sp fract
       * value satisfyies every comp bound in current list
       * so that the newly defined class is nested within existing class.
       * It can be identical to the existing class if the latter has fractional value
       */
      bool candBdIsWithinPreviousBound(true);
      bool candBdIsEqualToPreviousClassCpBd(false);

      if (spVarPt->first->type() == 'I') //|| (!curColClassHasFractionalValue))
      {
        ComponentSequence::const_iterator itcp = curClassCompBoundSet.begin();
        for (; itcp != curClassCompBoundSet.end(); itcp++)
        {
          if (printL(localPrintLevel))
            std::cout << "separateFrac:  compare  SP var "
                  << spVarPt->first->name() << " val  " << componentVal
                  << " with CP " << (*itcp);

          /// Require that new bound statisfy previous bound
          if (!(*itcp).satisfiedBy(spVarPt->first, componentVal))
          {
            if (printL(localPrintLevel))
              std::cout << "separateFrac:  cand CP on var "
                    << spVarPt->first->name() << " with val  "
              << componentVal << " does not statify CP " << (*itcp);

            candBdIsWithinPreviousBound = false;
          }

          /**
           * Require that new bound is different from previous bound;
           * if it was the same and lead to fract set it would have been
           * detected in explore in the case of binary var;
           * while if it is the same and leads to a column class
           * of frac cardinality (although not a frac orginal var val),
           * we must select another CP to further partition the class
           */
          if (*itcp == candCpBd)
          {
            if (printL(localPrintLevel))
              std::cout << "separateFrac:   cand CP is equal to previous bd  "
                    << (*itcp);

            candBdIsEqualToPreviousClassCpBd = true;
          }
        }
      }

      if (printL(localPrintLevel))
        std::cout << "separateFrac:  SP var" << spVarPt->first->name()
        << " index = " << index << " test componentVal = "
              << componentVal << " spVarFracPart = " << spVarFracPart
              << " orignalVarValIsFractional " << orignalVarValIsFractional
              << " nbOfSatisfyingFracCol " << nbOfSatisfyingFracCol
              << " nbOfNotSatisfyingFracCol " << nbOfNotSatisfyingFracCol
              << " totalCardinalityOfFractCol "
              << totalCardinalityOfFractCol
              << " curColClassHasFractionalValue "
              << curColClassHasFractionalValue
              << " totalCardinalityFractPart " << totalCardinalityFractPart
              << " candBdIsWithinPreviousBound "
              << candBdIsWithinPreviousBound
              << " candBdIsEqualToPreviousClassCpBd "
              << candBdIsEqualToPreviousClassCpBd << " candCpBd "
              << candCpBd << std::endl;

      /// Identification of branching set
      if (orignalVarValIsFractional && curColClassHasFractionalValue && candBdIsWithinPreviousBound)
      {
        if (printL(localPrintLevel))
          std::cout << "separateFrac: FOUND frac SP var"
                << spVarPt->first->name() << " with val = "
          << itemDissagregVect[index] << " spVarFracPart = "
                << spVarFracPart << " rounded sp var val = " << componentVal
                << " totalCardinalityFractPart = "
                << totalCardinalityFractPart << std::endl;

        updateGeneratedBrConstrGeneratorSet(curClassCompBoundSet, candCpBd,
                                            additionalNbOfCpBd + 1, true, 0.0,
                                            this,generatedBrConstrGeneratorSet);

        /// Or continue to test other potential fractional component
        foundCPforSepInThisCall2separate = true;
      }

      else if (candBdIsWithinPreviousBound
               && !candBdIsEqualToPreviousClassCpBd
               && (nbOfSatisfyingFracCol >= ((spVarPt->first->type() == 'I') ? 1 : 2)))
      {
        if (printL(localPrintLevel))
          std::cout << "separateFrac: candidateAsBestVar4FurtherSplit"
                << std::endl;

          if (keepOnlyBestVar4FurtherSplit)
          {
              bool candidateAsBestVar4FurtherSplit(true);
              /**
               * Need at least 2 columns in the binary case
               * and 1 in the integer case  in class for further splitting
               */
              Double sepIndex(Dfrac((double) (index + 1) / (double) (lastIndex + 1)));
              if (bestVar4FurtherSplit != NULL)
              {
                  if (spVarPt->first->priority() < bestVar4FurtherSplit->priority())
                      candidateAsBestVar4FurtherSplit = false;
                  else if ((spVarPt->first->priority()
                            <= bestVar4FurtherSplit->priority())
                           && (sepIndex < bestSepIndex))
                      candidateAsBestVar4FurtherSplit = false;
              }
              if (candidateAsBestVar4FurtherSplit)
              {
                  bestVar4FurtherSplit = spVarPt->first;
                  bestSepIndex = sepIndex;
                  bestCpBd4FurtherSplit = candCpBd;
                  bestNbOfSatisfyingFracCol = nbOfSatisfyingFracCol;
                  bestNbOfNotSatisfyingFracCol = nbOfNotSatisfyingFracCol;
                  if (printL(localPrintLevel))
                      std::cout << "separateFrac: bestVar4FurtherSplit = "
                      << spVarPt->first->name() << " bd = "
                      << candCpBd.val() << " sepIndex = " << sepIndex
                      << " nbOfSatisfyingFracCol = " << nbOfSatisfyingFracCol
                      << " nbOfNotSatisfyingFracCol = "
                      << nbOfNotSatisfyingFracCol << std::endl;
              }
          }
          else
          {
              furtherSplitCandidates.push_back(FurtherSplitCandidate(spVarPt->first, candCpBd,
                                                                     nbOfSatisfyingFracCol, nbOfNotSatisfyingFracCol));
              if (printL(localPrintLevel))
                  std::cout << "separateFrac: added to furtherSplitCandidates : "
                  << spVarPt->first->name() << " bd = "
                  << candCpBd.val() << " nbOfSatisfyingFracCol = " << nbOfSatisfyingFracCol
                  << " nbOfNotSatisfyingFracCol = "
                  << nbOfNotSatisfyingFracCol << std::endl;
          }
      } // else if (candBdIsWithinPreviousBound &&
    } // for index
  } // for spVarPt

  /// if selection strategy is highest priority, we need to check all possible branchings to
  /// select the highest priority one
  if ((!fullSearchForBestBranchingSet || (additionalNbOfCpBd >= 2))
      && (priorityRule() != SelectionStrategy::HighestPriority))
    {
      if (foundCPforSepInThisCall2separate)
        return;
    }

  /// Step e: Partition according to the component with highest branching priority
    if (keepOnlyBestVar4FurtherSplit && (bestVar4FurtherSplit != NULL))
        furtherSplitCandidates.push_back(FurtherSplitCandidate(bestVar4FurtherSplit, bestCpBd4FurtherSplit,
                                                               bestNbOfSatisfyingFracCol, bestNbOfNotSatisfyingFracCol));
  if (furtherSplitCandidates.empty())
  {
    if (printL(3))
      std::cout
            << "separateFrac: master solution within current column class is implicitly integer "
            << std::endl;

    /// Current solution is implicitly integer (not fractional)
    return;
  }

  for (std::list<FurtherSplitCandidate>::iterator candIt =
      furtherSplitCandidates.begin(); candIt != furtherSplitCandidates.end();
      ++candIt)
    {
      int nbFracColRequired(((candIt->varPtr->type() == 'I') ? 1 : 2));
      bool testFirstBranch(candIt->nbOfSatisfyingFracCol >= nbFracColRequired);
      bool testSecondBranch(candIt->nbOfNonSatisfyingFracCol >= nbFracColRequired);

      /**
       * This recursive call could be a dead end,
       * while a fractinal solution exists in another call;
       * hence, DO NOT require(testFirstBranch
       * || testSecondBranch || recordedBestBranchingSet.foundCbs(),
       * "separateFrac: one side of the partition at least should
       * be fit for further splitting");
       */
      if (testFirstBranch && testSecondBranch && (nbOfBranchToTest <= 1))
        {
          /// Go for larger fractional set
           if (2 * (int) fractColInvolvingSpVar[bestVar4FurtherSplit].size() >= (int) curListOfFractMastCol.size())
             testSecondBranch = false;
          else
            testFirstBranch = false;
        }

      /// Test set with bestCpBd4FurtherSplit
      if (testFirstBranch && (nbOfBranchToTest >= 1)) //(fractColInvolvingSpVar[bestVar4FurtherSplit].size() >= 2)
        {
          if (printL(localPrintLevel))
            std::cout
                << "separateFrac: TEST FIRST BRANCH with complementary component "
                << candIt->compBound << " additionalNbOfCpBd will be "
                << additionalNbOfCpBd + 1 << std::endl;

          /// Reset class of columns
          MasterColSolution newListOfFractMastCol;
          for (MasterColSolution::const_iterator colPt =
              curListOfFractMastCol.begin();
              colPt != curListOfFractMastCol.end(); colPt++)
            {
              /// Col satisfy cp bound
              if (candIt->compBound.satisfiedBy(colPt->first->spVarVal(candIt->compBound.varPtr())))
                {
                  if (printL(localPrintLevel))
                    std::cout << "separateFrac recursive call retains col "
                        << colPt->first->name() << std::endl;

                  newListOfFractMastCol.push_back(colPt->first, colPt->second);
                }
            }
          curClassCompBoundSet.push_back(candIt->compBound);
          curClassCompBoundSet.roundFracWeight();

          /// Binary case
          if (candIt->varPtr->curUb() == 1)
            {
              /// Remove it
              candVarForCompBoundBranching.erase(candIt->varPtr);

              separateFrac(newListOfFractMastCol, candVarForCompBoundBranching,
                           curClassCompBoundSet, totFracColVal, additionalNbOfCpBd + 1,
                           generatedBrConstrGeneratorSet);

              /// Reinsert it
              candVarForCompBoundBranching.insert(candIt->varPtr);
            }
          else
            {
              separateFrac(newListOfFractMastCol, candVarForCompBoundBranching,
                  curClassCompBoundSet, totFracColVal, additionalNbOfCpBd + 1,
                  //recordedBestBranchingSet,
                  generatedBrConstrGeneratorSet);
            }

          curClassCompBoundSet.pop_back();
          curClassCompBoundSet.roundFracWeight();
        }

      /**
       * Test set with complemented  bestCpBd4FurtherSplit
       * if its cardinality is at least 2
       */
      if (testSecondBranch && (nbOfBranchToTest >= 1))
        {
          candIt->compBound.complement();
          if (printL(localPrintLevel))
            std::cout
                << "separateFrac: TEST SECOND BRANCH with complementary component "
                << candIt->compBound << " additionalNbOfCpBd will be "
                << additionalNbOfCpBd + 1 << std::endl;

          /// Reset class of columns
          MasterColSolution newListOfFractMastCol;
          for (MasterColSolution::const_iterator colPt =
              curListOfFractMastCol.begin();
              colPt != curListOfFractMastCol.end(); colPt++)
            {
              if (candIt->compBound.satisfiedBy(colPt->first->spVarVal(candIt->compBound.varPtr())))
                {
                  if (printL(localPrintLevel))
                    std::cout << "separateFrac complementary  retains col "
                        << colPt->first->name() << std::endl;

                  newListOfFractMastCol.push_back(colPt->first, colPt->second);
                }
            }

          curClassCompBoundSet.push_back(candIt->compBound);
          curClassCompBoundSet.roundFracWeight();

          if (candIt->varPtr->ub() == 1)
            {
              /// Remove it
              candVarForCompBoundBranching.erase(candIt->varPtr);

              separateFrac(newListOfFractMastCol, candVarForCompBoundBranching,
                  curClassCompBoundSet, totFracColVal, additionalNbOfCpBd + 1,
                  generatedBrConstrGeneratorSet);

              /// Reinsert it
              candVarForCompBoundBranching.insert(candIt->varPtr);
            }
          else
            {
              separateFrac(newListOfFractMastCol, candVarForCompBoundBranching,
                  curClassCompBoundSet, totFracColVal, additionalNbOfCpBd + 1,
                  generatedBrConstrGeneratorSet);
            }

          curClassCompBoundSet.pop_back();
          curClassCompBoundSet.roundFracWeight();
        }
    }
  return;
}

bool CompBoundSetGenBranchConstr::genericMastColumnCount(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (printL(6))
    std::cout << "CompBoundSetGenBranchConstr::genericMastColumnCount : InstanciatedConstr "
              << icPtr->name() << std::endl;

  if (icPtr->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
  {
    CompSetInstMastBranchConstr * cpsIbcPtr = static_cast<CompSetInstMastBranchConstr *> (icPtr);
    if (cpsIbcPtr->compBoundSet().cgSpConfPtr() != colPtr->cgSpConfPtr())
      return false;
    return (cpsIbcPtr->CBsatisfiedBySol(colPtr->spSol()));
  }

  return false;
}

const LpCoef CompBoundSetGenBranchConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (printL(6))
    std::cout << "CompBoundSetGenBranchConstr::genericMastColumnCoef : InstanciatedConstr " << icPtr->name()
              << std::endl;

  if (icPtr->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
  {
    CompSetInstMastBranchConstr * cpsIbcPtr = static_cast<CompSetInstMastBranchConstr *> (icPtr);

    if ((cpsIbcPtr->compBoundSet().cgSpConfPtr() == colPtr->cgSpConfPtr())
        && (cpsIbcPtr->CBsatisfiedBySol(colPtr->spSol())))
      return LpCoef::UnitCoef;
    else
      return LpCoef::ZeroCoef;
  }

  return LpCoef::ZeroCoef;
}

bool CompBoundSetGenBranchConstr::genericCount(const InstanciatedConstr * const iconstrPtr,
                                               const InstanciatedVar * const ivarPtr) const
{
  return (false);
}

const LpCoef CompBoundSetGenBranchConstr::genericCoef(const InstanciatedConstr * const iconstrPtr,
                                                      const InstanciatedVar * const ivarPtr) const
{
  return LpCoef::ZeroCoef;
}

std::ostream& CompBoundSetGenBranchConstr::print(std::ostream& os) const
{
  return os << "CompBoundSetGenBranchConstr of genVar " << _genVarPtr->defaultName()
            << " with priorityLevel " << priorityLevel() << std::endl;
}

/*
 * Methods of class CompSetInstMastBranchConstr
 */

CompSetInstMastBranchConstr::CompSetInstMastBranchConstr(const ComponentSequence & cpList, 
							                             const IndexCell & id,
                                                         GenericConstr * genBrConstrPtr, 
							                             ProbConfig * probConfigPtr,
                                                         const std::string & name, 
							                             const Double & rhs,
							                             const char & kind) :
    InstMasterBranchingConstr(id,
			                  genBrConstrPtr,
			                  probConfigPtr,
			                  name,
			                  rhs,
                              (cpList.empty() && cpList.activeSense() == 'L' ? 'L' : 'G'),
			                  ' ',
                              (cpList.empty() ? 'I' : kind),
			                  'd'),
    _compBoundSet(cpList), _dirPredCSconstrPtr(NULL), _associatedPricingSPsolved(false),
    _pricingSpZetaVal(Bound::infDualBound(_modelPtr->objectiveSense())), _solPtr(NULL),
    _bestReducedCost(0, _modelPtr->objectiveSense())
{

  if ((probConfigPtr->param().mastInitMode().status() == MasterInitMode::localArtCol)
      || (probConfigPtr->param().StabilFuncOuterAngle() + probConfigPtr->param().StabilFuncInnerAngle() > 0))
    addLocalArtVar(probConfigPtr->modelPtr()->objectiveSense());

  return;
}

ColGenSpConf * CompSetInstMastBranchConstr::ColGenSpConfPtr() const
{
  return _compBoundSet.cgSpConfPtr();
}

void CompSetInstMastBranchConstr::setMembership()
{
  if (printL(5))
    std::cout << "CompSetInstMastBranchConstr::setMembership() brConstr " << name() << " of spConf "
              << (ColGenSpConfPtr() != NULL ? ColGenSpConfPtr()->name() : "undefined") << std::endl;

  if (!buildMembershipHasBeenPerformed())
  {
    genVarConstrPtr()->buildMembership(this);
    buildMembershipHasBeenPerformed(true);
  }

  ///  add Membership of Master columns
  bool noLowerCpBound(true);

  std::list<MastColumn *> mastColPtrConsideredForMembership;
  for (ComponentSequence::const_iterator cpPt = compBoundSet().begin(); cpPt != compBoundSet().end(); ++cpPt)
  {
    if (cpPt->sign() == 'G')
    {
      noLowerCpBound = false;

      if (!cpPt->varPtr()->isTypeOf(VcId::SubProbVariableMask))
        continue;

      SubProbVariable * spVarptr = static_cast<SubProbVariable *> (cpPt->varPtr());
      
      for (MapMastColumnPtr2Double::const_iterator mcPt = spVarptr->masterColumnMember2coefMap().begin();
           mcPt != spVarptr->masterColumnMember2coefMap().end(); ++mcPt)
      {
        if (cpPt->satisfiedBy(mcPt->second))
        {
          mastColPtrConsideredForMembership.push_back(mcPt->first);
        }
      }

    }
  }

  if (noLowerCpBound)
  {
    for (VarIndexManager::iterator varPtrIt = problemPtr()->probVarSet().begin(VcIndexStatus::Active, 'd');
            varPtrIt != problemPtr()->probVarSet().end(VcIndexStatus::Active, 'd'); varPtrIt++)
    {
      if((*varPtrIt)->isTypeOf(VcId::MastColumnMask) && 
              (*varPtrIt)->spSol()->probConfPtr() == compBoundSet().cgSpConfPtr())
      {
        mastColPtrConsideredForMembership.push_back(static_cast<MastColumn*>(*varPtrIt));
      }    
    }
    
    for (VarIndexManager::iterator varPtrIt = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
            varPtrIt != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); varPtrIt++)
    {
      if((*varPtrIt)->isTypeOf(VcId::MastColumnMask) && 
              (*varPtrIt)->spSol()->probConfPtr() == compBoundSet().cgSpConfPtr())
      {
        mastColPtrConsideredForMembership.push_back(static_cast<MastColumn*>(*varPtrIt));
      }    
    }
    
    for (VarIndexManager::iterator varPtrIt = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
            varPtrIt != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); varPtrIt++)
    {
      if((*varPtrIt)->isTypeOf(VcId::MastColumnMask) && 
              (*varPtrIt)->spSol()->probConfPtr() == compBoundSet().cgSpConfPtr())
      {
        mastColPtrConsideredForMembership.push_back(static_cast<MastColumn*>(*varPtrIt));
      }    
    }
  }
  bool cumulativeCoef(false);

  for (std::list<MastColumn *>::iterator colPt = mastColPtrConsideredForMembership.begin();
       colPt != mastColPtrConsideredForMembership.end(); ++colPt)
  {
    if (printL(6))
      std::cout << " considering mast colum " << (*colPt)->name() << std::endl;

    if (CBsatisfiedBySol((*colPt)->spSol()))
    {
      if (printL(6))
        std::cout << " include as member " << std::endl;
      includeMember(*colPt, 1, cumulativeCoef);
    }
  }

  Constraint::setMembership();

  return;
}

void CompSetInstMastBranchConstr::enumerativeSetMembership()
{
  if (presetMembership() == true)
  {
    std::cerr << "CompSetInstMastBranchConstr::enumerativeSetMembership(): should not have preset memebership"
              << std::endl;
    exit(1);
  }

  /**
   * The add other columns 
   * Even if presetMembership force constraint setmembership to include mast column 
   * Because MastColumn in Master constraints are not managed in buildMembership 
   */
  InstanciatedConstr::enumerativeSetMembership();

  return;
}

std::ostream & CompSetInstMastBranchConstr::print(std::ostream& os) const
{
  os << "CompSetInstMastBranchConstr" << std::endl;
  os << " name  = " << name() << std::endl;
  os << _compBoundSet;
  os << " rhs = " << curRhs() << std::endl;
  os << " margLvalue4DualBd = " << _margLvalue4DualBd << std::endl;
  os << " marginLvalue = " << _marginLvalue << std::endl;
  os << " depth = " << _setOfPredCSconstrPtr.size() << std::endl;

  /// Pointer to direct predecessor in tree of column classes
  if (_dirPredCSconstrPtr != NULL)
    os << " dirPredCSconstr = " << _dirPredCSconstrPtr->name() << std::endl;
  else
    os << " no dirPredCSconstr " << std::endl;

  os << " _associatedPricingSPsolved = " << _associatedPricingSPsolved
          << std::endl;
  os << " _sigma = " << _sigma << std::endl;
  InstMasterBranchingConstr::print(os);

  return (os);
}

std::vector<std::string> CompSetInstMastBranchConstr::forDotPrint() const
{
  list<string*> cmplines;
  string front = "";
  string back ="";
  
  front += (_compBoundSet._cgSpConfPtr != NULL ? _compBoundSet._cgSpConfPtr->name() : "undefined");
  
  string::size_type maxLineSize = 0;
  front += ":";
  
  if (!_compBoundSet.empty())
  {      
      int conditionCount = 0;      
      string cmpline = "";
      for (ComponentSequence::const_iterator it = _compBoundSet.begin(); it <= _compBoundSet.end(); ++it)
      {
        if (it != _compBoundSet.begin())
        {  
          if (++conditionCount == 2 || it == _compBoundSet.end()) 
          {            
            if(it != _compBoundSet.end()) 
            {
              cmpline += ",";
            }            
            cmplines.push_back(new string(cmpline));
            conditionCount = 0;
            maxLineSize = (std::max)(maxLineSize, cmpline.size());
            if(it != _compBoundSet.end()) 
            {
              cmpline = "";
            }
            else
            {
              break;
            }            
          }
          else 
          {  
            cmpline += ", ";
          }
        }   
        cmpline += it->_varPtr->name();
        switch (it->_sign) {
        case 'G':
          cmpline += " >= ";
          break;
        case 'L':
          cmpline += " <= ";
          break;
        }
        string sVal;
        stringstream strstream;
        strstream << (*it)._val;
        strstream >> sVal;
        cmpline += sVal;
      }
  }
  else
  {
    cmplines.push_back(new string("  "));
    maxLineSize = 2;
  }
  
  switch (_sense)
  {
    case 'G':
      back+= " >= ";
      break;
    case 'L':
      back+= " <= ";
      break;
    case 'E':
      back+= " == ";
      break;
    default:
      back+=" ?= ";
      break;
  }
  string sCostRhs;
  stringstream strstream;
  strstream << _costrhs;
  strstream >> sCostRhs;
  back += sCostRhs;
  
  vector<string> retLines;
  string leadingSpaces(front.size(), ' ');
  string tailingSpaces(back.size(), ' ');  
  for(list<string*>::iterator it = cmplines.begin() ; it!= cmplines.end() ; it++)
  {
    string exp;
    string midleSpaces(maxLineSize-(*it)->size(), ' ');
    if(it == cmplines.begin()) 
    {
      exp = front + '|' + (*(*it)) + midleSpaces + '|' + back;
    }
    else
    {
      exp = leadingSpaces + '|' + (*(*it)) + midleSpaces + '|' + tailingSpaces;
    }
    retLines.push_back(exp);
  }

  return retLines;
}

void CompSetInstMastBranchConstr::shortPrint(std::ostream& os) const
{
  os << "[ " << (_compBoundSet._cgSpConfPtr != NULL ? _compBoundSet._cgSpConfPtr->name() : "undefined");
  if (!_compBoundSet.empty())
    {
      os << " with ";
      for (ComponentSequence::const_iterator it = _compBoundSet.begin(); it != _compBoundSet.end(); ++it)
      {
        if (it != _compBoundSet.begin())
          os << ", ";
        os << it->_varPtr->name();
        switch (it->_sign) {
        case 'G':
          os << " >= ";
          break;
        case 'L':
          os << " <= ";
          break;
        }
        os << (*it)._val;
        //os << "(c=" << it->_cardinality << " cc=" << it->_complementCard << ")";
      }
  }
  os << " ]";
  switch (_sense)
  {
    case 'G':
      os << " >= ";
      break;
    case 'L':
      os << " <= ";
      break;
    case 'E':
      os << " == ";
      break;
    default:
      os << " ?= ";
      break;
  }
  os << _costrhs << " ";
}

/// Constraint enforces same class bounds that imlied by cpList

bool CompSetInstMastBranchConstr::enforces(const ComponentSequence & cpList) const
{
  if (cpList.size() != compBoundSet().size())
    return (false);

  ComponentSequence::const_iterator it1 = compBoundSet().begin();
  ComponentSequence::const_iterator it2 = cpList.begin();
  for (; it1 != compBoundSet().end(); ++it1, it2++)
    if (*it1 != *it2)
      return (false);

  return ((compBoundSet().classCardinality() == cpList.classCardinality()));
}

bool CompSetInstMastBranchConstr::incompatibleCsBd() const
{
  std::map<InstanciatedVar *, Double> spLb;
  std::map<InstanciatedVar *, Double> spUb;

  for (ComponentSequence::const_iterator it = compBoundSet().begin();
       it != compBoundSet().end(); ++it)
  {
    if (it->sign() == 'G')
    {
      if (!spLb.count(it->varPtr()))
        spLb[it->varPtr()] = it->val();
      else if (it->val() > spLb[it->varPtr()])
        spLb[it->varPtr()] = it->val();

      if (spUb.count(it->varPtr()))
        if (spLb[it->varPtr()] > spUb[it->varPtr()])
          return true;
    }

    if (it->sign() == 'L')
    {
      if (!spUb.count(it->varPtr()))
        spUb[it->varPtr()] = it->val();
      else if (it->val() < spUb[it->varPtr()])
        spUb[it->varPtr()] = it->val();
      if (spLb.count(it->varPtr()))
        if (spLb[it->varPtr()] > spUb[it->varPtr()])
          return true;
    }
  }

  return false;
}

bool CompSetInstMastBranchConstr::CBsatisfiedBySol(Solution * solPtr) const
{
  return compBoundSet().satisfiedBy(solPtr);
}

void CompSetInstMastBranchConstr::resetBounds(const Bound & dualBd,
                                              const Bound & redCost)
{
  _pricingSpZetaVal = dualBd;
  _bestReducedCost = redCost;

  return;
}

void CompSetInstMastBranchConstr::recSol(Solution * solPtr,
                                         const Bound & dualBd,
                                         const Bound & redCost)
{
  _solPtr = solPtr->clone(); /// local copy
  ///  sol of relaxed problem solves the more constrained problem
  if (_associatedPricingSPsolved == false)
  {
    _associatedPricingSPsolved = true;
    if (dirPredCSconstrPtr() != NULL)
    {

      dirPredCSconstrPtr()->_margLvalue4DualBd -= curRhs();
      if (printL(5))
        std::cout << " CompSetInstMastBranchConstr::recSol FOR " << name()  << std::endl
              << " whose rhs is " << curRhs() << std::endl
              << " the direct pred is " << dirPredCSconstrPtr()->name() << std::endl
              << " with rhs " << dirPredCSconstrPtr()->curRhs() << std::endl
              << " marLVal " << dirPredCSconstrPtr()->marginLvalue() << std::endl
              << " margLvalue4DualBd "
              << dirPredCSconstrPtr()->_margLvalue4DualBd << std::endl;
    }
  }

  resetBounds(dualBd, redCost);

  return;
}

void CompSetInstMastBranchConstr::reset()
{
    _marginLvalue = curRhs();
    _margLvalue4DualBd = curRhs();
    _setOfPredCSconstrPtr.clear();
    _dirPredCSconstrPtr = NULL;
    if (_solPtr != NULL)
        delete _solPtr;
    _solPtr = NULL;
    _associatedPricingSPsolved = false;
    _bestReducedCost = Bound(0, _modelPtr->objectiveSense());
    _sigma = 0;
    
    return;
}

void CompSetInstMastBranchConstr::recordInducedVarBounds()
{
  for (ComponentSequence::const_iterator cbPt = compBoundSet().begin(); cbPt != compBoundSet().end(); cbPt++)
  {
    if (printL(6))
      std::cout << name() << " recordInducedVarBounds() " << *cbPt;

    /// Set greatest value
    if (cbPt->sign() == 'G')
      cbPt->varPtr()->mapCompSetBrConstr2lowerBd()[this] = cbPt->val();
      ///  'L'
    else
      cbPt->varPtr()->mapCompSetBrConstr2upperBd()[this] = cbPt->val();
  }
}

bool operator==(const CompSetInstMastBranchConstr & a, const CompSetInstMastBranchConstr & b)
{
  if (a.compBoundSet().size() != b.compBoundSet().size())
    return (false);

  if (a.compBoundSet().activeSense() != b.compBoundSet().activeSense())
    return (false);

  if (a.compBoundSet().classCardinality() != b.compBoundSet().classCardinality())
    return (false);

  switch (compareCbS(a.compBoundSet(), b.compBoundSet())) {
  case ComponentSequence::identical:
    return (true);
  case ComponentSequence::different:
  case ComponentSequence::superclass:
  case ComponentSequence::subclass:
    ;
  }

  return (false);
}

bool CompSetInstMastBranchConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::CompSetInstMastBranchConstrMask, vcIdentifier);
}
