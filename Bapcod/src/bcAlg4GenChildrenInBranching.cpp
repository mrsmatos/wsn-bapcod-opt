/**
 *
 * This file bcAlg4GenChildrenInBranching.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4GenChildrenInBranching.hpp"
#include "bcAlg4EvalByColAndCutGen.hpp"
#include "bcAlg4EvalByLp.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcNetworkBasedBranching.hpp"
#include "bcNodeC.hpp"
#include "bcColGenSpConfC.hpp"


using namespace std;

const double CandidateBranchGroup::numberOfLeafs(const double & gap,
    const std::vector<double> & delta)
{
  double inf, mid, sup, exp;
  inf = 0.0;
  sup = 1e20;
  for (int it = 0; it < 100; it++)
    {
      mid = (inf + sup) / 2;
      if (sup - inf < sup / 1000000)
        break;
      exp = 0.0;
      for (std::vector<double>::const_iterator it = delta.begin();
          it != delta.end(); ++it)
        exp += pow(mid, -(*it) / gap);
      if (exp < 1.00)
        sup = mid;
      else
        inf = mid;
      if (mid > 0.999e20)
        return -1;
    }
  return mid;
}

CandidateBranchGroup::CandidateBranchGroup(const int & id_,
                                           BranchingConstrGenerator * generatorPtr,
                                           const BcObjStatus::MinMaxIntFloat & objStatus_,
                                           const bool & historyCandidate_) :
  objStatus(objStatus_), id(id_), productBranchingScore(-BapcodInfinity),
  linearBranchingScore(-BapcodInfinity), treeDepthBranchingScore(-BapcodInfinity), treeSize(BapcodInfinity),
  treeDepth(BapcodInfinity), genPtr(generatorPtr), pessimisticDualBound(Bound::infPrimalBound(objStatus)),
  optimisticDualBound(Bound::infDualBound(objStatus)), historyCandidate(historyCandidate_)
{
}

CandidateBranchGroup::~CandidateBranchGroup()
{
  delete genPtr;
  for (CandidateBranchGroup::iterator nodePtrIt = begin();
       nodePtrIt != end(); ++nodePtrIt)
    {
      delete *nodePtrIt; //notice that deleting a node also removes associations.
      *nodePtrIt = NULL;
    }
}

void CandidateBranchGroup
     ::computeBranchingScoreAsProduct(const Bound & parentFractionalDualBound,
                                      const Double & cutOffValue,
                                      const int & maxDescriptionLength,
                                      const int & currPhaseNumber,
                                      const long & elapsedTime,
                                      const int & candidatesEvaluated,
                                      double bestProductScores[3])
{
  std::stringstream sstream;
  sstream << std::fixed << std::setprecision(4);
  genPtr->nicePrint(sstream);
  int lengthDiff = maxDescriptionLength - sstream.str().length();

  double currentGap = cutOffValue - parentFractionalDualBound;
  productBranchingScore = 1;
  bool allBranchesAreAboveTheGap = true;
  std::vector<double> nodePrimalBounds;
  for (CandidateBranchGroup::const_iterator nodeIt = begin(); nodeIt != end(); ++nodeIt)
  {
    double curDelta = (*nodeIt)->nodeIncLpPrimalBound() - parentFractionalDualBound;
    if (curDelta < currentGap)
      allBranchesAreAboveTheGap = false;
    nodePrimalBounds.push_back((*nodeIt)->nodeIncLpPrimalBound());
  }
  /// we sort deltas in the non-decreasing order
  std::stable_sort(nodePrimalBounds.begin(), nodePrimalBounds.end());
  for (int boundId = 0; boundId < nodePrimalBounds.size(); ++boundId)
  {
    Bound curPrimalBound(nodePrimalBounds[boundId], objStatus);

    if (pessimisticDualBound > curPrimalBound)
      pessimisticDualBound = curPrimalBound;
    if (optimisticDualBound < curPrimalBound)
      optimisticDualBound = curPrimalBound;
    double curDelta = curPrimalBound - parentFractionalDualBound;
    if ((curDelta > currentGap) && (!allBranchesAreAboveTheGap || (boundId >= 2)))
      curDelta = currentGap + Double::precision * (curDelta - currentGap);
    if (objStatus == BcObjStatus::maxInt || objStatus == BcObjStatus::maxFloat)
      curDelta *= -1;
    if (curDelta < Double::precision)
      curDelta = Double::precision;
    if (boundId < 2)
      productBranchingScore *= curDelta;
    else
      productBranchingScore *= (curDelta / currentGap);
  }
  if (nodePrimalBounds.empty())
    productBranchingScore = currentGap * currentGap;
  else if (nodePrimalBounds.size() == 1)
    productBranchingScore *= currentGap;

  if (printL(-1))
  {
    bool showOutput = printL(0) || (currPhaseNumber > 1);
    if (!showOutput && printL(-1)) {
      /// we first verify whether the current score is one of the best 3
      if (productBranchingScore > bestProductScores[0])
      {
        showOutput = true;
        bestProductScores[2] = bestProductScores[1];
        bestProductScores[1] = bestProductScores[0];
        bestProductScores[0] = productBranchingScore;
      } else if (productBranchingScore > bestProductScores[1])
      {
        showOutput = true;
        bestProductScores[2] = bestProductScores[1];
        bestProductScores[1] = productBranchingScore;
      } else if (productBranchingScore > bestProductScores[2])
      {
        showOutput = true;
        bestProductScores[2] = productBranchingScore;
      }
    }
    if (showOutput)
    {
      std::cout << "SB phase " << currPhaseNumber << " cand. " << setw(2) << candidatesEvaluated << " branch on ";
      std::cout << std::fixed << std::setprecision(4);
      genPtr->nicePrint();
      std::cout << setw(lengthDiff) << "";
      std::cout << ": [";
      if (!nodePrimalBounds.empty())
      {
        std::cout << std::setw(10) << nodePrimalBounds[0];
        for (int boundId = 1; boundId < (int) nodePrimalBounds.size(); ++boundId)
          std::cout << ", " << std::setw(10) << nodePrimalBounds[boundId];
      }

      std::cout << std::setprecision(2) << std::setw(5) << "], score = " << productBranchingScore
                << ((historyCandidate) ? " (h)" : "    ");
      std::cout << std::setprecision(6);
      std::cout.unsetf(std::ios_base::floatfield);
      std::cout << "  <et=" << elapsedTime / (double) 100 << ">" << std::endl;
    }
  }

  if (printL(5))
    std::cout << "Node::computeBranchingScoreAsProduct() productBranchingScore = " << productBranchingScore
              << std::endl;
}

void CandidateBranchGroup
     ::computeBranchingTreeDepthScore(const Bound & parentFractionalDualBound,
                                      const Double & cutOffValue,
                                      const int & numEvaluatedChildren,
                                      const double notPromissingConservativeness)
{
  treeDepthBranchingScore = -BapcodInfinity;
  treeDepth = treeSize = BapcodInfinity;

  double currentGap = cutOffValue - parentFractionalDualBound;

  std::vector<double> delta;
  double maxRawDelta = 0;
  int nbZeroDeltas = 0;
  double sumDeltas = 0.0;
  for (CandidateBranchGroup::const_iterator nodeIt = begin(); nodeIt != end(); ++nodeIt)
    {
      double curDelta = (*nodeIt)->nodeIncLpPrimalBound() - parentFractionalDualBound;
      sumDeltas += curDelta;
      if (curDelta < 1e-6)
        nbZeroDeltas += 1;
      if (maxRawDelta < curDelta)
        maxRawDelta = curDelta;
      if (delta.size() >= numEvaluatedChildren)
        curDelta *= notPromissingConservativeness;
      if (curDelta < currentGap)
        delta.push_back(curDelta);
      else
        delta.push_back(currentGap);
    }

  if (nbZeroDeltas == size())
    return;

  if (nbZeroDeltas > 0)
  {
    for (std::vector<double>::iterator deltaIt = delta.begin(); deltaIt != delta.end(); ++deltaIt)
      if (*deltaIt < currentGap * 0.0001)
        *deltaIt = currentGap * 0.0001;
  }

  int deltaSize = delta.size();
  if (deltaSize == 0)
    {
      treeDepthBranchingScore = 0;
      treeSize = treeDepth = 0;
      return;
    }

  if (deltaSize == 1)
    {
      treeDepthBranchingScore = -treeDepthScoreNormConstant / delta[0];
      treeDepth = currentGap / delta[0];
      treeSize = treeDepth;
      return;
    }

  double numNormLeafs = numberOfLeafs(treeDepthScoreNormConstant, delta);

  if (numNormLeafs < 0)
    treeDepthBranchingScore = -BapcodInfinity;
  else
    treeDepthBranchingScore = -log(numNormLeafs) / log((double) deltaSize);

  double numRealLeafs = numberOfLeafs(currentGap, delta);
  if (numRealLeafs > 0)
    {
      treeSize = deltaSize * ((numRealLeafs - 1) / (deltaSize - 1));
      treeDepth = log(numRealLeafs) / log((double) deltaSize);
    }
}

/// for the case when master is solved by column generation only
void Alg4GenChildrenInBranching
     ::prepareNodeForTreatmentInStrBranchWithPhases(Node * nodePtr, int globalTreatOrder,
                                                    const StrongBranchingPhaseParameter & currPhaseParam)
{
  bool doFullSetup = (_firstNodeWasEvaluated && !_lastEvaluatedNodeWasLp)
                     || (!_firstNodeWasEvaluated && (nodePtr->probSetupInfoPtr()->treatOrderId != globalTreatOrder))
                     || nodePtr->probSetupInfoPtr()->fullSetupIsObligatory;
  if (doFullSetup)
    {
      nodePtr->setProblemSetupAlgorithm(new Alg4ProblemSetupFull(_masterCommons.masterCommons4ProblemSetup()));
    }
  else
    {
      Alg4ProblemSetupBranchingOnly* problemBranchingOnlySetupAlgPtr
        = new Alg4ProblemSetupBranchingOnly(_masterCommons.masterCommons4ProblemSetup());
      bool doSubproblemsSetup = (currPhaseParam.maxNumOfColGenIterations() > 0);
      problemBranchingOnlySetupAlgPtr->setOptionDoSubproblemsSetup(doSubproblemsSetup);
      nodePtr->setProblemSetupAlgorithm(problemBranchingOnlySetupAlgPtr);
    }

  if (currPhaseParam.maxNumOfColGenIterations() == 0)
    {
      /// solving only restricted master
      Alg4EvalByLp * evalAlgPtr = new Alg4EvalByLp(_masterCommons.masterCommons4ProblemSetup().problemList().front(),
                                                   _masterCommons.masterCommons4EvalAlg());
      /// dual bound should not be updated, as column generation is not performed
      evalAlgPtr->setOptionUpdateIncDualBound(false);
      nodePtr->setEvalAlg(evalAlgPtr);
      /// if column generation is not used, the problem cannot be changed,
      /// so we do not save the problem setup information
      nodePtr->setProblemSetDownAlgorithm(new Alg4ProblemSetDownOfNode(_masterCommons.masterCommons4ProblemSetup()));
    }
  else /// using column generation
    {
      /// preprocessing is used only if column generation is used
      if (param().ApplyPreprocessing())
        nodePtr->setPreprocessor(new Algorithm4PreprocessingAtNodeOtherThanRoot(
            _masterCommons.masterCommons4ProblemSetup().problemList()));

      Alg4EvalByColAndCutGen * evalAlgPtr
        = new Alg4EvalByColAndCutGen(_masterCommons.masterCommons4ProblemSetup().problemList().front(),
                                     _masterCommons.masterCommons4EvalAlg());
      evalAlgPtr->setOptionNeedBeConqueredIfSolIsInteger(currPhaseParam.exact());
      evalAlgPtr->setOptionMaxNbOfCgIterations(currPhaseParam.maxNumOfColGenIterations());
      evalAlgPtr->setOptionMaxNbOfCutRounds(currPhaseParam.maxNumCutRounds());
      if (currPhaseParam.exact())
      {
        evalAlgPtr->setOptionMinNbOfCutRounds(param().MinNumOfCutRoundsBeforeStopBySp());
      }
      else
      {
        evalAlgPtr->setOptionNonExactEvaluation(true);
        evalAlgPtr->setOptionMinNbOfCutRounds(currPhaseParam.minNumCutRounds());
      }
      if (currPhaseParam.minLevelOfSpRestriction() < param().MaxNbOfStagesInColGenProcedure())
        evalAlgPtr->setOptionMinLevelOfSpRestriction(currPhaseParam.minLevelOfSpRestriction());
      else
        evalAlgPtr->setOptionMinLevelOfSpRestriction(param().MaxNbOfStagesInColGenProcedure() - 1);
      if (currPhaseParam.doRedCostFixingAndEnumeration())
        evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(1);
      else
          evalAlgPtr->setOptionDoRedCostFixingAndEnumeration(0);
      evalAlgPtr->setOptionMaxNbOfPenaltyUpdates(param().ArtVarMaxNbOfPenaltyUpdates());
      if (currPhaseParam.exact())
      {
        int logFrequency = param().ColGenLogFrequency();
        if (!printL(0) && (logFrequency < 10))
          logFrequency = 10;
        evalAlgPtr->setOptionLogPrintFrequency(logFrequency);
      }
      else
      {
        evalAlgPtr->setOptionLogPrintFrequency(currPhaseParam.logPrintFrequency());
      }
      nodePtr->setEvalAlg(evalAlgPtr);
      nodePtr->setProblemSetDownAlgorithm(
        new ProblemFullSetDownAlgorithm(_masterCommons.masterCommons4ProblemSetup()));
    }
  _firstNodeWasEvaluated = true;
  _lastEvaluatedNodeWasLp = (currPhaseParam.maxNumOfColGenIterations() == 0);
}

/// these two structures are needed for the function Alg4GenChildrenInBranching::getConstrGeneratorsFromHistory()

struct GeneratorCurrentInfo
{
  BranchingConstrGenerator * generatorPtr;
  int maxPhaseEvaluated;
  double leftImprovement;
  double rightImprovement;
  double score;

  GeneratorCurrentInfo(BranchingConstrGenerator * generatorPtr_, const int & maxPhaseEvaluated_,
                       const double & leftImprovement_, double rightImprovement_) :
    generatorPtr(generatorPtr_), maxPhaseEvaluated(maxPhaseEvaluated_), leftImprovement(leftImprovement_),
    rightImprovement(rightImprovement_), score(rightImprovement_ * leftImprovement_)
  {
  }
};

struct GeneratorCurrentInfoComp
{
  bool operator()(const GeneratorCurrentInfo & infoA, const GeneratorCurrentInfo & infoB) const
  {
    if (infoA.maxPhaseEvaluated != infoB.maxPhaseEvaluated)
      return (infoA.maxPhaseEvaluated > infoB.maxPhaseEvaluated);
    if (infoA.score != infoB.score)
      return (infoA.score > infoB.score);
    return (*(infoA.generatorPtr) < *(infoB.generatorPtr));
  }
};

/// Ruslan : note that for the moment is valid only for the case with two child nodes
/// (the function needs to be adapted for the component set branching)
void Alg4GenChildrenInBranching
     ::getConstrGeneratorsFromHistory(std::vector<BranchingConstrGenerator *> & brConstrGeneratorsFromHistory,
                                      const int & maxNbOfCandidates, const Double & minPriorityLevel)
{
  Time startTime;

  std::set<GeneratorCurrentInfo, GeneratorCurrentInfoComp> generatorsSet;

  std::set<GenericBranchingConstr * , DynamicGenConstrSort>::const_iterator brGenConstrPtrIt;
  for (brGenConstrPtrIt = _masterCommons.candidateBranchingGenericConstr().begin();
       brGenConstrPtrIt != _masterCommons.candidateBranchingGenericConstr().end(); ++brGenConstrPtrIt)
    {
      GenericBranchingConstr * genBrConstrPtr = static_cast<GenericBranchingConstr *>(*brGenConstrPtrIt);
      if (genBrConstrPtr->priorityLevel() < minPriorityLevel)
        continue;

      for (BranchingGeneratorHistoryMap::iterator mapIt = genBrConstrPtr->branchingHistory().begin();
           mapIt != genBrConstrPtr->branchingHistory().end(); ++mapIt)
        {
          mapIt->first->computeLhs(_currentNodePtr->primalSol());
          if (mapIt->first->lFracPart().isZero())
            continue;
          double lFracPart = mapIt->first->lFracPart();
          int maxPhaseEvaluated = mapIt->second.numEvaluations.size();
          /// we search for the maximum phase for which all child nodes have beed evaluated at least once
          while ((maxPhaseEvaluated > 0) && !mapIt->second.numEvaluations[maxPhaseEvaluated - 1][0]
                 && !mapIt->second.numEvaluations[maxPhaseEvaluated - 1][1])
            maxPhaseEvaluated -= 1;
          if (maxPhaseEvaluated == 0)
            continue;
          std::vector<double> & pseudoCosts = mapIt->second.pseudoCosts[maxPhaseEvaluated - 1];
          double absoluteGap = _currentNodePtr->nodeIncIpPrimalBound() - _currentNodePtr->nodeIncLpDualBound();
          double leftImprovement = (std::min)(pseudoCosts[0] * ((mapIt->first->direction() == 'U')
                                                              ? (1 - lFracPart) : lFracPart), absoluteGap);
          double rightImprovement = (std::min)(pseudoCosts[1] * ((mapIt->first->direction() == 'U')
                                                               ? lFracPart : (1 - lFracPart)), absoluteGap);
          generatorsSet.insert(GeneratorCurrentInfo(mapIt->first, maxPhaseEvaluated,
                                                    leftImprovement, rightImprovement));
        }
    }

  int counter = 0;
  for (std::set<GeneratorCurrentInfo, GeneratorCurrentInfoComp>::iterator genInfoIt = generatorsSet.begin();
       genInfoIt != generatorsSet.end(); ++genInfoIt)
    {
      if (++counter > maxNbOfCandidates)
        break;
      brConstrGeneratorsFromHistory.push_back(genInfoIt->generatorPtr->clone());
    }

  bapcodInit().statistics().incrTimer("bcTimeSepFracSol", startTime.getElapsedTime_dbl());
}

void Alg4GenChildrenInBranching
     ::strongBranchingWithPhases(const double & curPriorityLevel, int & globalTreatOrder,
                                 std::vector<BranchingConstrGenerator *> & branchGeneratorsFromHistory,
                                 std::vector<BranchingConstrGenerator *> & newBranchGenerators)
{
  Time startTime;

  std::vector<BranchingConstrGenerator *>::iterator histBrGenPtrIt, newBrGenPtrIt;
  std::vector<CandidateBranchGroup *> candidateBranchingGroups;
  int candidatesEvaluatedOnThisPhase = 0;
  histBrGenPtrIt = branchGeneratorsFromHistory.begin();
  newBrGenPtrIt = newBranchGenerators.begin();
  while ((++candidatesEvaluatedOnThisPhase <= param().StrongBranchingPhaseOne().maxNumOfCandidates())
         && ( (newBrGenPtrIt != newBranchGenerators.end()) || (histBrGenPtrIt != branchGeneratorsFromHistory.end()) ) )
    {
      /// we will now alternate between candidates from the history and new candidates
      bool historyCandidate = (newBrGenPtrIt == newBranchGenerators.end())
                              || ((histBrGenPtrIt != branchGeneratorsFromHistory.end())
                                  && (candidatesEvaluatedOnThisPhase % 2));
      BranchingConstrGenerator * generatorPtr = ((historyCandidate) ? *histBrGenPtrIt : *newBrGenPtrIt);
      if (historyCandidate)
        ++histBrGenPtrIt;
      else
        ++newBrGenPtrIt;

      BranchingGeneratorHistory * generatorHistoryPtr = NULL;

      /// for the moment history is deactivated for component set branching constraints
      CompBdSetBranchConstrGenerator * compBdSetBranchConstrGenPtr
        = dynamic_cast<CompBdSetBranchConstrGenerator *>(generatorPtr);
      if (compBdSetBranchConstrGenPtr == NULL)
        {
          /// we first retrieve the branching history associated to this branching generator
          BranchingGeneratorHistoryMap & branchGenHistoryMap
            = generatorPtr->genericBrConstrPtr()->branchingHistory();
          BranchingGeneratorHistoryMap::iterator generatorHistIt = branchGenHistoryMap.find(generatorPtr);
          if (generatorHistIt != branchGenHistoryMap.end())
            generatorHistoryPtr = &(generatorHistIt->second);
          else
            {
              std::pair<BranchingGeneratorHistoryMap::iterator, bool> ret
                = branchGenHistoryMap.insert( std::make_pair(generatorPtr->clone(), BranchingGeneratorHistory()));
              generatorHistoryPtr = &(ret.first->second);
            }
          /// we add an entry (branching evaluation info) to this history
          /// pointer to this evaluation info will be passed to nodes
          generatorHistoryPtr->evaluationsInfo.push_back(new BranchingEvaluationInfo(_currentNodePtr->depth(),
                                                                                     _currentNodePtr->nodeIncLpDualBound(),
                                                                                     generatorPtr->lFracPart(),
                                                                                     generatorPtr->direction(), generatorHistoryPtr));
        }
      std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
      ConstrPtrSet dummyBranchingConstrSet;
      CandidateBranchGroup * newCandidatePtr = new CandidateBranchGroup(candidatesEvaluatedOnThisPhase, generatorPtr,
                                                                        _masterCommons.objStatus(), historyCandidate);
      int childNodeCounter = 0;
      while (generatorPtr->nextNodeBrConstr(_currentNodePtr,
                                            tmpLocalNodeBrConstrList, dummyBranchingConstrSet))
        {
          Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                            tmpLocalNodeBrConstrList);
          if (generatorHistoryPtr != NULL)
            newChildNodePtr->setBranchEvaluationInfo(generatorHistoryPtr->evaluationsInfo.back(),
                                                     ++childNodeCounter);
          std::list<BranchingConstrBaseType *>::const_iterator brConstrPtrIt;
          for (brConstrPtrIt = tmpLocalNodeBrConstrList.begin();
               brConstrPtrIt != tmpLocalNodeBrConstrList.end(); ++brConstrPtrIt)
            {
              if ((*brConstrPtrIt)->isTypeOf(VcId::InstMasterBranchingConstrMask))
                static_cast<InstMasterBranchingConstr *>(*brConstrPtrIt)->treatOrderId(_currentNodePtr->treatOrder());
              (*brConstrPtrIt)->depthWhenGenerated(newChildNodePtr->depth());
              (*brConstrPtrIt)->append2name(std::string() + "n" + newChildNodePtr->ref());
            }
          newCandidatePtr->push_back(newChildNodePtr);
        }
      candidateBranchingGroups.push_back(newCandidatePtr);
    }
  /// we delete remaining generators
  while (histBrGenPtrIt != branchGeneratorsFromHistory.end())
    {
      delete *histBrGenPtrIt;
      *histBrGenPtrIt = NULL;
    }
  branchGeneratorsFromHistory.clear();
  while (newBrGenPtrIt != newBranchGenerators.end())
    {
      delete *newBrGenPtrIt;
      newBranchGenerators.erase(newBrGenPtrIt++);
    }

  int maxPhaseNumber = ((param().StrongBranchingPhaseTwo().active())
                        ? ((param().StrongBranchingPhaseThree().active())
                           ? ((param().StrongBranchingPhaseFour().active()) ? 4 : 3) : 2) : 1);

  _firstNodeWasEvaluated = false;
  _lastEvaluatedNodeWasLp = false;
  Double bestTreeSize = BapcodInfinity;
  StrongBranchingPhaseParameter prevPhaseParam("");
  for (int currentPhaseNumber = 1; currentPhaseNumber <= maxPhaseNumber; ++currentPhaseNumber)
    {
      /// we get the parameters for the current phase
      StrongBranchingPhaseParameter currPhaseParam("");
      int numCandForNextPhase = 1;
      switch (currentPhaseNumber) {
        case 1:
          currPhaseParam = param().StrongBranchingPhaseOne();
          if (param().StrongBranchingPhaseTwo().active())
            numCandForNextPhase = param().StrongBranchingPhaseTwo().maxNumOfCandidates();
          break;
        case 2:
          currPhaseParam = param().StrongBranchingPhaseTwo();
          if (param().StrongBranchingPhaseThree().active())
            numCandForNextPhase = param().StrongBranchingPhaseThree().maxNumOfCandidates();
          break;
        case 3:
          currPhaseParam = param().StrongBranchingPhaseThree();
          if (param().StrongBranchingPhaseFour().active())
            numCandForNextPhase = param().StrongBranchingPhaseFour().maxNumOfCandidates();
          break;
        case 4:
          currPhaseParam = param().StrongBranchingPhaseFour();
          break;
        default:
          break;
        }

      /// we perform exact evaluation phase even if the number of candidates is one
      /// this is needed for the tree size estimation which is later used for
      /// stopping prematurely a strong branching phase (using treeSizeRatioToStop parameter)
      if ((candidateBranchingGroups.size() <= numCandForNextPhase) && !currPhaseParam.exact())
          continue;

      if (currentPhaseNumber == 2)
        startTime = Time();

      if (printL(-1))
        std::cout << "**** Strong branching phase "<< currentPhaseNumber << " is started *****" << std::endl;

      /// for nice printing, we compute the maximum description length
      int maxDescriptionLength = 0;
      std::vector<CandidateBranchGroup *>::iterator brGrPtrIt;
      for (brGrPtrIt = candidateBranchingGroups.begin(); brGrPtrIt != candidateBranchingGroups.end(); ++brGrPtrIt)
        {
          std::stringstream sstream;
          sstream << std::fixed << std::setprecision(4);
          (*brGrPtrIt)->genPtr->nicePrint(sstream);
          int length = sstream.str().length();
          if (maxDescriptionLength < length)
            maxDescriptionLength = length;
        }

      Double bestTreeSizeScore(-BapcodInfinity);
      bool aCandidateIsConquered = false;
      CandidateBranchGroup * bestCandPtr = NULL;
      candidatesEvaluatedOnThisPhase = 0;
      bool exactPhaseStopped = false;
      double bestProductScores[3] = {0.0, 0.0, 0.0}; /// to show only phase one candidate only if its score is one
                                                     /// on the three best (for phase 1 and for printlevel -1)
      for (brGrPtrIt = candidateBranchingGroups.begin();
           brGrPtrIt != candidateBranchingGroups.end(); ++brGrPtrIt)
        {
          if (!progStatus().doRun())
            break;
          Double treeSizeToTest = _currentNodePtr->estimatedSubtreeSize();
          if (treeSizeToTest > bestTreeSize)
            treeSizeToTest = bestTreeSize;
          if ( ( (currPhaseParam.treeSizeRatioToStop() > 0.0)
                 && (treeSizeToTest * currPhaseParam.treeSizeRatioToStop() < candidatesEvaluatedOnThisPhase) )
               || exactPhaseStopped )
            {
              (*brGrPtrIt)->resetAllScores();
              continue;
            }
          candidatesEvaluatedOnThisPhase += 1;

          if (currPhaseParam.maxNumOfColGenIterations() > 0)
            std::stable_sort((*brGrPtrIt)->begin(), (*brGrPtrIt)->end(), BestNodeLpValue());

          int numEvaluatedChildren = 0;
          CandidateBranchGroup::iterator nodePtrIt = (*brGrPtrIt)->begin();
          bool allNodesAreConquered = true;
          while (nodePtrIt != (*brGrPtrIt)->end())
            {
              if ((currPhaseParam.maxNumOfColGenIterations() > 0)
                  && ((currPhaseParam.logPrintFrequency() > 0) || currPhaseParam.exact()))
                {
                  if (printL(-1))
                  {
                    std::cout << "**** SB phase " << currentPhaseNumber << " evaluation of candidate "
                              << candidatesEvaluatedOnThisPhase << ", branch " << numEvaluatedChildren + 1 << " ( ";
                    for (std::list<BranchingConstrBaseType *>::const_iterator brConstrPtrIt =
                            (*nodePtrIt)->localNodeBrConstrList().begin();
                         brConstrPtrIt != (*nodePtrIt)->localNodeBrConstrList().end(); ++brConstrPtrIt)
                    {
                      if (brConstrPtrIt != (*nodePtrIt)->localNodeBrConstrList().begin())
                        std::cout << ", ";
                      (*brConstrPtrIt)->shortPrint();
                    }
                    std::cout << "), value = " << (*nodePtrIt)->nodeIncLpPrimalBound();
                    if ((*nodePtrIt)->debugSolutionAtThisNode())
                        std::cout << ", debug solution is here";
                    std::cout << std::endl;
                  }
                }

              if ((*nodePtrIt)->isConquered())
                {
                  if (printL(0))
                    std::cout << "Branch is already conquered!" << std::endl;
                  ++nodePtrIt;
                  ++numEvaluatedChildren;
                  continue;
                }

              if ((bestCandPtr != NULL) && currPhaseParam.exact()
                  && _masterCommons.candidateCutGenericConstr().empty())
                {
                  /// discard if not promising (only done if there are no cuts)
                  (*brGrPtrIt)->computeBranchingTreeDepthScore(_currentNodePtr->nodeIncLpDualBound(),
                                _currentNodePtr->probConfPtr()->cutOffValue(), numEvaluatedChildren,
                                param().StrongBranchingParameter4NotPromissingConservativeness());

                  if ((bestTreeSizeScore > -BapcodInfinity)
                      && (*brGrPtrIt)->treeDepthBranchingScore <= bestTreeSizeScore)
                    break;
                }

              /// if this is not phase one, we need to regenerate the node using problem setup information
              /// and evaluation information of the previous stage node
              if (prevPhaseParam.active())
                {
                  Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                                    (*nodePtrIt)->localNodeBrConstrList());
                  if (prevPhaseParam.maxNumOfColGenIterations() > 0)
                    {
                      newChildNodePtr->removeProblemSetupInfoAssociation();
                      newChildNodePtr->associateProblemSetupInfoPtr((*nodePtrIt)->probSetupInfoPtr());
                    }
                  /// we need to take the evaluation algorithm information from the previous phase,
                  /// as the simplex basis was changed there
                  newChildNodePtr->removeNodeEvalInfoAssociation();
                  newChildNodePtr->associateNodeEvalInfoPtr((*nodePtrIt)->nodeEvalInfoPtr());
                  newChildNodePtr->setBranchEvaluationInfoFromNode(*nodePtrIt);
                  delete *nodePtrIt;
                  *nodePtrIt = newChildNodePtr;
                }

              prepareNodeForTreatmentInStrBranchWithPhases(*nodePtrIt, globalTreatOrder, currPhaseParam);
              (*nodePtrIt)->evaluation(globalTreatOrder, _currentNodePtr->nodeIncIpPrimalBound());

              /// update the incumbent primal bound if a better one was found while solving this node
              if ((*nodePtrIt)->primalBoundIsUpdated())
                _currentNodePtr->updateNodeIncPrimalSolution((*nodePtrIt)->nodeIncIpPrimalSolPtr());

              if (!(*nodePtrIt)->isConquered())
                allNodesAreConquered = false;

              ++nodePtrIt;

              ++numEvaluatedChildren;
            }

          /// component set branching may produce no children as all of them may be redundant
          /// in this case, _currentNodePtr becomes conquered
          if ((*brGrPtrIt)->empty())
            {
              if (printL(-1))
                std::cout << " SB phase " << currentPhaseNumber << " candidate "
                          << candidatesEvaluatedOnThisPhase << " has no children !" << std::endl;
              _currentNodePtr->setDualBoundEqualToIncPrimalBound();
            }

          if (numEvaluatedChildren < (*brGrPtrIt)->size())
            {
              if (printL(-1))
                std::cout << " SB phase " << currentPhaseNumber << " evaluation of candidate "
                          << candidatesEvaluatedOnThisPhase << " is interrupted after "
                          << numEvaluatedChildren << " children" << std::endl;
              (*brGrPtrIt)->resetAllScores();
            }
          else if (allNodesAreConquered)
            {
              if (printL(-1))
                std::cout << " SB phase " << currentPhaseNumber << " candidate "
                          << candidatesEvaluatedOnThisPhase << " is conquered !" << std::endl;
              aCandidateIsConquered = true;
              bestCandPtr = *brGrPtrIt;
              break;
            }
          else if (currPhaseParam.exact())
            {
              (*brGrPtrIt)->computeBranchingTreeDepthScore(_currentNodePtr->nodeIncLpDualBound(),
                                                           _currentNodePtr->probConfPtr()->cutOffValue(),
                                                           numEvaluatedChildren, 1.0);
              if (printL(-1))
              {
                std::cout << "SB exact phase " << currentPhaseNumber << " branch on ";
                (*brGrPtrIt)->genPtr->nicePrint();
                std::cout << " : [";
                for (CandidateBranchGroup::const_iterator nodeIt = (*brGrPtrIt)->begin();
                     nodeIt != (*brGrPtrIt)->end(); ++nodeIt) {
                  if (nodeIt != (*brGrPtrIt)->begin())
                    std::cout << ", ";
                  std::cout << std::setprecision(8) << (*nodeIt)->nodeIncLpPrimalBound() ;
                }
                std::cout << std::setprecision(4) << "], tree depth = " << (*brGrPtrIt)->treeDepth
                          << ", tree size = " << (*brGrPtrIt)->treeSize << ", score = "
                          << (*brGrPtrIt)->treeDepthBranchingScore << (((*brGrPtrIt)->historyCandidate) ? " (h)" : "")
                          << std::setprecision(6) << std::endl;
              }
              if (bestTreeSizeScore < (*brGrPtrIt)->treeDepthBranchingScore)
                {
                  bestCandPtr = *brGrPtrIt;
                  bestTreeSizeScore = (*brGrPtrIt)->treeDepthBranchingScore;
                  bestTreeSize = (*brGrPtrIt)->treeSize;
                }
              /// make a parameter for stopping SB exact branching phase?
              if ((_currentNodePtr->estimatedSubtreeSize() < BapcodInfinity)
                  && ((*brGrPtrIt)->treeSize < _currentNodePtr->estimatedSubtreeSize())
                  && param().StrongBranchingExactPhaseInterrupt())
                {
                  if (printL(0))
                    std::cout << "SB exact phase " << currentPhaseNumber << " is stopped as the estimated tree size "
                              << " of the current candidate is smaller that one of the father ("
                              << _currentNodePtr->estimatedSubtreeSize() << ")" << std::endl;
                  exactPhaseStopped = true;
                }
            }
          else /// current phase is not exact, candidates are ranged using product score
            {
              (*brGrPtrIt)->computeBranchingScoreAsProduct(_currentNodePtr->nodeIncLpDualBound(),
                                                           _currentNodePtr->probConfPtr()->cutOffValue(),
                                                           maxDescriptionLength, currentPhaseNumber,
                                                           bapcodInit().startTime().getElapsedTime(),
                                                           candidatesEvaluatedOnThisPhase, bestProductScores);
            }
        }

      if (aCandidateIsConquered)
        {
          numCandForNextPhase = 1;
          /// this last candidate is conquered, we leave only it for the next phase
          /// by setting its scores to largest possible
          bestCandPtr->treeDepthBranchingScore = 0;
          bestCandPtr->productBranchingScore = BapcodInfinity;
        }

      /// we sort candidates by tree depth score (if exact phase) or by product score (non-exact phase)
      if (currPhaseParam.exact())
        std::stable_sort(candidateBranchingGroups.begin(), candidateBranchingGroups.end(),
                         CandidateSortByTreeDepthScore());
      else
        std::stable_sort(candidateBranchingGroups.begin(), candidateBranchingGroups.end(),
                         CandidateSortByProductScore());

      /// we cannot leave for the next phase candidates non evaluated on this phase
      if (numCandForNextPhase > candidatesEvaluatedOnThisPhase)
        numCandForNextPhase = candidatesEvaluatedOnThisPhase;

      /// this is change in BaPCod 052d
      /// we keep numCandForNextPhase best candidates according to the score
      /// if the phase is not exact, we delete not only bad candidates according to the score,
      /// but also candidates with similar bound increase to another one
      /// (we keep only one candidate with similar bound increase in all branches)
      int numKeptCandidates = 0;
      brGrPtrIt = candidateBranchingGroups.begin();
      std::vector<CandidateBranchGroup *>::iterator putPtrIt = candidateBranchingGroups.begin();
      while ((numKeptCandidates < numCandForNextPhase) && (brGrPtrIt != candidateBranchingGroups.end()))
      {
          std::vector<CandidateBranchGroup *>::iterator backPtrIt = brGrPtrIt;
          bool similarCandidateFound = false;
          if (!currPhaseParam.exact() && (backPtrIt != candidateBranchingGroups.begin()))
          {
              do
              {
                  --backPtrIt;

                  if (*backPtrIt == NULL)
                      continue;

                  if ((*backPtrIt)->productBranchingScore > (*brGrPtrIt)->productBranchingScore + Double::precision)
                      break;

                  bool similarCandidate = ((*brGrPtrIt)->size() == (*backPtrIt)->size());
                  CandidateBranchGroup::iterator candIt = (*brGrPtrIt)->begin();
                  CandidateBranchGroup::iterator backIt = (*backPtrIt)->begin();
                  while (similarCandidate && (candIt != (*brGrPtrIt)->end()) && (backIt != (*backPtrIt)->end()))
                  {
                      if (((*candIt)->nodeIncLpPrimalBound() < (*backIt)->nodeIncLpPrimalBound() - Double::precision)
                         || ((*candIt)->nodeIncLpPrimalBound() > (*backIt)->nodeIncLpPrimalBound() + Double::precision))
                          similarCandidate = false;
                      ++candIt;
                      ++backIt;
                  }
                  if (similarCandidate)
                  {
                      similarCandidateFound = true;
                      break;
                  }
              }
              while (backPtrIt != candidateBranchingGroups.begin());
          }
          if (similarCandidateFound)
          {
              delete *brGrPtrIt;
              *brGrPtrIt = NULL;
          }
          else
          {
              numKeptCandidates += 1;
              if (putPtrIt != brGrPtrIt)
              {
                  *putPtrIt = *brGrPtrIt;
                  *brGrPtrIt = NULL;
              }
              ++putPtrIt;
          }
          ++brGrPtrIt;
      }

      while (brGrPtrIt != candidateBranchingGroups.end())
      {
          delete *brGrPtrIt;
          *brGrPtrIt = NULL;
          ++brGrPtrIt;
      }

      candidateBranchingGroups.erase(putPtrIt, candidateBranchingGroups.end());

      prevPhaseParam = currPhaseParam;

      if (currentPhaseNumber == 1)
        bapcodInit().statistics().incrTimer("bcTimeSBphase1", startTime.getElapsedTime_dbl());
      else if (currentPhaseNumber == 2)
        bapcodInit().statistics().incrTimer("bcTimeSBphase2", startTime.getElapsedTime_dbl());
    }

  /// there should be only one candidate (it can be already conquered, but we need it for printing in the DOT file)
  if (candidateBranchingGroups.empty())
    {
      if (progStatus().doRun())
        std::cerr << "BaPCod WARNING : no candidates remain after strong branching with phases!" << std::endl;
      return;
    }
  if (candidateBranchingGroups.size() > 1)
    std::cerr << "BaPCod WARNING : several candidates remain after strong branching with phases!" << std::endl;

  /// copy nodes of the best candidate to to the current node as sons
  /// these nodes are then removed from the candidate (so that they will not be deleted after)
  for (CandidateBranchGroup::iterator nodePtrIt = candidateBranchingGroups.front()->begin();
       nodePtrIt != candidateBranchingGroups.front()->end(); nodePtrIt++)
    {
      /// if the last phase was not active (this can happen if only one candidate chosen in phase 0)
      /// or if the last phase was exact, then we do not need to recreate the node
      if (!prevPhaseParam.active() || prevPhaseParam.exact() || (*nodePtrIt)->isConquered())
        {
          (*nodePtrIt)->setEstimatedSubtreeSize(bestTreeSize); /// it is better to set a smaller value here
                                                               /// based on the DB improvement of this son
          (*nodePtrIt)->setBranchingPriorityLevel(curPriorityLevel);
          _currentNodePtr->sons().push_back(*nodePtrIt);
        }
      else
        {
          /// otherwise we need to evaluate nodes again
          Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                            (*nodePtrIt)->localNodeBrConstrList());
          newChildNodePtr->removeNodeEvalInfoAssociation();
          newChildNodePtr->associateNodeEvalInfoPtr((*nodePtrIt)->nodeEvalInfoPtr());
          /// the probelem info was saved in previous phase only if column generation was executed
          if (prevPhaseParam.maxNumOfColGenIterations() > 0)
            {
              newChildNodePtr->removeProblemSetupInfoAssociation();
              newChildNodePtr->associateProblemSetupInfoPtr((*nodePtrIt)->probSetupInfoPtr());
            }
          newChildNodePtr->setBranchEvaluationInfoFromNode(*nodePtrIt);
          newChildNodePtr->setBranchPhaseNumber(maxPhaseNumber + 1);
          delete *nodePtrIt;
          _currentNodePtr->sons().push_back(newChildNodePtr);
        }
    }
  /// we clear before deleting so that the associated nodes are not deleted
  candidateBranchingGroups.front()->clear();
  if (printL(0))
  {
    std::cout << "SB with phases chosed candidate ";
    candidateBranchingGroups.front()->genPtr->nicePrint();
    std::cout << std::endl;
  }
  delete candidateBranchingGroups.front();
}

/// replaces two functions: MasterConf::resetColGenSpListOfFractMastCol()
/// and MasterConf::ILOsortColGenSpListOfFractMastCol()
/// the copy of this function is in ColGenEvalAlg
/// TO DO : get rid of the two copies of the same function
void Alg4GenChildrenInBranching::resetAndILOsortColGenSpListOfFractMastCol(const MasterColSolution & listOfFractMastCol)
{
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    (*spcPt)->listOfFractMastColInColGenSp().clear();

  for (MasterColSolution::const_iterator colPt = listOfFractMastCol.begin();
       colPt != listOfFractMastCol.end(); colPt++)
    colPt->first->cgSpConfPtr()->listOfFractMastColInColGenSp().push_back(colPt->first, colPt->second);

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
  {
    ComponentSequence curClassCompBoundSet(*spcPt);
    MasterColSolution nonSortedListOfMastCol((*spcPt)->listOfFractMastColInColGenSp());
    (*spcPt)->listOfFractMastColInColGenSp().clear();
    CompBoundSetGenBranchConstr::ILOsortMastColumn(_currentNodePtr->treeOfColClasses(),
                                                   nonSortedListOfMastCol, curClassCompBoundSet,
                                                   (*spcPt)->listOfFractMastColInColGenSp());

    if (printL(6))
      std::cout << "ColGenSpConf " << (*spcPt)->name() << " has sorted list of frac col has size "
                << (*spcPt)->listOfFractMastColInColGenSp().size() << std::endl;
  }

}

void Alg4GenChildrenInBranching::findBranchingCandidatesWithColumns(const MasterColSolution & listOfFractMastCol,
                                                                    const int & maxNumOfCandidates,
                                                                    BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
  for (std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator
       dgcIt = _masterCommons.candidateBranchingGenericConstr().begin();
       dgcIt != _masterCommons.candidateBranchingGenericConstr().end(); ++dgcIt)
    {
      if (printL(5))
        std::cout << std::endl << "DynamicGenericConstr = " << (*dgcIt);
      (*dgcIt)->branchingSeparationFindCandidates(listOfFractMastCol, maxNumOfCandidates, generatedBrConstrGeneratorSet);
    }

  for (std::vector<ColGenSpConf *>::const_iterator spIt = _masterCommons.colGenSubProbConfPts().begin();
      spIt != _masterCommons.colGenSubProbConfPts().end(); ++spIt)
    {
      for (std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator
           dgcIt = (*spIt)->candidateBranchingGenericConstr().begin();
           dgcIt != (*spIt)->candidateBranchingGenericConstr().end(); ++dgcIt)
        {
          if (printL(5))
            std::cout << std::endl << "DynamicGenericConstr = " << (*dgcIt);
          (*dgcIt)->branchingSeparationFindCandidates(listOfFractMastCol, maxNumOfCandidates,
                                                      generatedBrConstrGeneratorSet);
        }
    }
}

void Alg4GenChildrenInBranching::performUsualBranching(MasterVarSolution & listOfFractPureMastVar,
                                                       MasterColSolution & listOfFractMastCol)
{
  Time startBcTimeSepFracSol;

  BranchGeneratorsSet generatedBrConstrGeneratorSet;
  if (param().BranchFirstOnPureMasterVariables)
    {
      if (listOfFractPureMastVar.size() >= 1)
        {

          if (printL(3))
            std::cout << "BranchingAlgorithm::run(): Branching in priority on PURE MAST VAR" << std::endl;

          for (std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator
                it = _masterCommons.candidateBranchingGenericConstr().begin();
                it != _masterCommons.candidateBranchingGenericConstr().end(); ++it)
            {
              if (printL(5))
                std::cout << std::endl << "DynamicGenericConstr = " << (*it);
              ///  Insert constraint prototype of the DynamicGenericConstr
              (*it)->branchingSeparationFindCandidates(listOfFractPureMastVar, listOfFractMastCol,
                                                       1, generatedBrConstrGeneratorSet);
            }
        }

      if (generatedBrConstrGeneratorSet.empty())
        {
          if (printL(3))
            std::cout << "BranchingAlgorithm::run(): Branching on projected SP VAR since no pure mast var is fractional"
                      << std::endl;
          findBranchingCandidatesWithColumns(listOfFractMastCol, 1, generatedBrConstrGeneratorSet);
        }
    }
  else if (param().BranchFirstOnSubProblemVariables)
    {
      if (listOfFractMastCol.size() >= 1)
        {
          if (printL(3))
            std::cout << "BranchingAlgorithm::run(): Branching in priority on projected SP VAR "
                      << "before testing if pure mast var are fractional" << std::endl;

          findBranchingCandidatesWithColumns(listOfFractMastCol, 1, generatedBrConstrGeneratorSet);

          if (printL(5))
            std::cout << std::endl << "BranchingAlgorithm::run(): generatedBrConstrGeneratorSet.size() = "
                      << generatedBrConstrGeneratorSet.size() << std::endl;
        }
      if (generatedBrConstrGeneratorSet.empty())
        {

          if (printL(3))
            std::cout << "BranchingAlgorithm::run(): Branching on PURE MAST VAR" << std::endl;

          for (std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator
               it = _masterCommons.candidateBranchingGenericConstr().begin();
               it != _masterCommons.candidateBranchingGenericConstr().end(); ++it)
            {
              if (printL(5))
                std::cout << std::endl << "DynamicGenericConstr = " << (*it);

              ///  Returns a constraint prototype of the DynamicGenericConstr
              (*it)->branchingSeparationFindCandidates(listOfFractPureMastVar, listOfFractMastCol,
                                                       1 ,generatedBrConstrGeneratorSet);
            }
        }
    }///  Mix all
  else
    {
      if (printL(3))
        std::cout << "Branching selection over PURE MAST VAR and masterColumn VAR" << std::endl;

      MasterVarSolution listOfVarInSol;

      for (MasterVarSolution::const_iterator colIt = listOfFractPureMastVar.begin();
           colIt != listOfFractPureMastVar.end(); ++colIt)
        listOfVarInSol.push_back(colIt->first, colIt->second);

      for (MasterColSolution::const_iterator colIt = listOfFractMastCol.begin();
           colIt != listOfFractMastCol.end(); ++colIt)
        listOfVarInSol.push_back(colIt->first, colIt->second);

      for (std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator
           it = _masterCommons.candidateBranchingGenericConstr().begin();
           it != _masterCommons.candidateBranchingGenericConstr().end(); ++it)
        {
          if (printL(5))
            std::cout << "DynamicGenericConstr = " << (*it) << std::endl;

          (*it)->branchingSeparationFindCandidates(listOfVarInSol, listOfFractMastCol, 1, generatedBrConstrGeneratorSet);
        }
    }

  bapcodInit().statistics().incrTimer("bcTimeSepFracSol", startBcTimeSepFracSol.getElapsedTime_dbl());
  
  if (generatedBrConstrGeneratorSet.empty())
    {
      std::cerr << "BaPCod WARNING: Alg4GenChildrenInBranching::performUsualBranching():"
                << " no candidates were found" << std::endl;
      if (printL(-1))
        std::cout << "BaPCod warning : solution is not integer but no branching candidates found" << std::endl;
      return;
    }
    
  for (BranchGeneratorsSet::const_iterator brGenPtrIt = ++generatedBrConstrGeneratorSet.begin();
       brGenPtrIt != generatedBrConstrGeneratorSet.end(); ++brGenPtrIt)
    delete (*brGenPtrIt);

  BranchingConstrGenerator * brConstrGenPtr = *(generatedBrConstrGeneratorSet.begin());
  brConstrGenPtr->computeLhs(_currentNodePtr->primalSol());

  std::list<BranchingConstrBaseType *> tmpLocalNodeBrConstrList;
  ConstrPtrSet dummyBranchingConstrSet;

  if (printL(0))
  {
    std::cout << "Chosen branch : ";
    brConstrGenPtr->nicePrint();
    std::cout << std::endl;
  }

  while (brConstrGenPtr->nextNodeBrConstr(_currentNodePtr, tmpLocalNodeBrConstrList, dummyBranchingConstrSet))
    {
      Node * newChildNodePtr = new Node(_masterCommons.getNodeCountAndIncreaseIt(), _currentNodePtr,
                                        tmpLocalNodeBrConstrList);
      for (std::list<BranchingConstrBaseType *>::const_iterator it = tmpLocalNodeBrConstrList.begin();
           it != tmpLocalNodeBrConstrList.end(); ++it)
        {
          (*it)->depthWhenGenerated(newChildNodePtr->depth()); /// added by Ruslan
          (*it)->append2name(std::string() + "n" + newChildNodePtr->ref());
        }
      _currentNodePtr->sons().push_back(newChildNodePtr);
    }

  delete brConstrGenPtr;
}

void Alg4GenChildrenInBranching::completeMasterVarSolution(const MasterColSolution & listOfFractMastCol,
                                                           const MasterVarSolution & listOfFractPureMastVar,
                                                           MasterVarSolution & listOfMastAndSubprobVar)
{
  for (MasterVarSolution::const_iterator listIt = listOfFractPureMastVar.begin();
       listIt != listOfFractPureMastVar.end(); ++listIt)
  {
    listOfMastAndSubprobVar.push_back(listIt->first, listIt->second);
  }

  /// project onto  OVF solution
  VarPtr2DoubleMap curAggregateMastSol;
  for (MasterColSolution::const_iterator listIt = listOfFractMastCol.begin(); listIt != listOfFractMastCol.end();
       ++listIt)
  {
    if (printL(5))
      std::cout << "consider Master Var " << listIt->first->name() << " with val "
                << listIt->second._value << std::endl;

    listIt->first->fillAggregateSol(curAggregateMastSol, listIt->second._value);
  }

  /// Record only frac spVar Value with Dfrac > param.BapCodCutViolationTolerance
  for (VarPtr2DoubleMap::iterator camIt = curAggregateMastSol.begin(); camIt != curAggregateMastSol.end(); camIt++)
  {
    if (printL(5))
      std::cout << "consider sp var " << (camIt->first)->name()
                << " with aggregate use " << camIt->second << std::endl;

    if (!camIt->first->isCandForBranching())
      continue;

    if (printL(5))
      std::cout << "check sp var = " << camIt->first->name() << " use " << camIt->second << std::endl;

    listOfMastAndSubprobVar.push_back(camIt->first, camIt->second);
  }
}

void Alg4GenChildrenInBranching::performStrongBranching(int & globalTreatOrder,
                                                        MasterVarSolution & listOfFractPureMastVar,
                                                        MasterColSolution & listOfFractMastCol)
{
  Time startBcTimeSepFracSol;

  std::set<double> priorityLevelsSet;
  std::vector<GenVarGenBranchConstr *> genVarGenBrGenConstrPts;
  std::vector<CompBoundSetGenBranchConstr *> compBoundSetBrGenConstrPts;
  std::vector<RyanAndFosterGenBranchConstr *> ryanAndFosterBrGenConstrPts;
#ifdef BCP_RCSP_IS_FOUND
  std::vector<PackSetResConsGenBranchConstr *> PackSetResConsBrGenConstrPts;
  std::vector<PackSetRyanFosterGenBranchConstr *> packSetRyanFosterBrGenConstrPts;
#endif //BCP_RCSP_IS_FOUND
  std::vector<GenericBranchingConstr *> otherBrGenConstrPts;
  std::set<GenericBranchingConstr *, DynamicGenConstrSort>::const_iterator genBrConstrPtrIt;
  for (genBrConstrPtrIt = _masterCommons.candidateBranchingGenericConstr().begin();
       genBrConstrPtrIt != _masterCommons.candidateBranchingGenericConstr().end(); ++genBrConstrPtrIt)
    {
      priorityLevelsSet.insert((*genBrConstrPtrIt)->priorityLevel());
      CompBoundSetGenBranchConstr * compBoundSetGenBrConstrPtr
        = dynamic_cast<CompBoundSetGenBranchConstr *>(*genBrConstrPtrIt);
      if (compBoundSetGenBrConstrPtr != NULL)
        {
          compBoundSetBrGenConstrPts.push_back(compBoundSetGenBrConstrPtr);
          continue;
        }
      GenVarGenBranchConstr * genVarGenBrConstrPtr = dynamic_cast<GenVarGenBranchConstr *>(*genBrConstrPtrIt);
      if (genVarGenBrConstrPtr != NULL)
        {
          genVarGenBrGenConstrPts.push_back(genVarGenBrConstrPtr);
          continue;
        }
#ifdef BCP_RCSP_IS_FOUND
      PackSetResConsGenBranchConstr * esrcGenBrConstrPtr
              = dynamic_cast<PackSetResConsGenBranchConstr *>(*genBrConstrPtrIt);
      if (esrcGenBrConstrPtr != NULL)
      {
        PackSetResConsBrGenConstrPts.push_back(esrcGenBrConstrPtr);
        continue;
      }
      PackSetRyanFosterGenBranchConstr * psrfGenBrConstrPtr
            = dynamic_cast<PackSetRyanFosterGenBranchConstr *>(*genBrConstrPtrIt);
      if (psrfGenBrConstrPtr != NULL)
      {
        packSetRyanFosterBrGenConstrPts.push_back(psrfGenBrConstrPtr);
        continue;
      }
#endif //BCP_RCSP_IS_FOUND
      GenericBranchingConstr * genBrConstrPtr = dynamic_cast<GenericBranchingConstr *>(*genBrConstrPtrIt);
      if (genBrConstrPtr != NULL)
        {
          otherBrGenConstrPts.push_back(genBrConstrPtr);
          continue;
        }
      else
        {
          std::cerr << "BaPCod error : candidateBranchingGenericConstr should be of type GenericBranchingConstr."
                    << std::endl;
          exit(1);
        }
    }
  std::vector<ColGenSpConf *>::const_iterator cgSpPtrIt;
  for (cgSpPtrIt = _masterCommons.colGenSubProbConfPts().begin();
       cgSpPtrIt != _masterCommons.colGenSubProbConfPts().end(); ++cgSpPtrIt)
    for (genBrConstrPtrIt = (*cgSpPtrIt)->candidateBranchingGenericConstr().begin();
         genBrConstrPtrIt != (*cgSpPtrIt)->candidateBranchingGenericConstr().end(); ++genBrConstrPtrIt)
      {
        priorityLevelsSet.insert((*genBrConstrPtrIt)->priorityLevel());
        RyanAndFosterGenBranchConstr * ryanAndFosterGenBrConstrPtr
          = dynamic_cast<RyanAndFosterGenBranchConstr *>(*genBrConstrPtrIt);
        if (ryanAndFosterGenBrConstrPtr != NULL)
          {
            ryanAndFosterBrGenConstrPts.push_back(ryanAndFosterGenBrConstrPtr);
            continue;
          }
        GenericBranchingConstr * genBrConstrPtr = dynamic_cast<GenericBranchingConstr *>(*genBrConstrPtrIt);
        if (genBrConstrPtr != NULL)
          {
            otherBrGenConstrPts.push_back(genBrConstrPtr);
            continue;
          }
        else
          {
            std::cerr << "BaPCod error : candidateBranchingGenericConstr should be of type GenericBranchingConstr."
                      << std::endl;
            exit(1);
          }
      }
    
  
  /// we get the estimated subtree size based on the history or (if it is not available),
  /// based on the gap increse of _currentNodePtr
  Double estimatedSubtreeSize(BapcodInfinity);
  bool averageSubreeSizeIsAvailable = false;
  if (param().StrongBranchingUseHistorySubtreeSize()
      && _masterCommons.getAverageSubtreeSize(_currentNodePtr->depth(), estimatedSubtreeSize))
    {
      if (printL(0))
        std::cout << "Average history subtree size rooted at depth " << _currentNodePtr->depth()
                  << " is " << estimatedSubtreeSize << std::endl;
      /// we correct estimated subtree size of a node, if historical data is available
      _currentNodePtr->setEstimatedSubtreeSize(estimatedSubtreeSize);
      averageSubreeSizeIsAvailable = true;
    }
  else
    {
      estimatedSubtreeSize = _currentNodePtr->estimatedSubtreeSize();
      if (printL(0) && (estimatedSubtreeSize < BapcodInfinity))
          std::cout << "Estimated subtree size of the father is " << estimatedSubtreeSize << std::endl;
    }
  int remNumOfCandidates = param().StrongBranchingPhaseOne().maxNumOfCandidates();
  /// if estimated subtree size is small enough, we decrease the number of candidates which will be
  /// evaluated at phase one
  double treeSizeNumOfCandidates = 1000000;
  if (param().StrongBranchingPhaseOne().treeSizeRatioToStop() > 0.0)
    treeSizeNumOfCandidates = ceil((double)estimatedSubtreeSize
                                   * param().StrongBranchingPhaseOne().treeSizeRatioToStop());
  if (treeSizeNumOfCandidates > 1000000)
    treeSizeNumOfCandidates = 1000000;
  if (remNumOfCandidates > treeSizeNumOfCandidates)
    remNumOfCandidates = treeSizeNumOfCandidates;

  std::vector<BranchGeneratorsSet> brGenSetVector(1, BranchGeneratorsSet());
  double minPriorityLevel = 0.0;
  double curPriorityLevel = 0.0;
  int numberOfBranchingCandidatesFound = 0;
  /// we try to genereate candidates in the non-increasing priority order
  for (std::set<double>::reverse_iterator setIt = priorityLevelsSet.rbegin();
       (numberOfBranchingCandidatesFound < remNumOfCandidates) && (setIt != priorityLevelsSet.rend()); ++setIt)
    {
      curPriorityLevel = *setIt;

      if ((curPriorityLevel < minPriorityLevel) && (numberOfBranchingCandidatesFound > 0))
        break;
      minPriorityLevel = floor(curPriorityLevel);

      if ((curPriorityLevel < _currentNodePtr->branchingPriorityLevel()) && !averageSubreeSizeIsAvailable)
      {
          remNumOfCandidates = param().StrongBranchingPhaseOne().maxNumOfCandidates();
          _currentNodePtr->setEstimatedSubtreeSize(BapcodInfinity);
      }

      /// we use the same set of branching constraints for GenVarGenBranchConstr corresponding to the same var name
      if (!brGenSetVector.back().empty())
        brGenSetVector.push_back(BranchGeneratorsSet());
      std::multiset<std::pair<std::string, GenVarGenBranchConstr *> > gvgPtrSet;
      for (std::vector<GenVarGenBranchConstr *>::iterator genBrConstrPtrIt = genVarGenBrGenConstrPts.begin();
           genBrConstrPtrIt != genVarGenBrGenConstrPts.end(); ++genBrConstrPtrIt)
        if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
          gvgPtrSet.insert(std::make_pair((*genBrConstrPtrIt)->genVarPtr()->defaultName(), *genBrConstrPtrIt));

      std::multiset<std::pair<std::string, GenVarGenBranchConstr *> >::iterator gvgSetIt;
      std::string varName("");
      for (gvgSetIt = gvgPtrSet.begin(); gvgSetIt != gvgPtrSet.end(); ++gvgSetIt)
        {
          if (gvgSetIt->first != varName)
            {
              varName = gvgSetIt->first;
              if (!brGenSetVector.back().empty())
                brGenSetVector.push_back(BranchGeneratorsSet());
            }
          gvgSetIt->second->branchingSeparationFindCandidates(listOfFractPureMastVar, listOfFractMastCol,
                                                              remNumOfCandidates, brGenSetVector.back());
          gvgSetIt->second->branchingSeparationFindCandidates(listOfFractMastCol, remNumOfCandidates,
                                                              brGenSetVector.back());
          numberOfBranchingCandidatesFound += brGenSetVector.back().size();
        }

      /// we use the same set of branching constraints for GenVarGenBranchConstr
      if (!brGenSetVector.back().empty())
        brGenSetVector.push_back(BranchGeneratorsSet());
      for (std::vector<CompBoundSetGenBranchConstr *>::iterator genBrConstrPtrIt = compBoundSetBrGenConstrPts.begin();
           genBrConstrPtrIt != compBoundSetBrGenConstrPts.end(); ++genBrConstrPtrIt)
        if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
          {
            (*genBrConstrPtrIt)->branchingSeparationFindCandidates(listOfFractMastCol, remNumOfCandidates,
                                                                   brGenSetVector.back());
            numberOfBranchingCandidatesFound += brGenSetVector.back().size();
          }
    
      /// we use the same set of branching constraints for RyanAndFosterGenBranchConstr
      if (!brGenSetVector.back().empty())
        brGenSetVector.push_back(BranchGeneratorsSet());
      for (std::vector<RyanAndFosterGenBranchConstr *>::iterator genBrConstrPtrIt = ryanAndFosterBrGenConstrPts.begin();
           genBrConstrPtrIt != ryanAndFosterBrGenConstrPts.end(); ++genBrConstrPtrIt)
        if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
          {
            (*genBrConstrPtrIt)->branchingSeparationFindCandidates(listOfFractMastCol, remNumOfCandidates,
                                                                   brGenSetVector.back());
            numberOfBranchingCandidatesFound += brGenSetVector.back().size();
          }

#ifdef BCP_RCSP_IS_FOUND
      if (!brGenSetVector.back().empty())
        brGenSetVector.push_back(BranchGeneratorsSet());
      for (std::vector<PackSetResConsGenBranchConstr *>::iterator genBrConstrPtrIt
              = PackSetResConsBrGenConstrPts.begin(); genBrConstrPtrIt != PackSetResConsBrGenConstrPts.end();
           ++genBrConstrPtrIt)
        if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
        {
          (*genBrConstrPtrIt)->branchingSeparationFindCandidates(listOfFractMastCol, remNumOfCandidates,
                                                                 brGenSetVector.back());
          numberOfBranchingCandidatesFound += brGenSetVector.back().size();
        }

    if (!brGenSetVector.back().empty())
      brGenSetVector.push_back(BranchGeneratorsSet());
    for (std::vector<PackSetRyanFosterGenBranchConstr *>::iterator genBrConstrPtrIt
            = packSetRyanFosterBrGenConstrPts.begin(); genBrConstrPtrIt != packSetRyanFosterBrGenConstrPts.end();
         ++genBrConstrPtrIt)
      if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
      {
        (*genBrConstrPtrIt)->branchingSeparationFindCandidates(listOfFractMastCol, remNumOfCandidates,
                                                               brGenSetVector.back());
        numberOfBranchingCandidatesFound += brGenSetVector.back().size();
      }
#endif //BCP_RCSP_IS_FOUND

      /// we use different set of branching constraints for other GenericBranchingConstr
      for (std::vector<GenericBranchingConstr *>::iterator genBrConstrPtrIt = otherBrGenConstrPts.begin();
           genBrConstrPtrIt != otherBrGenConstrPts.end(); ++genBrConstrPtrIt)
        if ((*genBrConstrPtrIt)->priorityLevel() == curPriorityLevel)
          {
            if (!brGenSetVector.back().empty())
              brGenSetVector.push_back(BranchGeneratorsSet());
            MasterVarSolution listOfMastAndSubprobVar;
            completeMasterVarSolution(listOfFractMastCol, listOfFractPureMastVar, listOfMastAndSubprobVar);

            (*genBrConstrPtrIt)->branchingSeparationFindCandidates(listOfMastAndSubprobVar, listOfFractMastCol,
                                                                   remNumOfCandidates,
                                                                   brGenSetVector.back());
            numberOfBranchingCandidatesFound += brGenSetVector.back().size();
          }
    }

  bapcodInit().statistics().incrTimer("bcTimeSepFracSol", startBcTimeSepFracSol.getElapsedTime_dbl());
  
  if (numberOfBranchingCandidatesFound == 0)
    {
      if (printL(5))
        std::cout << "Alg4GenChildrenInBranching::performStrongBranching(): no candidates were found" << std::endl;
      if (printL(-1))
        std::cout << "BaPCod warning : solution is not integer but no branching candidates found" << std::endl;
      return;
    }

  std::vector<BranchingConstrGenerator *> brConstrGeneratorsFromHistory;
  if (param().StrongBranchingUseHistory())
  {
      getConstrGeneratorsFromHistory(brConstrGeneratorsFromHistory, remNumOfCandidates / 2, curPriorityLevel);
      remNumOfCandidates -= brConstrGeneratorsFromHistory.size();
  }

  /// we copy generators from the history to the set to be able to recognise if
  /// a new generator coincides with one from the history or not
  std::set<BranchingConstrGenerator *, BranchingGeneratorPtrComp> generatorsHistorySet;
  std::vector<BranchingConstrGenerator *>::iterator histBrGenPtrIt;
  for (histBrGenPtrIt = brConstrGeneratorsFromHistory.begin();
       histBrGenPtrIt != brConstrGeneratorsFromHistory.end(); ++histBrGenPtrIt)
    generatorsHistorySet.insert(*histBrGenPtrIt);
    
  /// we now sort sets of branching generators by their size (in non-decreasing order)
  std::vector<std::pair<int, int> > brGenSetSizes;
  for (int brGenSetId = 0; brGenSetId < (int)brGenSetVector.size(); ++brGenSetId)
    if (!brGenSetVector[brGenSetId].empty())
      brGenSetSizes.push_back(std::make_pair(brGenSetVector[brGenSetId].size(), brGenSetId));

  std::stable_sort(brGenSetSizes.begin(), brGenSetSizes.end());
    
  std::vector<BranchingConstrGenerator *> newBranchingGenerators;
  int numSetsOfGenerators = (int)brGenSetSizes.size();
  for (std::vector<std::pair<int, int> >::iterator pairIt = brGenSetSizes.begin();
       pairIt != brGenSetSizes.end(); ++pairIt)
    {
      BranchGeneratorsSet & brGeneratorsSet = brGenSetVector[pairIt->second];
      int remNumGeneratorsFromThisSet = remNumOfCandidates / numSetsOfGenerators;
      BranchGeneratorsSet::iterator brConstrGenPtrIt;
      for (brConstrGenPtrIt = brGeneratorsSet.begin(); brConstrGenPtrIt != brGeneratorsSet.end(); ++brConstrGenPtrIt)
        {
          if ((remNumGeneratorsFromThisSet > 0) && !generatorsHistorySet.count(*brConstrGenPtrIt))
            {
              (*brConstrGenPtrIt)->computeLhs(_currentNodePtr->primalSol());
              newBranchingGenerators.push_back(*brConstrGenPtrIt);
              remNumGeneratorsFromThisSet -= 1;
              remNumOfCandidates -= 1;
            }
          else
            delete *brConstrGenPtrIt;
        }
      brGeneratorsSet.clear();
      numSetsOfGenerators -= 1;
    }

  strongBranchingWithPhases(curPriorityLevel, globalTreatOrder, brConstrGeneratorsFromHistory, newBranchingGenerators);
}

void Alg4GenChildrenInBranching::run(int & globalTreatOrder)
{
  MasterVarSolution listOfFractPureMastVar;
  MasterColSolution listOfFractMastCol;

  for (SolutionVarInfoPtrList::const_iterator varInfoPtrIt = _currentNodePtr->primalSol().begin();
       varInfoPtrIt != _currentNodePtr->primalSol().end(); ++varInfoPtrIt)
    {
        ValueRecord rec((*varInfoPtrIt)->value, param().BapCodIntegralityTolerance());
        double veryHighPrecision = Double::precision * 1e-4; /// used to be VERYHIGHPRECISION parameter with value 1e-10
        if (rec._dfracValue > veryHighPrecision)
        {
            Variable * varPtr = (*varInfoPtrIt)->varPtr;
            if (varPtr->isTypeOf(VcId::MastColumnMask))
                listOfFractMastCol.push_back(varPtr, rec);
            else if ((varPtr->type() == 'B') || (varPtr->type() == 'I'))
                listOfFractPureMastVar.push_back(varPtr, rec);
        }
    }

  if (listOfFractPureMastVar.empty() && listOfFractMastCol.empty())
    {
      if (printL(3))
        std::cout << "BranchingAlgorithm::run(): no more fract VAR" << std::endl;

      bapcodInit().check(1, "BranchingAlgorithm::run() the master solution should have been detected as integer",
                         ProgStatus::run);
      return;
    }

  /// for branchingSeparationFindCandidates() functions,
  /// we need to put col gen subproblems trees of col classes
  /// TO DO: to use directly _currentNodePtr->cgSpConfTreeOfColClassesMap() in branchingSeparationFindCandidates()
  for (std::vector<ColGenSpConf *>::const_iterator spcPt = _masterCommons.colGenSubProbConfPts().begin();
       spcPt != _masterCommons.colGenSubProbConfPts().end(); ++spcPt)
    if (_currentNodePtr->cgSpConfTreeOfColClassesMap().count(*spcPt))
      (*spcPt)->treeOfColClasses() = _currentNodePtr->cgSpConfTreeOfColClassesMap()[*spcPt];

  resetAndILOsortColGenSpListOfFractMastCol(listOfFractMastCol);

  if (param().StrongBranchingPhaseOne().active())
    performStrongBranching(globalTreatOrder, listOfFractPureMastVar, listOfFractMastCol);
  else
    performUsualBranching(listOfFractPureMastVar, listOfFractMastCol);

}

Alg4GenChildrenInBranching::Alg4GenChildrenInBranching(MasterCommons4GenChildNodesAlgorithm & masterCommons) :
  Alg4GenChildrenOfNode(masterCommons), _firstNodeWasEvaluated(false), _lastEvaluatedNodeWasLp(false)
{
}

Alg4GenChildrenInBranching::~Alg4GenChildrenInBranching()
{
}

