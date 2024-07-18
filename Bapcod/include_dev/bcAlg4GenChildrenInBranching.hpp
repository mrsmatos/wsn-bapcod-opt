/**
 *
 * This file bcAlg4GenChildrenInBranching.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCBRANCHING_HPP__
#define BCBRANCHING_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4GenChildrenOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"


class Alg4EvalOfNode;
class Algorithm4PreprocessingAtNodeOtherThanRoot;
class Alg4ProblemSetDownOfNode;
class ProblemFullSetDownAlgorithm;
class Alg4ProblemSetupFull;

class CandidateBranchGroup : public std::vector<Node *>
{
protected:

  BcObjStatus::MinMaxIntFloat objStatus;
  static const double treeDepthScoreNormConstant;

  const double numberOfLeafs(const double & gap, const std::vector<double> & delta);

public:
  int id;
  double productBranchingScore; /// double is important for applications like bin packing 
  Double linearBranchingScore;
  Double treeDepthBranchingScore;
  Double treeSize;
  Double treeDepth;
  BranchingConstrGenerator * genPtr;
  Bound pessimisticDualBound;
  Bound optimisticDualBound;
  bool historyCandidate;

  CandidateBranchGroup(const int & id_,
                       BranchingConstrGenerator * generatorPtr,
                       const BcObjStatus::MinMaxIntFloat & objStatus_,
                       const bool & historyCandidate_ = false);

  virtual ~CandidateBranchGroup();

  void computeBranchingScoreAsProduct(const Bound & parentFractionalDualBound,
                                      const Double & cutOffValue,
                                      const int & maxDescriptionLength,
                                      const int & currPhaseNumber,
                                      const long & elapsedTime,
                                      const int & candidatesEvaluated,
                                      double bestProductScores[3]);

  void computeBranchingTreeDepthScore(const Bound & parentFractionalDualBound,
                                      const Double & cutOffValue,
                                      const int & numEvaluatedChildren,
                                      const double notPromissingConservativeness);
  inline void resetAllScores()
  {
    productBranchingScore = -BapcodInfinity;
    linearBranchingScore = -BapcodInfinity;
    treeDepthBranchingScore = -BapcodInfinity;
    treeSize = treeDepth = BapcodInfinity;
  }

} ;

struct CandidateSortByProductScore
{
  bool operator()(CandidateBranchGroup * a, CandidateBranchGroup * b) const
  {
    if (a->productBranchingScore == b->productBranchingScore)
    {
      return a->id < b->id;
    }
    return (a->productBranchingScore > b->productBranchingScore);
  }
};

struct CandidateSortByTreeDepthScore
{
  bool operator()(CandidateBranchGroup * a, CandidateBranchGroup * b) const
  {
    if (a->treeDepthBranchingScore == b->treeDepthBranchingScore)
    {
      return a->front()->ref() < b ->front()->ref();
    }
    return (a->treeDepthBranchingScore > b->treeDepthBranchingScore);
  }
};

class Alg4GenChildrenInBranching : public Alg4GenChildrenOfNode
{
  bool _firstNodeWasEvaluated;
  bool _lastEvaluatedNodeWasLp;

  void resetAndILOsortColGenSpListOfFractMastCol(const MasterColSolution & listOfFractMastCol);

  void findBranchingCandidatesWithColumns(const MasterColSolution & listOfFractMastCol,
                                          const int & maxNumOfCandidates,
                                          BranchGeneratorsSet & generatedBrConstrGeneratorSet);
  void prepareNodeForTreatmentInStrBranchWithPhases(Node * nodePtr, int globalTreatOrder,
                                                    const StrongBranchingPhaseParameter & currPhaseParam);
  virtual void strongBranchingWithPhases(const double & curPriorityLevel, int & globalTreatOrder,
                                         std::vector<BranchingConstrGenerator *> & branchGeneratorsFromHistory,
                                         std::vector<BranchingConstrGenerator *> & newBranchGenerators);
  virtual void getConstrGeneratorsFromHistory(std::vector<BranchingConstrGenerator *> & brConstrGeneratorsFromHistory,
                                              const int & maxNbOfCandidates, const Double & minPriorityLevel);
  void performUsualBranching(MasterVarSolution & listOfFractPureMastVar,
                             MasterColSolution & listOfFractMastCol);

  void completeMasterVarSolution(const MasterColSolution & listOfFractMastCol,
                                 const MasterVarSolution & listOfFractPureMastVar,
                                 MasterVarSolution & listOfMastAndSubprobVar);

  void performStrongBranching(int & globalTreatOrder, MasterVarSolution & listOfFractPureMastVar,
                              MasterColSolution & listOfFractMastCol);
public:

  Alg4GenChildrenInBranching(MasterCommons4GenChildNodesAlgorithm & masterCommons);

  virtual ~Alg4GenChildrenInBranching();

  virtual bool setupAlgo(Node * nodePtr)
  {
    return Alg4GenChildrenOfNode::setupAlgo(nodePtr);
  }

  virtual void run(int & globalTreatOrder);

  virtual void setDownAlgo()
  {
    Alg4GenChildrenOfNode::setDownAlgo();
  }

} ;


#endif /* BCBRANCHING_HPP_ */
