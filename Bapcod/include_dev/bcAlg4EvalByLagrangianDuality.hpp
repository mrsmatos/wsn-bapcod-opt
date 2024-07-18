/**
 *
 * This file bcAlg4EvalByLagrangianDuality.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCLAGREVALALG_HPP_
#define BCLAGREVALALG_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalOfNode.hpp"

class ColGenStabilization;

class PricingStrategy : public Parameter<SelectionStrategy>
{
public:
  enum PricingStrategyEnum
  {
    AllSubproblems = 0,
    Cycling = 1,
    Intensification = 2,
    Diversification = 3
  };
private:
  PricingStrategyEnum _selectedRule;
public:
  PricingStrategy(const int & stat = -1)
  {
    set(stat);
  }
  PricingStrategy(const PricingStrategyEnum & stat): _selectedRule(stat) {}
  virtual ~PricingStrategy()
  {
  }
  const PricingStrategyEnum &
  selectedRule() const
  {
    return _selectedRule;
  }
  bool set(const int & stat) {
    switch (stat)
    {
      case 0:
        _selectedRule = AllSubproblems;
        break;
      case 1:
        _selectedRule = Cycling;
        break;
      case 2:
        _selectedRule = Intensification;
        break;
      case 3:
        _selectedRule = Diversification;
      default:
        return false;
    }
    return true;

  }
};


class Alg4EvalByLagrangianDuality : public Alg4EvalOfNode
{
protected:
  long _savedMasterLPTime;
  long _savedPricingTime;
  long _savedNbColumns;

  /// Eval options
  int _maxNbOfPenaltyUpdates;
  int _maxNbOfCgIterations;
  int _minLevelOfSpRestriction;
  int _minNbOfCutRounds;
  int _maxNbOfCutRounds;
  int _doRedCostFixingAndEnumeration;
  int _logPrintFrequency;
  bool _nonExactEvaluation;

  std::vector<ColGenSpConf *>::const_iterator _lastSpcPt;
  bool _isLastSpcPtInitialized;

  ColGenStabilization * _colGenStabilizationPtr;

  //std::list<Variable *> _nonStabArtVarPtrList;
  std::list< MastColumn *> _candidatColWithNegRedCost;

  bool _need2reintroduceArtVarInMast;
  bool _currentlyPerformingPhase1;

  PricingStrategy _pricingStrat;

  /// Indices of suproblems which when last tested found a col with neg reduced cost
  std::list<int> _promisingOrOpenSpIndices;

  /// Indices of suproblems which when last tested did not find a col with neg reduced cost
  std::list<int> _unpromisingSpIndices;

  void printIntermediateStatistics(std::ostream & os, const int & curMaxLevelOfSubProbRestriction,
                                   const int & nbNumColumns,
                                   const int & nbCgIterations, long & elapsedTime,
                                   const bool & doNotCheckLogFrequency = false,
                                   const bool & printIncumbentDual = false);

  bool earlyCGtermType1();

  int searchNegRedCostInactiveCol();

  virtual void compMastDualBoundContrib(Bound & mastDualBoundContrib);

  void cleanupRestrictedMastColumns(const int & nbCgIterations);
  void cleanupRestrictedMastCuts();

  void recordColInForm();

  int genNewColumns(int & maxLevelOfSubProbRestriction, int & allAddedColumns, int & addedNegRedCostColumns);

  void updatePricingSolverCutsMessageId();

  void updateLagrangianDualBound(bool updateDualBnd);

  bool phase1isCompleted();

  virtual bool setPurePhaseI()
  {
    _currentlyPerformingPhase1 = true;
    return false;
  }

  virtual bool unsetPurePhaseI()
  {
    _currentlyPerformingPhase1 = false;
    return false;
  }

public:

  Alg4EvalByLagrangianDuality(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons);

  virtual ~Alg4EvalByLagrangianDuality();

  virtual bool setupAlgo(Node * nodePtr) = 0;

  virtual bool eval() = 0;

  virtual NodeEvalInfo * recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr = NULL) = 0;

  virtual void setDownAlgo() = 0;

  virtual void setOptionLogPrintFrequency(const int value);
  virtual void setOptionMaxNbOfCgIterations(const int value);
  virtual void setOptionMaxNbOfPenaltyUpdates(const int value);
  virtual void setOptionMinLevelOfSpRestriction(const int value);
  virtual void setOptionMinNbOfCutRounds(const int value);
  virtual void setOptionMaxNbOfCutRounds(const int value);
  virtual void setOptionDoRedCostFixingAndEnumeration(const int value);
  virtual void setOptionNonExactEvaluation(const bool value);
};
#endif /* BCLAGREVALALG_HPP_ */
