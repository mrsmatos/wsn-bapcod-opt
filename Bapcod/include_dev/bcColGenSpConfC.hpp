/**
 *
 * This file bcColGenSpConfC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCCOLGENSPCONFC_H_
#define BCCOLGENSPCONFC_H_

#include "bcProbConfigC.hpp"
#include "bcUsefulHeadFil.hpp"

class ColGenSpConf: public ProbConfig
{
 private:
  MasterConf * _mastConfPtr;
  BcObjStatus::MinMaxIntFloat _objStatus;

  /// Belongs to MasterConf
  MissingColumn * _misColPtr;
  std::list<InstanciatedConstr *> _tempMastConstrPtrList4Insertion;
  std::list<MastColumn *> _tempColPtrList4Insertion; /// to store generated column after test for insertion and until
                                                     /// next round of genNewColum (reseted in updateConf)

  MasterColSolution _listOfFractMastColInColGenSp; //(Issam)
 protected:

  /**
   * Use to know wheter CSbranching constraint induces tigher bounds on SP variables
   */
  bool _spConfHasClassInducedSpVarBounds;

  /**
   * That are memorised in _defaultClassLb and _defaultClassUb Belongs to ColGenSpConf
   */
  bool _implicitlyFixCardinality; /// if true it means that _upperBound and _lowerBound are implicitly equal, 
                                  /// i.e. one of the convexity constraint has not been introduced in the model
                                  /// because it is equal to the other and not binding

  Double * _upperBoundPtr;

  /// Belongs to ColGenSpConf
  Double * _lowerBoundPtr;
  Double _defaultDualVal4UbConstr;
  Double _defaultDualVal4LbConstr;
  InstMastConvexityConstr * _lowerBoundMastConstrPtr;
  InstMastConvexityConstr * _upperBoundMastConstrPtr;
  Double _fixedCost;
  Double _spRootReducedCost;
  bool _toSplit;
  bool _rollbackPointSaved;

  /**
   * fixedDualCost is the dual cost associated to convexity constraints
   */
  Double _fixedDualCost;
  /**
   * Target is the min sp value required
   * in order to obtained a negative reduced
   * cost column (= infty if no requirement)
   * if (minCost >= target) no attempt is made
   * generate new columns
   */
  Bound _target;

  /**
   * Temporary record of spcof solution multiplicit
   */
  Double _mult;

  /// added by Ruslan, needed for prioriterizing columns in the diving heuristic
  Double _priorityLevel;

  /**
   * Reset min Cost to min Sp value and target according du dual master value
   */
  void updateTarget(const bool currentlyPerformingPhaseI);
  void correctDualboundContribAndBestSol(const Bound & rootDualBound, Bound & bestRedCost, Solution * & bestSolPtr);
  void clearSolutions();

public:

   virtual MasterColSolution & listOfFractMastColInColGenSp() {return _listOfFractMastColInColGenSp;}

  ColGenSpConf(const std::string genericName,
	           const IndexCell & id,
	           MasterConf * mastConfPtr,
	           const Double & fixedCost,
	           const bool & implicitlyFixCardinality,
               Double * upperBoundPtr, 
	           Double * lowerBoundPtr,
	           const Double & defaultDualVal4UbConstr,
               const Double & defaultDualVal4LbConstr, 
	           Problem * problemPtr);

  virtual ~ColGenSpConf();

  /**
   * Return true if a feasible primal solution yielding a negative reduced cost col was found, false otherwise
   */
  virtual MastColumn * checkColumn4Insertion(MastColumn * colPtr, bool inPurePhaseOne, const int & insertionLevel = 2);
  virtual InstanciatedConstr * checkConstraint4Insertion(Constraint * constrPtr, const int & insertionLevel = 2);
  virtual InstanciatedConstr * checkConstraint4Insertion(InstanciatedConstr * iconstrPtr,
                                                         const int & insertionLevel = 2);
  virtual Constraint * castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately = false);
  virtual InstanciatedConstr * castAndAddConstraint(InstanciatedConstr * iconstrPtr,
                                                    const bool & insertImmediately = false);
  virtual Variable * castAndAddVariable(Variable * varPtr, const bool & insertImmediately = false);
  virtual InstanciatedVar * castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately = false);

  virtual const bool & spConfHasClassInducedSpVarBounds() const;
  virtual void spConfHasClassInducedSpVarBounds(const bool & f);
  int genNewColumns(bool currentlyPerformingPhaseI, int & maxLevelOfSubProbRestriction, int & allAddedColumns,
                    int & addedNegRedCostColumns);
  bool performSpRelaxationLigntening(const bool & masterConverged, const int & callMode);
  virtual bool cannotGenerateAnyMoreCol();
  virtual bool needNotGenerateAnyMoreCol();
  virtual bool cardinalityIsFixed();


  virtual MastColumn * recordSubproblemSolution(Solution * spSolPtr,
                                                bool inPurePhaseOne,
                                                const int & insertionLevel = 2,
                                                Solution * masterSolPtr = NULL,
                                                bool changeEnumeratedFlag = false);

  virtual int insertConstraintsInMaster();
  virtual void insertColumnsInMaster(int & allAddedColumns, int & addedNegRedCostColumns);
  virtual int insertAllColumnsInMaster();
  void clearColPtrList4Insertion();
  virtual void prepareProbConfig();
  virtual Double * upperBoundPtr() const ;
  virtual Double * lowerBoundPtr() const ;
  virtual void upperBoundPtr(Double * ubPtr);
  virtual void lowerBoundPtr(Double * lbPtr);
  virtual MissingColumn * misColPtr() const;
  virtual MasterConf * mastConfPtr() const;
  virtual void misColPtr(MissingColumn * misCptr);
  virtual void lowerBoundMastConstrPtr(InstMastConvexityConstr * mccPtr);
  virtual void upperBoundMastConstrPtr(InstMastConvexityConstr * mccPtr);
  virtual InstMastConvexityConstr * lowerBoundMastConstrPtr()  const;
  virtual InstMastConvexityConstr * upperBoundMastConstrPtr()  const;
  virtual const Double & defaultDualVal4UbConstr() const;
  virtual const Double & defaultDualVal4LbConstr() const;
  virtual Solution * solvePC(int & maxLevelOfSubProbRestriction, bool inPurePhaseOne);
  virtual void getSolPC(bool inPurePhaseOne);

  /**
   * Retrieve implicit solution \delta,
   * compute mult and dualBoundContrib
   */
  void computeSpDualBoundContrib();
  virtual const Double & fixedCost()  const;
  virtual void fixedCost(const Double & b);
  virtual const Double & fixedDualCost()  const;
  virtual void fixedDualCost(const Double & b);
  virtual const Double & spRootReducedCost() const;

  virtual const Bound & target()  const;
  virtual void target(const Bound & t);
  virtual const Double & mult() const;

  /**
   * Set var cost, target, minCost,
   */
  bool updateConf(bool inPurePhaseOne);
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  virtual void toSplit() { _toSplit = true; } /// added by Ruslan
  
  void setRollbackPointSavedStatus(const bool & status) {_rollbackPointSaved = status;}
  
  const bool & getRollbackPointSavedStatus() {return _rollbackPointSaved;}
  
  virtual void priorityLevel(const Double & prLevel)
  {
    _priorityLevel = prLevel;
  }

  virtual const Double & priorityLevel() const
  {
    return _priorityLevel;
  }
  
private:
  MastColumn * _curBestMastColumnPtr;
  MastColumn * _incBestMastColumnPtr;

public:
  MastColumn * curBestMastColumnPtr() /// added by Artur
  {
    return _curBestMastColumnPtr;
  }
  /// start added by Ruslan
  MastColumn * incBestMastColumnPtr() /// added by Artur
  {
    return _incBestMastColumnPtr;
  }
  void saveBestMastColumnPtrAsIncumbentBest()
  {
    _incBestMastColumnPtr = _curBestMastColumnPtr;
  }

  virtual bool isTypeOf(const PcId::PcIdentifier & pcIdentifier) const
  {
    return compareIdentifier(PcId::ColGenSpConfMask, pcIdentifier);
  }

};

inline std::ostream& operator<<(std::ostream& os, const ColGenSpConf & that)
{
  return that.print(os);
}


#endif /* BCCOLGENSPCONFC_H_ */
