/**
 *
 * This file bcMastColumnC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef MastColumnC_h
#define MastColumnC_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp" 

class Node;
class ProbConfig;
class MasterConf;
class ColGenSpConf;

/**
 * A base class for Master Variables 
 * associated to colGen subproblem solutions 
 */
class MastColumn: public AggregateVariable, public Variable
{
  long _mcref;
  int _treatOrderId; /// treat order of the node where the column has been generated (needed for problem setup)
  int _nodeTreatOrder; /// treat order of the node on which this column has been generated

  /**
   * Where one can find amongst other  thing upper and lower bound on variable
   */
  ColGenSpConf * _cgSpConfPtr;

 public:
  
    
  MastColumn(MasterConf * masterConfPtr, 
	         ColGenSpConf * cgSpConfPtr,
	         Solution * spSol = NULL,
	         const std::string & name = "MC");
  
  MastColumn(const MastColumn & that);
  virtual ~MastColumn();
  virtual void fillAggregateSol(VarPtr2DoubleMap & curAggregateMastSol, const Double & val) const;
  virtual void fillMapOfIntSpVar(VarPtr2DoubleMap & curSolMap, const Double & val ) const;
  virtual void fillSpIndicatorMap(std::map< Variable *,  
					              std::multiset< std::pair< MastColumn *, ValueRecord > ,
							      SortMastColPerDecreasingSpVal>, Variable::PrioritySort > & curMap,
				                  const ValueRecord & valRec,
                                  const VarPtrSet & candVarForCompBoundBranching);
  virtual void fillSpIndicatorColMultiset(std::multiset< std::pair < MastColumn *, ValueRecord > ,
                                          SortMastColPerDecreasingSpVal > & curMultiSet,
					                      Variable * spVarPtr,
					                      const ValueRecord & valRec);

  /**
   * Measures the variable contribution to the satisfaction of member constraints (for use in greedy heuristic)
   */
  virtual const Double contrib(const Double & use = 1);
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void recordInMembershipOfSubProbVar();
  virtual bool defaultComputeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef defaultComputeCoef(ConstVarConstrConstPtr vcPtr);
  virtual bool suitableToFixValue(const Double & value); /// added by Ruslan
  int maxValueInCurrentMasterProblem(); /// added by Ruslan
  bool suitableForResidualProb(const Double & useLevel = 1.0);

  virtual Solution * spSol() const
  {
    return _spSol;
  }

  void treatOrderId(const int value )
  {
    _treatOrderId = value;
  }

  const int treatOrderId() const
  {
    return _treatOrderId;
  }
  
  virtual const long & mcref() const
  {
    return _mcref;
  }

  virtual ColGenSpConf * cgSpConfPtr() const;

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual std::ostream& printColVector(std::ostream& os = std::cout) const;

  /**
   * Also do behind the scene membership for MasterConstr and SubProblemVariables
   */
  virtual void setMembership();
  
  /**
   * Assumes spSol var have been set before
   */
  virtual void agvSetMembership(Variable * varPtr);
  virtual const Double & curCost() const;
  virtual const Double & costrhs() const;
  virtual void costrhs(const Double & newCostrhs) {return Variable::costrhs(newCostrhs);}

  const Double spVarVal(Variable * varPtr) const;
  bool spVarCount(Variable * varPtr) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;

  virtual bool operator<(const VarConstr& b) const;
};

struct SortMastColumnPerNonDecreasingRedCost
{

  bool operator()(const MastColumn * a, const MastColumn * b) const
  {
    if (a->reducedCost() < b->reducedCost())
      return true;
    if (a->reducedCost() > b->reducedCost())
      return false;
    return (a->mcref() < b->mcref());
  }
};


/**
 * A base class for Aritificial Master Variables  replacing  a missing subproblem solutions
 */
class MissingColumn: public MastColumn, public ArtificialVar
{
 public:
  MissingColumn(MasterConf * masterConfPtr, ColGenSpConf * cgSpConfPtr, const Double & artVarCost);
  virtual ~MissingColumn();
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void addMember(VarConstr * vcPtr);

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

inline std::ostream& operator<<(std::ostream& os, MastColumn & that)
{
  return that.print(os);
}

struct LexicographicMastColValSorting
{
  bool operator()(const std::pair < MastColumn *, ValueRecord > & a, const std::pair < MastColumn *, ValueRecord > & b) const
  { return LexicographicallyST(a.first , b.first);}

};

struct MasterColSolution : public  std::list< std::pair < MastColumn *, ValueRecord > >
{
  void push_back(Variable * varPtr, const ValueRecord  & rec);
};


#endif // MastColumnC_h
