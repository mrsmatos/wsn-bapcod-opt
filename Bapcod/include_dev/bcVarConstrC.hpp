/**
 *
 * This file bcVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef VarConstrClasses_h
#define VarConstrClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcIndexC.hpp"
#include "bcPrintC.hpp"
#include "bcBapcodInit.hpp"
#include "bcVcIndexInProbC.hpp"
#include "bcVcIdentifierC.hpp"

#if(BOOST_VERSION >= 104000)
#include <boost/unordered_map.hpp>
#else
//#include <ext/hash_map>
#include <unordered_map>
#endif //BOOST_LIB_VERSION
class BaseFormulation;
class Constraint;
class VarConstrStabInfo;
class Variable;
class ValueRecord;
class VarConstr;
class Model;
class Problem;
class OvfVar;
class OvfConstr;

struct VarConstrSort /// equivalent to VcRefSorting
{
  bool operator()(const VarConstr * a, const VarConstr * b) const;
};

class GenericVarConstr;
class GenericVar;
class GenericConstr;
class InstanciatedVar;
class InstMasterConstr;
class InstMasterVar;
class MastColumn;
class InstanciatedConstr;
class NonLinearInstConstr;
class MasterConstr;
class SubProbVariable;
class Solution;
class CompSetInstMastBranchConstr;
class AggregateVariable;
class MastColumn;

typedef VarConstr * VarConstrPointer;
typedef Variable* VarPointer;
typedef VarConstr * ConstVarConstrConstPtr;

typedef VarConstr * VCPTR;

class Cpr3
{
public:

  bool operator()(const VCPTR &key1, const VCPTR &key2) const
  {
    /// Equality of two keys
    return (key1 == key2);
  }

  size_t operator()(const VCPTR &key) const
  {
    /// Explanation: http://stackoverflow.com/q/20953390/1200528
    static const size_t shift = (size_t)log2(1 + sizeof(VCPTR));
    /// Hash of a key
    return (size_t)(key) >> shift;
  }
};


class Cpr2
{
public:

  bool operator()(const VCPTR &key1, const VCPTR &key2) const
  {
    /// Equality of two keys
    return (key1 == key2);
  }

  size_t operator()(const VCPTR &key) const
  {
    unsigned long z1 = (unsigned long) key;

    /// Hash of a key
    return z1;
  }
};

#if (BOOST_VERSION >= 104000) //#if(0) //BOOST_VERSION
typedef boost::unordered_map<VCPTR, Double, Cpr3> ConstVarConstrPtr2Double;
#else
#ifndef _MSC_VER
// Only available on Unix
typedef __gnu_cxx::hash_map<VCPTR, Double, Cpr2> ConstVarConstrPtr2Double;
#else
// Only available on Windows
typedef stdext::hash_map<VCPTR, Double> ConstVarConstrPtr2Double;
#endif //_MSC_VER
//
//#else
//typedef std::map < VarConstr *, Double, VarConstrSort > ConstVarConstrPtr2Double;
#endif

typedef std::map<MastColumn *, Double, VarConstrSort> MapMastColumnPtr2Double;
typedef std::map<SubProbVariable *, Double, VarConstrSort> MapSubProbVariablePtr2Double;

typedef std::set<VarConstr *> PlainVarConstrPointerSet;

typedef std::set<Constraint *, VarConstrSort> ConstrPtrSet;
typedef std::set<InstMasterVar *, VarConstrSort> InstMasterVarPtrSet;
typedef std::set<InstMasterConstr *, VarConstrSort> InstMasterConstrPtrSet;
typedef std::set<VarConstr *, VarConstrSort> VarConstrPtrSet; // not well defined since var and constr numbering are independent
typedef std::set<Variable *, VarConstrSort> VarPtrSet;
typedef std::set<CompSetInstMastBranchConstr *, VarConstrSort> CompSetInstBrConstrPtrSet;

typedef VarConstrIndexManager<Constraint> ConstrIndexManager;
typedef VarConstrIndexManager<Variable> VarIndexManager;

typedef std::list<Variable*> VarPtrList;
typedef std::list<Constraint*> ConstrPtrList;

typedef std::map<Constraint *, Double, VarConstrSort> ConstrPtr2DoubleMap;
typedef std::map<CompSetInstMastBranchConstr *, Double, VarConstrSort> CompSetConstrPtr2DoubleMap;
typedef std::map<Variable *, Double, VarConstrSort> VarPtr2DoubleMap;

class VarConstr : public VcIndexInProb
{
private:

  int _VCref;

  std::string _name;

  /// 'U' or 'D'
  char _directive;

  /**
   * A higher priority means that variable
   * is selected first for branching or diving
   */
  Double _priority;

protected:
  Model * _modelPtr;

  /**
   * Cost for a variable,
   * rhs for a constraint
   */
  Double _costrhs;

  /**
   * Variables:
   * sense : 'P' = positive
   * sense : 'N' = negative
   * sense : 'F' = free
   */
  /**
   * Constraints: 'G', 'L', or 'E';
   */
  char _sense;
  int _sign;

  /** _type
   * For Variables:
   * 'C' = continuous,
   * 'B' = binary,
   * or 'I' = integer;
   * For Constraints:
   * type = 'C' for core (required for the IP formulation),
   * type = 'F' for facultative (only helpfull for the LP formulation),
   * type = 'S' for constraints defining a subsystem in column generation for extended formulation approach
   * type = 'M' for constraints defining a pure master constraint
   * type = 'X' for constraints defining a subproblem convexity constraint in the master
   */
  char _type;

  /** _kind
   * For Constraints:
   * 'E' = explicit,
   * 'I' = implicit,
   * 'R' = delayed constraints
   * Added only after solving the
   * LP root node without them
   */
  char _kind;

  /** _flag
   * 's' (by  default) for static VarConstr belonging
   * to the problem (and erased when the problem is erased);
   * 'd' for generated dynamic VarConstr
   *  not belonging to the problem
   * 'a' for artificial VarConstr.
   */
  char _flag;

  bool _inCurProb;
  bool _inCurForm;

  /**
   * Index number of entity in formulation
   */
  int _index;
  MultiIndex _multiIndex;

  /**
   * Primal Value for a variable,
   * dual value a constraint;
   */
  Double _val;

  /**
   * Temprarity recorded primal Value for a variable
   * (dual value a constraint) used for functions returning fract part
   */
  Double _tmpVal;
  /**
   * Temprarity recorded primal Value for a variable in rounding or fixing heuristic
   * (dual value a constraint) used for functions returning fract part
   */
  Double _challengerRoundedValue;

  /**
   * Temprarity recorded primal Value for a variable
   * (dual value a constraint) used for functions of Solution
   */
  Double _solVal;

  /**
   * incumbent value
   */
  Double _incumbentVal;

  /**
   * To represent global upper bound
   * on variable primal value
   * or on constraint dual value
   */
  Double _upperBound;

  /**
   * To represent global lower bound
   * on  variable primal value
   * or on constraint dual value
   */
  Double _lowerBound;

  /**
   * To represent global upper bound
   * on sp variable primal value
   */
  Double _globalUb;

  /**
   * To represent global lower bound
   * on sp variable primal value
   */
  Double _globalLb;
  Double _memorisedCurUb;
  Double _memorisedCurLb;
  Double _curUb;
  Double _curLb;
  Double _mult;

  bool _presetMembership;
  bool _membershipUpToDate;
  bool _buildMembershipHasBeenPerformed;
  ConstVarConstrPtr2Double _member2coefMap;
  VarConstrPtrSet _nonMemberSet;
  Problem * _problemPtr;

  bool _infoIsUpdated; /// added by Ruslan, needed for VarConstrResetInfo
  bool _inPreprocessedList; /// added by Ruslan, needed for preprocessing
  
  int _participation; //added by Issam, to know when we can delete a dynamic variable.
  /**
   * Temporary record updated by function
   * virtual const Double & violation(const Double & lhs);
   */
  Double _violation;

  Double _reducedCost;

  /// to hold Info Need for stabilisation of constraint in Col Gen approach
  VarConstrStabInfo * _stabInfoPtr;

  std::list<OvfVar *> _ovfVarList;
  std::list<OvfConstr *> _ovfConstrList;

  /// added by Boris: cache of membership in constraints
  std::vector<std::pair<VCPTR,double> > _cachedMemberShip;

  VarConstr() = delete;
  VarConstr(const VarConstr & vc) = delete;

public:

  VarConstr(const VarConstr & vc, const int & ref);

  /**
   * @param  costrhs Cost for a variable, rhs for a cosntraint
   * @param   sense  Constraints: 'G', 'L',
   * or 'E'; Variables 'P' = positive, 'N' = negative, 'F' = free;
   * @param  type For Variables: 'C' = continuous, 'B' = binary, or 'I' = integer;
   * for a constraint 'C' for core (required for the IP formulation,
   * 'F' for facultative (only helpfull for the LP formulation).
   * @param  kind  For Constraints: 'E' = explicit, 'I' = implicit,
   * 'R' = delayed constraints added only after solving the LP root node without them
   * @param flag 's' (by  default) for static VarConstr belonging to the problem
   * (and erased when the problem is erased);
   * 'd' for generated dynamic VarConstr not belonging to the problem
   * @param val Indinital primal Value for a variable, dual value a constraint;
   * @param upperBound Represents a global upper bound on variable primal value
   * or on constraint dual value
   * @param  _lowerBound representsglobal lower bound
   * on variable primal value or on constraint dual value
   * @param globalUb represents a global upper bound on sp variable primal value
   * @param globalLb represents a global lower bound on sp variable primal value
   * @param directive for branching
   * @param priority for branching or separation (the higher the most mpriority)
   **/
  VarConstr(Model * modelPtr, const long & ref, const std::string & name,
            const Double & costrhs, const char & sense, const char & type = 'I',
            const char & kind = 'E', const char & flag = 's', const int & index = -1,
            const Double & val = 0, const Double & upperBound = BapcodInfinity,
            const Double & lowerBound = -BapcodInfinity, const Double & globalUb =
            BapcodInfinity, const Double & globalLb = -BapcodInfinity,
            const char & directive = 'U', const Double & priority = 1.0,
            const bool & presetMembership = true);

  /// Includes clearMembership()
  virtual ~VarConstr();

  Model * modelPtr() const
  {
    return _modelPtr;
  }

  virtual const MultiIndex & multiIndex() const
  {
    return _multiIndex;
  }

  const int & ref() const
  {
    return _VCref;
  } 

  const int & VCref() const
  {
    return _VCref;
  }
  const std::string & name() const;
  void name(const std::string & cname);
  void name(char * cname);

  virtual const std::string & genericName() {return name();}
  virtual const char & sense() const;
  virtual const int & sign() const;
  virtual const char & sense(const char & sense);
  virtual const char & type() const;
  virtual const char & kind() const;
  virtual const char & type(const char & s);
  virtual const char & kind(const char & s);
  virtual const char & flag() const;
  virtual void flag(const char & cflag);

  virtual const int & index() const;
  virtual const int & index(const Double & newIndex);
  virtual const Double & priority() const;
  virtual const char & directive() const;
  virtual void directive(const char & dir);
  virtual void priority(const Double & prior);
  virtual const Double & costrhs() const;
  virtual const Double & rhs() const;
  virtual const Double cost() const;
  virtual void costrhs(const Double & newCostrhs);
  virtual void resetCostFromDefaultCost(const Double & factor = 0.2);
  virtual const Double curRhs() const;

  /**
   * Active varConstr are in the implicit formulation
   * and coonsider by preprocessor as part of the current formulation
   */
  virtual const bool inCurProb();

  /**
   * inCurForm varConstr are those set for the explicit
   * formulation and later included in the current formulation
   */
  virtual const bool & inCurForm() const;
  virtual void activate();
  virtual void desactivate();
  virtual void setInForm();
  virtual void unsetInForm();
  virtual const Double & val() const;
  virtual void val(const Double & newVal);
  virtual const Double & tmpVal() const;
  virtual void tmpVal(const Double & newVal);
  virtual const Double & solVal() const;
  virtual void solVal(const Double & newVal);
  virtual const Double fracPart() const;
  virtual const Double fracPartRelativeToTarget() const;
  virtual const Double lFracPart() const;
  virtual const Double uFracPart() const;
  virtual const Double & globalUb() const;
  virtual void globalUb(const Double & ub);
  virtual const Double & globalLb() const;
  virtual void globalLb(const Double & lb);
  virtual const Double & ub() const;
  virtual void ub(const Double & ub);
  virtual const Double & lb() const;
  virtual void lb(const Double & lb);

  virtual const Double & curUb() const
  {
    return _curUb;
  }

  virtual const Double & curLb() const
  {
    return _curLb;
  }

  /**
   * Differs from local bound only for SubproblemVariable
   */
  virtual const Double & curGlobUb() const
  {
    return _curUb;
  }

  virtual const Double & curGlobLb() const
  {
    return _curLb;
  }

  virtual void curUb(const Double & ub)
  {
    _curUb = ub;
  }

  virtual void curLb(const Double & lb)
  {
    _curLb = lb;
  }

  /**
   * Measures the variable min contribution
   * to the satisfaction of a specific  constraint
   */
  virtual const Double lhsMinContrib(ConstVarConstrConstPtr vcPtr);

  /**
   * Measures the variable min contribution
   * to the satisfaction of a specific  constraint
   */
  virtual const Double lhsMaxContrib(ConstVarConstrConstPtr vcPtr);
  virtual const Double & mult() const;
  virtual const Double & incumbentVal() const;
  virtual void incumbentVal(const Double & newVal);
  virtual void mult(const Double & m);

  virtual void problemPtr(Problem * pPtr)
  {
    _problemPtr = pPtr;
  }

  virtual Problem * problemPtr() const
  {
    return _problemPtr;
  }
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr) = 0;

  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr) = 0;

  virtual ProbConfig * probConfPtr() const
  {
    return NULL;
  }

  virtual ConstVarConstrPtr2Double & member2coefMap();
  virtual const ConstVarConstrPtr2Double & member2coefMap() const;

  // virtual const VarConstrPtrSet & member() const;
  virtual const Double & includeAsMember(VarConstr * vcPtr, 
					 const Double & coef, 
					 const bool & cumulativeCoef);
  virtual void eraseAsMember(VarConstr * vcPtr);

  /**
   * Check whether this has non zero coef in vcPtr,
   * if so, set insert this in vcPtr.member and vice versa
   */
  virtual void addMember(VarConstr * vcPtr);

  /// Insert this in vcPtr.member and vice versa
  virtual const Double & includeMember(VarConstr * vcPtr, 
				       const Double & coef, 
				       const bool & cumulativeCoef);

  /// Remove this from vcPtr.member and vice versa
  virtual void recordNonMember(VarConstr * vcPtr);
  virtual const Double & upToDateMembCoef(ConstVarConstrConstPtr vcPtr) const;
  virtual void clearMembership();
  virtual void deleteFromProb();
  virtual void addToProb(Problem * probPtr);

  virtual const VarConstrPtrSet & nonMemberSet() const
  {
    return _nonMemberSet;
  }

  /**
   * Also do behind the schene membership for MasterConstr and SubProblemVariables
   */
  virtual void setMembership();

  virtual const bool & presetMembership() const
  {
    return _presetMembership;
  }

  virtual void presetMembership(const bool & flag)
  {
    _presetMembership = flag;
  }

  inline const bool & membershipUpToDate() const //Issam this method is not redefined so i removed the virtual to allow inline
  {
    return _membershipUpToDate;
  }

  virtual const bool & buildMembershipHasBeenPerformed() const
  {
    return _buildMembershipHasBeenPerformed;
  }

  virtual void buildMembershipHasBeenPerformed(const bool & flag)
  {
    _buildMembershipHasBeenPerformed = flag;
  }

  void incrParticipation(int where = -1);
  
  void decrParticipation(int where = -1);
  
  int participation()
  {
    return _participation;
  }

  /// start added by Ruslan
  
  /// for inBasis mark we use _infoIsUpdated for the moment to use less memory
  void isInBasis(const bool value)
  {
    _infoIsUpdated = value;
  }
    
  const bool isInBasis() const
  {
    return _infoIsUpdated;
  }

  void infoIsUpdated(const bool value)
  {
    _infoIsUpdated = value;
  }

  const bool infoIsUpdated() const
  {
    return _infoIsUpdated;
  }

  const bool inPreprocessedList() const
  {
    return _inPreprocessedList;
  }

  void removeFromPreprocessedList()
  {
    _inPreprocessedList = false;
  }

  virtual void addToPreprocessedList()
  {
  }

  /// end added by Ruslan
  const std::list<OvfVar *> & ovfVarList() const
  {
    return _ovfVarList;
  }

  const std::list<OvfConstr *> & ovfConstrList() const
  {
    return _ovfConstrList;
  }
  /**
   * Temporary record updated by function
   * virtual const Double & violation(const Double & lhs);
   */
  const Double & violation() const;
  void violation(const Double & viol);

  virtual const Double & reducedCost() const;
  virtual const Double & computeReducedCost();

  /**
   * Measures the variable contribution
   * to the satisfaction of member constraints
   * (for use in greedy heuristic)
   */
  virtual const Double contrib(const Double & use = 1)
  {
    return 1;
  }

  virtual const Double greedyCost();

  BapcodInit & bapcodInit() const;
  const ControlParameters & param() const;
  ControlParameters & param();

  virtual void reducedCost(const Double & newRedCost);

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual bool operator<(const VarConstr & that) const;
  virtual bool operator==(const VarConstr & that) const;
  virtual bool operator!=(const VarConstr & that) const;

  virtual GenericVarConstr * genVarConstrPtr() const
  {
    return NULL;
  }


  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;

  VarConstrStabInfo * stabInfoPtr() const
  {
    return _stabInfoPtr;
  }
  
  void deleteStabInfoPtr();
  
  const Double & valOrSepPointVal() const;

};

inline std::ostream& operator<<(std::ostream& os, const VarConstr * that)
{
  return that->print(os);
}

inline std::string getDebugInfo(const VarConstr * const a)
{
  std::stringstream def;
  def << "rhs " << a->rhs() << " - ";
  ConstVarConstrPtr2Double::const_iterator it1 = a->member2coefMap().begin();
  for (; (it1 != a->member2coefMap().end()); ++it1)
    {
      def << "(" << it1->first->VCref() << "," << it1->second << ") - ";
    }
  return def.str();
}


inline bool LexicographicallyST(const VarConstr * const a,
                                const VarConstr  * const b) 
{
  return (a->operator<(*b));
}

inline bool VcRefST(const VarConstr * const a, const VarConstr * const b)
{
  return (a->VCref() < b->VCref());
}

struct LexicographicSorting
{
  bool operator()(const VarConstr * const a, const VarConstr * const b) const
  { return LexicographicallyST(a,b);}

};

class Variable : public VarConstr
{
protected:
  Double _memorisedCurCost;
protected:
  bool _isCandForBranching; 

  /**
   * Holds variables of type AggregateVar to the SpSol of which spVar belongs
   */
  VarPtrSet _setOfAggregateVarToWhichItBelongs;
  CompSetConstrPtr2DoubleMap _mapCompSetBrConstr2upperBd;
  CompSetConstrPtr2DoubleMap _mapCompSetBrConstr2lowerBd;

  Variable() = delete;

public:

  struct PrioritySort
  {
    bool operator()(const VarPointer & a, const VarPointer & b) const;
  };

  /**
   * @param  name string defining a name for the varaible
   * @param  costrhs Cost for the variable
   * @param   sense  'P' = positive, 'N' = negative, 'F' = free;
   * @param  type  'C' = continuous, 'B' = binary, or 'I' = integer;
   * @param  kind  unsed for Variables
   * @param flag 's' (by  default) for  Var belonging to the problem
   * (and erased when the problem is erased);
   * 'd' for generated  Var not belonging to the problem
   * @param val Indinital primal Value for a variable, dual value a constraint;
   * @param upperBound Represents a global upper bound on variable primal value
   * or on constraint dual value
   * @param  lowerBound represents global lower bound
   * on variable primal value or on constraint dual value
   * @param globalUb represents a global upper bound on sp variable primal value
   * @param globalLb represents a global lower bound on sp variable primal value
   * @param directive for branching
   * @param priority for branching or separation Priorities must be integer within [0,1000];
   * lower priority index means chosen first for branching
   *
   * @return
   */
  Variable(Model * modelPtr, 
	       const std::string & name,
	       const Double & costrhs = 0,
           const char & sense = 'P', 
	       const char & type = 'I',
	       const char & kind = 'E',
           const Double & upperBound = BapcodInfinity,
	       const Double & lowerBound = 0,
	       const char & flag = 's',
	       const char & directive = 'U',
           const Double & priority = 1.0, 
	       const Double & val = 0,
           const Double & globalUb = BapcodInfinity,
	       const Double & globalLb = 0,
           const bool & presetMembership = true, 
	       const int & index = -1);
  Variable(const Variable & that);
  virtual ~Variable();

  /**
   * Virtual copy constructor of type needs
   * to be defined for derived classes
   */
  bool isCandForBranching() const {return _isCandForBranching;}
  void isCandForBranching(bool flag) {_isCandForBranching = flag;}

  /**
   * Dummy fonctions to avoid upcasting
   */
  virtual ColGenSpConf * cgSpConfPtr() const {return NULL;}
  virtual GenericVar * genVarPtr() const  {return NULL;}
  virtual Solution * spSol() const  {return NULL;}
  virtual void fillAggregateSol(VarPtr2DoubleMap & curAggregateMastSol, const Double & val) const{return;}

  /**
   * Measures the variable contribution to the satisfaction
   * of member constraints (for use in greedy heuristic)
   */
  virtual const Double contrib(const Double & use = 1);

  virtual const Double cost() const;
  virtual const Double & curCost() const;
  virtual const Double & costrhs() const;
  virtual void costrhs(const Double & newCostRhs);
  virtual void resetCost(const bool & inPurePhaseOne);
  virtual void resetCurCostByValue(const Double & newcost);

  /// start added by Ruslan for new var constr reset and preprocessing

  // note that for master vars global and local bounds are the same
  // The following methods are thus overwritten in bcSpVarConstrC
  virtual void globalCurUb(const Double & ub)
  {    
    _curUb = _memorisedCurUb = ub;
  }

  virtual void globalCurLb(const Double & lb)
  {
    _curLb = _memorisedCurLb = lb;
  }

  virtual const Double & globalCurUb() const
  {
    return _memorisedCurUb;
  }

  virtual const Double & globalCurLb() const
  {
    return _memorisedCurLb;
  }

  virtual void globalUb(const Double & ub)
  {
    _upperBound = _curUb = _memorisedCurUb = ub;
  }

  virtual void  globalLb(const Double & lb)
  {
    _lowerBound = _curLb = _memorisedCurLb = lb;
  }

  virtual const Double & globalUb() const
  {
    return _upperBound;
  }

  virtual const Double & globalLb() const
  {
    return _lowerBound;
  }

  virtual void addToPreprocessedList();

  virtual bool boundsAreDifferentFromDefaultInForm();

  virtual void resetBoundsAndCostToDefaults();
    
  virtual bool activateVariable(bool toSetInForm);

  virtual bool desactivateVariable(const VcIndexStatus::VcStatus & status, bool toUnsetInForm = true);

  /// end added by Ruslan for new var constr reset and preprocessing

  /**
   * Sp Var Default Component curUb in SP:
   * specific def only for SP var
   */
  virtual const Double & defaultSpVarUb() const
  {
    return curUb();
  }

  /**
   * Sp Var Default Component curLb  in SP:
   * specific def only for SP var
   */
  virtual const Double & defaultSpVarLb() const
  {
    return curLb();
  }

  /**
   * Also do behind the schene membership
   * for MasterConstr and SubProblemVariables
   */
  virtual void setMembership();

  /**
   * Put in _member2coefMap the subset
   * of constrSet to which Var belongs and
   * Includes Var in those var member set
   */
  virtual void setMembership(const ConstrIndexManager & constrSet);

  virtual void clearMembership();

  /**
   * Set membership by enumerating problem constraints
   */
  virtual void enumerativeSetMembership();

  virtual void includeAggregateVarAsMember(Variable * mvPtr)
  {
    _setOfAggregateVarToWhichItBelongs.insert(mvPtr);
  }

  virtual void eraseAggregateVarAsMember(Variable * mvPtr)
  {
    _setOfAggregateVarToWhichItBelongs.erase(mvPtr);
  }

  virtual const VarPtrSet & setOfAggregateVarToWhichItBelongs()
  {
    return _setOfAggregateVarToWhichItBelongs;
  }

  virtual void recallMemorisedBounds();

  virtual bool suitableToFixValue(const Double & value);
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);

  /**
   * Check whether this has non zero coef in vcPtr,
   * if so, set insert this in vcPtr.member and vice versa
   */
  virtual void addMember(VarConstr * vcPtr);

  /**
   * Ub < Lb ?
   */
  virtual bool infeasible() const;

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  virtual bool greedySmallerThan(Variable * that)
  {
    return (greedyCost() < that->greedyCost());
  }

  virtual void includeOvfPtr(OvfVar * ovfPtr);

  /**
   * To memorize Branching constraints that imply a bound on variable
   */
  CompSetConstrPtr2DoubleMap & mapCompSetBrConstr2upperBd()
   {
     return _mapCompSetBrConstr2upperBd;
   }

   CompSetConstrPtr2DoubleMap & mapCompSetBrConstr2lowerBd()
   {
     return _mapCompSetBrConstr2lowerBd;
   }

  /**
   * Test if Variable is candidtate for behing chosen next to varPtr when branching on component bound set
   */
  virtual bool consecutive2varWhenBrOnCBS(InstanciatedVar * varPtr);

  /**
    * Test if the current class is type of vcIdentifier.
    * @param vcIdentifier
    * @return true if current class is type of vcIdentifier, false otherwise.
    */
   virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

struct SolutionVarInfo
{
public:
  Variable * varPtr;
  Double value;
  Double reducedCost;
  Double incumbentValue;
  Double cost;
  Double priorityLevel;
  bool canRoundDown;
  bool canRoundUp;

  SolutionVarInfo(Variable * varPtrV);
  SolutionVarInfo(const SolutionVarInfo & that);

  virtual ~SolutionVarInfo();
};

typedef std::list<SolutionVarInfo *> SolutionVarInfoPtrList;

/**
 * Generic functions required for generic subroutines  
 */
struct GreedyOrderingOfVar
{

  virtual bool operator()(Variable * a, Variable * b)
  {
    return (a->greedySmallerThan(b));
  }

  virtual ~GreedyOrderingOfVar() { }
};

struct SortVariablePtrByIncreasingRedCost
{

  bool operator()(const Variable * a, const Variable * b)
  {
    return (a->reducedCost() < b->reducedCost());
  }
};

struct SortColumnsForLocalSearchHeurSol
{

  bool operator()(const Variable * a, const Variable * b)
  {
    if (a->tmpVal() < b->tmpVal())
      return true;
    if (a->tmpVal() > b->tmpVal())
      return false;
    return (a->ref() < b->ref());
  }
};

/// To model a side behavior of Variable

class AggregateVariable
{
protected:
  /**
   * Solution belongs to AggregateVariable
   * and is deleted by AggregateVariable
   */
  Solution * _spSol;
  Variable * _varPtr;

public:
  AggregateVariable(Solution * spSoly);

  /**
   * Need to transmit the own copy of solution
   * that shall belong to AggregateVariable
   */
  AggregateVariable(const AggregateVariable & that);
  virtual ~AggregateVariable();
  virtual bool operator<(const VarConstr & that) const;



  virtual void setAggregateVariable(Variable * varPtr);
  virtual void unsetAggregateVariable();

  /// Assumes spSol var have been set before
  virtual void agvSetMembership(Variable * varPtr);
  virtual bool agvComputeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef agvComputeCoef(ConstVarConstrConstPtr vcPtr);
  virtual std::ostream & printAgv(std::ostream& os) const;
  void includeVar(Variable * aggVarPtr, Variable * varPtr, const Double & val, const bool & cumulativeVal); /// added by Ruslan

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class ArtificialVar
{
protected:
  Double _defaultCost;
public:
  ArtificialVar(const Double & defaultCost = 1);
  virtual ~ArtificialVar();
  virtual const Double & defaultCost();
  virtual void setCost(const Double & cost);

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class GlobalArtificialVar : public Variable, public ArtificialVar
{
  char _senseType; /// 'G' for positive coefficient in 'G' (>=) constraints,
  /// 'L' for  for positive coefficient in 'L' (<=) constraints
public:
  GlobalArtificialVar(Model * modelPtr, 
		              const Double & cost,
                      const char & senseType = 'G', 
		              const BcObjStatus::MinMaxIntFloat & minmax = BcObjStatus::minInt,
		              const std::string & name = "globArtVar");
  virtual ~GlobalArtificialVar();
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void resetCostFromDefaultCost(const Double & factor = 0.2);
  virtual const Double & curCost() const;
  virtual const Double & costrhs() const;

  virtual void costrhs(const Double & newCostrhs)
  {
    return Variable::costrhs(newCostrhs);
  }

  const char & senseType() const
  {
    return _senseType;
  }

  /**
   * Put in _member2coefMap the subset
   * of constrSet to which Var belongs and
   */
  //virtual void setMembership(const ConstrPtrSet & constrSet);
  virtual void setMembership(const ConstrIndexManager & constrSet);

  virtual void setMembership(const ConstrIndexManager & constrSet, const VcIndexStatus::VcStatus& status, char flag);

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class LocalArtificialVar : public Variable, public ArtificialVar
{
public:
  enum LocalArtClassId ///added by Ruslan
  {
    PosLocalId,
    NegLocalId,
    PosOuterId,
    NegOuterId,
    PosInnerId,
    NegInnerId
  };

private:
  Constraint * _constrPtr;
  Double _tpVal;
  LocalArtClassId _locClassId; ///added by Ruslan

public:
  LocalArtificialVar(Constraint * constrPtr, 
                     const LocalArtClassId locClassId, ///added by Ruslan
		     const BcObjStatus::MinMaxIntFloat & minmax = BcObjStatus::minInt, 
		     const std::string & name = "la", 
		     const Double & defCost = 1.0,
		     const Double & ub = BapcodInfinity);
  virtual ~LocalArtificialVar();
  virtual void setConstraintPtr(Constraint * constrPtr);
  inline bool canBeDeleted() const
  {
    return (_constrPtr == NULL);
  }
  /// to delete, only for debugging
  inline Constraint * constrPtr() const
  {
    return _constrPtr;
  }
  virtual const Double & defaultCost();
  virtual const Double & curCost() const;
  virtual const Double & costrhs() const;

  virtual void costrhs(const Double & newCostrhs)
  {
    return Variable::costrhs(newCostrhs);
  }
  virtual void resetCostFromDefaultCost(const Double & factor = 0.2);
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);

  /**
   * Also do behind the schene membership
   * for MasterConstr and SubProblemVariables
   */
  virtual void setMembership();

  const LocalArtClassId & localClassId() const; ///added by Ruslan

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class Constraint : public VarConstr
{
protected:
  Double _memorisedCurRhs;
  bool _lhsComputed;
  bool _toBeUsedInPreprocessing;
  bool _considerAsEqualityInPreprocessing;
  LocalArtificialVar * _posLocalArtVarPtr;
  LocalArtificialVar * _negLocalArtVarPtr;

  /// start added by Ruslan: needed for the new preprocessing
  Double _minSlack;
  Double _maxSlack;
  Double _curMinSlack;
  Double _curMaxSlack;
  bool _inPropagateList;
  /// end added by Ruslan

  Constraint() = delete;
public:
  /**
   * @param  costrhs = rhs for a cosntraint
   * @param   sense  Constraints: 'G', 'L', or 'E';
   * @param  type
   * type = 'C' for core (required for the IP formulation,
   * type = 'F' for facultative (only helpfull for the LP formulation),
   * type = 'S' for constraints defining a subsystem in column generation for extended formulation approach
   * type = 'M' for constraints defining a pure master constraint
   * type = 'X' for constraints defining a subproblem convexity constraint in the master
   * @param  kind  For Constraints: 'E' = explicit, 'I' = implicit,
   * 'R' = delayed constraints added only after solving the LP root node without them
   * @param flag 's' (by  default) for static VarConstr belonging to the problem
   * (and erased when the problem is erased);
   * 'd' for generated dynamic VarConstr not belonging to the problem
   * @param val Inital dual value a constraint;
   * @param upperBound Represents a global upper bound on on  constraint dual value
   * @param  lowerBound represents global lower bound on constraint dual value
   * @param directive for branching
   * @param priority for branching or separation (the higher the most mpriority)
   **/
  Constraint(Model * modelPtr, 
	         const std::string & name,
	         const Double & costrhs,
             const char & sense = 'E', 
	         const char & type = 'C',
	         const char & kind = 'E',
             const char & flag = 's', 
	         const int & index = -1,
	         const Double & val = 0,
             const Double & upperBound = BapcodInfinity,
	         const Double & lowerBound = -BapcodInfinity,
	         const char & directive = 'U',
             const Double & priority = 1.0, 
	         const bool & presetMembership = true,
             const bool & toBeUsedInPreprocessing = true,
             const bool & considerAsEqualityInPreprocessing = false);
  Constraint(const Constraint & that);
  virtual ~Constraint();
  virtual bool operator<(const VarConstr & that) const;

  /**
   * Is there  only to  handle
   * the case of artificial variables
   */
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);

  /**
   * Default count() return 1 for global art var
   * in 'G' constraint or for localartVar
   * Is there only to handle the case of artificial variables
   */
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);

  virtual const Double curRhs() const;
  virtual void curRhs(const Double & rhs);

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual Double computeLhs(const std::list<std::pair<Variable *, Double> > & curSol);
  virtual Double computeLhs(const VarPtr2DoubleMap & solVarValMap);
  virtual Double computeLhs(const VarPtrSet & curSol);
  virtual Double computeLhs(const std::list<Variable *> & curSol);
  virtual const Double & computeViolation(const Double & lhs);

  /// start added by Ruslan
  const Double & maxSlack() const
  {
    return _maxSlack;
  }

  const Double & minSlack() const
  {
    return _minSlack;
  }

  void maxSlack(const Double & slack)
  {
    _maxSlack = _curMaxSlack = slack;
  }

  void minSlack(const Double & slack)
  {
    _minSlack = _curMinSlack = slack;
  }

  const bool & inPropagateList() const
  {
    return _inPropagateList;
  }

  void addToPropagateList(ConstrPtrList & constrsListToPropagate)
  {
    if (!_inPropagateList)
      {
        constrsListToPropagate.push_back(this);
        _inPropagateList = true;
      }
  }

  void removeFromPropagateList()
  {
    _inPropagateList = false;
  }

  const Double & curMaxSlack() const
  {
    return _curMaxSlack;
  }

  void curMaxSlack(const Double & value)
  {
    _curMaxSlack = value;
  }

  void incrCurMaxSlack(const Double & delta)
  {
    _curMaxSlack += delta;
  }

  const Double & curMinSlack() const
  {
    return _curMinSlack;
  }

  void incrCurMinSlack(const Double & delta)
  {
    _curMinSlack += delta;
  }

  void curMinSlack(const Double & value)
  {
    _curMinSlack = value;
  }

  virtual void addToPreprocessedList();

  void resetSlacksAndRhsToDefaults();
    
  virtual bool activateConstraint(bool toSetInForm);

  virtual bool desactivateConstraint(const VcIndexStatus::VcStatus & status,
                                     bool toUnsetInForm = true);

  /// end added by Ruslan

  /**
   * Also do behind the schene membership
   * for MasterConstr and SubProblemVariables
   */
  virtual void setMembership();
  virtual void setMembership(const VarIndexManager & varSet, const VcIndexStatus::VcStatus& status, char flag);

  /**
   * Set membership by enumerating problem variables
   */
  virtual void enumerativeSetMembership();

  /**
   * Put in  _member2coefMap the subset of varSet
   * to which constrs belongs and
   */
  virtual void setMembership(const VarPtrSet & varSet);
  virtual void setMembership(const VarIndexManager & varSet);

  virtual bool violated(Variable * varPtr, const Double & useLevel = 1);
  virtual void resetRhs();

  virtual void toBeUsedInPreprocessing(const bool value)
  {
    _toBeUsedInPreprocessing = value;
  }
  virtual bool toBeUsedInPreprocessing() const
  {
    return _toBeUsedInPreprocessing;
  }
  virtual bool considerAsEqualityInPreprocessing() const
  {
    return _considerAsEqualityInPreprocessing;
  }

  virtual GenericConstr * genConstrPtr() const
  {
    return NULL;
  }

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;

  LocalArtificialVar * posLocalArtVarPtr() const;
  LocalArtificialVar * negLocalArtVarPtr() const;
  void posLocalArtVarPtr(LocalArtificialVar * ptr);
  void negLocalArtVarPtr(LocalArtificialVar * ptr);

  void includeOvfPtr(OvfConstr * ovfPtr);
};

/**
 * @brief to recognize non linear instanciated Csontraints and imply Specific behaviour of routines 
 */
class Base4NonLinearConstraint
{
public:

  virtual ~Base4NonLinearConstraint() { }

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	  return compareIdentifier(VcId::Base4NonLinearConstraintMask,vcIdentifier);
  }
};

struct ValueRecord
{
  Double _value;
  Double _dfracValue;
  Double _lfracValue;
  bool _isFractional;
  ValueRecord(const Double & val): 
  _value(val) , _dfracValue(0), _lfracValue(0), _isFractional(false) {} 
  ValueRecord(const Double & val, const double & precision): 
    _value(val) , _dfracValue(0), _lfracValue(Lfrac(val)), _isFractional(false) 
  { 
    _dfracValue = Dmin(Ufrac(val), _lfracValue);
    if (precision < _dfracValue) _isFractional = true; 
  } 
};

struct MasterVarSolution : public  std::list< std::pair < Variable *, ValueRecord > >
{
  void push_back(Variable * varPtr, const ValueRecord  & rec);
};

#endif // VarConstrClasses_h

