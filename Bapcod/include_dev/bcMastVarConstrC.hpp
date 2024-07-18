/**
 *
 * This file bcMastVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef MastVarConstrClasses_h
#define MastVarConstrClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"

class Node;
class ProbConfig;
class SubProbVariable;
class MasterConf;
class ColGenSpConf;
class InstMasterConstr;
class MastColumn;
class CompSetInstMastBranchConstr;
class Base4NonLinearConstraint;

struct SortMastColPerDecreasingSpVal
{
  SortMastColPerDecreasingSpVal(Variable * spVarPtr) : _spVarPtr(spVarPtr){}
  
  SortMastColPerDecreasingSpVal(){}
  Variable * _spVarPtr;
  bool operator()(const std::pair < MastColumn *, ValueRecord > & a, const std::pair < MastColumn *, ValueRecord > & b) const;
};

/**
 * A base class for pure Master Variables
 */
class InstMasterVar: public InstanciatedVar
{
 public:
  InstMasterVar(InstanciatedVar * iv);
  /** 
   * @param priority Lower priority index means chosen first for branching
   */
  InstMasterVar(const IndexCell & id, 
		        GenericVar * genVarPtr,
		        MasterConf * masterConfPtr,
		        const std::string & name,
		        const Double & costrhs = 0,
		        const char & sense = 'P',
		        const char & type = 'I',
		        const char & kind = 'E',
		        const Double & upperBound = BapcodInfinity,
		        const Double & lowerBound = 0,
		        const char & flag = 's',
		        const char & directive = 'U',
		        const Double & priority = 0.0,
		        const Double & val = 0,
		        const Double & globalUb = BapcodInfinity,
		        const Double & globalLb = 0,
		        const bool & presetMembership = false);
  
  virtual ~InstMasterVar();
  virtual void enumerativeSetMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);
};

/**
 * A Marker  Master Constraints: manages MastColumn count and coef and _subProbVarMember2coefMap
 */
class MasterConstr
{
 protected:
  MasterConf * _masterConfPtr;
  MapSubProbVariablePtr2Double _subProbVarMember2coefMap;
  ConstVarConstrPtr2Double _pureMastVarMember2coefMap; 
 public:
  MasterConstr(MasterConf * masterConfPtr);
  virtual ~MasterConstr();

	virtual const MapSubProbVariablePtr2Double & subProbVarMember2coefMap() const
	{
		return _subProbVarMember2coefMap;
	}

  virtual MapSubProbVariablePtr2Double & subProbVarMember2coefMap()
  {
    return _subProbVarMember2coefMap;
  }
  
  virtual const Double & upToDateSubProbVarMembCoef(SubProbVariable * vcPtr) const;
  virtual const Double & includeSubProbVarAsMember(SubProbVariable * spvPtr, const Double & coef);
  virtual void eraseSubProbVarAsMember(SubProbVariable * spvPtr);
  virtual void clearSubProbVarMember();
  virtual const Double & includePureMastVarAsMember(InstMasterVar * pmvPtr, const Double & coef);
  virtual void erasePureMastVarAsMember(VarConstr * spvPtr);
  void clearPureMastVarMember();
  virtual bool countColumn(ConstVarConstrConstPtr vcPtr);
  virtual Double coefColumn(ConstVarConstrConstPtr vcPtr);
  MasterConf * masterConfPtr() const {return _masterConfPtr;}

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};


/**
 * A base class for Master Constraints
 */
class InstMasterConstr: public MasterConstr, public InstanciatedConstr
{
 
 int _treatOrderId;
 
 public:
   
  InstMasterConstr(InstanciatedConstr* ic);
  /** 
   * @param directive For branching on constraints
   */
  InstMasterConstr(const IndexCell& id,
                   GenericConstr* genConstrPtr,
                   ProbConfig* probConfigPtr,
                   const std::string& name,
                   const Double& costrhs,
                   const char& sense,
                   const char& type = ' ',
                   const char& kind = 'E',
                   const char& flag = 's',
                   const Double& val = 0,
                   const Double& upperBound = BapcodInfinity,
                   const Double& lowerBound = - BapcodInfinity,
		   const char & directive = 'U',
		   const Double & priority = 1.0, 
		   const bool & presetMembership = true, 
		   const bool & toBeUsedInPreprocessing = true,
		   const bool & considerAsEqualityInPreprocessing = false);
  
  InstMasterConstr(InstMasterConstr * imcPtr,  
		           GenericConstr * genBrConstrPtr,
		           const std::string & name,
		           const Double & rhs,
		           const char & sense,
		           const char & type,
		           const char & kind,
		           const char & flag);

  InstMasterConstr(const InstMasterConstr& ic);
  
  virtual ~InstMasterConstr();
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);

  /**
   * Insert this in vcPtr.member and vice versa
   */
  virtual const Double & includeMember(VarConstr * vcPtr, const Double & coef, const bool & cumulativeCoef);
  virtual void setMembership();
  
  /**
   * Set membership by enumerating problem variables
   */
  virtual void enumerativeSetMembership();
  virtual void clearMembership();
  virtual bool addLocalArtVar(const BcObjStatus::MinMaxIntFloat & minmax = BcObjStatus::minInt);
  void createStabInfo(const BcObjStatus::MinMaxIntFloat & minmax = BcObjStatus::minInt);
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
  
  void treatOrderId(const int value )
  {
    _treatOrderId = value;
  }

  const int treatOrderId() const
  {
    return _treatOrderId;
  }

};

inline std::ostream& operator<<(std::ostream& os, const InstMasterConstr & that)
{
  return that.print(os);
}

/**
 * A base class for non linear Master Constraints 
 */
class NonLinearInstMastConstr: public InstMasterConstr, public Base4NonLinearConstraint
{
 public:
  NonLinearInstMastConstr(InstanciatedConstr * ic);

  /** 
   * @param directive For branching on constraints
   */
  NonLinearInstMastConstr(const IndexCell& id,
			              GenericConstr* genConstrPtr,
			              ProbConfig* probConfigPtr,
			              const std::string& name,
			              const Double& costrhs,
			              const char& sense,
			              const char& type = ' ',
			              const char& kind = 'E',
			              const char& flag = 's',
			              const Double& val = 0,
			              const Double& upperBound = BapcodInfinity,
			              const Double& lowerBound = - BapcodInfinity,
			              const char & directive = 'U',
			              const Double & priority = 1.0);

  NonLinearInstMastConstr(const NonLinearInstMastConstr& ic);

  virtual ~NonLinearInstMastConstr();
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

inline std::ostream& operator<<(std::ostream& os, const NonLinearInstMastConstr & that)
{
  return that.print(os);
}

/**
 * A base class for  Master Constraints 
 * implementing lower and upper bound 
 * on the number of subproblme solutions 
 * that can be selected
 */
class InstMastConvexityConstr: public InstMasterConstr
{
  /// spConf
  ColGenSpConf * _cgSpConfPtr;
  CompSetInstMastBranchConstr * _memorisedCsBcConstr;
  bool _locallyValidRhsDefined;
  Double _locallyValidRhs;

 public:
  enum BoundsWhanInactive {lowerBoundWhenInactive = 0, upperBoundWhenInactive = 1000000};
  InstMastConvexityConstr(const IndexCell & id,  
			  GenericConstr * genConstrPtr, 
			  MasterConf * masterConfPtr, 
			  ColGenSpConf * cgSpConfPtr,
			  const std::string & name, 
			  const Double & costrhs, 
			  const char & sense,
			  const char & type = 'X', 
			  const char & kind = 'E', 
			  const char & flag = 's', 
			  const int & index = -1, 
			  const Double & val = 0,
			  const Double & upperBound = BapcodInfinity,
			  const Double & lowerBound = 0);
  
  virtual ~InstMastConvexityConstr();
  virtual ColGenSpConf * cgSpConfPtr() const;
  virtual void memorisedCsBcConstr(CompSetInstMastBranchConstr * ptr) 
  {
    _memorisedCsBcConstr = ptr; 
  }
  virtual void defineLocalRhs(const Double & localRhs);
  void newLocalRhs(const Double & localRhs);
  const Double & newLocalRhs() const;
  virtual void resetLocalRhs();
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void addMember(VarConstr * vcPtr);
  virtual void setMembership();
  virtual void enumerativeSetMembership();
  virtual void clearMembership();
  virtual bool addLocalArtVar(const BcObjStatus::MinMaxIntFloat & minmax = BcObjStatus::minInt);
  virtual const Double & costrhs() const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

inline std::ostream& operator<<(std::ostream& os, const InstMastConvexityConstr & that)
{
  return that.print(os);
}

/**
 * A base class for  Master branching Constraints
 * Base class for all branching constraints 
 * that are enforced explicitly or implicitly 
 * in the master
 */
class InstMasterBranchingConstr:  public InstMasterConstr, public BranchingConstrBaseType
{
  /* /// spConf where the branching constraint induces modifications */
  /*   std::map<IndexCell, ColGenSpConf *> _spConfMap; */
 public:
  InstMasterBranchingConstr(const IndexCell & id, 
			    GenericConstr * genBrConstrPtr, 
			    ProbConfig * probConfigPtr,
                            const std::string & name, 
			    const Double & rhs, 
			    const char & sense, 
			    const char & type = ' ',
                            const char & kind = 'I', 
			    const char & flag = 'd', 
			    const Double& val = 0, 
			    const Double& upperBound = BapcodInfinity,
			    const Double& lowerBound = - BapcodInfinity,
			    const char & directive = 'U', 
			    const Double & priority = 1.0);
  
  InstMasterBranchingConstr(InstMasterConstr * imcPtr,  
			    GenericConstr * genBrConstrPtr,
			    const std::string & name, 
			    const Double & rhs, 
			    const char & sense,
			    const char & type, 
			    const char & kind, 
			    const char & flag);
  
  /// Cast operator
  InstMasterBranchingConstr(InstanciatedConstr * icPtr);
  InstMasterBranchingConstr(const InstMasterBranchingConstr & that);
  virtual ~InstMasterBranchingConstr();

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  void shortPrint(std::ostream& os = std::cout) const {}

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

inline std::ostream& operator<<(std::ostream& os, const InstMasterBranchingConstr & that)
{
  return that.print(os);
}

/**
 * A class for  branching on the number 
 * of columns with a specific master variable  
 * or the aggregate value of a SP variable
 */
class GenVarInstMastBranchConstr: public InstMasterBranchingConstr
{
  InstanciatedVar * _instVarPtr;
 public:
  GenVarInstMastBranchConstr(const IndexCell & id, 
			     GenericConstr * genBrConstrPtr, 
			     ProbConfig * probConfigPtr, 
			     InstanciatedVar * instVarPtr,
			     const std::string & name, 
			     const Double & rhs, 
			     const char & sense, 
			     const char & type = ' ',
			     const char & kind = 'I', 
			     const char & flag = 'd', 
			     const Double& val = 0, 
			     const Double& upperBound = BapcodInfinity,
			     const Double& lowerBound = - BapcodInfinity,
			     const char & directive = 'U', 
			     const Double & priority = 1.0);
  
  virtual ~GenVarInstMastBranchConstr();
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  void shortPrint(std::ostream& os = std::cout) const;
  std::vector<std::string> forDotPrint() const;


  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

/**
 * A class for  branching on the number 
 * of columns with a specific component 
 * of a master variable,  or the aggregate 
 * value of a SP variable, or an instantiated 
 * constraint 
 */
class BasicConstrInstMastBranchingConstr: public InstMasterBranchingConstr 
{
  InstanciatedConstr * _instConstrPtr;
  std::string _description;
 public:
  BasicConstrInstMastBranchingConstr(const IndexCell & id, 
				     GenericConstr * genBrConstrPtr, 
				     ProbConfig * probConfigPtr, 
				     InstanciatedConstr * instConstrPtr,
                     const std::string & description,
				     const std::string & name, 
				     const Double & rhs, 
				     const char & sense, 
				     const char & type = ' ',
				     const char & kind = 'I', 
				     const char & flag = 'd');
  
  BasicConstrInstMastBranchingConstr(InstMasterConstr * imcPtr,  
				     GenericConstr * genBrConstrPtr,
                     const std::string & description,
				     const std::string & name, 
				     const Double & rhs,
				     const char & sense,
				     const char & type, 
				     const char & kind, 
				     const char & flag);
  
  virtual ~BasicConstrInstMastBranchingConstr();
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void shortPrint(std::ostream& os = std::cout) const;
  virtual std::vector<std::string> forDotPrint() const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

/**
 * A class for  branching on a set of SP variable components
 */
class CompSetInstMastBranchConstr: public InstMasterBranchingConstr, public Base4NonLinearConstraint
{
 public:
  struct BCDepthWhenGeneratedComparator
  {
    bool operator()(CompSetInstMastBranchConstr * a, CompSetInstMastBranchConstr * b)
    {
      return (a->_depthWhenGenerated > b->_depthWhenGenerated);
    }
  };
  struct CSbrConstDepthComparator
  {
    bool operator()(CompSetInstMastBranchConstr * a, CompSetInstMastBranchConstr * b)
    {
      int deptha(a->_setOfPredCSconstrPtr.size());
      int depthb(b->_setOfPredCSconstrPtr.size());
      if (deptha < depthb) 
	return true;
      
      return false;
    }
  };
  struct CSbrGreedyComparator
  {
    bool operator()(CompSetInstMastBranchConstr * a, CompSetInstMastBranchConstr * b)
    {
      int deptha(a->_setOfPredCSconstrPtr.size());
      int depthb(b->_setOfPredCSconstrPtr.size());
      if (deptha == depthb)
        return ((a->_marginLvalue * a->_sigma)
                > (b->_marginLvalue * b->_sigma));
      return (deptha < depthb);
    }
  };
 private:
  ComponentSequence _compBoundSet;
  
  /// Pointer to  predecessors in tree of column classes
  CompSetInstBrConstrPtrSet _setOfPredCSconstrPtr;

  /// Pointer to direct predecessor in tree of column classes
  CompSetInstMastBranchConstr * _dirPredCSconstrPtr;
  bool _associatedPricingSPsolved;

  /// sp solution value
  Bound _pricingSpZetaVal;
  
  /**
   * Total dual value = \sigma_k + \sum_{j \in pred(k)} \sigma_j
   */
  Double _sigma;
  
  /// \overline{L}_k
  Double _marginLvalue;
  Double _margLvalue4DualBd;
  Solution * _solPtr;
  Bound _bestReducedCost;
 public:
  CompSetInstMastBranchConstr(const ComponentSequence & cpList, 
			      const IndexCell & id, 
			      GenericConstr * genBrConstrPtr, 
			      ProbConfig * probConfigPtr,
			      const std::string & name, 
			      const Double & rhs, 
			      const char & kind = 'E');

  const ComponentSequence & compBoundSet() const
  {
    return _compBoundSet;
  }
  
  /// Virtual const Double & costrhs() const;
  bool CBsatisfiedBySol(Solution * solPtr) const;
  virtual ~CompSetInstMastBranchConstr()
  {
    if (_solPtr != NULL) delete _solPtr;
  }
  Solution * solPtr() const {return  _solPtr;}
  void solPtr(Solution * ptr)
  {
    if (_solPtr != NULL) delete _solPtr;
    _solPtr = ptr;
  }
  const bool & associatedPricingSPsolved() const {return _associatedPricingSPsolved;}
  void associatedPricingSPsolved(const bool & flag)
  {
    if (!flag) _margLvalue4DualBd = curRhs();
    _associatedPricingSPsolved = flag;
  }
  const Bound & bestReducedCost() const {return _bestReducedCost;}
  void bestReducedCost(const Bound & bd) {_bestReducedCost = bd;}
  const CompSetInstBrConstrPtrSet & setOfPredCSconstrPtr() const {return _setOfPredCSconstrPtr;}
  CompSetInstBrConstrPtrSet & setOfPredCSconstrPtr() { return  _setOfPredCSconstrPtr;}

  virtual bool enforces(const ComponentSequence & cpList) const;
  virtual bool incompatibleCsBd() const;
  
  /// Constraint enforces same class bounds that imlied by cpList
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  void shortPrint(std::ostream& os = std::cout) const;
  virtual std::vector<std::string> forDotPrint() const;
  void addToSigma(const Double & s) {_sigma += s;}
  
  /// sp solution value
  const Bound & pricingSpZetaVal() const 
  {
    return _pricingSpZetaVal;
  }

  CompSetInstMastBranchConstr * dirPredCSconstrPtr() const 
  {
    return _dirPredCSconstrPtr; 
  }
  
  void setDirPredCSconstrPtr(CompSetInstMastBranchConstr* dirPredCSconstrPtr)
  {
    _dirPredCSconstrPtr = dirPredCSconstrPtr;
  }

  void recSol(Solution * solPtr, const Bound & dualBd, const Bound & redCost);
  
  void resetBounds(const Bound & dualBd, const Bound & redCost);

  ColGenSpConf * ColGenSpConfPtr() const;
  
  const Bound & reducedCost() const 
  {
    return  _bestReducedCost;
  }
  
  /** 
   * @return total dual value = sigma_i +  sum(j in pred(i)) sigma_j
   */
  const Double & sigma() const 
  {
    return _sigma;
  }
  
  /**
   * Total dual value = sigma_i +  sum(j in pred(i)) sigma_j
   */
  void setSigma(const Double & s = 0.0)  
  {
    _sigma = s;
  }
  
  void reset();
  virtual void recordInducedVarBounds();
  virtual void setMembership(); 

  /**
   * Set membership by enumerating problem variables
   */
  virtual void enumerativeSetMembership();
  
  /// \overline{L}_k
  const   Double & margLvalue4DualBd() const 
  {
    return  _margLvalue4DualBd;
  }
  const   Double & marginLvalue() const 
  {
    return  _marginLvalue;
  }
  
  Double & marginLvalue()
  {
    return  _marginLvalue;
  }

  friend bool operator==(const CompSetInstMastBranchConstr & a, const CompSetInstMastBranchConstr & b);

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

inline std::ostream& operator<<(std::ostream& os, const  CompSetInstMastBranchConstr & that)
{
  return that.print(os);
}

bool operator==(const CompSetInstMastBranchConstr & a, const CompSetInstMastBranchConstr & b);

struct BcCustomNonLinearCutInfo;
class GenericCustomNonLinearCutConstr;

class CustomNonLinearCut: public InstMasterConstr, Base4NonLinearConstraint
{
  friend class BcCustomNonLinearCut;
  friend class GenericCustomNonLinearCutConstr;
  const BcCustomNonLinearCutInfo * _cutInfoPtr;
  GenericCustomNonLinearCutConstr * _genUserNonLinearCutConstrPtr;
public:
  CustomNonLinearCut(const IndexCell& id,
                     GenericCustomNonLinearCutConstr * genConstrPtr,
                     ProbConfig* probConfigPtr,
                     const std::string & name,
                     const BcCustomNonLinearCutInfo * cutInfoPtr);
  virtual ~CustomNonLinearCut();
  virtual void setMembership();
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class GenericSoftConflictsCutConstr;

class SoftConflictsCut: public InstMasterConstr, Base4NonLinearConstraint
{
	friend class GenericSoftConflictsCutConstr;
	int _cutType; /// type = 0 : coefficient of the column in the cut is equal to the number of conflicts in the column
	             ///            (in the subproblem, we need an indication variable per conflict)
	             /// type = 1 : coefficient of the column in the cut is 1 if there are conflicts in the column,
	             ///            0 otherwise (in the subproblem, we need an indicator variable per cut)
	GenericSoftConflictsCutConstr * _genSoftConflictConstrPtr;
	const std::vector<std::pair<SubProbVariable *, SubProbVariable *> > _conflicts;
public:
	SoftConflictsCut(const IndexCell& id,
					 GenericSoftConflictsCutConstr * genConstrPtr,
					 ProbConfig* probConfigPtr,
					 const std::string & name,
					 const Double & rhs,
					 const int cutType,
					 const std::vector<std::pair<SubProbVariable *, SubProbVariable *> > & conflicts);
	virtual ~SoftConflictsCut();
	virtual void setMembership();
	virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
    const std::vector<std::pair<SubProbVariable *, SubProbVariable *> > & conflicts() const;
    const int cutType() const {return _cutType;}
};

#endif // MastVarConstrClasses_h
