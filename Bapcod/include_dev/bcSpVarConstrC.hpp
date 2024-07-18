/**
 *
 * This file bcSpVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef SpVarConstrClasses_h
#define SpVarConstrClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcVarConstrC.hpp"

class CompSetInstMastBranchConstr;

/**
 * Keeps a link to the master to compute the subproblem cost 
 * which involves dual master information 
 * and variable bounds implied by master constraints
 */
class SubProbVariable: public InstanciatedVar
{
  MasterConf * _masterConfPtr;
  ColGenSpConf * _cgSpConfPtr;

  ConstVarConstrPtr2Double _masterConstrMember2coefMap;
  MapMastColumnPtr2Double _masterColumnMember2coefMap;

  Double _memorisedCurGlobUb;
  Double _memorisedCurGlobLb;

 public:

  /** 
   * @param priority Lower priority index means chosen first for branching
   * 
   * @return 
   */
  SubProbVariable(MasterConf * masterConfPtr, 
		          const IndexCell & id,
		          GenericVar * genVarPtr,
		          ProbConfig * probConfigPtr,
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
		          const bool & presetMembership = true);
  
  SubProbVariable(InstanciatedVar * iv, MasterConf * masterConfPtr);
  virtual ~SubProbVariable();
  
  /**
   * Virtual copy constructor
   */
  virtual ConstVarConstrPtr2Double & masterConstrMember2coefMap()
  {
    return _masterConstrMember2coefMap;
  }

    virtual const ConstVarConstrPtr2Double & masterConstrMember2coefMap() const
    {
        return _masterConstrMember2coefMap;
    }

    virtual MapMastColumnPtr2Double & masterColumnMember2coefMap()
  {
    return _masterColumnMember2coefMap;
  }

  virtual const Double & includeMasterConstrAsMember(VarConstr * mcPtr, const Double & coef);
  virtual const Double & includeMasterColAsMember(MastColumn * colPtr, const Double & coef);
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);
  virtual const Double & upToDateMasterConstrMembCoef(ConstVarConstrConstPtr vcPtr) const;
  virtual void eraseMasterConstrAsMember(VarConstr * mcPtr);

  /**
   * Insert this in vcPtr.member and vice versa
   */
  virtual const Double & includeMember(VarConstr * vcPtr, const Double & coef, const bool & cumulativeCoef);
  virtual void setMembership();
  virtual void clearMembership();
  
  
  /**
   * Set membership by enumerating problem constraints
   */
  virtual void enumerativeSetMembership();
  virtual const Double & curCost() const;
  virtual const Double & costrhs() const;
  virtual void costrhs(const Double & newCostrhs)
  {
      return Variable::costrhs(newCostrhs);
  }
  void updateCurCostWithConstraint(const Constraint * constrPtr);
  void calculateScaledCurCostForSafeDualBound(const bool & inPurePhaseOne);
  virtual void resetCost(const bool & inPurePhaseOne);
  virtual void recallMemorisedBounds();

  virtual bool infeasible() const;
  
  /**
   * Measures the variable min contribution 
   * to the satisfaction of a specific  constraint
   */
  virtual const Double lhsMinContrib(ConstVarConstrConstPtr  vcPtr);

  /**
   * Measures the variable min contribution 
   * to the satisfaction of a specific  constraint
   */
  virtual const Double lhsMaxContrib(ConstVarConstrConstPtr vcPtr); 
  virtual const Double & curGlobUb() const 
  { 
    return _memorisedCurGlobUb;
  }
  
  virtual const Double & curGlobLb() const 
  { 
    return _memorisedCurGlobLb;
  }

  /// start added by Ruslan for new var constr reset and preprocessing

  virtual void globalCurUb(const Double & ub)
  {
    _memorisedCurGlobUb = ub;
  }

  virtual void  globalCurLb(const Double & lb)
  {
    _memorisedCurGlobLb = lb;
  }

  virtual const Double & globalCurUb() const
  {
    return _memorisedCurGlobUb;
  }

  virtual const Double & globalCurLb() const
  {
    return _memorisedCurGlobLb;
  }

  virtual void localCurUb(const Double & ub)
  {
    _curUb = _memorisedCurUb = ub;
  }

  virtual void  localCurLb(const Double & lb)
  {
    _curLb = _memorisedCurLb = lb;
  }

  virtual const Double & localCurUb() const
  {
    return _memorisedCurUb;
  }

  virtual const Double & localCurLb() const
  {
    return _memorisedCurLb;
  }

  virtual void globalUb(const Double & ub)
  {
    _globalUb = _memorisedCurGlobUb = ub;
  }

  virtual void  globalLb(const Double & lb)
  {
    _globalLb = _memorisedCurGlobLb = lb;
  }

  virtual const Double & globalUb() const
  {
    return _globalUb;
  }

  virtual const Double & globalLb() const
  {
    return _globalLb;
  }

  virtual void localUb(const Double & ub)
  {
    _upperBound = _curUb = _memorisedCurUb = ub;
  }

  virtual void  localLb(const Double & lb)
  {
    _lowerBound = _curLb = _memorisedCurLb = lb;
  }

  virtual const Double & localUb() const
  {
    return _upperBound;
  }

  virtual const Double & localLb() const
  {
    return _lowerBound;
  }

  virtual void resetBoundsAndCostToDefaults();

  /// end added by Ruslan for new var constr reset and preprocessing

  virtual const Double minGlobCurLb() const;
  virtual const Double maxGlobCurUb() const;

  virtual MasterConf * masterConfPtr() const;
  virtual ColGenSpConf * cgSpConfPtr() const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

/**
 * InstSubProbBranchingConstr 
 */
class InstSubProbBranchingConstr:  public InstanciatedConstr, public BranchingConstrBaseType
{
 public:
  /** 
   * @param directive For branching on constraints
   * @param priority For branching on constraints
   */
  InstSubProbBranchingConstr(const IndexCell& id,  
			                 GenericConstr * genConstrPtr,
			                 ProbConfig * probConfigPtr,
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
  
  virtual ~InstSubProbBranchingConstr(){}
  
  virtual bool violated(Variable * varPtr, const Double & useLevel = 1);
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

/**
 * RyanAndFosterInstSubProbBranchConstr
 */
class RyanAndFosterInstSubProbBranchConstr:  public InstSubProbBranchingConstr
{
  friend class BcRyanAndFosterBranchConstr;

  /// Specific instanciated first var on which we branch
  InstanciatedVar * _variPtr;
  
  /// Specific instanciated second var on which we branch
  InstanciatedVar * _varjPtr;
 public:
  /** 
   * 
   * 
   * @param directive For branching on constraints
   * @param piority  For branching on constraints 
   * 
   */
  RyanAndFosterInstSubProbBranchConstr(InstanciatedVar * variPtr,
                                       InstanciatedVar * varjPtr,
				                       GenericConstr * genConstrPtr,
				                       ProbConfig * probConfigPtr,
				                       const std::string & name,
				                       const Double & costrhs,
				                       const char & sense,
				                       const char & type = ' ',
				                       const char & kind = 'E',
				                       const char & flag = 'd',
				                       const Double & val = 0,
				                       const Double & upperBound = BapcodInfinity,
				                       const Double & lowerBound = - BapcodInfinity,
				                       const char & directive = 'U',
				                       const Double & priority = 1.0);
  
  virtual ~RyanAndFosterInstSubProbBranchConstr()
  {
  }
  
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void shortPrint(std::ostream& os) const;
  virtual std::vector<std::string> forDotPrint() const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

#endif // SpVarConstrClasses_h
