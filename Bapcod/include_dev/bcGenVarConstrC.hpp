/**
 *
 * This file bcGenVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef GenVarConstrClasses_h
#define GenVarConstrClasses_h

#include "bcVcIdentifierC.hpp"


#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"

class GenericVarConstr;
class GenericVar;
class GenericConstr;
class NonLinearGenericConstr;
class InstanciatedVarConstr;
class InstanciatedVar;
class InstanciatedConstr;
class MastColumn;
class Node;
class ProbConfig;
class OvfConf;
class MasterConf;
class Model;
class Base4NonLinearConstraint;
class DynamicGenericConstr;
class ColGenSpConf;
struct MasterColSolution;
class BcModel;
class BcConstr;
class BcAddConstrFunctor;
class BcSolverOracleFunctor;
class BcCutSeparationFunctor;
class CustomNonLinearCut;
class SoftConflictsCut;

typedef std::map < IndexCell, InstanciatedVar * , IndexCellSort > IndexCell2InstancVarPtrMap;
typedef std::map < IndexCell, InstanciatedConstr * , IndexCellSort > IndexCell2InstancConstrPtrMap;

/**
 * by Artur: Sorting procedure to select 
 * among branching variables: 
 * by genericVarConstr _priorityLevel 
 * than by constraint PriorityRule 
 * and priority of fract val.
 */
class BrVarPriorityCalcAndComp
{
public:
  virtual ~BrVarPriorityCalcAndComp()
  {
  }
  virtual bool operator()(Variable * a, Variable * b) = 0;
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir) = 0;
};

class BrVarPriorityComparator
{
  BrVarPriorityCalcAndComp * _compPtr;
public:
  BrVarPriorityComparator(BrVarPriorityCalcAndComp* compPtr) :
      _compPtr(compPtr)
  {}
  bool operator()(Variable * a, Variable * b)
  {
    return (*_compPtr)(a,b);
  }
};

class BrVarPriorityCalcAndComp_FirstFound : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_HighestPriority : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_FracWeightedPriority : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_MostFractional : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_LeastFractional : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_Closest2RoundUp : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

class BrVarPriorityCalcAndComp_Closest2RoundDown : public BrVarPriorityCalcAndComp
{
public:
  virtual bool operator()(Variable * a, Variable * b);
  virtual void calculateVarViolAndDir(Variable* var, Double & viol, char & dir);
};

/**
 * Sorting procedure to select among cuts: 
 * by genericVarConstr _priorityLevel than 
< * by constraint PriorityRule and priority 
 * of violation or fract val.
 */
struct CutSeparationPriorityComp
{
  bool operator()(InstanciatedConstr *  a, InstanciatedConstr * b) const;
};

/**
 * BaseType for Generic Var and Constr
 */
class GenericVarConstr
{
private:
  BcVarConstrType::BcVcType _gvcType;
  ProbConfig * _probConfigPtr;
  int _ref;
  /// Set of formulations  (eihter Master or ColGenSpConf)  in which the global genVarConstr has duplicate that  is included
  std::set<ProbConfig *> _vcProbConfigPtrSet;
  
  /// Belongs to a Model
  Model * _model;
  
  /**
   * Class priority rule used with 
   * this class for separation of 
   * dynamic constraints or branching
   */
  SelectionStrategy _priorityRule;
  Double _priorityLevel;
  bool _toBeUsedInPreprocessing;
  bool _considerAsEqualityInPreprocessing;

  std::string _defaultName;
  char _defaultType;
  char _defaultKind;
  char _defaultFlag;
  char _defaultSense;
  char _defaultDirective;
  Double _defaultCostRhs;
  Double _defaultUb;
  Double _defaultLb;
  Double _defaultGlobalUb;
  Double _defaultGlobalLb;
  Double _defaultVal;
  MultiIndexNames _multiIndexNames;
  MapValByMultiIndex _curSolVal;
  int _dimension; //added by Issam to fix the dimension of arrays. 
                  //(not allowing indices with different number of elements)
                  //the first instanciated var/constr will fix the dimension.

 public:

  GenericVarConstr(Model * modelPtr,
                   const BcVarConstrType::BcVcType &  vctype,
                   ProbConfig * probConfPtr,
                   const std::string & name,
                   const MultiIndexNames & multiIndexNames = MultiIndexNames(),
                   const SelectionStrategy & priorityRule = SelectionStrategy::NotConsideredForSelection,
                   const Double & priorityLevel = -1,
                   const bool & toBeUsedInPreprocessing = true);

  virtual ~GenericVarConstr();
  virtual const BcModel & bcModel() const;
  virtual BcModel & bcModel();
  const BcVarConstrType::BcVcType & gvcType() const {return _gvcType;}
  int ref() const {return _ref;}

  virtual const std::set<ProbConfig *> & vcProbConfigPtrSet() const 
  {
    return _vcProbConfigPtrSet;
  }
  virtual ProbConfig * probConfPtr() const {return _probConfigPtr;}
  virtual void probConfPtr(ProbConfig * pcPtr) {_probConfigPtr = pcPtr;}

  virtual void recordInstanciation(InstanciatedVar * iPtr) 
  {
  }
  
  virtual void recordInstanciation(InstanciatedConstr * iPtr) 
  {
  }
  
  virtual void model(Model * instPtr) 
  {
    _model = instPtr;
  }

  virtual Model * modelPtr()
  {
    return _model;
  }

  virtual Model * modelPtr() const
  {
    return _model;
  }
  
  virtual int dimension() const
  {
    return _dimension;
  }
  
  virtual void dimension( int dim)
  {
    _dimension = dim;
  }

  virtual const Double & genericCost(const InstanciatedVar * const ivarconstrPtr) const;
  virtual const Double & genericRhs(const InstanciatedConstr * const ivarconstrPtr) const;
  virtual bool genericCount(const InstanciatedConstr * const icPtr,
                            const InstanciatedVar * const ivPtr) const = 0;
  
  virtual const LpCoef genericCoef(const InstanciatedConstr * const icPtr,
                                   const InstanciatedVar * const ivPtr) const = 0;
  
  virtual void buildMembership(InstanciatedVar * ivPtr) 
  {
    if (printL(0)) 
      std::cout << "GenericVarConstr::buildMembership should not be called" << std::endl;
  }
  
  virtual void buildMembership(InstanciatedConstr * iconstrPtr) 
  {
    if (printL(0)) 
      std::cout << "GenericVarConstr::buildMembership should not be called" << std::endl;
  }
  
  virtual bool consecutiveVarWhenBrOnCBS(InstanciatedVar * var1Ptr, InstanciatedVar * var2Ptr) const  = 0;
  virtual const SelectionStrategy & prioritySelectionRule() const 
  {
    return _priorityRule;
  }
  
  virtual const SelectionStrategy::PriorityEnum & priorityRule() const 
  {
    return _priorityRule.selectedRule();
  }
  
  /// The higher priority dynamic constraint will be generated first
  virtual const Double & priorityLevel() const 
  {
    return _priorityLevel;
  }
  
  virtual void priorityLevel(const Double & pl);
  
  void priorityRule(const SelectionStrategy & sp) 
  {
    _priorityRule = sp;
  }
  
  virtual bool toBeUsedInPreprocessing() const 
  {
    return _toBeUsedInPreprocessing;
  }
  void toBeUsedInPreprocessing(const bool & flag) 
  {
    _toBeUsedInPreprocessing = flag;
  }
  virtual bool considerAsEqualityInPreprocessing() const
  {
    return _considerAsEqualityInPreprocessing;
  }
  void considerAsEqualityInPreprocessing(const bool & flag)
  {
    _considerAsEqualityInPreprocessing = flag;
  }

  virtual void defaultName(const std::string & name) {_defaultName = name;}
  virtual void defaultType(const char & type) {_defaultType = type;}
  virtual void defaultKind(const char & kind) {_defaultKind = kind;}
  virtual void defaultFlag(const char & flag) {_defaultFlag = flag;}
  virtual void defaultSense(const char & sense) {_defaultSense = sense;}
  virtual void defaultDirective(const char & directive) {_defaultDirective = directive;}
  virtual void defaultCostRhs(const Double & costrhs) {_defaultCostRhs = costrhs;}
  virtual void defaultUb(const Double & ub) {_defaultUb = ub;}
  virtual void defaultLb(const Double & lb) {_defaultLb = lb;}
  virtual void defaultGlobalUb(const Double & ub) {_defaultGlobalUb = ub;}
  virtual void defaultGlobalLb(const Double & lb) {_defaultGlobalLb = lb;}
  virtual void defaultVal(const Double & val) {_defaultVal = val;}
  virtual const std::string & defaultName() const {return _defaultName;}
  virtual const char & defaultType()  const {return _defaultType;}
  virtual const char & defaultKind()  const {return _defaultKind;}
  virtual const char & defaultFlag()  const {return _defaultFlag;}
  virtual const char & defaultSense()  const {return _defaultSense;}
  virtual const char & defaultDirective()  const {return _defaultDirective;}
  virtual const Double & defaultCostRhs()  const {return _defaultCostRhs;}
  virtual const Double & defaultUb()  const {return _defaultUb;}
  virtual const Double & defaultLb()  const {return _defaultLb;}
  virtual const Double & defaultGlobalUb()  const {return _defaultGlobalUb;}
  virtual const Double & defaultGlobalLb()  const {return _defaultGlobalLb;}
  virtual const Double & defaultVal()  const {return _defaultVal;}
  virtual const MultiIndexNames multiIndexNames()  const {return _multiIndexNames;}
  virtual void  multiIndexNames(const MultiIndexNames & mi)  {_multiIndexNames = mi;}

  virtual std::ostream & print(std::ostream& os = std::cout) const;


  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;
  const ControlParameters& param() const;
  ControlParameters& param();
  
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	  return compareIdentifier(VcId::GenVarConstrMask, vcIdentifier);
  }

};

inline std::ostream& operator<<(std::ostream& os, const GenericVarConstr * that)
{
  return that->print(os);
}

/**
 * Type Base4NonLinearGenericConstr recognizes  generic constraint 
 * where maste column have coefficient have implicitcefficient that 
 * are non linear function of subproblem variables: genericCount 
 * and genericCoef are defined for MastColumn *  
 */
class Base4NonLinearGenericConstr
{
 protected: 
  NonLinearGenericConstr * _nlGenericConstrPtr;
 
 public:
 Base4NonLinearGenericConstr(NonLinearGenericConstr * nlGenericConstrPtr): 
  _nlGenericConstrPtr(nlGenericConstrPtr) 
  {
  }
  
  virtual bool genericMastColumnCount(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
  virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
  virtual ~Base4NonLinearGenericConstr()
    {
    }
};

class BrVarPriorityCalcAndComp;
class GenVarGenBranchConstr;

/**
 * Generic Template for a Variable 
 */
class GenericVar: public GenericVarConstr
{
 protected:
  Double _genericBranchingOnAggregateVarPL;
  Double _compBoundSetBranchingPL;
  Double _ryanFosterBranchingPL;
  BrVarPriorityCalcAndComp * _brCmpPtr;
  Double _target4MostFractionalBranchingPriority;
  
  IndexCell2InstancVarPtrMap _indexCell2InstancVarPtrMap;
  MapInstVarPtrByMultiIndex _multiIndex2VarPtrMap;
  
  std::vector<InstanciatedVar *> _varPtrVec1D;
  std::vector<std::vector<InstanciatedVar *> > _varPtrVec2D;
  std::vector<std::vector<std::vector<InstanciatedVar *> > > _varPtrVec3D;
  bool _isDirectAccessForInstanciatedVars;
  
  std::list<GenVarGenBranchConstr *> _locVarBranchingConstrPtrList;
  CompBoundSetGenBranchConstr * _compBoundSetGenBranchConstrPtr;

 public:
  /** 
   * @param modelPtr
   * @param priorityRule 
   * @param genericBranchingOnAggregateVarPL The highest priority branching constraint 
   * will be generated first,
   * @param first 
   * @param compBoundSetBranchingPL A negative priority means that branching on aggregate 
   * value of this variable will not be considered, the highest priority branching constraint 
   * will be generated first,
   * @param genBranchingConstrPtr A negative priority means that branching on sets of component 
   * bounds on this variable will not be considered 
   * SelectionStrategy
   * @return 
   */
  GenericVar(Model * modelPtr,
  	         const BcVarConstrType::BcVcType &  gvctype,
	         ProbConfig * probConfPtr,
	         const std::string & name,
	         const MultiIndexNames & multiIndexNames = MultiIndexNames(),
	         const char & type = 'I',
	         const Double & cost = 0,
	         const Double & ub = BapcodInfinity,
	         const SelectionStrategy & branchingPriorityRule = SelectionStrategy::NotConsideredForSelection,
	         const Double & genericBranchingOnAggregateVarPL = -1.0,
	         const Double & compBoundSetBranchingPL = -1.0,
	         const char & flag = 's',
             int firstIndexMax = -1, int secondIndexMax = -1, int thirdIndexMax = -1);
  virtual ~GenericVar(); 

  const Double & target4MostFractionalBranchingPriority() const {return _target4MostFractionalBranchingPriority;}
  CompBoundSetGenBranchConstr * compBoundSetGenBranchConstrPtr() {return _compBoundSetGenBranchConstrPtr;}
  void target4MostFractionalBranchingPriority(const Double & target);

  virtual void addVarPtr2MultiIndexMap(const MultiIndex & id, InstanciatedVar * ivarPtr);
  virtual InstanciatedVar * getVarPtr(const MultiIndex & id) const;
  virtual MapInstVarPtrByMultiIndex multiIndex2VarPtrMap() const {return _multiIndex2VarPtrMap;}

  virtual InstanciatedVar * newInstanciation(const IndexCell & id,                    
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

  InstanciatedVar * newInstanciation(const IndexCell & id,                    
				     GenericVarConstr * genericVarConstrPtr,  
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
				     const bool & presetMembership = true)
  {return newInstanciation(id,
			   probConfigPtr, 
			   name,
			   costrhs,
			   sense, 
			   type,                 
			   kind,                 
			   upperBound , 
			   lowerBound,           
			   flag,                 
			   directive,            
			   priority,           
			   val,                  
			   globalUb,   
			   globalLb,            
			   presetMembership);
  } 
  virtual InstanciatedVar * checkIfInstanciationAlreadyExist(const IndexCell & id);
                 
  virtual void recordInstanciation(InstanciatedVar * iPtr);
  virtual void deleteInstanciation(InstanciatedVar * iPtr);

  const IndexCell2InstancVarPtrMap & indexCell2InstancVarPtrMap() const 
  {
    return _indexCell2InstancVarPtrMap;
  }
  virtual void setupGenericBranchingConstr();
  const Double & genericBranchingOnAggregateVarPL() const
  { 
    return _genericBranchingOnAggregateVarPL;
  }
  const Double &  compBoundSetBranchingPL() const
  { 
    return _compBoundSetBranchingPL;
  }
  virtual void genericBranchingOnAggregateVarPL(const Double & genericBranchingOnAggregateVarPLVal)
  { 
    _genericBranchingOnAggregateVarPL = genericBranchingOnAggregateVarPLVal;
  }
  virtual void compBoundSetBranchingPL(const Double & compBoundSetBranchingPLVal)
  { 
    _compBoundSetBranchingPL = compBoundSetBranchingPLVal;
  }
  virtual void ryanFosterBranchingPL(const Double & ryanFosterBranchingPLVal)
  {
    _ryanFosterBranchingPL = ryanFosterBranchingPLVal;
  }
  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr,
                            const InstanciatedVar * const ivarPtr) const;

  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr ,
                                   const InstanciatedVar * const ivarPtr) const;
  
  virtual void buildMembership(InstanciatedVar * ivPtr);
  virtual bool consecutiveVarWhenBrOnCBS(InstanciatedVar * var1Ptr, InstanciatedVar * var2Ptr) const;
  BrVarPriorityCalcAndComp * brCmpPtr() const {return _brCmpPtr;}

  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const GenericVar * that)
{
  return that->print(os);
}

/**
 * Generic Template for a Constraint
 */
class GenericConstr: public GenericVarConstr
{
 protected:
    BcAddConstrFunctor * _addConstrFunctorPtr;
    IndexCell2InstancConstrPtrMap _indexCell2InstancConstrPtrMap;
    MapInstConstrPtrByMultiIndex _multiIndex2ConstrPtrMap;

 public:

  GenericConstr(Model * modelPtr, 
		        const BcVarConstrType::BcVcType &  gvcType,
		        ProbConfig * probConfPtr,
		        const std::string & name,
		        const MultiIndexNames & multiIndexNames = MultiIndexNames(),
		        const char & sense = 'E',
		        const Double & rhs = 1,
		        const SelectionStrategy & priority = SelectionStrategy::NotConsideredForSelection,
		        const Double & priorityLevel = -1,
		        const bool & toBeUsedInPreprocessing = true,
		        const char & flag = 's');

  virtual ~GenericConstr();

  virtual void addConstrFunctorPtr(BcAddConstrFunctor * addConstrRoutinePtr);

  virtual BcAddConstrFunctor * addConstrFunctorPtr() const {return _addConstrFunctorPtr;}

  virtual InstanciatedConstr * newInstanciation(const IndexCell & id,
                                                ProbConfig* probConfigPtr,
                                                const std::string& name,
                                                const Double& rhs = 0,
                                                const char& sense = 'E',
                                                const char& type = 'C',
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

  virtual InstanciatedConstr * newInstanciation(const IndexCell& id,
                                                GenericVarConstr* genericVarConstrPtr,
                                                ProbConfig* probConfigPtr,
                                                const std::string& name,
                                                const Double& rhs = 0,
                                                const char& sense = 'E',
                                                const char& type = 'C',
                                                const char& kind = 'E',
                                                const char& flag = 's',
                                                const Double& val = 0,
                                                const Double& upperBound = BapcodInfinity,
                                                const Double& lowerBound = - BapcodInfinity,
                                                const char & directive = 'U',
                                                const Double & priority = 1.0,
                                                const bool & presetMembership = true,
                                                const bool & toBeUsedInPreprocessing = true,
                                                const bool & considerAsEqualityInPreprocessing = false)
  { 
    return newInstanciation(id, probConfigPtr, name, rhs, sense, type, kind, flag, val, upperBound, lowerBound,
                            directive, priority, presetMembership, toBeUsedInPreprocessing,
                            considerAsEqualityInPreprocessing);
  }
  virtual InstanciatedConstr * checkIfInstanciationAlreadyExist(const IndexCell & id);
  virtual void recordInstanciation(InstanciatedConstr * iPtr);
  virtual void deleteInstanciation(InstanciatedConstr * iPtr);

  const IndexCell2InstancConstrPtrMap & indexCell2InstancConstrPtrMap() const 
    {
      return _indexCell2InstancConstrPtrMap;
    }

  virtual void addConstrPtr2MultiIndexMap(const MultiIndex & id, InstanciatedConstr * iconstrPtr);
  virtual InstanciatedConstr * getConstrPtr(const MultiIndex & id) const;
  virtual MapInstConstrPtrByMultiIndex multiIndex2ConstrPtrMap() const {return _multiIndex2ConstrPtrMap;}

  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr,
                            const InstanciatedVar * const ivarPtr) const;
  
  virtual bool consecutiveVarWhenBrOnCBS(InstanciatedVar * var1Ptr,
                                         InstanciatedVar * var2Ptr) const
  { 
    bapcodInit().check(false, "GenericConstr::consecutiveVarWhenBrOnCBS should not be called");
    return true;
  }
  
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr,
                                   const InstanciatedVar * const ivarPtr) const
  { 
    return Double::staticZero;
  }
  
  virtual void buildMembership(InstanciatedConstr * iconstrPtr);
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrintAllConstraints(std::ostream& os = std::cout) const;

  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	  return compareIdentifier(VcId::GenericConstrMask, vcIdentifier);
  }
};
inline std::ostream& operator<<(std::ostream& os, const GenericConstr * that)
{
  return that->print(os);
}

/**
 * Generic Template for a Non-Linear Constraint
 */
class NonLinearGenericConstr: public GenericConstr
{
 public:
  virtual ~NonLinearGenericConstr(){}
  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr, 
			                const VarPtr2DoubleMap & varValMap) const = 0;

  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr,
				                   const InstanciatedVar * const ivarPtr) const;
  
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr , 
				                   const VarPtr2DoubleMap & varValMap) const = 0;
  
  virtual InstanciatedConstr * newInstanciation(const IndexCell & id,
                                                ProbConfig* probConfigPtr,
                                                const std::string& name,
						                        const Double& rhs,
						                        const char& sense,
						                        const char& type = 'C',
						                        const char& kind = 'E',
						                        const char& flag = 's',
						                        const Double& val = 0,
						                        const Double& upperBound = BapcodInfinity,
						                        const Double& lowerBound = - BapcodInfinity,
						                        const char & directive = 'U',
						                        const Double & priority = 1.0,
						                        const bool & presetMembership = true,
						                        const bool & toBeUsedInPreprocessing = true);

  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const NonLinearGenericConstr * that)
{
  return that->print(os);
}

/**
 * BaseType to recognize var 
 * and constraint that are generated 
 * dynamically in the course of the solution 
 * algorithm 
 */
class DynamicVarConstr
{
 public:
  virtual ~DynamicVarConstr(){}
  
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	  return compareIdentifier(VcId::DynamicVarConstrMask, vcIdentifier);
  }
};

/// Sort by priority Level
struct DynamicGenConstrSort
{
  bool operator()(DynamicGenericConstr * a, DynamicGenericConstr * b) const;
};

/// Sort by priority Level at root
struct DynamicGenConstrSortAtRoot
{
  bool operator()(DynamicGenericConstr * a, DynamicGenericConstr * b) const;
};

/**
 * Generatic template for constraint 
 * that are generated dynamically 
 * in the course of the solution algorithm 
 */
class DynamicGenericConstr: public DynamicVarConstr, public GenericConstr 
{
 protected:
  bool _oracleDefined;
  char _type;
  std::list<InstanciatedConstr *> _constrPrototypes;
  Double _rootPriorityLevel;
 public:
  /**
   * @param modelPtr = Model to which it belongs
   * @param type = 'C' for core constraint required for a correct IP formulation, 
   * 'F' for facultatif constraint use to improve the LP formulation
   * @param priorityRule strategy for prioritysing 
   * separation and sorting instanciation, the higher priority dynamic 
   * constraint will be generated first
   * @param priorityLevel the level of priority compared 
   * to other DynamicGenericConstr (the higher the most priority)
   * @paramy toBeUsedInPreprocessing == true if constraint can be used 
   * for preprocesing purposes by bapcod, == false otherwise
   **/
  DynamicGenericConstr(Model * modelPtr, 
                       ProbConfig * probConfPtr,
                       const std::string & name,
		               const char & type = 'E',
		               const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                       const Double & nonRootPriorityLevel = 1.0,
                       const Double & rootPriorityLevel = 1.0,
		               const bool & toBeUsedInPreprocessing = true);
  
  virtual ~DynamicGenericConstr();
  const bool & oracleDefined() const {return _oracleDefined;}
  void oracleDefined(const bool & flag)  {_oracleDefined = flag;}
  const char & type() const 
  {
    return _type;
  }
  
  void insertConstrPrototypes(InstanciatedConstr * cpPtr)
  {
    _constrPrototypes.push_back(cpPtr);
  }
   
  virtual bool prepareSeparation() = 0;
  
  /// thise three virtual functions are needed for the case when cut separation procedure works with phases
  /// (each subsequent phase is more expensive than preceding), in case of tailing off, we pass to the next phase
  virtual void resetSeparationPhase(){}
  virtual void increaseSeparationPhase(){}
  virtual bool separationPhaseIsMaximum() const
  {
    return true;
  }

  void rootPriorityLevel(const double & value);

  const Double & rootPriorityLevel() const
  {
    return _rootPriorityLevel;
  }

  const Double & dynamicPriorityLevel(bool rootNode) const
  {
    if (rootNode)
      return rootPriorityLevel();
    return priorityLevel();
  }

	virtual void nicePrintAllConstraints(std::ostream& os = std::cout) const;

  /**
   * Default Separation routine: 
   * need to be redefined in derived class 
   * for a more sophisticated separation
   */
  
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	  return compareIdentifier(VcId::DynamicGenericConstrMask, vcIdentifier);
  }
  
  friend struct DynamicGenConstrSort;
};

inline std::ostream& operator<<(std::ostream& os, const DynamicGenericConstr * that)
{
  return that->print(os);
}

/**
 * Generatic template for constraint 
 * that are generated as cuts in the master program
 */
class GenericCutConstr: public DynamicGenericConstr
{
protected:
  BcCutSeparationFunctor * _cutSeparationFunctorPtr;
public:
  /**
   * @param modelPtr Model to which it belongs
   * @param masterConfPtr the problem configuration to which it belongs
   * @param type = 'C' for core constraint required for a correct IP formulation, 
   * 'F' for facultatif constraint use to improve the LP formulation
   * @param priorityRule strategy for prioritysing separation and 
   * sorting instanciation
   * @param priorityLevel the level of priority copared to other 
   * DynamicGenericConstr (the higher the most priority). 
   * The higher priority dynamic constraint will be generated first 
   * @paramy toBeUsedInPreprocessing == true if constraint can be used 
   * for preprocesing purposes by bapcod, == false otherwise
   **/
  GenericCutConstr(Model * modelPtr, 
		           ProbConfig * probConfPtr,
		           const std::string & name,
		           const char& type = 'C',
		           const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
		           const Double & nonRootPriorityLevel = 1.0,
		           const Double & rootPriorityLevel = 1.0,
		           const bool & toBeUsedInPreprocessing = true);
  
  virtual ~GenericCutConstr();

  BcCutSeparationFunctor * cutSeparationFunctorPtr() const
  {
    return  _cutSeparationFunctorPtr;
  }

  void cutSeparationFunctorPtr(BcCutSeparationFunctor * cutSepfunctPtr);

  virtual bool prepareSeparation();
 
  /**
   * Default Separation routine: check violation of each constr Prototype
   * Default Separation routine: need to be redefined in derived class for a more sophisticated separation
   */
  virtual void cutSeparationRoutine(const VarPtrSet & curSol, 
				                    std::multiset < InstanciatedConstr * ,
                                                    CutSeparationPriorityComp > & generatedCutConstrSet);

  virtual void cutSeparationBasedOnFixedSol(const VarPtr2DoubleMap & oldPartialSol, const VarPtr2DoubleMap & fixedSol,
									        std::multiset < InstanciatedConstr * ,
											                CutSeparationPriorityComp > & generatedCutConstrSet);

	virtual void buildMembership(InstanciatedConstr * iconstrPtr);
};


inline std::ostream& operator<<(std::ostream& os, const GenericCutConstr * that)
{
  return that->print(os);
}

class BcCustomNonLinearCutArrayFunctor;

class GenericCustomNonLinearCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
  BcCustomNonLinearCutArrayFunctor * _userFunctorPtr;

public:
  GenericCustomNonLinearCutConstr(Model * modelPtr,
                                  ProbConfig * probConfPtr,
                                  const std::string & name,
                                  const char & type,
                                  const SelectionStrategy & priorityRule,
                                  const Double & nonRootPriorityLevel,
                                  const Double & rootPriorityLevel,
                                  BcCustomNonLinearCutArrayFunctor * userFunctorPtr);
  virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
  virtual void buildMembership(InstanciatedConstr * iconstrPtr);
  virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                    std::multiset < InstanciatedConstr * ,
                                                    CutSeparationPriorityComp > & generatedCutConstrSet);
  const LpCoef getMastColumnCoeff(CustomNonLinearCut * cutPtr, MastColumn * colPtr) const;
    
  virtual ~GenericCustomNonLinearCutConstr();
};

class BcSoftConflictsCutArrayFunctor;

class GenericSoftConflictsCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
	BcSoftConflictsCutArrayFunctor * _softConflictsCutSepFunctorPtr;
	std::map<ColGenSpConf *, GenericVar *> _genIndicVarPts;

	void updateSubprobemsWithIndicatorVarAndConstr(const std::list<BcConstr> & cutList);
public:
	GenericSoftConflictsCutConstr(Model * modelPtr,
								  ProbConfig * probConfPtr,
								  const std::string & name,
								  const char & type,
								  const SelectionStrategy & priorityRule,
								  const Double & nonRootPriorityLevel,
								  const Double & rootPriorityLevel,
								  BcSoftConflictsCutArrayFunctor * softConflictsCutSepFunctorPtr);
	virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
	virtual bool prepareSeparation();
	virtual void buildMembership(InstanciatedConstr * iconstrPtr);
	virtual void cutSeparationBasedOnFixedSol(const VarPtr2DoubleMap & oldPartialSol, const VarPtr2DoubleMap & fixedSol,
											  std::multiset < InstanciatedConstr * ,
													  CutSeparationPriorityComp > & generatedCutConstrSet);
	virtual void cutSeparationRoutine(const VarPtrSet & curSol,
									  std::multiset < InstanciatedConstr * ,
											  CutSeparationPriorityComp > & generatedCutConstrSet);
	const LpCoef getMastColumnCoeff(SoftConflictsCut * cutPtr, MastColumn * colPtr) const;

	virtual ~GenericSoftConflictsCutConstr();
};

#endif // GenVarConstrClasses_h


