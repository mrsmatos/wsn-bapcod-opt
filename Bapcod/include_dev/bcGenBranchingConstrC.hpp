/**
 *
 * This file bcGenBranchingConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef GenBranchingConstrClasses_h
#define GenBranchingConstrClasses_h

#include "bcVcIdentifierC.hpp"

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcGenVarConstrC.hpp"

class BranchingConstrBaseType;
class CompSetInstMastBranchConstr;
class Base4NonLinearConstraint;
class BranchingConstrGenerator;
struct BranchingSeparationPriorityComp;
class BcDisjunctiveBranchingConstrSeparationFunctor;
class BcVarBranchingPriorityFunctor;

typedef std::set<BranchingConstrGenerator *, BranchingSeparationPriorityComp> BranchGeneratorsSet;

typedef std::vector<CompSetInstMastBranchConstr *> ColClassesVector;
typedef std::map<ProbConfig *, ColClassesVector > ColClassesVectorsMap;

typedef std::map < IndexCell, InstanciatedVar * , IndexCellSort > IndexCell2InstancVarPtrMap;
typedef std::map < IndexCell, InstanciatedConstr * , IndexCellSort > IndexCell2InstancConstrPtrMap;

/**
 * Sorting procedure to select 
 * among branching constraint: 
 * by genericVarConstr _priorityLevel 
 * than by constraint PriorityRule 
 * and priority of violation or fract val.
 */
struct BranchingSeparationPriorityComp
{
  bool operator()(BranchingConstrGenerator * a, BranchingConstrGenerator * b) const;
};

/**
 * Sorting branching generators by score
 */
struct BranchingSeparationScoreComp
{
  bool operator()(BranchingConstrGenerator * a, BranchingConstrGenerator * b) const;
};


/// Virtual Generator of child node branching constraints
class BranchingConstrGenerator
{
 protected:
  GenericBranchingConstr * _genericBrConstrPtr;
  int _ref;
  std::string _description;
  
  /// Records branching direction (up or down): 'U' or 'D'
  char _direction;

  /**
   * Specific constraint that serves  as a prototype for branching constraints
   */
  InstanciatedConstr * _prototypeInstConstrPtr;

  /**
   * Intermetidate cmputation attributes
   * Fractional rhs to be rounded up or down
   * Internal counter for function nextNodeBrConstr()
   */
  Double _candidateLhs;
  
  int _childNbCounter;

  void instanciateBrConstr(const int & parentNodeNb,
                           const int & childNb,
                           const Double & rhs,
                           const char & sense,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);
 public:
  BranchingConstrGenerator(GenericBranchingConstr * gbcPtr,
			               const char & priorityDir = 'U',
			               const Double & candidateLhs = 0,
			               InstanciatedConstr * prototypeInstConstrPtr = NULL,
                           const std::string & description = "");
  int ref() const {return _ref;}

  BranchingConstrGenerator(const BranchingConstrGenerator & that);
    
  virtual ~BranchingConstrGenerator();
  
  virtual BranchingConstrGenerator * clone() const
  {
    return new BranchingConstrGenerator(*this);
  }

  virtual GenericBranchingConstr * genericBrConstrPtr() const
  {
    return _genericBrConstrPtr;
  }
    
  /**
   * Assuming  branching on a fractional 
   * aggregate auxiliary variable
   */
  virtual const Double violation() const
  {
    return fracPart();
  }
  
  virtual const Double priority() const
  {
    return violation();
  }
  
  virtual const char direction() const
  {
    return _direction;
  }
  
  virtual const Double fracPart() const
  {
    return Dfrac(_candidateLhs);
  }
  
  virtual const Double lFracPart() const
  {
    return Lfrac(_candidateLhs);
  }

  virtual const Double fracPartRelativeToTarget() const
  {
    return Lfrac(_candidateLhs);
  }
  
  virtual const Double uFracPart() const
  {
    return Ufrac(_candidateLhs);
  }

  virtual const Double defaultComparisonFactor() const
  {
    return ref();
  }

  virtual void prototypeInstConstrPtr(InstanciatedConstr * icPtr) 
  {
    _prototypeInstConstrPtr = icPtr;
  }
  
  virtual InstanciatedConstr * prototypeInstConstrPtr() const 
  {
    return _prototypeInstConstrPtr;
  }
  
  virtual void reset()
  {
    _childNbCounter = 0;
  }
  virtual void computeLhs(const SolutionVarInfoPtrList & curSol);
  virtual const int numberOfChilds() const
  {
    return _childNbCounter;
  }
  
  virtual bool nextNodeBrConstr(Node * parentNodePtr, 
                                std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
                                const ConstrPtrSet & existingMasterBranchingConstr);

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;
  BapcodInit & bapcodInit() const;

  /// this comparison operator needed to retrieve the branching history of the generator
  /// we use description string for comparison
  inline bool operator<(const BranchingConstrGenerator & otherGenerator) const
  {
    return (_description < otherGenerator._description);
  }
};

inline std::ostream& operator<<(std::ostream& os, const BranchingConstrGenerator & that)
{
  return that.print(os);
}

struct BranchingPhaseEvoluationInfo
{
  std::vector<bool> exists;
  std::vector<double> lpValueImprovement;
  std::vector<double> evaluationTime;
  
  inline unsigned int size() const
  {
    return exists.size();
  }
  
  inline void resize(unsigned int n)
  {
    exists.resize(n, false);
    lpValueImprovement.resize(n, 0.0);
    evaluationTime.resize(n, 0.0);
  }
  
  inline void setEvalEntry(unsigned int n, const double & lpValueImprovement_,
                           const double & evaluationTime_)
  {
    exists[n] = true;
    lpValueImprovement[n] = lpValueImprovement_;
    evaluationTime[n] = evaluationTime_;
  }
    
};

struct BranchingGeneratorHistory;

struct BranchingEvaluationInfo
{
  int depth;
  double dualBound;
  double lhsFracPart;
  char direction;
  std::vector<BranchingPhaseEvoluationInfo> phaseInfo;
  BranchingGeneratorHistory * historyPtr;
  
  BranchingEvaluationInfo(const int & depth_, const Double & dualBound_,
                          const Double & lhsFracPart_, const char & direction_,
                          BranchingGeneratorHistory * const historyPtr_);

  void addEvalEntry(const int & phaseNumber, const int & childNumber,
                    const double & lpValue, const double & evalTime);
};

struct BranchingGeneratorHistory
{
  std::vector<BranchingEvaluationInfo *> evaluationsInfo;
  std::vector< std::vector<double> > pseudoCosts;
  std::vector< std::vector<int> > numEvaluations;
};

struct BranchingGeneratorPtrComp
{
  bool operator()(const BranchingConstrGenerator * const genAPtr,
                  const BranchingConstrGenerator * const genBPtr) const
  {
    return (*genAPtr < *genBPtr);
  }
};

typedef std::map<BranchingConstrGenerator *, BranchingGeneratorHistory,
                 BranchingGeneratorPtrComp> BranchingGeneratorHistoryMap;


/**
 * Generator of child node branching 
 * disjunctive constraints on Generic 
 * Var instanciated components
 */
class GenVarBranchConstrGenerator: public BranchingConstrGenerator
{
  /// Specific intantianted var on which we branch
  Variable * _varPtr;

  void instanciateBrConstr(const int & parentNodeNb,
                           const int & childNb,
                           const Double & rhs,
                           const char & sense,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);

public:
 GenVarBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
			                 Variable * varPtr,
                             const char & priorityDir,
			                 const Double & candidateLhs);
 GenVarBranchConstrGenerator(const GenVarBranchConstrGenerator & that);
  
 virtual ~GenVarBranchConstrGenerator();
  
 virtual BranchingConstrGenerator * clone() const
 {
   return new GenVarBranchConstrGenerator(*this);
 }

  virtual bool nextNodeBrConstr(Node * parentNodePtr,
				                std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
				                const ConstrPtrSet & existingMasterBranchingConstr);

  virtual void computeLhs(const SolutionVarInfoPtrList & curSol);

  virtual const Double violation() const
  {
    return(Dfrac(_candidateLhs));
  }
  
  virtual const Double priority() const
  {
    return _varPtr->priority();
  }
  
  virtual const Double fracPart() const
  {
    return Dfrac(_candidateLhs);
  }
  
  virtual const Double lFracPart() const
  {
    return Lfrac(_candidateLhs);
  }

  virtual const Double fracPartRelativeToTarget() const
  {
    return _varPtr->fracPartRelativeToTarget();
  }
  
  virtual const Double uFracPart() const
  {
    return Ufrac(_candidateLhs);
  }
  virtual const Double defaultComparisonFactor() const
  {
    return priority();
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const GenVarBranchConstrGenerator & that)
{
  return that.print(os);
}

/**
 * Generator of child node branching disjunctive 
 * constraints for branching on a pari od SP 
 * instanciatedVar components
 */
class RyanAndFosterBranchConstrGenerator: public BranchingConstrGenerator
{
  /// Specific instanciated first var on which we branch
  InstanciatedVar * _variPtr;

  /// Specific instanciated second var on which we branch
  InstanciatedVar * _varjPtr;

  void instanciateBrConstr(const int & parentNodeNb,
                           const int & childNb,
                           const Double & rhs,
                           const char & sense,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);

public:
  RyanAndFosterBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
                                     InstanciatedVar * variPtr,
                                     InstanciatedVar * varjPtr,
				                     const Double & candidateLhs,
				                     const char & priorityDir = 'U');
  RyanAndFosterBranchConstrGenerator(const RyanAndFosterBranchConstrGenerator & that);
 
  virtual ~RyanAndFosterBranchConstrGenerator();

  virtual BranchingConstrGenerator * clone() const
  {
    return new RyanAndFosterBranchConstrGenerator(*this);
  }

  virtual bool nextNodeBrConstr(Node * parentNodePtr, 
				                std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
			 	                const ConstrPtrSet & existingMasterBranchingConstr);

  virtual void computeLhs(const SolutionVarInfoPtrList & curSol);

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const RyanAndFosterBranchConstrGenerator & that)
{
  return that.print(os);
}


/**
 * Generator of child node branching constraints 
 * for branching on component sets
 */
class CompBdSetBranchConstrGenerator: public BranchingConstrGenerator
{
  ComponentSequence _compBoundSet;
  ComponentSequence _currentCompBoundSet;

  void instanciateBrConstr(const int & parentNodeNb,
                           const int & childNb,
                           const ComponentSequence & compBdSeq,
                           std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);

public:
  CompBdSetBranchConstrGenerator(GenericBranchingConstr * gbcPtr,
                                 const ComponentSequence & compBoundSet,
				                 const char & priorityDir = 'U');
  CompBdSetBranchConstrGenerator(const CompBdSetBranchConstrGenerator & that);
  virtual ~CompBdSetBranchConstrGenerator();

  virtual BranchingConstrGenerator * clone() const
  {
    return new CompBdSetBranchConstrGenerator(*this);
  }

  virtual bool nextNodeBrConstr(Node * parentNodePtr,
				                std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
				                const ConstrPtrSet & existingMasterBranchingConstr);

  virtual void computeLhs(const SolutionVarInfoPtrList & curSol);
  virtual void reset()
  {
    BranchingConstrGenerator::reset();
    _currentCompBoundSet = _compBoundSet;
  }

  virtual const Double violation() const
  {
    return(Dfrac(_compBoundSet.fracWeight()));
  }
  
  virtual const Double priority() const;

  virtual const Double fracPart() const
  {
    return Dfrac(_compBoundSet.fracWeight());
  }
  
  virtual const Double lFracPart() const
  {
    return (_compBoundSet.fracWeight() - _compBoundSet.fracWeightRD());
  }

  virtual const Double fracPartRelativeToTarget() const
  {
    return uFracPart();
  }
  
  virtual const Double uFracPart() const
  {
    return (_compBoundSet.fracWeightRU() - _compBoundSet.fracWeight());
  }

  virtual const Double defaultComparisonFactor() const
  {
    return _compBoundSet.additionalNbOfCpBd();
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const CompBdSetBranchConstrGenerator & that)
{
  return that.print(os);
}

/**
 * Generatic template for constraint that are generated when branching
 */
class GenericBranchingConstr: public DynamicGenericConstr
{
 protected:
  BranchingGeneratorHistoryMap _branchingHistory;

  BcDisjunctiveBranchingConstrSeparationFunctor  * _branchingSeparationFunctorPtr;

public:
  GenericBranchingConstr(Model * modelPtr,
			             ProbConfig * probConfPtr,
			             const std::string & name,
			             const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                         const Double & nonRootPriorityLevel = 1.0,
                         const Double & rootPriorityLevel = 1.0,
			             const bool & toBeUsedInPreprocessing = true);

  virtual ~GenericBranchingConstr();

  BranchingGeneratorHistoryMap & branchingHistory()
  {
    return _branchingHistory;
  }

  BcDisjunctiveBranchingConstrSeparationFunctor * branchingSeparationFunctorPtr() const
  {
    return  _branchingSeparationFunctorPtr;
  }

  void branchingSeparationFunctorPtr(BcDisjunctiveBranchingConstrSeparationFunctor * branchSepfunctPtr) 
  {
    _branchingSeparationFunctorPtr = branchSepfunctPtr;
    oracleDefined(true);
  }

  double computeLhs(const MasterVarSolution & curListOfMastVarInSol, const InstanciatedConstr * instConstrPtr);

  virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  virtual void branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                 const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  virtual void cutSeparationRoutine(const VarPtrSet & curSol,
				    std::multiset < InstanciatedConstr * , CutSeparationPriorityComp > & generatedCutConstrSet) 
  {
    bapcodInit().check(1, "GenericBranchingConstr::cutSeparationRoutine(const VarPtrSet & curSol) not defined");
    return;
  }
  
  virtual bool prepareSeparation();
  
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
	return compareIdentifier(VcId::GenericBranchingConstrMask, vcIdentifier);
  }
};

inline std::ostream& operator<<(std::ostream& os, const GenericBranchingConstr * that)
{
  return that->print(os);
}

class GenAggrSubProbVarBranchingConstr: public GenericBranchingConstr
{
  std::string _genVarName;
  double _targetFraction;
  int _numIgnoredIndices;
  int _numGeneratedBrConstrs;
  BcVarBranchingPriorityFunctor * _branchingVarPriorityFunctorPtr;

public:
  GenAggrSubProbVarBranchingConstr(Model * modelPtr,
			                       ProbConfig * probConfPtr,
			                       const std::string & name,
                                   const std::string & genVarName,
                                   const double & targetFraction = 0.5,
                                   const int & numIgnoredIndices = 0,
                                   const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                                   const Double & nonRootPriorityLevel = 1.0,
                                   const Double & rootPriorityLevel = 1.0,
			                       const bool & toBeUsedInPreprocessing = true);

  virtual ~GenAggrSubProbVarBranchingConstr();

  const BcVarBranchingPriorityFunctor * branchingVarPriorityFunctorPtr() const
  {
    return _branchingVarPriorityFunctorPtr;
  }

  void branchingVarPriorityFunctorPtr(BcVarBranchingPriorityFunctor * branchVarPriorityFuncPtr);

  virtual void branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                 const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);
};

/**
 * Generatic template for constraint 
 * that are generated when branching 
 * on variables instanciated from genericVar
 */
class GenVarGenBranchConstr: public GenericBranchingConstr
{
 protected:
  /// Pointer to the associated generic Var
  GenericVar * _genVarPtr;
 public:
  GenVarGenBranchConstr(Model * modelPtr, 
                        ProbConfig * probConfPtr,
                        GenericVar * genVarPtr,
                        const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
			            const Double & priorityLevel = 1.0);
  
  virtual ~GenVarGenBranchConstr();

  GenericVar * genVarPtr() const 
  {
    return _genVarPtr;
  }

  BrVarPriorityCalcAndComp * brCmpPtr() const
  {
    return _genVarPtr->brCmpPtr();
  }

  virtual const Double & genericCostRhs(const InstanciatedVarConstr * const ivarconstrPtr) const;
  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr, const InstanciatedVar * const ivarPtr) const;
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr,
                                   const InstanciatedVar * const ivarPtr) const;
  virtual void buildMembership(InstanciatedConstr * iconstrPtr);

  virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
						                         BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  virtual void branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                 const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
                                                 BranchGeneratorsSet & generatedBrConstrGeneratorSet);
 
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
  {
    return compareIdentifier(VcId::GenVarGenBranchConstrMask, vcIdentifier);
  }
};

inline std::ostream& operator<<(std::ostream& os, const GenVarGenBranchConstr * that)
{
  return that->print(os);
}

struct FurtherSplitCandidate
{
  Variable * varPtr;
  ComponentBound compBound;
  int nbOfSatisfyingFracCol;
  int nbOfNonSatisfyingFracCol;
  FurtherSplitCandidate(Variable * vPtr, const ComponentBound & cBound, const int & satCol, const int & nonSatCol) :
      varPtr(vPtr), compBound(cBound), nbOfSatisfyingFracCol(satCol), nbOfNonSatisfyingFracCol(nonSatCol)
  {
  }
};

/**
 * Generatic template for constraint 
 * that are generated when branching 
 * on component sets
 */
class CompBoundSetGenBranchConstr: public GenVarGenBranchConstr, public Base4NonLinearGenericConstr
{
 public:
  CompBoundSetGenBranchConstr(Model * modelPtr, 
			      GenericVar * genVarPtr, 
			      const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
			      const Double & priorityLevel = 1.0);
  
  virtual ~CompBoundSetGenBranchConstr();

  virtual const Double & genericCostRhs(const InstanciatedVarConstr * const ivarconstrPtr) const
  {
    bapcodInit().check(true,"GenVarGenBranchConstr::genericCostRhs(): error should not be called");
    return Double::staticZero;
  }

  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr, const InstanciatedVar * const ivarPtr) const;
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr ,
                                   const InstanciatedVar * const ivarPtr) const;
  virtual void buildMembership(InstanciatedConstr * iconstrPtr);
  virtual bool genericMastColumnCount(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
  virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;

  virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  static void ILOsortMastColumn(const ColClassesVector & treeOfColClasses,
                                const MasterColSolution & curListOfMasterCol,
                                ComponentSequence & curClassCompBoundSet,
                                MasterColSolution & ILOsortedListOfMastCol);
  
  /**
   * Takes columns of curlistofofmastcol 
   * and returns them in a sorted list 
   * ILOsortedlistofmastcol in the lexicographic 
   * order dictated by
   * The class partition of current branching constaints 
   * csbrconstrlist, whose component bound already treated 
   * 
   * @param treeOfColClasses C
   * @param curListOfFractMastCol F 
   * @param candVarForCompBoundBranching I
   * @param curClassCompBoundSet S, p is implicitly given by the size of S
   * @param totFracColVal Records the overall value of the fractional columns in the current master s
   * @param additionalNbOfCpBd incremented  in each recursive call to separate
   * @param recordedBestBranchingSet Record of the selected component set for branching along with charateristics

   */
  virtual void exploreFrac(const ColClassesVector & treeOfColClasses,
                           const MasterColSolution & curListOfFractMastCol,
                           VarPtrSet & candVarForCompBoundBranching, 
                           ComponentSequence & curClassCompBoundSet,
                           Double & totFracColVal, int additionalNbOfCpBd,
			               BranchGeneratorsSet & generatedBrConstrGeneratorSet);
  /**
   * @param curListOfFractMastCol F
   * @param candVarForCompBoundBranching I 
   * @param curClassCompBoundSet S, p is implicitly given by the size of S
   * @param totFracColVal Records the overall value of the fractional 
   * columns in the current master solution
   * @param additionalNbOfCpBd Incremented  in each recursive call to separate
   * @param recordedBestBranchingSet Record of the selected component set for branching along with charateristics
   */
  virtual void separateFrac(const MasterColSolution & curListOfFractMastCol,
			                VarPtrSet & candVarForCompBoundBranching,
			                ComponentSequence & curClassCompBoundSet,
                            Double & totFracColVal,
			                int additionalNbOfCpBd,
        			        BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  void separateCardinality(const MasterColSolution & curListOfFractMastCol,
                           ComponentSequence & curClassCompBoundSet,
                           BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  void updateGeneratedBrConstrGeneratorSet(const ComponentSequence & curCPSet,
					   const ComponentBound & lastCP, 
					   const int & additionalNbOfCpBd,
					   const bool & pushAdditionalCP,
                       const Double & cardinality,
					   CompBoundSetGenBranchConstr * curCompBoundSetGenBranchConstr, 
					   BranchGeneratorsSet & generatedBrConstrGeneratorSet);

  
  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const CompBoundSetGenBranchConstr * that)
{
  return that->print(os);
}

class RyanAndFosterGenBranchConstr: public GenVarGenBranchConstr
{
  public:
  RyanAndFosterGenBranchConstr(Model * modelPtr,
			                   GenericVar * genVarPtr,
			                   const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                               const Double & priorityLevel = 1.0);

  virtual ~RyanAndFosterGenBranchConstr();
  
  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr,
                            const InstanciatedVar * const ivarPtr) const;
  
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr , 
				                   const InstanciatedVar * const ivarPtr) const;

  virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);
  
  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const RyanAndFosterGenBranchConstr * that)
{
  return that->print(os);
}

#endif // GenVarConstrClasses_h

