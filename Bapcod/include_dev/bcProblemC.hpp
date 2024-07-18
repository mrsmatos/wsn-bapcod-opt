/**
 *
 * This file bcProblemC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef ProblemClasses_h
#define ProblemClasses_h
#include "bcUsefulHeadFil.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"

class LPform;
class ProbConfig;
class BcSolution;
class BcDualSolution;
class Node;
class BcSolverOracleFunctor;
class BcFracSolBasedHeuristicFunctor;
class BcMasterHeuristicFunctor;
class BcDivingFixingFunctor;
class BcEnumSolBasedHeuristicFunctor;
class BcSolverOracleInfo;
class BcRyanAndFosterBranchConstr;
class BcSoftConflictsCut;

struct VarPtr_MpFormIndexStatus
{
  Variable * _varPtr;
  int _statusInBasicSol;

  VarPtr_MpFormIndexStatus(Variable * varPtr, const int & statusInBasicSol) :
    _varPtr(varPtr), _statusInBasicSol(statusInBasicSol)
  {
  }
} ;

struct ConstrPtr_MpFormIndexStatus
{
  Constraint * _constrPtr;
  int _statusInBasicSol;

  ConstrPtr_MpFormIndexStatus(Constraint * constrPtr, const int & statusInBasicSol) :
    _constrPtr(constrPtr), _statusInBasicSol(statusInBasicSol)
  {
  }

} ;

struct LpBasisRecord
{
  std::string _name;
  std::vector< VarPtr_MpFormIndexStatus > _varInBasis;
  std::vector< ConstrPtr_MpFormIndexStatus > _constrInBasis;

  LpBasisRecord() :
    _name("basis"), _varInBasis(), _constrInBasis()
  {
  }

  LpBasisRecord(const  std::string & name) :
    _name(name), _varInBasis(), _constrInBasis()
  {
  }

  LpBasisRecord(const LpBasisRecord & that) :
    _name(that._name),
    _varInBasis(that._varInBasis),
    _constrInBasis(that._constrInBasis)
  {
  }

  LpBasisRecord(const LpBasisRecord & that, const  std::string & name) :
    _name(name),
    _varInBasis(that._varInBasis),
    _constrInBasis(that._constrInBasis)
  {
  }

  virtual ~LpBasisRecord()
  {
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
	  os << "LpBasisRecord " << _name;
      os << " #v = " << _varInBasis.size() << " : ";
      for (std::vector<VarPtr_MpFormIndexStatus>::const_iterator it = _varInBasis.begin();
           it != _varInBasis.end(); ++it)
          os << it->_varPtr->name() << "(" << it->_statusInBasicSol << "), ";
      os << std::endl;
      os << " #c = " << _constrInBasis.size() << " : ";
      for (std::vector<ConstrPtr_MpFormIndexStatus>::const_iterator it = _constrInBasis.begin();
           it != _constrInBasis.end(); ++it)
          os << it->_constrPtr->name() << "(" << it->_statusInBasicSol << "), ";
      os << std::endl;
    return os;
  }
    
  void clear(const bool removeMarksInVars = false, const bool removeMarksInConstrs = false);
} ;

inline std::ostream& operator<<(std::ostream& os, const LpBasisRecord & that)
{
  return that.print(os);
}

/// needed for partial solution
struct VariableSolInfo
{
public:
  Variable * const varPtr;
  Double value;

  VariableSolInfo(Variable * varPtrV, const Double & val) :
      varPtr(varPtrV), value(val)
  {
  }

  virtual ~VariableSolInfo()
  {
  }

  virtual void applyVarInfo() const;
};

/// Problem is an intermediate layer between ProbConfig and Form
/// In principle, Problem may be merged with ProbConfig
class Problem
{
  friend class MipProblem;
  friend class BcSolution;

public:

  /// Methods to retrieve and load advanced basis
  void retrieveBasis(LpBasisRecord * curBasisPtr, bool markInVars = false, bool markInConstrs = false);
  void reloadMemorizedBasis(LpBasisRecord * basisPtr);

protected:
  int _ref;
  std::string _name;
  ProbConfig * _probConfPtr;
  bool _solverOracleFunctorIsDefined;
  BcSolverOracleFunctor * _solverOracleFunctorPtr;
  BcFracSolBasedHeuristicFunctor * _fracSolBasedHeuristicFunctorPtr;
  BcMasterHeuristicFunctor * _masterHeuristicFunctorPtr;
  BcDivingFixingFunctor * _divingFixingFunctorPtr;
  BcEnumSolBasedHeuristicFunctor * _enumSolBasedHeuristicFunctor;

  BcObjStatus::MinMaxIntFloat _objStatus;
  bool _probInfeasibleFlag;

  SolutionMethod _solMode;
  Double _minCost;
  Double _maxCost;
private:
  bool _probIsBuilt;
  BaseFormulation * _formulationPtr;
  LPform * _primalFormulationPtr;
  GlobalArtificialVar * _posGlobalArtVarPtr;
  GlobalArtificialVar * _negGlobalArtVarPtr;

protected:
  Double _objVal;
  Double _primalBound;
  Double _dualBound;

  ConstrIndexManager _probConstrManager;
  VarIndexManager _probVarManager;

  /// Holds non linear constraints in the current problem
  ConstrPtrSet _probNonLinearConstrSet;

  /// Information needed to pass on problem solution
  Solution * _primalSolPtr;
  DualSolution * _dualSolPtr;

  VarPtrSet _inPrimalSol;
  VarPtrSet _nonZeroRedCostVars;
  ConstrPtrSet _inDualSol;
  Double _partialSolutionValue;
  VarPtr2DoubleMap _partialSolution;
  std::list <Solution *> _recordedSolPtr;
  SolutionStatus _requiredStatus;
  SolutionStatus _probStatus;
  bool _LPpreprocessorOn;
  bool _LPprobingOn;
  char _LPsolverSelect;
  
  double _rightHandSideZeroTol;
  double _reducedCostTolerance;

  /// start added by Ruslan: needed for new preprocessing
  ConstrPtrList _preprocessedConstrsList;
  VarPtrList _preprocessedVarsList;
  const Node * _curNodePtr;
  /// end added by Ruslan
    
  /// added by Issam for more efficiency and to fix bug #
  /// after columns are cleaned we can t ask for red costs
  /// before the MIPSolver solves the master again.
  /// It is put to true in retrieveRedCosts()
  /// It is put to false in resetSolution()
  bool _isRetrievedRedCosts;

public:

  int ref() const {return _ref;}
  const std::string & name() const {return _name;}

  virtual void setMIPRequiredStatus(const SolutionStatus & newStatus);

  Problem(const Problem & that) = delete;

  Problem(const int & ref,
          const double & rightHandSideZeroTol,
          const double & reducedCostTolerance,
          const BcObjStatus::MinMaxIntFloat & minmaxStatus = BcObjStatus::minInt,
          const SolutionMethod & solMode = SolutionMethod::lpSolver,
          const std::string & name = " ",
          const SolutionStatus & requiredSolStat = SolutionStatus::Optimum,
          const bool & LPpreprocessorOn = true,
          const bool & LPprobingOn = true,
          const char & LPsolverSelect = 'd');

  void solverOracleFunctorPtr(BcSolverOracleFunctor * solverOracleFunctPointer);
  const BcSolverOracleFunctor * solverOracleFunctorPtr();
  bool prepareSolverOracleFunctor();
  bool solverOracleFunctorDefined();
  bool getEnumeratedStatus();
  void masterHeuristicFunctorPtr(BcMasterHeuristicFunctor * masterHeuristicFunctPointer);
  void fracSolBasedHeuristicFunctorPtr(BcFracSolBasedHeuristicFunctor * heurFunctorPtr);
  void divingFixingFunctorPtr(BcDivingFixingFunctor * divingFixingFunctPointer);
  void enumSolBasedHeuristicFunctorPtr(BcEnumSolBasedHeuristicFunctor * enumSolBasedHeuristicFunctor);
  bool enumSolBasedHeuristicFunctorDefined();
  bool runEnumSolBasedHeuristicFunctor(const std::vector<Solution *> & enumSolPts, const Bound & incPrimalIpBound);

  bool runMasterHeuristicFunctor();
  bool runFracSolBasedHeuristicFunctor(const std::list<VariableSolInfo> & partialSol,
                                       const SolutionVarInfoPtrList & primalSol);
  bool divingFixingFunctorDefined();
  void runDivingColCutGenTerminatedFunctor(const SolutionVarInfoPtrList & primalSol, const VarPtrSet & tabuVarSet,
                                           std::vector<MastColumn *> & columnsToAdd);
  Solution * runDivingFixingFunctor(const std::list<VariableSolInfo> & partialSol,
                                    const SolutionVarInfoPtrList & primalSol, const VarPtrSet & tabuVarSet);

  virtual ~Problem();

  const BcObjStatus::MinMaxIntFloat & objStatus() const
  {
    return  _objStatus;
  }

  const bool & probInfeasibleFlag() const
  {
    return  _probInfeasibleFlag;
  }

  const ProgStatus & progStatus() const;
  
  const ControlParameters & param() const;

  virtual void defineFormulation();

  ProbConfig * probConfPtr() const
  {
    return _probConfPtr;
  }

  virtual void probConfPtr(ProbConfig * ptr);

  virtual void recordPosGlobalArtVar(GlobalArtificialVar * globalArtVarPtr);
  virtual void recordNegGlobalArtVar(GlobalArtificialVar * globalArtVarPtr);

  virtual GlobalArtificialVar * posGlobalArtVarPtr() const
  {
    return _posGlobalArtVarPtr;
  }

  virtual GlobalArtificialVar * negGlobalArtVarPtr() const
  {
    return _negGlobalArtVarPtr;
  }

  void addConstraintToProblem(Constraint * constrPtr);

  void insertActiveConstr(Constraint * constrPtr, const int & updateForm);
  void insertConstr(Constraint * constrPtr, const VcIndexStatus::VcStatus & status);

  virtual void updateLocArtVarList(Constraint * constrPtr, std::list<Variable *> & varPtrList);

  template<typename Container>
  void addConstrSet(const Container & newConstrSet, const int & flag, const int & updateForm);
  virtual void addConstr(Constraint * constrPt, const int & flag = 1, const int & updateForm = 0);

  template<typename Container>
  void addVarSet(const Container & newVarSet, const int & flag,const int & updateForm);
  virtual int addVar(Variable * varPtr, const int & flag = 1, const int & updateForm = 0);

  template<typename Container>
  void delConstrSet(const Container & oldConstrSet, const int & flag, const int & updateForm);
  virtual void delConstr(Constraint * constrPtr, const int & flag = 3, const int & updateForm = 0);

  template<typename Container>
  void delVarSet(const Container &  oldVarSet, const int & flag, const int & updateForm);
  virtual void delVar(Variable * varPtr, const int & flag = 3, const int & updateForm = 0);

  virtual void deleteForm();

  void deactivateAndRemoveAllVarsAndConstrsFromMemory(); /// to be used from the destructor
  void removeUnusedDynamicVarsFromMemory(const bool & checkArtVars = false);
  void removeUnusedDynamicConstrsFromMemory();

  const VarPtrSet & nonZeroRedCostVars() const
  {
    return _nonZeroRedCostVars;
  }

  const VarPtrSet & inPrimalLpSol() const
  {
    return _inPrimalSol;
  }

  const ConstrPtrSet & inDualSol() const
  {
    return _inDualSol;
  }

  const Double & minCost()  const
  {
    return _minCost;
  }

  const Double & maxCost()  const
  {
    return _maxCost;
  }

  const bool & preprocessorOn()  const
  {
    return  _LPpreprocessorOn;
  }

  const bool & probingOn()  const
  {
    return  _LPprobingOn;
  }
  
  void addVarsSimplyInForm(std::list<Variable *> & varsToAddToForm);
  void addVarsSimplyInForm(std::list<InstanciatedVar *> & varsToAddToForm);
  void addConstrsSimplyInForm(std::list<Constraint *> & constrsToAddToForm);
  void delVarsSimplyInForm(std::list<Variable *> & varsToRemovefromForm);
  void delConstrsSimplyInForm(std::list<Constraint *> & constrsToRemovefromForm);
  
protected:
  bool checkIfConstrTempValIsFeasible();
  void recordCurRhsInConstrTempVal();

  virtual void removeActiveConstr(Constraint * constrPtr, const int & updateForm = 0);
  virtual void removeVar(Variable * varPtr, const int & updateForm = 0);
public:
  virtual void setConstr2Form(Constraint * constrPtr);
  virtual void setVar2Form(Variable * varPtr);
  virtual void unsetVar2Form(Variable * varPtr);
  virtual void unsetConstr2Form(Constraint * constrPtr);
  void addConstrInForm();
  void addVarInForm();
  void delConstrInForm();
  void delVarInForm();

  virtual const std::list <Solution *> & recordedSolList() const
  {
    return _recordedSolPtr;
  }

  void storeDualSolution(ConstrPtr2DoubleMap & dualSolution) const;
  void resetDualSolution(const ConstrPtr2DoubleMap & dualSolution);
  virtual void resetSolution(const char & flag = 'p');

  virtual BaseFormulation * formulationPtr() const
  {
    return _formulationPtr;
  }

  virtual LPform * primalFormulationPtr() const
  {
    return _primalFormulationPtr;
  }

protected:
  virtual void addConstr2Prob(Constraint * constrPtr);
  virtual void addVar2Prob(Variable * varPtr);
  virtual void delVarFromProb(Variable * varPtr);
  virtual void delConstrFromProb(Constraint * constrPtr);
public:
  virtual const Double & objVal() const ;
  virtual const Double & primalBound() const ;
  void setPrimalLpBound(const Double & bd);
  virtual const Double & dualBound() const ;

  void setDualBound(const Double & bd);
  virtual const Double & partialSolutionValue() const;
  virtual void updatePartialSolution(Variable * varPtr, const Double & value = 1);
  virtual const VarPtr2DoubleMap & partialSolution() const;

  virtual const SolutionStatus & probStatus() const
  {
    return _probStatus;
  }
  virtual void setProbStatus(const int & stat);
  virtual void setProbStatus(const SolutionStatus & stat);
  virtual void setStatusAfterSol();

  /**
   *  Prepare matrix and generate formulation
   */
  virtual void buildProblem();

  void updateProbVarCostWithSoftConflictCuts();
  bool rankOneCutsArePresent();

   /**
   * Compute upper and lower bound on obj value; return true if  problem infeasible
   */
  virtual bool updateProbForColGen(bool inPurePhaseOne);

  virtual bool updateProblem();
  virtual void hardResetObjective(char flag);
  virtual void hardResetConstraintRHS(char flag);

  virtual bool updateProbConstr(char flag);
  virtual bool updateProbVar(bool inPurePhaseOne, int printlevel, char flag);

  virtual std::ostream& printActiveDynamicConstraints(std::ostream& os = std::cout) const;
  virtual std::ostream& printDetailedPrimalSol(std::ostream& os = std::cout) const;
  virtual std::ostream& printSol(std::ostream& os = std::cout) const;
  virtual std::ostream& printSolVal(std::ostream& os = std::cout) const;
  virtual std::ostream& printDualSol(std::ostream& os = std::cout, bool compact = false) const;
  virtual std::ostream& printProb(std::ostream& os = std::cout) const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual std::ostream & printForm(std::ostream& os = std::cout) const;
  virtual std::ostream& printPartialSolution(std::ostream& os = std::cout) const;
  virtual int solveProb(int & maxLevelOfRestriction, const char & flag = ' ', const bool & ifPrint = false);
  virtual bool solveProbLP(const char & flag = ' ', const bool & ifPrint = false);
  virtual bool solveProbMIP(const char & flag = ' ', const bool & ifPrint = false);
  virtual void retrieveRedCosts();

  /**
   * Return true if a feasible solution  to the problem was found, false otherwise
   */
  const VarIndexManager & probVarSet() const;
  VarIndexManager & probVarSet();
  const ConstrIndexManager & probConstrSet() const;
  ConstrIndexManager & probConstrSet();
  virtual const ConstrPtrSet & probNonLinearConstrSet() const;

  virtual void updateInDualSol();

  virtual void updateObjCoeffsInForm(const std::list<Variable *> & varPtrList);
  virtual void updateBoundsInForm(const std::list<Variable *> & varPtrList);
  virtual void updateConstrRhsInForm(const std::list<Constraint *> & constrPtrList);

  void clearPreprocessingLists();

  VarPtrList & preprocessedVarsList()
  {
    return _preprocessedVarsList;
  }

  ConstrPtrList & preprocessedConstrsList()
  {
    return _preprocessedConstrsList;
  }

  const Node * curNodePtr() const
  {
    return _curNodePtr;
  }

  void curNodePtr(const Node * nodePtr)
  {
    _curNodePtr = nodePtr;
  }

  virtual void resetPartialSolution();
  bool solIsInt();
  bool solIsInt(const VarPtr2DoubleMap & curPrimalSol);

  virtual bool primalSolIsFeasible();
  virtual bool primalSolIsFeasible(const VarPtr2DoubleMap & curPrimalSol);

  virtual void  clearRecordedSol();
  virtual Solution * retrieveCurPrimalLpSol(const bool & recordPartialSolOnly = false); // const;
  virtual Solution * recordSolution(Solution * solPtr);
  virtual Solution * incumbentSol() const;
  virtual Solution * extractIncumbent();

  const SolutionMethod & solMode() const;

  /**
   * @param flag unused
   * @param ifPrint indicate whether the routine is expedted 
   * to print details on the procedure to the standard output 
   * @param requiredSolStat provides a set 
   * of admissible solution status
   * @param objVal return value for the objective 
   * value of the best solution found
   * @param primalBound return value for 
   * the best primal bound on the problem solution
   * @param dualBound return value for the best dual bound on the problem solution
   * @param PrimalSol return set of pair <variables, value> involved 
   * in the best primalsolution (if any has been found)
   * @param DualSol return set of of pair <constraint, value>  involved 
   * in the best dual solution (if any has been found)
   *
   * @return // True if a feasible solution to the problem was found, false otherwise
   */

  virtual bool customizedSolver(int & maxLevelOfRestriction,
                                const SolutionStatus & requiredSolStat,
                                Double & objVal,
                                Double & primalBound,
                                Double & dualBound,
                                VarPtrSet & inPrimalSol,
                                ConstrPtrSet & inDualSol,
                                BcSolution & primalSolPtr,
                                BcDualSolution & dualSolPtr);

  virtual bool setupCustomizedSolver(const BcSolverOracleInfo * infoPtr);
  virtual BcSolverOracleInfo * recordSolverOracleInfo();
  void getActiveRyanAndFosterBranchingConstraints(std::list<BcRyanAndFosterBranchConstr>
                                                  & ryanAndFosterBranchConstrList) const;

  /// returns -1 if could not enumerate, otherwise returns the number of enumerated solutions
  /// adds enumerated columns to the solution chain defined by solPtr
  int enumerateAllColumns(Solution * solPtr);

  void getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions);

  /// returns -1 if not enumerated, otherwise returns the number of enumerated solutions
  virtual void reducedCostFixingAndEnumeration(const int & enumerationMode, const Double & threshold);

  virtual void checkEnumeratedSolutions(const std::vector<Solution *> & solPts, std::vector<bool> & solIsEnumerated);
  virtual void getEnumeratedSolutions(const int & maxNumOfSolutions, Solution * solPtr, std::vector<double> & redCosts);
  virtual void getDebugSolution(Solution * solPtr);
  virtual bool setDebugSolution(const std::vector<std::vector<int> > & orderedSolutions, bool vertexBased);

  virtual bool isProperSolution(Solution * solPtr);
  virtual bool solSatisfiesCurrentSpRelaxation(Solution * solPtr);
  virtual void callColGenTerminationCallBack(bool afterRedCostFixing,  int nodeOrder, int nodeDepth, int cutSepRound,
                                             double dualBound, double elapsedTime, bool masterConverged);
  virtual bool improveCurrentSpRelaxation(std::vector<MastColumn *> & colsInMasterSolution,
                                          const bool & masterConverged);
  virtual bool lightenCurrentSpRelaxation(const int & masterConverged, const int & callMode);
  virtual bool drawPrimalSolutionToDotFile(std::vector<MastColumn *> & colsInMasterSolution,
                                           const std::string & filename) const;

  virtual int getMessageIdToCutGeneration() const;
  virtual int getNumberOfEnumeratedSolutions() const;

  void printDynamicVarConstrStats(std::ostream & os, const bool & completePrint = false);

  void updateInNonZeroRedCostVarsSet(Variable * varPtr);
  void removeVarsNotInProblemFromNonZeroRedCostVars();
  void removeVarsNotInProblemFromPrimalSolution();

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr();
} ;

inline std::ostream& operator<<(std::ostream& os, const Problem & that)
{
  return that.print(os);
}

class MipProblem : public Problem
{
protected:
  SolutionStatus _mipRequiredStatus;
  SolutionStatus _mipProbStatus;
  bool _MIPpreprocessorOn;
  bool _MIPprobingOn;
  bool _MIPautomaticCuttingPlanesOn;
  char _MIPsolverSelect;

public:

  MipProblem(const int & ref,
             const double & rightHAndSideZeroTol,
             const double & reducedCostTolerance,
             const BcObjStatus::MinMaxIntFloat & minmaxStatus = BcObjStatus::minInt,
             const SolutionMethod & solMode = SolutionMethod::mipSolver,
             const std::string & name = "sol",
             const SolutionStatus & LPrequiredSolStat = SolutionStatus::Optimum,
             const bool & LPpreprocessorOn = true,
             const bool & LPprobingOn = true,
             const SolutionStatus & MIPrequiredSolStat = SolutionStatus::Optimum,
             const bool & MIPpreprocessorOn = true,
             const bool & MIPprobingOn = true,
             const bool & MIPautomaticCuttingPlanesOn = true,
             const char & MIPsolverSelect = 'd');

  virtual ~MipProblem()
  {
  }

  virtual bool primalSolIsFeasible();
  virtual bool primalSolIsFeasible(const VarPtr2DoubleMap & curPrimalSol);
  const bool & preprocessorOn()  const
  {
    return  _MIPpreprocessorOn;
  }

  const bool & probingOn()  const
  {
    return  _MIPprobingOn;
  }

  const bool & automaticCuttingPlanesOn()  const
  {
    return  _MIPautomaticCuttingPlanesOn;
  }

  virtual const SolutionStatus & probStatus() const
  {
    return _mipProbStatus;
  }

  virtual int solveProb(int & maxLevelOfRestriction, const char & flag = ' ', const bool & ifPrint = false);
  virtual bool solveProbMIP(const char & flag = ' ', const bool & ifPrint = false);

  virtual void resetSolution(const char & flag = 'p');

  virtual void setProbStatus(const int & stat);
  virtual void setProbStatus(const SolutionStatus & stat);
  virtual void setStatusAfterSol();
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  // added by Boris (modified by Guillaume)
  virtual void setMIPRequiredStatus(const SolutionStatus &newStatus)
  {
	  _mipRequiredStatus = newStatus;
  }
  // added by Boris
  void getMIPRequiredStatus(SolutionStatus &currentStatus)
  {
	  currentStatus = _mipRequiredStatus;
  }
};

inline std::ostream& operator<<(std::ostream& os, const MipProblem & that)
{
  return that.print(os);
}

#endif // ProblemClasses_h
