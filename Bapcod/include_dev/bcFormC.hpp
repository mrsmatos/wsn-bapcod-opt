/**
 *
 * This file bcFormC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef FormulationClasses_h
#define FormulationClasses_h
#include "bcUsefulHeadFil.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcErrorC.hpp"
#include "bcMathProgSolverInterfaceC.hpp" //Activated in order to define a typrdef (hsen)
#include "bcVarConstrC.hpp"

//class MathProgSolverInterface; //Commented out because the header included above
struct LpBasisRecord;

/**
   @brief Formulation represents a linear program with mixed varaibles
   @details it serves as a common interface to several mip solvers such as Cplex, Gurobi or Clp or Express-mp
*/

class BaseFormulation
{
  Problem * _problemPtr;
 protected:
  /// 1 if minimization, -1 if maximisation;
  int _minmax;
  SolutionStatus _status;
 private:
  BaseFormulation();
 public:
  BaseFormulation(Problem * problemPtr);
  Problem * problemPtr() const {return _problemPtr;}

  /// Problem interface
  virtual void setConstr2Form(Constraint * constrPtr, const bool & fillConstrMatrix = true) = 0;

  virtual void unsetConstr2Form(Constraint * constrPtr) = 0;
  virtual void setVar2Form(Variable * varPtr) = 0;
  virtual void unsetVar2Form(Variable * varPtr) = 0;
  virtual void addConstr2Formulation() = 0;
  virtual void addVar2Formulation() = 0;
  virtual void delConstrFromFormulation() = 0;
  virtual void delVarFromFormulation() = 0;
  virtual void resetRhs(Constraint * constrPtr) = 0;
  virtual void resetObjCoef(Variable * varPtr) = 0;
  virtual void resetBounds(Variable * varPtr) = 0;
  virtual void updateConstrRhsInFormulation() = 0;
  virtual void updateObjectiveInFormulation() = 0;
  virtual void updateBoundsInFormulation() = 0;
  virtual void buildFormulation() = 0;
  virtual void printForm(std::ostream& os = std::cout) = 0;
  virtual void clearColFormulationDataStruct() = 0;
  virtual void clearRowFormulationDataStruct() = 0;
  virtual const SolutionStatus & status() const
  {
    return _status;
  }

  virtual std::ostream& printMatrix(std::ostream& os = std::cout) const = 0;

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;
  const ControlParameters& param() const;
  ControlParameters& param();

  virtual void retrieveBasis(LpBasisRecord & basis, const bool markInVars, const bool markInConstrs){}
  virtual void loadBasis(const LpBasisRecord & basis){}
};

class LPform: public BaseFormulation
{
protected:
  MathProgSolverInterface * _interfacePtr;

  /// Counter setup when filling data structure
  int _probRowCnt;

  /// Counter setup when filling data structure
  int _probColCnt;
  ProbCoefContainer _objectRow;
  ProbRowMatrixContainer _rowMatrix;
  ProbColMatrixContainer _colMatrix;
  ProbRhsContainer _rhsv;
  ProbBoundContainer _bounds;
  std::set<ProbType> _varTypes;
  std::set<ProbIntC> _varBrDirective;
  std::set<ProbSetCoef> _soSets;

  /// For delCols
  std::set<int> _indexSetOfCol2Delete;

  /// For delRows
  std::set<int> _indexSetOfRow2Delete;
  std::map < int , std::string > _mapSeqnb2Cname;
  std::map < int , std::string > _mapSeqnb2Rname;

  /**
   * Information needed to retrieve problem solution from Form Class
   */
#define MAPCONSTR 1
#if MAPCONSTR
  std::map < int , Constraint * > _mapConstrSeqNb2ConstrPtr;
  std::map < int , Variable * > _mapVarSeqNb2VarPtr;
#else
  std::deque<Constraint *> _vectConstrPtr;
  std::deque<Variable *> _vectVarPtr;
#endif

  Double _objScalFact;

 protected:
  virtual void checkFormulation();
  void retrieveSol(const char & flag,
		           const bool & ifPrint,
                   VarPtrSet  & inPrimalSol,
                   ConstrPtrSet  & inDualSol);

  virtual void scaleObjectif();
 public:
  LPform(Problem * probPtr, const bool & defineInterface = true);
  
  MathProgSolverInterface * interfacePtr()
  {
    return _interfacePtr;
  }

  virtual ~LPform();
  virtual void chgObjCoef(const ProbCoefContainer & newObjCoef);
  void getSol();
  void retrieveRedCosts(const bool & ifPrint, VarPtrSet & posRedCostVars); /// fast reduced costs preparation (Ruslan)
  void retrievePrimalSol(const std::map<int, Double> & primalSolVect, VarPtrSet  & inPrimalSol) const;
  void fillDataStruct(Constraint * constrPtr, int & sense, double & rhs, std::vector<int> & cutind,
                      std::vector<double> & cutval) const;

  virtual bool solve(const double & BarrierConvergenceTolerance,
                     const double & rightHAndSideZeroTol,
                     const double & reducedCostTolerance,
                     const char & flag, const bool & ifPrint,
                     const  SolutionStatus  & requiredStatus,
                     Double & objVal, Double & primalBound, Double & dualBound,
                     VarPtrSet  & inPrimalSol,
                     ConstrPtrSet  & inDualSol,
                     const bool & preprocessorOn = true,
                     const bool & probingOn = true,
                     const bool & MIPautomaticCuttingPlanesOn = false,
		             const char & solverSelection = 'd');

  virtual void setBounds(Double & objVal, Double & primalBound, Double & dualBound);

  /// Problem Interface
  virtual void setConstr2Form(Constraint * constrPtr, const bool & fillConstrMatrix = true);
  virtual void unsetConstr2Form(Constraint * constrPtr);
  virtual void setVar2Form(Variable * varPtr);
  virtual void unsetVar2Form(Variable * varPtr);
  virtual void fillDataStruct(Constraint * constrPtr, const bool & fillConstrMatrix = true);
  virtual void fillDataStruct(Variable * varPtr);
  virtual void addConstr2Formulation();
  virtual void addVar2Formulation();
  virtual void delConstrFromFormulation();
  virtual void delVarFromFormulation();
  virtual void resetRhs(Constraint * constrPtr);
  virtual void resetObjCoef(Variable * varPtr);
  virtual void replaceObjCoef(Variable * varPtr);
  virtual void resetBounds(Variable * varPtr);
  virtual void updateConstrRhsInFormulation();
  virtual void updateObjectiveInFormulation();
  virtual void updateBoundsInFormulation();
  virtual void buildFormulation();
  virtual void printForm(std::ostream& os = std::cout);
  virtual void clearColFormulationDataStruct();
  virtual void clearRowFormulationDataStruct();
  virtual std::ostream& printMatrix(std::ostream& os = std::cout) const;

  virtual void retrieveBasis(LpBasisRecord & basis, const bool markInVars,
                             const bool markInConstrs);
  virtual void loadBasis(const LpBasisRecord & basis);
};

class MIPform: public LPform
{
 protected:
  SolutionStatus _mipStatus;

  virtual void checkFormulation();
 private:
  MIPform();
 protected:
  virtual std::ostream& printMatrix(std::ostream& os = std::cout) const;
 public:
  MIPform(Problem * probPtr);

  virtual ~MIPform(){}

  virtual bool solve(const double & BarrierConvergenceTolerance,
                     const double & rightHAndSideZeroTol,
                     const double & reducedCostTolerance,
                     const char & flag,
                     const bool & ifPrint,
                     const  SolutionStatus & requiredStatus,
                     Double & objVal,
		             Double & primalBound,
                     Double & dualBound,
                     VarPtrSet & inPrimalSol,
                     ConstrPtrSet & inDualSol,
                     const bool & MIPpreprocessorOn = true,
                     const bool & MIPprobingOn = true,
                     const bool & MIPautomaticCuttingPlanesOn = true,
		             const char & solverSelection = 'd');

  virtual const SolutionStatus & status() const
  {
    return _mipStatus;
  }

  virtual void setBounds(Double & objVal, Double & primalBound, Double & dualBound);

  virtual void fillDataStruct(Variable * varPtr);

  //added by Issam
  void resetMIPpartOfFormulation(const Double & cutOffValue, const bool exactSolution = true);

  void setLazyConstraintsCallback(MasterConf * masterConfPtr);

  void resetAfterMIP();

  void setTimeLimit(const double timeLimit);

  virtual void addVar2Formulation();

  /// Problem Interface
  virtual void buildFormulation();
  virtual void clearColFormulationDataStruct();
};

#endif // FormulationClasses_h
