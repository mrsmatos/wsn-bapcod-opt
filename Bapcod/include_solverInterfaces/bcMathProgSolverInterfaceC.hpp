/**
 *
 * This file bcMathProgSolverInterfaceC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMathProgSolverInterfaceC_H_
#define BCMathProgSolverInterfaceC_H_

#include "bcDoubleC.hpp"

class ProbBound;
class ProbCoef;
class ProbIntC;
class BapcodInit;
class ControlParameters;
class ProbCoefColSmallerThan;
class ProbType;
class ProbCoefRowSmallerThan;
class ProbSetCoef;
class SolutionStatus;

//Added by (hsen): Define one and only one of the following: { PROBCONTAINERSET, PROBCONTAINERVECTOR, PROBCONTAINERDEQUE }
//#define PROBCONTAINERDEQUE
#define PROBCONTAINERSET
/* PROBCONTAINER
 * SET was the implementation in place at the time of the change and had considerable amount of overhead (both in DW and BD)
 * VECTOR is the fastest (i.e., least overhead) but it may lead to heap fragmentation in "huge" applications
 * DEQUE is the best of both worlds (i.e., allocation in blocks of equal sizes, which are chained together)
 * However, as long as it is possible, one should use VECTOR. If there is an application which has an unsatisfactory
 * performance due to fragmanted heap memory (i.e., slowed calls to new) then DEQUE should be tried before SET.
 */

#ifdef PROBCONTAINERSET
typedef std::set<ProbCoef>  ProbCoefContainer;
typedef std::set<ProbBound> ProbRhsContainer;
typedef std::set<ProbBound> ProbBoundContainer;
#endif

#ifdef PROBCONTAINERVECTOR
typedef std::vector<ProbCoef>  ProbCoefContainer;
typedef std::vector<ProbBound> ProbRhsContainer;
typedef std::vector<ProbBound> ProbBoundContainer;
#endif

#ifdef PROBCONTAINERDEQUE
typedef std::deque<ProbCoef>  ProbCoefContainer;
typedef std::deque<ProbBound> ProbRhsContainer;
typedef std::deque<ProbBound> ProbBoundContainer;
#endif


//Added by (hsen): Define one and only one of the following: { PROBMATRIXCONTAINERSET, PROBMATRIXCONTAINERVECTOR, PROBMATRIXCONTAINERDEQUE }
//#define PROBMATRIXCONTAINERVECTOR
#define PROBMATRIXCONTAINERSET
/* PROBMATRIXCONTAINER
 * SET was the implementation in place at the time of the change and had considerable amount of overhead (both in DW and BD)
 * However, it was sorted with ProbCoefColSmallerThan for ProbColMatrixContainer and
 * with ProbCoefRowSmallerThan for ProbRowMatrixContainer.
 * In DEQUE and VECTOR implementations, the actual containers are sorted whenever they need to be.
 * SET implementation is kept in case the other container types lead to a bug. For example, LpCplexInterface::addCols leads to
 * differrent dual solutions (from Cplex), even though the printed models are completely the same.
 * For the advantages of DEQUE and VECTOR see the above explanation for PROBCONTAINER
 */

#ifdef PROBMATRIXCONTAINERSET
typedef std::set<ProbCoef, ProbCoefColSmallerThan> ProbColMatrixContainer; //LpCplexInterface::loadFormulation assumes colMatrix is sorted! LpCplexInterface::addCols behaves more predictable when colMatrix is sorted!
typedef std::set<ProbCoef, ProbCoefRowSmallerThan> ProbRowMatrixContainer; //
#endif

#ifdef PROBMATRIXCONTAINERVECTOR
typedef std::vector<ProbCoef> ProbColMatrixContainer; //LpCplexInterface::loadFormulation assumes colMatrix is sorted! LpCplexInterface::addCols behaves more predictable when colMatrix is sorted!
typedef std::vector<ProbCoef> ProbRowMatrixContainer; //ProbRowMatrixContainer is not sorted!
#endif

#ifdef PROBMATRIXCONTAINERDEQUE
typedef std::deque<ProbCoef> ProbColMatrixContainer; //LpCplexInterface::loadFormulation assumes colMatrix is sorted! LpCplexInterface::addCols behaves more predictable when colMatrix is sorted!
typedef std::deque<ProbCoef> ProbRowMatrixContainer; //ProbRowMatrixContainer is not sorted!
#endif

class MasterConf;

class MathProgSolverInterface
{
    friend struct LazyCallbackStruct;
public:
  enum BasisStatus    // by Artur
  {
    UndefinedStatus = -1,
    AtLowerVarConstr = 0,
    BasicVarConstr = 1,
    AtUpperVarConstr = 2,
    FreeSuperVarConstr = 3
  };

protected:
  BapcodInit* _bapcodInit;

  int _ref;
  bool _formCurrentlyLoaded;
  bool _formCurrentlySaved;
  bool _pureLP;

  /// Counter representing current value
  long _ncol;

  /// Counter representing current value
  long _nrow;

 public:
  
   
  // THESE METHODS WHERE INITIALLY IN MathProgSolverInit 
   
  virtual void setPresolve(const bool & flag = true) = 0;
  virtual void setTimeLimit(const double & seconds = 3600) = 0;
  virtual void setMultiThread(const int & maxNbThread = 0) = 0; /// automatic use of max nb of threads
  virtual void setRelativeMipGapLimit(const double & relativeGap = 0.0001) = 0;
  virtual void setSearchPriority(const int & flag = 0) = 0;
  virtual void setWorkingMemorySpace(const double & sizeInMB = 4000) = 0;
  virtual void setLPoptimalityTolerance(const double& tolerance = 1e-6) = 0;
  virtual void setLPfeasibilityTolerance(const double& tolerance = 1e-6) =0; 
  virtual void setScreenOutput(const bool & flag) {}
  virtual void setMaxBBNodes(const int & maxBBNodes) {}  
  
  // THESE METHODS WHERE INITIALLY IN MathProgSolverInit 
   
  MathProgSolverInterface(BapcodInit* bapcodInit);

  MathProgSolverInterface(BapcodInit* bapcodInit, const int & refconst, const std::string & name);
  virtual ~MathProgSolverInterface();
  virtual const int & ref() const;

  virtual void loadFormulation(const std::string & name,
                               const int & minmaxStatus,
                               const int &  ncol,
                               const int &  nrow,
                               const std::map < int, std::string > & mapSeqnb2Cname,
                               const std::map < int, std::string > & mapSeqnb2Rname,
                               const ProbCoefContainer & objectRow,
                               ProbColMatrixContainer & colMatrix,
                               const ProbRhsContainer & rhsv,
                               const ProbBoundContainer & bounds,
                               const std::set<ProbType> & types,
                               const std::set<ProbIntC> & directs,
                               const std::set<ProbSetCoef> & sets) = 0;

  virtual void unLoadFormulation() = 0;

  virtual void addCols(const ProbCoefContainer & objectiveRow,
                       ProbColMatrixContainer & colMatrix,
                       const ProbBoundContainer & bounds,
                       const std::map < int, std::string > & mapSeqnb2Cname) = 0;

  virtual void delCols(const std::set<int> & indexSetOfCol2Delete) = 0;

  virtual void addRows(const ProbRhsContainer & rhsv,
                       const ProbRowMatrixContainer & _rowMatrix,
                       const std::map < int, std::string > & mapSeqnb2Rname) = 0;

  virtual void delRows(const std::set<int> & indexSetOfRow2Delete) = 0;
  virtual void addDirectives(const std::set<ProbIntC> & varBrDirection) = 0;
  virtual void getObjCoef(ProbCoef & pc) = 0;
  virtual void chgObjCoef(const ProbCoef & pc) = 0;
  virtual void chgBds(const ProbBoundContainer & newBounds) = 0;
  virtual void chgMatCoef(const ProbCoef & pc) = 0;
  virtual void chgRhs(const ProbBound & pb) = 0;
  virtual void chgRhs(const ProbRhsContainer & newRhs) = 0;
  virtual void chgColType(const std::set<ProbType> & newTypes) = 0;
  virtual void getObjVal(Double & objV) = 0;
  virtual void getObjValLp(Double & objV) = 0;
  virtual void getDualBound(Double & val) = 0;
  virtual void getPrimalBound(Double & val) = 0;
  virtual bool getOptimStatus(SolutionStatus & lpStatus, SolutionStatus & mipStatus) = 0;
  virtual void getSol(std::map<int, Double> & primSol, const bool & ifPrint = false) = 0;
  virtual void getSol(std::map<int, Double> & primSol, std::map<int, Double> & dualSol,
                      const int & minmaxStatus, const bool & ifPrint = false) = 0;
          
  virtual void getSolPool(std::vector< std::map<int, Double> > & primSolPool, const bool & ifPrint = false) = 0;
  
  virtual void getReducedCost(std::map<int, Double> & redCost, const bool & ifPrint = false) = 0;

  //added by Ruslan
  virtual void setLazyConstraintsCallback(MasterConf * mastConfPtr);
  virtual void removeLazyConstraintsCallback();
  virtual void resetUpperCutOffValue() = 0;
  virtual void setUpperCutOffValue(const Double & cutOff) = 0;
  virtual void setSolveFromScratch(const bool exactSolution = true) = 0;
  virtual void resetSolveFromScratch() = 0;

  // by Artur: TODO: implement for solvers other than CPLEX
  virtual void getBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus) = 0;
  virtual void setBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus) = 0;

  // by Aurélien
  // Makes the solver stops once the primal bound computed by the solver becomes equal to the given value)
  virtual void setAbortValue(const Double& abortValue) = 0;

  //  bool preProcessorChangedForm(int * ptrProbRowCnt, const bool & checkVarOnly = false);
  virtual void LPwrite(const int & minmaxStatus, std::ostream& os = std::cout) = 0;
  virtual void MPSwrite() = 0;
  virtual void printRow(int irow,
			            const int & readNcol,
                        int  *readMclind,
			            double *readDmatval,
			            char  *readCnames) = 0;

  virtual void optimise(const int & minmaxStatus,
            const double & BarrierConvergenceTolerance,
            const double & rightHAndSideZeroTol,
            const double & reducedCostTolerance,
			const bool & preprocessorOn = true,
			const bool & probingOn = true,
			const bool & automaticCuttingPlanesOn = false,
			const char & solverSelection = 'p') = 0;

  virtual void optimiseLp(const int & minmaxStatus,
            const double & BarrierConvergenceTolerance,
            const double & rightHAndSideZeroTol,
            const double & reducedCostTolerance,
  			const bool & preprocessorOn = true,
  			const bool & probingOn = true,
  			const bool & automaticCuttingPlanesOn = false,
  			const char & solverSelection = 'p') = 0;

  virtual void reset() = 0;
  virtual void makeSpaceForLoadingForm() = 0;
  virtual void saveCopyOfCurForm() = 0;
  virtual void load() = 0;
  virtual void unload(const bool & deleteMat = false) = 0;

  const long & ncol() const;
  const long & nrow() const;

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;
  const ControlParameters& param() const;
  ControlParameters& param();
};

#endif /* BCMathProgSolverInterfaceC_H_ */
