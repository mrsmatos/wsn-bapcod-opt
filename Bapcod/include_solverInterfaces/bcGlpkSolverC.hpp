/**
 *
 * This file bcGlpkSolverC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef bcGlpkInterfaceC_h
#define bcGlpkInterfaceC_h

#if _GLPK_FOUND

#include "bcMathProgSolverInterfaceC.hpp"

#include <glpk.h>

class MipGlpkInit : public MathProgSolverInit
{
 public:
  /// Contructor
  MipGlpkInit(BapcodInit* bapcodInit);
  ~MipGlpkInit();

  void setPresolve(const bool & flag = true);
  void setTimeLimit(const double & seconds = 3600);
  void setMultiThread(const int & maxNbThread = 0); /// automatic use of max nb of threads
  void setRelativeMipGapLimit(const double & relativeGap = 0.0001);
  void setSearchPriority(const int & flag = 0);
  void setWorkingMemorySpace(const double & sizeInMB = 4000);
};

class LpGlpkInterface : public MathProgSolverInterface
{
 public:
  static int GLPKLpCounter;
  static int GLPKMipCounter;
  glp_prob * glpkProbPtr;
  glp_smcp glpkParamPtr;

 public:
  LpGlpkInterface(BapcodInit* bapcodInit, const int & refconst, const std::string & name);
  virtual ~LpGlpkInterface();
  const int & ref() const
  {
    return _ref;
  }
  virtual void loadFormulation(const std::string & name,
                               const int & minmaxStatus,
                               const int &  ncol,
                               const int &  nrow,
                               const std::map < int, std::string > & mapSeqnb2Cname,
                               const std::map < int, std::string > & mapSeqnb2Rname,
                               const std::set<ProbCoef> & objectRow,
                               const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
                               const std::set<ProbBound> & rhsv,
                               const std::set<ProbBound> & bounds,
                               const std::set<ProbType> & types,
                               const std::set<ProbIntC> & directs,
                               const std::set<ProbSetCoef> & sets);

  virtual void unLoadFormulation();

  virtual void addCols(const std::set<ProbCoef> & objectiveRow,
                       const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
                       const std::set<ProbBound> & bounds,
                       const std::map < int, std::string > & mapSeqnb2Cname);

  virtual void delCols(const std::set<int> & indexSetOfCol2Delete);

  virtual void addRows(const std::set<ProbBound> & rhsv,
                       const std::set<ProbCoef, ProbCoefRowSmallerThan> & _rowMatrix,
                       const std::map < int, std::string > & mapSeqnb2Rname);

  virtual void delRows(const std::set<int> & indexSetOfRow2Delete);
  virtual void addDirectives(const std::set<ProbIntC> & varBrDirection);
  virtual void getObjCoef(ProbCoef & pc);
  virtual void chgObjCoef(const ProbCoef & pc);
  virtual void chgBds(const std::set<ProbBound> & newBounds);
  virtual void chgMatCoef(const ProbCoef & pc);
  virtual void chgRhs(const ProbBound & pb);
  virtual void chgRhs(const std::set<ProbBound> & newRhs);
  virtual void chgColType(const std::set<ProbType> & newTypes);
  virtual void getObjVal(Double & objV);
  inline virtual void getObjValLp(Double & objV) { LpGlpkInterface::getObjVal(objV); }
  virtual void getDualBound(Double & val);
  virtual void getPrimalBound(Double & val);
  virtual bool getOptimStatus(SolutionStatus & lpStatus, SolutionStatus & mipStatus);
  virtual void getSol(std::map<int, Double> & primSol, const bool & ifPrint = false);
  virtual void getSol(std::map<int, Double> & primSol,
		      std::map<int, Double> & dualSol,
		      const int & minmaxStatus,
		      const bool & ifPrint = false);
  virtual void getReducedCost(std::map<int, Double> & redCost, const bool & ifPrint = false);

  //  bool preProcessorChangedForm(int * ptrProbRowCnt, const bool & checkVarOnly = false);
  virtual void LPwrite(const int & minmaxStatus, std::ostream& os = std::cout);
  virtual void MPSwrite();
  virtual void printRow(int irow,
			const int & readNcol,
			int  *readMclind,
			double *readDmatval,
			char  *readCnames);

  virtual void optimise(const int & minmaxStatus,
			const bool & preprocessorOn = true,
			const bool & probingOn = true,
			const bool & automaticCuttingPlanesOn = false,
			const char & solverSelection = 'p');

  virtual void optimiseLp(const int & minmaxStatus,
			  const bool & preprocessorOn = true,
			  const bool & probingOn = true,
			  const bool & automaticCuttingPlanesOn = false,
			  const char & solverSelection = 'p')
  {
    LpGlpkInterface::optimise(minmaxStatus, preprocessorOn, probingOn, automaticCuttingPlanesOn, solverSelection);
  }
  virtual void reset();
  virtual void makeSpaceForLoadingForm();
  virtual void saveCopyOfCurForm();
  virtual void load();
  virtual void unload(const bool & deleteMat = false);
};

class MipGlpkInterface: public LpGlpkInterface
{
 public:
  glp_iocp glpkParamMipPtr;

 public:
  MipGlpkInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name);
  ~MipGlpkInterface() {}

  void loadFormulation(const std::string & name,
		       const int & minmaxStatus,
		       const int &  ncol,
		       const int &  nrow,
		       const std::map < int, std::string > & mapSeqnb2Cname,
		       const std::map < int, std::string > & mapSeqnb2Rname,
		       const std::set<ProbCoef> & objectRow,
		       const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
		       const std::set<ProbBound> & rhsv,
		       const std::set<ProbBound> & bounds,
		       const std::set<ProbType> & types,
		       const std::set<ProbIntC> & directs,
		       const std::set<ProbSetCoef> & sets);

  void addDirectives(const std::set<ProbIntC> & varBrDirection);
  void chgColType(const std::set<ProbType> & newTypes);
  void getObjVal(Double & objV);
  void getDualBound(Double & val);
  void getPrimalBound(Double & val);
  bool getOptimStatus(SolutionStatus & lpStatus, SolutionStatus & mipStatus);
  void optimise(const int & minmaxStatus,
		const bool & preprocessorOn = true,
                const bool & probingOn = true,
                const bool & automaticCuttingPlanesOn = true,
                const char & solverSelection = 'p');

  void getSol(std::map<int, Double> & primSol, const bool & ifPrint = false);
  void reset();
};


#endif // _Glpk_FOUND
#endif // bcGlpkInterfaceC_h
