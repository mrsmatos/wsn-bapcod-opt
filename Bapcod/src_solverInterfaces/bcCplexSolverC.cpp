/**
 *
 * This file bcCplexSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _CPLEX_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcCplexSolverC.hpp"
#include "bcPrintC.hpp"
#include "bcModelC.hpp"
#include "bcMasterConfC.hpp"
#include "bcFormC.hpp"

//#include "cplex.h"
#include <ilcplex/cplex.h>

//#define SAVE_ALL_MIPS

/// Fixed length of names given to sub-solvers
const int wordSize = 16;
const char endOfWord = '\0';
const std::string filling = "________________";

/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

CPXENVptr cplexEnv = NULL;
FILE* log_cplex = NULL;

using namespace std;


static int CPXPUBLIC lazycallback(CPXCENVptr env, void *cbdata, int wherefrom,
                                  void *cbhandle, int *useraction_p)
{
    *useraction_p = CPX_CALLBACK_DEFAULT;

    LazyCallbackStruct * lcStructPtr = (LazyCallbackStruct *)cbhandle;
    MasterConf * masterPtr = lcStructPtr->mastConfPtr;


#ifdef _MSC_VER
    std::vector<double> x(lcStructPtr->numCols);
#else
    double x[lcStructPtr->numCols];
#endif
    int status = CPXgetcallbacknodex(env, cbdata, wherefrom, &x[0], 0,lcStructPtr->numCols - 1);
    if (status != 0)
        return status;

    std::map<int, Double> primalSolVect;
    for (int colId = 0; colId < lcStructPtr->numCols; ++colId)
    {
        if (x[colId] > Double::precision)
            primalSolVect[colId] = x[colId];
    }

    masterPtr->probPtr()->resetSolution();
    VarPtrSet inPrimalSol;
    masterPtr->probPtr()->primalFormulationPtr()->retrievePrimalSol(primalSolVect, inPrimalSol);

    std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> generatedCutConstrSet;
    for (auto genCutConstrPtr : lcStructPtr->mastConfPtr->candidateCutGenericConstr())
        if (genCutConstrPtr->type() == 'C')
        {
            genCutConstrPtr->cutSeparationRoutine(inPrimalSol, generatedCutConstrSet);
        }

    std::vector<int> cutind;
    std::vector<double> cutval;
    for (auto constrPtr : generatedCutConstrSet)
    {
        InstanciatedConstr * instConstrPtr = masterPtr->castAndAddConstraint(constrPtr, false);
        instConstrPtr->toBeUsedInPreprocessing(false);
        /// add the cut to the master problem, its local artificial variables will be also added
        instConstrPtr->addToProb(masterPtr->probPtr());

        int sense;
        double rhs;
        cutind.clear();
        cutval.clear();
        masterPtr->probPtr()->primalFormulationPtr()->fillDataStruct(instConstrPtr, sense, rhs, cutind, cutval);

        int vectorSize = (int)cutind.size();
        status = CPXcutcallbackadd(env, cbdata, wherefrom, vectorSize, rhs, sense, &cutind[0],
                                   &cutval[0], CPX_USECUT_FORCE);
        if (status != 0)
            return status;
    }

    /// we should reset solution
    for (auto varPtr : inPrimalSol)
        varPtr->val(0);
    inPrimalSol.clear();

    return 0;
}

/// New functionalities to add to mip inteface: begin
void LpCplexInterface::setPresolve(const bool & flag)
{
  if (flag)
    bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_PREIND, CPX_ON),
                       "Failure to turn on presolve ");
  else
    bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_PREIND, CPX_OFF),
                       "Failure to turn off presolve ");

  return;
}

void LpCplexInterface::setTimeLimit(const double & seconds)
{

  if (seconds < 0.1)
    bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_TILIM,0.1),
                       "Failure to set time limit ");
  else
    bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_TILIM, seconds),
                       "Failure to set time limit ");

  return;
}

void LpCplexInterface::setMultiThread(const int & maxNbThread)
{
  bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_THREADS, maxNbThread),
                     "Failure to set number of threads ");

  return;
}


void LpCplexInterface::setSearchPriority(const int & flag)
{
    switch (flag)
    {
        case 0 :
        {
            bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BALANCED),
                               "Failure to set mipemphasis indicator ");
            break;
        }
        case 1 :
        {
            bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_FEASIBILITY),
                               "Failure to set mipemphasis indicator ");
            break;
        }
        case 2 :
        {
            bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY),
                               "Failure to set mipemphasis indicator ");
            break;
        }
        case 3 :
        {
            bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BESTBOUND),
                               "Failure to set mipemphasis indicator ");
            break;
        }
        case 4 :
        {
            bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_HIDDENFEAS),
                               "Failure to set mipemphasis indicator ");
            break;
        }
    }

    return;
}

void LpCplexInterface::setMaxBBNodes(const int & maxBBNodes)
{
  bapcodInit().check(CPXsetlongparam (cplexEnv, CPX_PARAM_NODELIM, maxBBNodes),
                     "Failure to turn on screen indicator");
}


void LpCplexInterface::setScreenOutput(const bool & flag)
{
  if (flag)
    {
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_SCRIND, CPX_ON),
                         "Failure to turn on screen indicator");
    }
  else
    {
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_SCRIND, CPX_OFF),
                         "Failure to turn off screen indicator");
    }
}

void LpCplexInterface::setRelativeMipGapLimit(const double & relativeGap)
{
  bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_EPGAP, relativeGap),
                     "Failure to set relativeGap ");

  return;
}

void LpCplexInterface::setWorkingMemorySpace(const double & sizeInMB)
{
  bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_WORKMEM, sizeInMB),
                     "Failure to reset work memory ");

  return;
}

void LpCplexInterface::setLPoptimalityTolerance(const double& tolerance)
{
  bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_EPOPT, tolerance),
                     "Failure to set LP optimality tolerance ");
  bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_EPGAP, tolerance),
                     "Failure to set MIP optimality tolerance ");
}

void LpCplexInterface::setLPfeasibilityTolerance(const double& tolerance)
{
  bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_EPRHS, tolerance),
                     "Failure to set LP feasibility tolerance ");
}

LpCplexInterface::LpCplexInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  MathProgSolverInterface(bapcodInit, ref, name),
  cplexProbPtr(NULL),
  cplexPrecision(10e-8) /// used to be parameter HIGHPRECISION
{
  int status(0);
  cplexEnv = CPXopenCPLEX(&status);
  if ( cplexEnv == NULL || status != 0 )
    {
      std::cout << "CPXopenCPLEX status = " << status <<endl;
      std::cout << "Could not open CPLEX environment.\n";
      char errmsg[1024];
      CPXgeterrorstring (cplexEnv, status, errmsg);
      std::cout << errmsg << std::endl;
      exit (1);
    }

  if (printL(1))
    {
      std::cout << "CPXopenCPLEXdevelop status = " << status <<endl;
      std::cout << "CPLEX version is " << CPXversion(cplexEnv) <<endl;
      bapcodInit->check(CPXsetintparam (cplexEnv, CPX_PARAM_SCRIND, CPX_ON),
                        "Failure to turn on screen indicator");
    }
  else
    {
      bapcodInit->check(CPXsetintparam (cplexEnv, CPX_PARAM_SCRIND, CPX_OFF),
                        "Failure to turn off screen indicator");
    }
  
  char *probname = new char[name.size() + 1];
  snprintf(probname, 100, "%s", name.c_str());
  probname[name.size()] = endOfWord;
  cplexProbPtr = CPXcreateprob(cplexEnv, &status, probname);
  bapcodInit->check(cplexProbPtr == NULL,"CPXcreateprob: Failed to create LP");

  bapcodInit->check(CPXchgprobtype(cplexEnv, cplexProbPtr, CPXPROB_LP),
                    " CPXchgprobtype: Failed to change typ\e to LP.");

  delete [] probname; probname = NULL;
  
  return;
}

LpCplexInterface::~LpCplexInterface()
{
    bapcodInit().check(CPXfreeprob(cplexEnv, &cplexProbPtr),"CPXfreeprob failed",
                       ProgStatus::terminate);
  
    bapcodInit().check(CPXcloseCPLEX(&cplexEnv),"could not close cplex and free memory",
                       ProgStatus::terminate);

  return;
}


void LpCplexInterface::loadFormulation(const std::string & name,
				                       const int & minmaxStatus,
				                       const int & ncol,
				                       const int & nrow,
				                       const std::map<int, std::string> & mapSeqnb2Cname,
				                       const std::map<int, std::string> & mapSeqnb2Rname,
				                       const ProbCoefContainer & objectRow,
				                       ProbColMatrixContainer & colMatrix,
				                       const ProbRhsContainer & rhsv,
				                       const ProbBoundContainer & bounds,
				                       const std::set<ProbType> & types,
				                       const std::set<ProbIntC> & directs,
	                			       const std::set<ProbSetCoef> & sets)
{
  if (printL(3))
    std::cout << "LpCplexInterface::loadFormulation() ncol = " << ncol << " , nrow = " << nrow << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  /// Assumes sorted colMatrix
#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
  std::sort(colMatrix.begin(), colMatrix.end());
#endif
  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);

  for (ProbCoefContainer::const_iterator oPtr = objectRow.begin(); oPtr != objectRow.end(); oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  double *rhs = new double[nrow];
  fill_n(rhs, nrow, 0);
  char *rsense = new char[nrow];
  fill_n(rsense, nrow, ' ');
  for (ProbRhsContainer::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
    {
      rhs[rvPtr->ref] = zero(rvPtr->bound);
      rsense[rvPtr->ref] = rvPtr->sense;
    }

  int *matbeg = new int[ncol + 1];
  fill_n(matbeg, ncol + 1, 0);
  int *matcnt = new int[ncol + 1];
  fill_n(matcnt, ncol + 1, 0);
  int nnze = colMatrix.size();
  int *matind = new int[nnze];
  fill_n(matind, nnze, -1);
  double *matval = new double[nnze];
  fill_n(matval, nnze, 0);
  int cntNbNz(0);
  ProbColMatrixContainer::const_iterator mPtr = colMatrix.begin();
  int curCol(0);
  for (;curCol<ncol;curCol++)
  {
    matbeg[curCol] = cntNbNz;
    matcnt[curCol] = 0;

    while ((mPtr != colMatrix.end()) && (mPtr->colRef==curCol))
    {
      matind[cntNbNz] = mPtr->rowRef;
      matval[cntNbNz] = zero(mPtr->coef);
      cntNbNz++;
      matcnt[curCol]++;
      mPtr++;
    }
  }
  /// For dummy  element ncol + 1
  matbeg[curCol] = cntNbNz;
  matcnt[curCol] = 0;


  double *dlb = new double[ncol];
  fill_n(dlb, ncol, 0);
  double *dub = new double[ncol];
  fill_n(dub, ncol, BapcodInfinity);
  for (ProbBoundContainer::const_iterator bPtr = bounds.begin(); bPtr != bounds.end(); bPtr++)
    {
      if (bPtr->sense == 'U')
	    dub[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'L')
	    dlb[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'F')
        {
          dlb[bPtr->ref] = zero(bPtr->bound);
          dub[bPtr->ref] = dlb[bPtr->ref];
        }
    }

  double *range = NULL;
  std::string curName;
  char *colnamestore = new char[ncol * placeSize];
  fill_n(colnamestore, ncol * placeSize, endOfWord);
  char **colnames = new char*[ncol];
  fill_n(colnames, ncol, static_cast<char*>(NULL));
  for(int ic = 0; ic < ncol; ic++)
    {
      if (mapSeqnb2Cname.count(ic))
	    curName = mapSeqnb2Cname.at(ic) + filling;
      else
	    curName = filling;
      strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
      colnamestore[(ic +1) * placeSize - 1] = endOfWord;
      colnames[ic]= &(colnamestore[ic*placeSize]);
    }

  char *rownamestore = new char[nrow * placeSize];
  fill_n(rownamestore, nrow * placeSize, endOfWord);
  char **rownames = new char*[nrow];
  fill_n(rownames, nrow, static_cast<char*>(NULL));
  for(int ir = 0; ir < nrow; ir++)
    {
      if (mapSeqnb2Rname.count(ir))
	    curName = mapSeqnb2Rname.at(ir) + filling;
      else
	    curName = filling;

      strncpy(&(rownamestore[ir * placeSize]), curName.c_str(), wordSize);
      rownamestore[(ir +1) * placeSize - 1] = endOfWord;
      rownames[ir]= &(rownamestore[ir*placeSize]);
    }

  bapcodInit().check(cplexProbPtr == NULL, "LpCplexInterface.load(): problem must have been  created");
  bapcodInit().check(CPXcopylpwnames(cplexEnv, cplexProbPtr, ncol, nrow, minmaxStatus, obj, rhs, rsense, matbeg,
                                     matcnt, matind, matval, dlb, dub, range, colnames, rownames),
                     "CPXcopylpwnames: Failed to copy problem data.");

  delete [] obj;
  obj = NULL;
  delete [] rhs;
  rhs = NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] matbeg;
  matbeg = NULL;
  delete [] matcnt;
  matcnt = NULL;
  delete [] matind;
  matind = NULL;
  delete [] matval;
  matval = NULL;
  delete [] dlb;
  dlb = NULL;
  delete [] dub;
  dub = NULL;
  delete [] range;
  range = NULL;
  delete [] colnames;
  colnames = NULL;
  delete [] rownames;
  rownames = NULL;
  delete [] colnamestore;
  colnamestore = NULL;
  delete [] rownamestore;
  rownamestore = NULL;
  return;
}

void LpCplexInterface::unLoadFormulation()
{
  return;
  _ncol = 0;
  _nrow = 0;
  return;
}

void LpCplexInterface::addCols(const ProbCoefContainer & objectiveRow,
                               ProbColMatrixContainer & colMatrix,
			                   const ProbBoundContainer & bounds,
			                   const std::map<int, std::string> & mapSeqnb2Cname)
{
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  int readNcol(0);
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNcol == _ncol,"LpCplexInterface::addCols: readNcol != _ncol");

  if (printL(7))
    std::cout << "_ncol = " <<  _ncol << "newcol = " << newcol << std::endl;

  int newnzInCols = colMatrix.size();

  double *eobjx = new double[newcol];
  fill_n(eobjx, newcol, 0);

  for (ProbCoefContainer::const_iterator oPtr = objectiveRow.begin(); oPtr != objectiveRow.end(); oPtr++)
    {
      if (printL(7))
	    std::cout << "eobjx mPtr->colRef  = " << oPtr->colRef << ", oPtr->colRef - _ncol = " << oPtr->colRef - _ncol
		          << ", oPtr->coef= " << oPtr->coef << std::endl;
      eobjx[oPtr->colRef - _ncol] = zero(oPtr->coef);
    }

  int *ematbeg = new int[newcol + 1];
  fill_n(ematbeg, newcol + 1, 0);

  int *ematind = new int[newnzInCols];
  fill_n(ematind, newnzInCols, 0);

  double *ematval = new double[newnzInCols];
  fill_n(ematval, newnzInCols, 0);
#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
  std::sort(colMatrix.begin(), colMatrix.end());
#endif
  int cnt(0);
  ProbColMatrixContainer::const_iterator mPtr = colMatrix.begin();
  int newColRef(0);
  for (newColRef = 0; newColRef < newcol; newColRef++)
    {
      ematbeg[newColRef] = cnt;

      if (printL(7))
	    std::cout << " newColRef = " << newColRef << " cnt = " << cnt << " ematbeg[newColRef] = " << ematbeg[newColRef]
		          << std::endl;

      if (mPtr == colMatrix.end())
	    continue;

      if (newColRef < (mPtr->colRef - _ncol))
	    continue;

      if (printL(7))
	    std::cout << "ColMatrix = " << *mPtr << " newColRef = " << newColRef << " mPtr->colRef = " << mPtr->colRef
		          << " _ncol = " <<  _ncol << " newcol = " << newcol
                  << " mPtr->colRef - _ncol = " << mPtr->colRef - _ncol << std::endl;

      while (newColRef == (mPtr->colRef - _ncol))
        {
          ematind[cnt] = mPtr->rowRef;
          ematval[cnt] = zero(mPtr->coef);

         if (printL(7))
	        std::cout << " cnt = " << cnt << " ematind[cnt] = " << ematind[cnt]
                      << " ematval[cnt] = " << ematval[cnt] << std::endl;

          cnt++;
          mPtr++;
          if (mPtr == colMatrix.end())
	        break;
        }
    }
  /// For dummy element newcol + 1
  ematbeg[newColRef] = cnt;

  double * ebdl = new double[newcol];
  fill_n(ebdl, newcol, 0);

  double * ebdu = new double[newcol];
  fill_n(ebdu, newcol, BapcodInfinity);

  if (bounds.size() > 0)
    {
      for (ProbBoundContainer::const_iterator bPtr = bounds.begin(); bPtr != bounds.end(); bPtr++)
        {
          if (bPtr->sense == 'U')
	        ebdu[bPtr->ref - _ncol] = zero(bPtr->bound);
          else if (bPtr->sense == 'L')
	        ebdl[bPtr->ref - _ncol] = zero(bPtr->bound);
          else if (bPtr->sense == 'F')
            {
              ebdl[bPtr->ref - _ncol] = zero(bPtr->bound);
              ebdu[bPtr->ref - _ncol] = zero(bPtr->bound);
            }
        }
    }

  std::string curName;
  char *colnamestore = new char[newcol * placeSize];
  fill_n(colnamestore, newcol * placeSize, endOfWord);

  char **colnames = new char*[newcol];
  fill_n(colnames, newcol, static_cast<char*>(NULL));

  for (int ic = 0; ic < newcol; ic++)
    {
      int colRefNb = ic + _ncol;
      if (mapSeqnb2Cname.count(colRefNb))
	    curName = (mapSeqnb2Cname.at(colRefNb)) + filling;
      else
	    curName = filling;

      strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
      colnamestore[(ic +1) * placeSize - 1] = endOfWord;
      colnames[ic]= &(colnamestore[ic*placeSize]);
    }

  bapcodInit().check(CPXaddcols(cplexEnv, cplexProbPtr, newcol, newnzInCols, eobjx, ematbeg, ematind, ematval,
                                ebdl, ebdu, colnames),
                     "CPXaddcols: could not add new cols");

  delete[] colnamestore;
  colnamestore = NULL;
  delete[] colnames;
  colnames = NULL;
  delete[] eobjx;
  eobjx = NULL;
  delete[] ematbeg;
  ematbeg = NULL;
  delete[] ematind;
  ematind = NULL;
  delete[] ematval;
  ematval = NULL;
  delete[] ebdl;
  ebdl = NULL;
  delete[] ebdu;
  ebdu = NULL;
  _ncol += newcol;

  return;
}

void LpCplexInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  int readNcol;

  int begIndex, endIndex;

  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNcol <= _ncol, "LpCplexInterface::delCols: readNcol > _ncol");
  bapcodInit().require(nbCol2Delete <= readNcol,"LpCplexInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin(); iPtr != indexSetOfCol2Delete.end(); iPtr++)
    {
      dindex[cnt++] = *iPtr;
    }


  begIndex = dindex[cnt-1];
  endIndex = dindex[cnt-1];
  for(int cpt = cnt-2; cpt>=0; cpt--)
    {
      if( dindex[cpt] == (begIndex - 1) )
	    begIndex = dindex[cpt];
      else
        {
	    bapcodInit().check(CPXdelcols(cplexEnv, cplexProbPtr, begIndex, endIndex),
                           "CPXdelcols: could not delete cols");
#ifndef _MSC_VER
	  /**
	   * Only available on Unix, 
	   * due to an ILOG portability issue, 
	   * it doesn't work on Windows.
	   */
	  if ( log_cplex != NULL )
	    fflush(static_cast<FILE *>(log_cplex));
#endif
	  begIndex = dindex[cpt];
	  endIndex = dindex[cpt];
	}
    }
  bapcodInit().check(CPXdelcols(cplexEnv, cplexProbPtr, begIndex, endIndex),
                     "CPXdelcols: could not delete cols");


#ifndef _MSC_VER
  /**
   * Only available on Unix, 
   * due to an ILOG portability issue, 
   * it doesn't work on Windows.
   */
  if ( log_cplex != NULL )
    fflush(static_cast<FILE *>(log_cplex));
#endif

  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}

void LpCplexInterface::addRows( const ProbRhsContainer & rhsv, const ProbRowMatrixContainer & rowMatrix,
                                const std::map<int, std::string> & mapSeqnb2Rname)
{
  int newrows = rhsv.size();
 
  if (newrows <= 0)
    return;
 
  int newnzInNewRows = rowMatrix.size();
  int readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  if (printL(7))
    std::cout << "LpCplexInterface::addRows() newrows = " << newrows
	          << "  newnzInNewRows = " << newnzInNewRows << "  readNrow = " << readNrow << "  _nrow = " << _nrow
	          << std::endl;

 bapcodInit().require(readNrow == _nrow, "LpCplexInterface::addRows: readNrow != _nrow");

  double *erhs = new double[newrows];
  fill_n(erhs, newrows, 0);

  char *esense = new char[newrows];
  fill_n(esense, newrows, ' ');

  for (ProbRhsContainer::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
    {
      if (printL(7))
	    std::cout << "rhsv = " << *rvPtr;
      esense[rvPtr->ref - _nrow] = rvPtr->sense;
      erhs[rvPtr->ref - _nrow] = zero(rvPtr->bound);
    }

  int *ematbeg = new int[newrows + 1];
  fill_n(ematbeg, newrows + 1, 0);

  int *ematind = new int[newnzInNewRows];
  fill_n(ematind, newnzInNewRows, 0);

  double *ematval = new double[newnzInNewRows];
  fill_n(ematval, newnzInNewRows, 0);

  int cnt(0);
  ProbRowMatrixContainer::const_iterator mPtr = rowMatrix.begin();
  int newRowRef(0);
  for (newRowRef = 0; newRowRef < newrows; newRowRef++)
    {
      ematbeg[newRowRef] = cnt;
      if (printL(7))
          std::cout << "rowMatrix = " << *mPtr << " newRowRef = " << newRowRef << " mPtr->rowRef - _nrow = "
                    << mPtr->rowRef - _nrow << " ematbeg[newRowRef] = " << ematbeg[newRowRef] << std::endl;

      if (mPtr == rowMatrix.end())
    	continue;

      if (newRowRef < (mPtr->rowRef - _nrow))
	    continue;

      while (newRowRef == (mPtr->rowRef - _nrow))
        {
          ematind[cnt] = mPtr->colRef;
          ematval[cnt] = zero(mPtr->coef);
          if (printL(7))
	        std::cout << "rowMatrix = " << *mPtr << " cnt = " << cnt << " ematind[cnt] = " << ematind[cnt]
		              << " ematval[cnt] = " << ematval[cnt] << std::endl;
          cnt++;
          mPtr++;
          if (mPtr == rowMatrix.end()) break;
        }
    }

  /// For dummy element newrows + 1
  ematbeg[newRowRef] = cnt;

  std::string curName;

  char * emptyName = new char[newrows * placeSize];
  fill_n(emptyName, newrows * placeSize, endOfWord);

  char** rownames = new char*[newrows];
  fill_n(rownames, newrows, static_cast<char*>(NULL));

  for(int ir = 0; ir < newrows; ir++)
    {
      int newRowRef = ir + _nrow;
      if (mapSeqnb2Rname.count(newRowRef))
	    curName = mapSeqnb2Rname.at(newRowRef) + filling;
      else
	    curName = filling;
      strncpy(&(emptyName[ir * placeSize]), curName.c_str(), wordSize);
      emptyName[(ir +1) * placeSize - 1] = endOfWord;
      rownames[ir]= &(emptyName[ir*placeSize]);
    }


  int newlyAddedColsAlongSideNewRows(0);
  char ** namesOfNewlyAddedCols = NULL;

  bapcodInit().check(CPXaddrows(cplexEnv, cplexProbPtr, newlyAddedColsAlongSideNewRows, newrows, newnzInNewRows,
                                erhs, esense, ematbeg, ematind, ematval, namesOfNewlyAddedCols, rownames),
                     "CPXaddrows: could not add new rows");

  delete[] rownames;
  rownames = NULL;
  delete[] emptyName;
  emptyName = NULL;
  delete[] esense;
  esense = NULL;
  delete[] erhs;
  erhs = NULL;
  delete[] ematbeg;
  ematbeg = NULL;
  delete[] ematind;
  ematind = NULL;
  delete[] ematval;
  ematval = NULL;
  _nrow += newrows;

  return;
}

void LpCplexInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  int readNrow;

  int begIndex, endIndex;

  readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNrow <= _nrow, "LpCplexInterface::delRowss: readNrow > _nrow");
  bapcodInit().require(nbRow2Delete <= readNrow, "LpCplexInterface::delRows: nbRow2Delete > readNrow");

  int *dindex = new int[nbRow2Delete + 1];
  fill_n(dindex, nbRow2Delete, -1);
  int cnt(0);

  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin(); iPtr != indexSetOfRow2Delete.end(); iPtr++)
    {
      dindex[cnt++] = *iPtr;
    }

  begIndex = dindex[cnt-1];
  endIndex = dindex[cnt-1];
  for(int cpt = cnt-2; cpt>=0; cpt--)
    {
      if( dindex[cpt] == (begIndex - 1) )
	    begIndex = dindex[cpt];
      else
        {
	        bapcodInit().check(CPXdelrows(cplexEnv, cplexProbPtr, begIndex, endIndex),
                               "CPXdelrows: could not delete cols");

#ifndef _MSC_VER
	  /**
	   * Only available on Unix, 
	   * due to an ILOG portability issue, 
	   * it doesn't work on Windows.
	   */
	  if ( log_cplex != NULL )
	    fflush(static_cast<FILE *>(log_cplex));
#endif


	  begIndex = dindex[cpt];
	  endIndex = dindex[cpt];
	}
    }
  bapcodInit().check(CPXdelrows(cplexEnv, cplexProbPtr, begIndex, endIndex), "CPXdelcols: could not delete cols");

#ifndef _MSC_VER
  /**
   * Only available on Unix, due to an ILOG portability issue, it doesn't work on Windows.
   */
  if ( log_cplex != NULL )
    fflush(static_cast<FILE *>(log_cplex));
#endif

  delete[] dindex;
  dindex = NULL;
  _nrow -= nbRow2Delete;

  return;
}

void LpCplexInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  bapcodInit().check(1, "LP can not have directives");
  return;
}

void LpCplexInterface::getObjCoef(ProbCoef & pc)
{
  double *coef = new double[1];
  bapcodInit().check(CPXgetobj(cplexEnv, cplexProbPtr, coef, pc.colRef, pc.colRef), "could not getobj");
  pc.coef = *coef;
  delete [] coef;
  return;
}

void LpCplexInterface::chgObjCoef(const ProbCoef & pc)
{
  int *ref = new int[1];
  *ref = pc.colRef;
  double *coef = new double[1];
  *coef = pc.coef;
  bapcodInit().check(CPXchgobj(cplexEnv, cplexProbPtr, 1, ref, coef), "could not chgobj");

  delete [] ref;
  delete [] coef;
  return;
}

void LpCplexInterface::chgMatCoef(const ProbCoef & pc)
{
  bapcodInit().check(CPXchgcoef(cplexEnv, cplexProbPtr, pc.rowRef, pc.colRef, double(pc.coef)), "could not chgcof");
  return;
}

void LpCplexInterface::chgRhs(const ProbBound & pb)
{
  int *index = new int[1];
  index[0] = pb.ref;
  double *rhs = new double[1];
  rhs[0] = zero(pb.bound);
  bapcodInit().check(CPXchgrhs(cplexEnv, cplexProbPtr, 1, index, rhs), "could not chgths");
  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}

void LpCplexInterface::chgRhs(const ProbRhsContainer & newRhs)
{
  if (newRhs.empty()) return;
  int nrhs(0);
  int *index = new int[newRhs.size()];
  double *rhs = new double[newRhs.size()];

  for (ProbRhsContainer::const_iterator bPtr = newRhs.begin(); bPtr != newRhs.end(); bPtr++)
    {
      index[nrhs] = bPtr->ref;
      rhs[nrhs] = bPtr->bound;
      nrhs++;
    }

  bapcodInit().check(CPXchgrhs(cplexEnv, cplexProbPtr, nrhs, index, rhs), "could not chgths");

  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}

void LpCplexInterface::chgBds(const ProbBoundContainer & newBounds)
{
  if (newBounds.empty())
    return;
  int Lnbnds(0);
  int * Lmindex = new int[newBounds.size()];
  fill_n(Lmindex, newBounds.size(), 0);

  char *Lqbtype = new char[newBounds.size()];
  fill_n(Lqbtype, newBounds.size(), ' ');

  double *Lbnd = new double[newBounds.size()];
  fill_n(Lbnd, newBounds.size(), 0);

  for (ProbBoundContainer::const_iterator bPtr = newBounds.begin(); bPtr != newBounds.end(); bPtr++)
    {
      if (printL(7))
	    std::cout << "bound = " << *bPtr;

      Lmindex[Lnbnds] = bPtr->ref;
      Lbnd[Lnbnds]= zero(bPtr->bound);
      if (bPtr->sense == 'L')
        {
          Lqbtype[Lnbnds] = 'L';
        }
      else if ((bPtr->sense == 'U') || (bPtr->sense == 'I') || (bPtr->sense == 'B'))
        {
          Lqbtype[Lnbnds] = 'U';
        }
      else if (bPtr->sense == 'F')
        {
          Lqbtype[Lnbnds] = 'B';
        }
      Lnbnds++;
    }
  bapcodInit().check(CPXchgbds(cplexEnv, cplexProbPtr, Lnbnds, Lmindex, Lqbtype, Lbnd), "could not chgbds");

  delete[] Lmindex;
  Lmindex = NULL;
  delete[] Lqbtype;
  Lqbtype = NULL;
  delete[] Lbnd;
  Lbnd = NULL;
  return;
}

void LpCplexInterface::chgColType(const std::set<ProbType> & newTypes)
{
  bapcodInit().check(1,"LP Form cannot have ColType");

  return;
}

void LpCplexInterface::getObjVal(Double & objV)
{
  double objval;
  bapcodInit().check(CPXgetobjval(cplexEnv, cplexProbPtr, &objval), "could not get objective value");
  objV = zero(objval, cplexPrecision);

  return;
}


void LpCplexInterface::getDualBound(Double & val)
{
  double dualval = 0.0;
  
  int nrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  double * dualsol = new double[nrow];
  bapcodInit().check(CPXgetpi(cplexEnv, cplexProbPtr, dualsol, 0, nrow-1), " could not get dual solution");

  double * rhs = new double[nrow];
  bapcodInit().check(CPXgetrhs(cplexEnv, cplexProbPtr, rhs, 0, nrow-1), " could not get rhs");

  int ncol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  
  double * reducedcost = new double[ncol];
  bapcodInit().check(CPXgetdj (cplexEnv, cplexProbPtr, reducedcost, 0, ncol-1), " could not get reducedcosts");

  double * x = new double[ncol];
  bapcodInit().check(CPXgetx(cplexEnv, cplexProbPtr, x, 0, _ncol-1), " could not get primal solution");

  for(int i = 0; i < ncol; i++)
    dualval +=  x[i] * reducedcost[i];

  for(int i = 0; i < nrow; i++)
    dualval += dualsol[i] * rhs[i];

  val = zero(dualval, cplexPrecision);

  return;
}


void LpCplexInterface::getPrimalBound(Double & val)
{
  double primalval = 0.0;
  bapcodInit().check(CPXgetobjval(cplexEnv, cplexProbPtr, &primalval), "could not get primal bound");
  val = zero(primalval, cplexPrecision);
  return;
}

bool LpCplexInterface::getOptimStatus(SolutionStatus & lpStatus, SolutionStatus & mipStatus)
{
  int status;
  status = CPXgetstat(cplexEnv, cplexProbPtr);

  if (status == CPX_STAT_OPTIMAL)
    {
      lpStatus = SolutionStatus::Optimum;
      mipStatus = SolutionStatus::Optimum;

      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: LP optimal" << std::endl;

      return(true);
    }

  if (printL(4))
    MPSwrite();

  if ((status == CPX_STAT_INFEASIBLE) || (status == CPX_STAT_INForUNBD))
    {
      lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: LP infeasible" << std::endl;

      return(false);
    }

  if (status == CPX_STAT_UNBOUNDED)
    {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
	     std::cout << "LpCplexInterface::getOptimStatus: LP unbounded" << std::endl;

      if (printL(-1))
         std::cout << "BaPCod WARNING : status of Cplex is 'Unbounded'" << std::endl;

      std::cerr << "BaPCod WARNING : status of Cplex is 'Unbounded'" << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_OBJ_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(0))
        {
	      std::cout << "LpCplexInterface::getOptimStatus: objective limit exceeded in phase II";
          double lowerLimit, upperLimit;
          CPXgetdblparam(cplexEnv, CPX_PARAM_OBJLLIM, &lowerLimit);
          CPXgetdblparam(cplexEnv, CPX_PARAM_OBJULIM, &upperLimit);
          std::cout << " lower limit " << lowerLimit
                    << " upper limit " << upperLimit << std::endl;
        }
      return(false);
    }

  if (status == CPX_STAT_ABORT_IT_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: iteration limit exceeded" << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_TIME_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: time limit exceeded" << std::endl;

      return(false);
    }

  if (status == CPX_STAT_NUM_BEST)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: problem non-optimal" << std::endl;

      return(false);
    }

    if (status == CPX_STAT_OPTIMAL_INFEAS)
    {
      mipStatus = SolutionStatus::OptimumUnscalInfeas;
      lpStatus = SolutionStatus::OptimumUnscalInfeas;

      if (printL(0))
        std::cout << "BaPCod info : status of Cplex is 'Optimal solution found, unscaled infeasibilities'" << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_USER)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "LpCplexInterface::getOptimStatus: Aborded" << std::endl;

      return(false);
    }
  std::cout << "LpCplexInterface::getOptimStatus: undefined status" << std::endl;

  return(false);
}


void LpCplexInterface::getSol(std::map<int, Double> & primSol, const bool & ifPrint)
{
  primSol.clear();
  int readNcol = 0;
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNcol <= _ncol, "LpCplexInterface::getSol: readNcol > _ncol");

  double *x = new double[_ncol];
  bapcodInit().check(CPXgetx(cplexEnv, cplexProbPtr, x, 0, _ncol-1), " could not get primal solution");

  if (printL(6))
    std::cout << "readNcol = " << readNcol << "  _ncol = " << _ncol << std::endl;

  for (int jcol = 0; jcol < _ncol; jcol++)
    if (zero(x[jcol], cplexPrecision) != 0)
    {
      if (printL(6))
        std::cout << "primSol[" << jcol << "] = " << x[jcol] << std::endl;
      
      primSol[jcol] = x[jcol];
    }
  
  delete[] x;
  x = NULL;
  return;
}


void LpCplexInterface::getSol(std::map<int, Double> & primSol,
			      std::map<int, Double> & dualSol,
			      const int & minmaxStatus,
			      const bool & ifPrint)
{
  getSol(primSol, ifPrint);
  dualSol.clear();
  int readNrow(0);
  readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNrow <= _nrow, "LpCplexInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[_nrow];
  fill_n(dual, _nrow, 0);

  char *rsense = new char[_nrow];
  fill_n(rsense, _nrow, ' ');

  bapcodInit().check(CPXgetpi(cplexEnv, cplexProbPtr, dual, 0, _nrow-1), " could not get dual solution");
  
  bapcodInit().check(CPXgetsense(cplexEnv, cplexProbPtr, rsense, 0, _nrow-1), " could not get row type");

  for (int irow = 0; irow < _nrow; irow++)
    if (zero(dual[irow], cplexPrecision) != 0)
      {
          if (rsense[irow] == 'L')
              dualSol[irow] = minmaxStatus * flipsign * dual[irow];
          else if (rsense[irow] == 'G')
              dualSol[irow] = minmaxStatus * flipsign * dual[irow];
          else
              dualSol[irow] = flipsign * dual[irow];

          if (printL(6))
          {
              std::cout << "dual[" << irow << "] = " << dual[irow] << std::endl;
              std::cout << " dualSol[" << irow << "] = " << dualSol[irow] << std::endl;
          }
      }
  delete[] dual;
  dual = NULL;
  delete[] rsense;
  rsense = NULL;
  return;
}


void LpCplexInterface::getSolPool(std::vector< std::map<int, Double> > & primSolPool, const bool & ifPrint)
{
  bapcodInit().check(1, "LpCplexInterface::getSolPool() should not be called since LP's do not have a solution pool");
  return;
}

void LpCplexInterface::getReducedCost(std::map<int, Double> & redCost, const bool & ifPrint)
{
  redCost.clear();
  int readNcol = 0;
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNcol <= _ncol, "LpCplexInterface::getSol: readNcol > _ncol");

  double *dj = new double[_ncol];
  fill_n(dj, _ncol, 0);
  bapcodInit().check(CPXgetdj(cplexEnv, cplexProbPtr, dj, 0, _ncol-1), " could not get reduced cost");

  if (printL(6))
    std::cout << "readNcol = " << readNcol << "  _ncol = " << _ncol << std::endl;

  for (int jcol = 0; jcol < _ncol; jcol++)
    if (zero(dj[jcol], cplexPrecision) != 0)
      {
          if (printL(6))
              std::cout << "redCost[" << jcol << "] = " << dj[jcol] << std::endl;
          redCost[jcol] = dj[jcol];
      }

  delete[] dj;
  dj = NULL;
  return;
}

void LpCplexInterface::resetUpperCutOffValue()
{
}

void LpCplexInterface::setUpperCutOffValue(const Double & cutOff)
{
}

void LpCplexInterface::setSolveFromScratch(const bool exactSolution)
{
}

void LpCplexInterface::resetSolveFromScratch()
{
}

void LpCplexInterface::setAbortValue(const Double& abortValue)
{
 // DO NOTHING
}

void LpCplexInterface::getBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
  colStatus.resize(CPXgetnumcols(cplexEnv, cplexProbPtr));
  rowStatus.resize(CPXgetnumrows(cplexEnv, cplexProbPtr));
  if (CPXgetbase(cplexEnv, cplexProbPtr, &colStatus[0], &rowStatus[0]))
  {
      if (printL(2))
        std::cout << "BaPCod info: could not get basis" << std::endl;
  }
  // don't need to convert status because constants are defined as CPLEX does
}

void LpCplexInterface::setBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
  // don't need to convert status because constants are defined as CPLEX does

  if (!rowStatus.empty())
  {
    bapcodInit().check(CPXcopybase(cplexEnv, cplexProbPtr, &colStatus[0],&rowStatus[0]),
                       " could not set basis");
  }
}

void LpCplexInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  int irow, jcol, iel, nchar;
  int readNrow, readNcol, idone, readNglents;
  char coltype;
  double readBdl, readBdu;
  char probname[4 * wordSize];
  int pnsurplus(0);

  bapcodInit().check(CPXgetprobname(cplexEnv, cplexProbPtr, probname, 4 * wordSize, &pnsurplus),
                     " could not getprobname");

  bapcodInit().require(pnsurplus >= 0,"CPXgetprobname did not have enough space ");

  readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);

  int cnsurplus(0);
  char *colnamestore = new char[readNcol * placeSize];
  fill_n(colnamestore, readNcol * placeSize, endOfWord);
  char **colnames = new char*[readNcol];
  fill_n(colnames, readNcol, static_cast<char*>(NULL));

  for(int ic = 0; ic < readNcol; ic++)
    colnames[ic] = &(colnamestore[ic*placeSize]);

  if (param().MipSolverRecordNamesInFormulation)
    bapcodInit().check(CPXgetcolname(cplexEnv, cplexProbPtr, colnames, colnamestore,
                                     readNcol * placeSize, &cnsurplus, 0,readNcol-1),
                       " could not getcolname");

  bapcodInit().require(cnsurplus >= 0,"CPXgetcolname did not have enough space ");
  printf("PROBLEM: %s whose formulation ref number is %d\n", probname, ref());

  /// Output Obj
  if ( minmaxStatus ==1 )
      printf("Minimize\n");
  else
      printf("Maximize\n");

  int *readMclind = new int[readNcol];
  fill_n(readMclind, readNcol, 0);
  double *readDmatval = new double[readNcol];
  fill_n(readDmatval, readNcol, 0);
  bapcodInit().check(CPXgetobj(cplexEnv, cplexProbPtr, readDmatval,0,readNcol-1),
                     " could not getobj");

  nchar = printf(" Objectif:");
  for (iel = 0; iel < readNcol; iel++)
    if (Double(readDmatval[iel]) != 0)
      {
          if(param().MipSolverRecordNamesInFormulation)
          {
              nchar += printf(" %+.9g %s", readDmatval[iel], &(colnamestore[placeSize * iel]));
          } else {
              nchar += printf(" %+.9g v_%d", readDmatval[iel], iel);
          }

          if (nchar >= 80-9-13)
          {
              printf("\n");
              nchar = 0;
          }
      }

  /// Output Matrix
  printf("\nSubject To\n");
  for(irow = 0; irow < readNrow; irow++)
    printRow(irow, readNcol, readMclind, readDmatval, colnamestore);
  printf("Bounds\n");

  for (jcol = 0; jcol < readNcol; jcol++)
    {
      bapcodInit().check(CPXgetlb(cplexEnv, cplexProbPtr, &readBdl, jcol, jcol)," could not getlb");

      bapcodInit().check(CPXgetub(cplexEnv, cplexProbPtr, &readBdu, jcol, jcol)," could not getub");

      if (param().MipSolverRecordNamesInFormulation)
      {
        if (readBdu < CPX_INFBOUND / 10)
	        printf(" %.9g <= %s <= %.9g\n", readBdl, &colnamestore[placeSize * jcol], readBdu);
        else if (readBdl != 0.0)
	        printf("      %s >= %.9g\n", &colnamestore[placeSize * jcol], readBdl);
      } 
      else {
        if (readBdu < CPX_INFBOUND / 10)
	        printf(" %.9g <= v_%d <= %.9g\n", readBdl, jcol, readBdu);
        else if (readBdl != 0.0)
	        printf("      v_%d >= %.9g\n", jcol, readBdl);
      }
    }

  /// Always Pure for LPSolve yet
  if (_pureLP)
    printf("End\n");
  else
    {
      readNglents = CPXgetnumint(cplexEnv, cplexProbPtr) + CPXgetnumbin(cplexEnv, cplexProbPtr);
      CPXgetnumsos(cplexEnv, cplexProbPtr);
      if ( readNglents )
        {
          printf("Integers\n");
          for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
            {
              bapcodInit().check(CPXgetctype(cplexEnv, cplexProbPtr, &coltype, jcol, jcol),
                                 " could not getctype");

              if (coltype != 'C')
              {
                  if(param().MipSolverRecordNamesInFormulation)
                  {
                      printf(" %s ", &colnamestore[placeSize * jcol]);
                  } else {
                      printf(" v_%d ", jcol);
                  }
                  idone++;
              }
              if (idone >= wordSize) {
		        printf("\n");
                idone = 0;
              }
            }
          if (idone != 0)
	        printf("\n");
        }

      printf("End\n");

    }

  delete[] readMclind;
  readMclind = NULL;
  delete[] readDmatval;
  readDmatval = NULL;
  delete[] colnamestore;
  colnamestore = NULL;
  delete[] colnames;
  colnames = NULL;
  return;
}


void LpCplexInterface::MPSwrite()
{
  bapcodInit().check(CPXwriteprob(cplexEnv, cplexProbPtr,"curprob.lp","LP"),
                     " could not write prob",
                     ProgStatus::terminate);

  bapcodInit().check(CPXwriteprob(cplexEnv, cplexProbPtr,"curprob.mps","MPS"),
                     " could not write prob",
                     ProgStatus::terminate);

  bapcodInit().check(CPXwriteprob(cplexEnv, cplexProbPtr, "curprob.sav","SAV"),
                     " could not write prob",
                     ProgStatus::terminate);

  bapcodInit().check(CPXmbasewrite(cplexEnv, cplexProbPtr, "curbasis.bas"),
                     " could not write basis information",
                     ProgStatus::run);

  bapcodInit().check(CPXwriteparam(cplexEnv,"curparams.prm"),
                     " could not write the current parameters",
                     ProgStatus::terminate);

  return;
}


void LpCplexInterface::printRow(int irow, const int & readNcol, int *readMclind, double *readDmatval, char *readCnames)
{
  int iel, readNels(0), nchar, nbRows(1);
  double readDrhs;
  char readRtype, readSense[3];
  int readMstart[2];
  int rsurplus(0);
  bapcodInit().check(CPXgetrows(cplexEnv, cplexProbPtr, &readNels, readMstart, readMclind, readDmatval,
                                readNcol, &rsurplus,irow, irow),
                     "could not getrows");

  bapcodInit().check(CPXgetrhs(cplexEnv, cplexProbPtr, &readDrhs, irow, irow),"could not getrhs");

  bapcodInit().check(CPXgetsense(cplexEnv, cplexProbPtr, &readRtype, irow, irow),"could not getrowtype");

  int rnsurplus(0);
  char *rownamestore = new char[placeSize * nbRows];
  fill_n(rownamestore, placeSize * nbRows, endOfWord);

  char **rownames = new char*[nbRows];
  fill_n(rownames, nbRows, static_cast<char*>(NULL));

  for(int ir = 0; ir < nbRows; ir++)
    rownames[ir] = &(rownamestore[ir*placeSize]);

  if (param().MipSolverRecordNamesInFormulation)
    bapcodInit().check(CPXgetrowname(cplexEnv, cplexProbPtr, rownames, rownamestore, placeSize * nbRows,
                                     &rnsurplus, irow, irow),
                       "could not getnames");

  bapcodInit().require(rnsurplus >= 0,"CPXgetrowname did not get enough space");

  strcpy(readSense,readRtype=='O' ? "$" : (readRtype=='E' ?"=" :(readRtype=='L' ?"<=" :">=" ) ));

  if (param().MipSolverRecordNamesInFormulation) {
    nchar = printf(" %s:", rownamestore);
  } else {
    nchar = printf(" c_%d:", irow);
  }      
  
  if (readNels > 0)
    for (iel = 0; iel < readNels; iel++)
      {
        if(param().MipSolverRecordNamesInFormulation)
        {
          nchar += printf(" %+.9g %s", readDmatval[iel], &readCnames[placeSize * readMclind[iel]]);
        } else {
          nchar += printf(" %+.9g v_%d", readDmatval[iel], readMclind[iel]);
        }
        if (nchar >= 80-9-13) {
            printf("\n");
            nchar = 0;
        }
      }

  printf( " %s %.9g", readSense, readDrhs);
  printf("\n");

  //added by Issam to fix a memory leak
  delete[] rownamestore;
  delete[] rownames;
  return;
}

void LpCplexInterface::optimise(const int & minmaxStatus,
                                const double & BarrierConvergenceTolerance,
                                const double & rightHAndSideZeroTol,
                                const double & reducedCostTolerance,
				                const bool & preprocessorOn,
				                const bool & probingOn,
				                const bool & automaticCuttingPlanesOn,
				                const char & solverSelection)
{
  if (printL(1))
    std::cout << "LpCplexInterface::optimise setting param rightHAndSideZeroTol to " << rightHAndSideZeroTol
              << std::endl;

  bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_EPRHS, rightHAndSideZeroTol),
                     "Failure to turn on rightHAndSideZeroTol ");

  if (printL(1))
    std::cout << "LpCplexInterface::optimise setting param reducedCostTolerance to " << reducedCostTolerance
              << std::endl;

  bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_EPOPT, reducedCostTolerance),
                     "Failure to turn on reducedCostTolerance ");

  if (preprocessorOn)
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_PREIND,CPX_ON),
                       "Failure to turn on presolve ");
  else
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_PREIND,CPX_OFF),
                       "Failure to turn off presolve ");

  bapcodInit().require(_formCurrentlyLoaded, "Form not Currently Loaded", ProgStatus::quit, 3);
  if (printL(8))
    MPSwrite();

  bapcodInit().check(CPXchgprobtype(cplexEnv, cplexProbPtr, CPXPROB_LP), "could not change type to LP");

  switch (solverSelection)
  {
    case 'p' :
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL),
                         "Failure to change the LP method");
      break;
    case 'd' :
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL),
                         "Failure to change the LP method");
      break;
    case 'n' :
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_NET),
                         "Failure to change the LP method");
      break;
    case 'b' :
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER),
                         "Failure to change the LP method to Barrier");

      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_BARCROSSALG, -1),
                         "Failure to turn off cross-over");

      bapcodInit().check(CPXsetdblparam(cplexEnv, CPXPARAM_Barrier_ConvergeTol, BarrierConvergenceTolerance),
                         "Failure to change CPXPARAM_Barrier_ConvergeTol");

      bapcodInit().check(CPXsetintparam(cplexEnv, CPXPARAM_Barrier_StartAlg, 4),
                         "Failure to change CPXPARAM_Barrier_StartAlg");
      break;
    case 'c' :
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER),
                         "Failure to change the LP method");
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_BARCROSSALG, 0),
                         "Failure to change crossover mode in LP method");
      break;
    case 'a' :
    default:
      bapcodInit().check(CPXsetintparam (cplexEnv, CPX_PARAM_LPMETHOD, CPX_ALG_AUTOMATIC),
                         "Failure to change the LP method");
      break;
  }

  /// Ruslan : sometimes these parameters are changed during optimization by Cplex (do not know why)
  CPXsetdblparam(cplexEnv, CPX_PARAM_OBJLLIM, -1e+75);
  CPXsetdblparam(cplexEnv, CPX_PARAM_OBJULIM, 1e+75);

  if (printL(1))
    bapcodInit().check(CPXwriteprob(cplexEnv, cplexProbPtr, "curprob.sav","SAV"),
                       " could not write prob", ProgStatus::run);

  bapcodInit().check(CPXlpopt(cplexEnv,cplexProbPtr),"could not solve LP by CPXlpopt");

  return;
}



void LpCplexInterface::makeSpaceForLoadingForm()
{
  _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << " _formCurrentlyLoaded = " << _formCurrentlyLoaded << std::endl;
  return;
}

void LpCplexInterface::saveCopyOfCurForm()
{
  if (printL(6))
    std::cout << "LpCplexInterface::saveCopyOfCurForm(): formCurrentlyLoaded = " << _formCurrentlyLoaded << std::endl;
  return;
}

void LpCplexInterface::load()
{
  if (printL(6))
    std::cout << "LpCplexInterface::load(): formCurrentlyLoaded = " << _formCurrentlyLoaded << std::endl;
  return;
}

void LpCplexInterface::unload(const bool & deleteMat)
{
  return;
}


void LpCplexInterface::reset()
{
  return;
}

// ************************************
// ***** Class MipCplexInterface ******
// ************************************

MipCplexInterface::MipCplexInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name)
  : LpCplexInterface(bapcodInit, ref, name), _lazyCallbackStruct(NULL, 0)
{
  _pureLP = false;
  bapcodInit->check(cplexProbPtr == NULL,"MipSolverInterface: problem must have been  created");

  bapcodInit->check(CPXchgprobtype(cplexEnv, cplexProbPtr, CPXPROB_MILP),
                    " CPXchgprobtype: Failed to change typ\e to MILP.");
  return;
}

void MipCplexInterface::loadFormulation(const std::string & name,
					                    const int & minmaxStatus,
					                    const int & ncol,
					                    const int & nrow,
					                    const std::map<int, std::string> & mapSeqnb2Cname,
					                    const std::map<int, std::string> & mapSeqnb2Rname,
					                    const ProbCoefContainer & objectRow,
                                        ProbColMatrixContainer & colMatrix,
					                    const ProbRhsContainer & rhsv,
					                    const ProbBoundContainer & bounds,
					                    const std::set<ProbType> & types,
					                    const std::set<ProbIntC> & directs,
					                    const std::set<ProbSetCoef> & sets)
{
  if (printL(5))
    std::cout << "MipCplexInterface::loadFormulation() ncol = " << ncol << " , nrow = " << nrow << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  _pureLP = false;
  /// Assumes sorted colMatrix
#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
  std::sort(colMatrix.begin(), colMatrix.end());
#endif
  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);
  for (ProbCoefContainer::const_iterator oPtr = objectRow.begin(); oPtr != objectRow.end(); oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  double *rhs = new double[nrow];
  fill_n(rhs, nrow, 0);

  char *rsense = new char[nrow];
  fill_n(rsense, nrow, ' ');

  for (ProbRhsContainer::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
    {
      rhs[rvPtr->ref] = zero(rvPtr->bound);
      rsense[rvPtr->ref] = rvPtr->sense;
    }

  int *matbeg = new int[ncol + 1];
  fill_n(matbeg, ncol + 1, 0);

  int *matcnt = new int[ncol + 1];
  fill_n(matcnt, ncol + 1, 0);

  int nnze = colMatrix.size();

  int *matind = new int[nnze];
  fill_n(matind, nnze, -1);

  double *matval = new double[nnze];
  fill_n(matval, nnze, 0);

  int cntNbNz(0);
  ProbColMatrixContainer::const_iterator mPtr = colMatrix.begin();
  int curCol(0);
  for (;curCol<ncol;curCol++)
    {
      matbeg[curCol] = cntNbNz;
      matcnt[curCol] = 0;

      while ((mPtr != colMatrix.end()) && (mPtr->colRef==curCol))
      {
          matind[cntNbNz] = mPtr->rowRef;
          matval[cntNbNz] = zero(mPtr->coef);
          cntNbNz++;
          matcnt[curCol]++;
          mPtr++;
      }
    }
  /// For dummy  element ncol + 1
  matbeg[curCol] = cntNbNz;
  matcnt[curCol] = 0;


  double *dlb = new double[ncol];
  fill_n(dlb, ncol, 0);

  double *dub = new double[ncol];
  fill_n(dub, ncol, BapcodInfinity);

  for (ProbBoundContainer::const_iterator bPtr = bounds.begin(); bPtr != bounds.end(); bPtr++)
    {
        if (bPtr->sense == 'U')
            dub[bPtr->ref] = zero(bPtr->bound);
        else if (bPtr->sense == 'L')
            dlb[bPtr->ref] = zero(bPtr->bound);
        else if (bPtr->sense == 'F')
        {
            dlb[bPtr->ref] = zero(bPtr->bound);
            dub[bPtr->ref] = dlb[bPtr->ref];
        }
    }

  double *range = NULL;
  std::string curName;
  char *colnamestore(NULL);
  char **colnames = NULL;
  char *rownamestore(NULL);
  char **rownames = NULL;

  if (param().MipSolverRecordNamesInFormulation)
    {
      colnamestore = new char[ncol * placeSize];
      fill_n(colnamestore, ncol * placeSize, endOfWord);

      colnames = new char*[ncol];
      fill_n(colnames, ncol, static_cast<char*>(NULL));

      rownamestore = new char[nrow * placeSize];
      fill_n(rownamestore, nrow * placeSize, endOfWord);

      rownames = new char*[nrow];
      fill_n(rownames, nrow, static_cast<char*>(NULL));

      for(int ic = 0; ic < ncol; ic++)
      {
          if (mapSeqnb2Cname.count(ic))
              curName = mapSeqnb2Cname.at(ic) + filling;
          else
              curName = filling;

          strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
          colnamestore[(ic +1) * placeSize - 1] = endOfWord;
          colnames[ic]= &(colnamestore[ic*placeSize]);
      }

      for (int ir = 0; ir < nrow; ir++)
      {
          if (mapSeqnb2Rname.count(ir))
              curName = mapSeqnb2Rname.at(ir) + filling;
          else
              curName = filling;

          strncpy(&(rownamestore[ir * placeSize]), curName.c_str(), wordSize);
          rownamestore[(ir +1) * placeSize - 1] = endOfWord;
          rownames[ir]= &(rownamestore[ir*placeSize]);
      }

      /// Not required since it is done in CPXcopyctype
      bapcodInit().check(cplexProbPtr == NULL,
                         "MipCplexInterface.load(): problem must have been  created");

      bapcodInit().check(CPXcopylpwnames(cplexEnv, cplexProbPtr, ncol, nrow, minmaxStatus, obj, rhs, rsense,
                                         matbeg, matcnt, matind, matval, dlb, dub, range, colnames, rownames),
                         "CPXcopylpwnames: Failed to copy problem data.");
    }
  else
    {
      bapcodInit().check(cplexProbPtr == NULL, "MipCplexInterface.load(): problem must have been  created");

      bapcodInit().check(CPXcopylp(cplexEnv, cplexProbPtr, ncol, nrow, minmaxStatus, obj, rhs, rsense,
                                   matbeg, matcnt, matind, matval, dlb, dub, range),
                         "CPXcopylp: Failed to copy problem data.");
    }

  char *ctype = new char[ncol];
  fill_n(ctype, ncol, 'C');
  for (set<ProbType>::const_iterator typePt = types.begin(); typePt != types.end(); typePt++)
    {
      if (typePt->type == 'I')
	    ctype[typePt->ref] = 'I';

      if (typePt->type == 'B')
	    ctype[typePt->ref] = 'B';

      if (typePt->type == 'S')
	    ctype[typePt->ref] = 'S';

      if (typePt->type == 'N')
	    ctype[typePt->ref] = 'N';
    }

  bapcodInit().check(CPXcopyctype(cplexEnv, cplexProbPtr, ctype ), "CPXcopyctype: Failed to add integer type info");

  int ndir = 0;
  int *pricind = NULL;
  int *prival = NULL;
  int *pridir = NULL;

  if (!directs.empty())
    {
      pricind = new int[directs.size()];
      fill_n(pricind, directs.size(), -1);

      prival = new int[directs.size()];
      fill_n(prival, directs.size(), 0);

      pridir = new int[directs.size()];
      fill_n(pridir, directs.size(), 0);

      for (set<ProbIntC>::const_iterator dPtr = directs.begin(); dPtr != directs.end(); dPtr++)
	    /// Neg value refer to sets.
        if (dPtr->ref >= 0)
          {
            pricind[ndir] = dPtr->ref;
            prival[ndir] = (int) dPtr->val;
            if (dPtr->type == 'D')
                pridir[ndir] = CPX_BRANCH_DOWN;
            else if (dPtr->type == 'U')
	            pridir[ndir] = CPX_BRANCH_UP;
            else
	            pridir[ndir] = CPX_BRANCH_GLOBAL;
            ndir++;
          }
    }

  /// Higher priority index means chosen first for branching
  if (ndir > 0)
    {
      bapcodInit().check(CPXcopyorder(cplexEnv, cplexProbPtr, ndir, pricind, prival, pridir ),
                         "CPXcopyorder: Failed to set branching priorities info",
                         ProgStatus::terminate);
    }

  delete [] obj;
  obj = NULL;
  delete [] rhs;
  rhs = NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] matbeg;
  matbeg = NULL;
  delete [] matcnt;
  matcnt = NULL;
  delete [] matind;
  matind = NULL;
  delete [] matval;
  matval = NULL;
  delete [] dlb;
  dlb = NULL;
  delete [] dub;
  dub = NULL;
  delete [] range;
  range = NULL;
  delete [] colnames;
  colnames = NULL;
  delete [] rownames;
  rownames = NULL;
  delete [] colnamestore;
  colnamestore = NULL;
  delete [] rownamestore;
  rownamestore = NULL;
  delete [] ctype;
  ctype = NULL;
  delete [] pricind;
  pricind = NULL;
  delete [] prival;
  prival = NULL;
  delete [] pridir;
  pridir = NULL;

  if (!sets.empty())
    {
      int numsos = 0;
      int numsosnz = 0;
      char *sostype = NULL;
      int *sospri = NULL;
      int *sosbeg = NULL;
      int *sosind = NULL;
      double *sosref = NULL;

      std::set< int > setNbs;
      {
        for (set<ProbSetCoef>::const_iterator sPtr = sets.begin(); sPtr != sets.end(); sPtr++)
          setNbs.insert(sPtr->setRef);
      }

      numsos = setNbs.size();
      sostype = new char[numsos];
      fill_n(sostype, numsos, ' ');
      sospri = new int[numsos];
      fill_n(sospri, numsos, 1);
      sosbeg = new int[numsos];
      fill_n(sosbeg, numsos, 0);
      sosind = new int[sets.size()];
      fill_n(sosind, sets.size(), -1);
      sosref = new double[sets.size()];
      fill_n(sosref, sets.size(), 0);

      set<ProbSetCoef>::const_iterator sPtr = sets.begin();
      while (sPtr != sets.end())
        {
          int curSet = sPtr->setRef;
          sostype[curSet] = sPtr->setType;
          sospri[curSet] = (int) sPtr->setPriorityIndex;
          sosbeg[curSet] = numsosnz;
          do
            {
              sosind[numsosnz] = sPtr->colRef;
              sosref[numsosnz] = sPtr->coef;
              numsosnz++;
              sPtr++;
            }
          while ((sPtr != sets.end()) && (sPtr->setRef == curSet));
        }

      bapcodInit().check(CPXcopysos(cplexEnv, cplexProbPtr, numsos, numsosnz, sostype, sosbeg, sosind, sosref, NULL),
                         "CPXcopysos: Failed to add sos info");

      delete [] sostype;
      sostype = NULL;

      delete [] sosbeg;
      sosbeg = NULL;

      delete [] sosind;
      sosind = NULL;

      delete [] sosref;
      sosref = NULL;
    }

  return;
}


void MipCplexInterface::chgColType(const std::set<ProbType> & newTypes)
{
  int Lnels = 0;
  int *Lmindex = new int[newTypes.size()];
  fill_n(Lmindex, newTypes.size(), 0);
  char *Lqctype = new char[newTypes.size()];
  fill_n(Lqctype, newTypes.size(), ' ');

  for (set<ProbType>::const_iterator bPtr = newTypes.begin(); bPtr != newTypes.end(); bPtr++)
    {
      /// Set data structure for chgcoltype
      Lmindex[Lnels] = bPtr->ref;
      Lqctype[Lnels] = bPtr->type;
      Lnels++;
    }
  bapcodInit().check(CPXchgctype(cplexEnv, cplexProbPtr, Lnels, Lmindex, Lqctype), "could not chgcoltype");

  /// Delete extra arrays
  delete[] Lmindex;
  Lmindex = NULL;
  delete[] Lqctype;
  Lqctype = NULL;

  return;
}

void MipCplexInterface::resetUpperCutOffValue()
{
    bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_CUTUP, 1e+75), "Failure to set CutUp parameter");
}

void MipCplexInterface::setUpperCutOffValue(const Double & cutOffValue)
{
  if (cutOffValue < BapcodInfinity)
    bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_CUTUP, cutOffValue.val()), "Failure to set CutUp parameter");
}

void MipCplexInterface::setLazyConstraintsCallback(MasterConf * mastConfPtr)
{
    if (mastConfPtr != NULL)
    {
        bapcodInit().check(CPXsetintparam(cplexEnv, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF));
        bapcodInit().check(CPXsetintparam(cplexEnv, CPXPARAM_Preprocessing_Linear, 0));

        _lazyCallbackStruct = LazyCallbackStruct(mastConfPtr, _ncol);
        bapcodInit().check(CPXsetlazyconstraintcallbackfunc(cplexEnv, lazycallback, &_lazyCallbackStruct),
                           "Failure to set lazy callback");
    }
}

void MipCplexInterface::removeLazyConstraintsCallback()
{
    bapcodInit().check(CPXsetlazyconstraintcallbackfunc(cplexEnv, NULL, NULL),
                       "Failure to remove lazy callback");
    bapcodInit().check(CPXfreelazyconstraints(cplexEnv, cplexProbPtr), "Failure to remove lazy constraints");
    bapcodInit().check(CPXsetintparam(cplexEnv, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_ON));
    bapcodInit().check(CPXsetintparam(cplexEnv, CPXPARAM_Preprocessing_Linear, 1));
}

void MipCplexInterface::setSolveFromScratch(const bool exactSolution)
{
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_ADVIND, 0), "Failure to set AdvInt parameter");

    if (!exactSolution)
    {
        int MIPemphasis = param().MIPemphasisInRestrictedMasterIpHeur();
        if ((MIPemphasis >= 0) && (MIPemphasis <= 4))
            bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_MIPEMPHASIS, MIPemphasis),
                               "Failure to set MIPEmphasis parameter");


        if (param().PolishingAfterTimeInRestrictedMasterIpHeur() > 0)
            bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_POLISHAFTERTIME,
                                              param().PolishingAfterTimeInRestrictedMasterIpHeur()),
                               "Failure to set PolishAfterTime parameter");
    }

}

void MipCplexInterface::resetSolveFromScratch()
{
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_ADVIND, 1),
                       "Failure to set AdvInt parameter");

    bapcodInit().check(CPXsetdblparam(cplexEnv, CPX_PARAM_POLISHAFTERTIME, 1e+75),
                       "Failure to set PolishAfterTime parameter");

    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BALANCED),
                       "Failure to set MIPEmphasis parameter");

}

void MipCplexInterface::setTimeLimit(const double & seconds)
{
    if (seconds < 0.1)
        bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_TILIM, 0.1),
                           "Failure to set time limit ");
    else
        bapcodInit().check(CPXsetdblparam (cplexEnv, CPX_PARAM_TILIM, seconds),
                           "Failure to set time limit ");
}

void MipCplexInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
     if(newDirects.empty()) return;
     int Lndir = 0;

     int *Lmcols = new int[newDirects.size()]; 
     fill_n(Lmcols, newDirects.size(), -1);

     int *Lmpri = new int[newDirects.size()]; 
     fill_n(Lmpri, newDirects.size(), 0);

     int *Lqbr = new int[newDirects.size()]; 
     fill_n(Lqbr, newDirects.size(), 0);

      for (set<ProbIntC>::const_iterator dPtr = newDirects.begin(); dPtr != newDirects.end(); dPtr++)
        {
          if (printL(6)) 
   			std::cout << *dPtr;

          Lmcols[Lndir] = dPtr->ref;
          Lmpri[Lndir] = (int) dPtr->val;

		  if (dPtr->type == 'D')
			Lqbr[Lndir] = CPX_BRANCH_DOWN;
          else if (dPtr->type == 'U')
			Lqbr[Lndir] = CPX_BRANCH_UP;
          else
			Lqbr[Lndir] = CPX_BRANCH_GLOBAL;
          Lndir++;
        }

	if (Lndir > 0)
    {
      bapcodInit().check(CPXcopyorder(cplexEnv, cplexProbPtr, Lndir, Lmcols, Lmpri, Lqbr ),
                         "CPXcopyorder: Failed to set branching priorities info",
                         ProgStatus::terminate);
    }

      if (printL(6))
        {
          std::cout << "Directives\n";
          for (int iel = 0; iel < Lndir; iel++)
   	        cout <<  "col[" << Lmcols[iel] << "] " << ", priority = " << Lmpri[iel] << ", sense = " << Lqbr[iel]
   	             << std::endl;
        }

      // delete extra arrays
      delete[] Lmcols; 
      Lmcols = NULL;

      delete[] Lmpri; 
      Lmpri = NULL;

      delete[] Lqbr; 
      Lqbr = NULL;

  return;
}

void MipCplexInterface::getObjVal(Double & objV)
{
  double objval = 0.0;
  bapcodInit().check(CPXgetobjval(cplexEnv, cplexProbPtr, &objval), "could not get objective value");

  objV = zero(objval, cplexPrecision);

  return;
}

void MipCplexInterface::getDualBound(Double & val)
{
  double dualval = 0.0;
  int status;
  status = CPXgetstat(cplexEnv, cplexProbPtr);
  if (status == CPXMIP_OPTIMAL)
    {
      bapcodInit().require(1, "MIP optimal");
      bapcodInit().check(CPXgetobjval(cplexEnv, cplexProbPtr, &dualval), "could not get dual bound");
    }
  else
  {
    bapcodInit().check(CPXgetbestobjval(cplexEnv, cplexProbPtr, &dualval), "could not get dual bound");
  }

  val = zero(dualval, cplexPrecision);

  return;
}

void MipCplexInterface::getPrimalBound(Double & val)
{
  double primalval = 0.0;
  bapcodInit().check(CPXgetobjval(cplexEnv, cplexProbPtr, &primalval), "could not get primal bound");

  val = zero(primalval, cplexPrecision);
  return;
}


bool MipCplexInterface::getOptimStatus(SolutionStatus & lpStatus, SolutionStatus & mipStatus)
{
  int status;
  status = CPXgetstat(cplexEnv, cplexProbPtr);
  if (status == CPXMIP_OPTIMAL)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: global search complete and an integer solution has been found"
		          << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }

  if (printL(4))
    MPSwrite();

    if (status == CPX_STAT_OPTIMAL_INFEAS)
    {
      mipStatus = SolutionStatus::OptimumUnscalInfeas;
      lpStatus = SolutionStatus::OptimumUnscalInfeas;
      if (printL(0))
        std::cout << "BaPCod info: Cplex status is 'Optimal solution found, unscaled infeasibilities'" << std::endl;

      return(false);
    }

  if (status == CPXMIP_OPTIMAL_TOL)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP optimal solution found within tolerance" << std::endl;

      // Optimal solution found
      return(true);
    }

  if ((status == CPXMIP_INFEASIBLE) || (status == CPXMIP_INForUNBD))
      {
      mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	     std::cout << "MipCplexInterface::getOptimStatus: MIP infeasible" << std::endl;

      return(false);
    }

  if (status == CPXMIP_SOL_LIM)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	     std::cout << "MipCplexInterface::getOptimStatus: MIP solutions limit exceeded" << std::endl;

      return(false);
    }

  if (status == CPXMIP_NODE_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
        std::cout << "MipCplexInterface::getOptimStatus: MIP node limit exceeded, integer solution exists" << std::endl;

      return(true);
    }

  if (status == CPXMIP_NODE_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP node limit exceeded, no integer solution exists"
                  << std::endl;

      return(false);
    }

  if (status == CPXMIP_TIME_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP time limit exceeded, integer solution exists" << std::endl;

      return(true);
    }

  if (status == CPXMIP_TIME_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP time limit exceeded, no integer solution exists"
                  << std::endl;

      return(false);
    }

  if (status == CPXMIP_FAIL_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP error termination, integer solution exists" << std::endl;

      return(true);
    }

  if (status == CPXMIP_FAIL_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP error termination, no integer solution exists"
                  << std::endl;

      return(false);
    }

  if (status == CPXMIP_MEM_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP Treememory limit exceeded, integer solution exists"
		          << std::endl;

      return(true);
    }

  if (status == CPXMIP_MEM_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	     std::cout << "MipCplexInterface::getOptimStatus: MIP Treememory limit  exceeded, no integer solution exists"
		           << std::endl;

      return(false);
    }

  if (status == CPXMIP_ABORT_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP aborded, integer solution exists" << std::endl;

      return(true);
    }

  if (status == CPXMIP_ABORT_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	     std::cout << "MipCplexInterface::getOptimStatus: MIP aborded, no integer solution exists" << std::endl;

      return(false);
    }

  if (status == CPXMIP_OPTIMAL_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP optimal with unscaled infeasibilities" << std::endl;

      return(false);
    }

  if (status == CPXMIP_FAIL_FEAS_NO_TREE)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP out of memory, no tree, integer solution exists"
                  << std::endl;

      return(true);
    }

  if (status == CPXMIP_FAIL_INFEAS_NO_TREE)
    {
      mipStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP out of memory, no tree, no integer solution exists"
		          << std::endl;

      return(false);
    }

  if (status == CPXMIP_NODE_LIM_FEAS)
    {
      mipStatus = lpStatus = SolutionStatus::PrimalFeasSolFound;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP node file size limit exceeded, integer solution exists"
		          << std::endl;

      return(true);
    }

  if (status == CPXMIP_NODE_LIM_INFEAS)
    {
      mipStatus = lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	    std::cout << "MipCplexInterface::getOptimStatus: MIP node file size limit  exceeded, no integer solution exists"
		          << std::endl;
      return(false);
    }

  if (status == CPXMIP_TIME_LIM_FEAS )
      {
          mipStatus = lpStatus = SolutionStatus::PrimalFeasSolFound;
          if (printL(3))
              std::cout << "LpCplexInterface::getOptimStatus: time limit exceeded" << std::endl;

          return(false);
      }

    if (status == CPXMIP_TIME_LIM_INFEAS )
    {
        mipStatus = lpStatus = SolutionStatus::UnSolved;
        if (printL(3))
            std::cout << "LpCplexInterface::getOptimStatus: time limit exceeded" << std::endl;

        return(false);
    }

    if (printL(1)) std::cout << "MipCplexInterface::getOptimStatus: undefined MIP status = " << status << std::endl;

    return LpCplexInterface::getOptimStatus(lpStatus, mipStatus);
}



void MipCplexInterface::getSol(std::map<int, Double> & primSol, const bool & ifPrint)
{
  primSol.clear();
  LpCplexInterface::getSol(primSol, ifPrint);
  return;
}

void MipCplexInterface::getSolPool(std::vector< std::map<int, Double> > & primSolPool, const bool & ifPrint)
{
  int readNcol = 0;
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNcol <= _ncol, "LpCplexInterface::getSol: readNcol > _ncol");

  const int nSol = CPXgetsolnpoolnumsolns(cplexEnv, cplexProbPtr);

  primSolPool.clear();

  primSolPool.reserve(nSol);

  bool noSolnIsNew = true;

  for (int soln = 0; soln < nSol; ++soln)
  {
    double objVal     = -1
    , sumVarVal       = -1
    , sumIntInfeas    = -1
    , sumPrimInfeas   = -1
    , maxVarVal       = -1
    , maxIntInfeas    = -1
    , maxPrimInfeas   = -1;


      noSolnIsNew = false;

      primSolPool.push_back(std::map<int, Double>());

      double *x = new double[_ncol];

      bapcodInit().check(CPXgetsolnpoolx(cplexEnv, cplexProbPtr, soln, x, 0, _ncol-1),
                         " could not get primal solution");

      if (printL(6))
        std::cout << "readNcol = " << readNcol << "  _ncol = " << _ncol << std::endl;

      for (int jcol = 0; jcol < _ncol; jcol++)
        if (zero(x[jcol], cplexPrecision) != 0)
        {
          if (printL(6))
            std::cout << "primSol[" << jcol << "] = " << x[jcol] << std::endl;

          primSolPool[primSolPool.size()-1][jcol] = x[jcol];
        }

      delete[] x;
      x = NULL;
    }


  primSolPool.push_back(std::map<int, Double>());
    
  getSol(primSolPool[0], ifPrint); // Get the incumbent anyways

  return;
}

void MipCplexInterface::reset()
{
  return;
}

struct CplexCurrentBounds
{
  double startTime;
  double primalBound;
  double dualBound;
  bool print;
};
 
 int mycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
 {
   CplexCurrentBounds * curBoundsPtr = (CplexCurrentBounds *) cbhandle;
   
   double incPB;   
   double incDB;   
   int status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &incPB);
   if ( status ) 
   {
     std::cerr << "Failure to turn on CPXgetcallbackinfo " << std::endl;
     exit(1);
   }
   
   status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &incDB);
   if ( status ) 
   {
     std::cerr << "Failure to turn on CPXgetcallbackinfo " << std::endl;
     exit(1);
   }
   
   if (incPB < curBoundsPtr->primalBound || incDB > curBoundsPtr->dualBound)
   {
     double timenow;
     status = CPXgettime (env, &timenow);
     if ( status ) 
     {
       std::cerr << "Failure to turn on CPXgettime " << std::endl;
       exit(1);
     }
     
     if (curBoundsPtr->print)
        std::cout << "<CPLEX_MIP> <t=" << timenow - curBoundsPtr->startTime  << "> <DB=" << incDB << ">"
                  << "<PB=" << incPB << ">" << std::endl;
     
     curBoundsPtr->primalBound = incPB;
     curBoundsPtr->dualBound = incDB;
   }
   
   return status;
 }
 

  

void MipCplexInterface::optimise(const int & minmaxStatus,
                                 const double & BarrierConvergenceTolerance,
                                 const double & rightHAndSideZeroTol,
                                 const double & reducedCostTolerance,
				                 const bool & preprocessorOn,
				                 const bool & probingOn,
				                 const bool & automaticCuttingPlanesOn,
				                 const char & solverSelection)
{
  bapcodInit().require(_formCurrentlyLoaded, "Form not Currently Loaded/2", ProgStatus::quit,3);
  if (printL(8))
    MPSwrite();

  CPXINT currentParamSetting = 0;

  CPXgetintparam(cplexEnv, CPX_PARAM_PREIND, &currentParamSetting);
  if (currentParamSetting != preprocessorOn)
  {
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_PREIND, preprocessorOn),
                       "Failure to change presolve ");
  }

  CPXgetintparam(cplexEnv, CPX_PARAM_PROBE, &currentParamSetting);
  if ( currentParamSetting != ((probingOn)?0:-1) )
  {
    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_PROBE,((probingOn)?0:-1)),
                       "Failure to change probing ");
  }

    CPXINT currentStartAlgo = 0;

    CPXgetintparam(cplexEnv,  CPX_PARAM_STARTALG, &currentStartAlgo);
    if (printL(5))
        cout << "MIP Start algo = " << currentStartAlgo << endl;

    bapcodInit().check(CPXsetintparam(cplexEnv, CPX_PARAM_THREADS,
                                      bapcodInit().param().MipSolverMultiThread() ),
                       "Failure to set number of threads ");

    //> Aurélien (changed the value from 1 to 2 to disable printing lp while still having the cplex output in the console)
    if (printL(2))
        CPXwriteprob (cplexEnv, cplexProbPtr, "myprob.lp",NULL);


    /// Ruslan : sometimes these parameters are changed during optimization by Cplex (do not know why)
    CPXsetdblparam(cplexEnv, CPX_PARAM_OBJLLIM, -1e+75);
    CPXsetdblparam(cplexEnv, CPX_PARAM_OBJULIM, 1e+75);

    bapcodInit().check(CPXmipopt(cplexEnv, cplexProbPtr),"could not solve MIP");

    return;
}


struct AbortObject
{
     double abortValue;
     volatile int * terminatePtr;
};

#ifndef _MSC_VER
#define __stdcall
#endif

int (__stdcall abortCallback)(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{
    AbortObject * obj = (AbortObject *) cbhandle;

    double incPB;
    int status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &incPB);
    if ( status )
    {
        std::cerr << "Failure to turn on CPXgetcallbackinfo " << std::endl;
        exit(1);
    }
    cout << "ABORT CALLBACK " << incPB << " vs " << obj->abortValue << endl;
    if (incPB == obj->abortValue)
    {
        std::cout<<"Stop due to the abort value"<<endl;
        (*(obj->terminatePtr))=1;
    }
    return status;
}


void MipCplexInterface::setAbortValue(const Double& abortValue)
{
    volatile int terminate_init = 0;
    CPXsetterminate(cplexEnv,&terminate_init);
    AbortObject * object = new AbortObject();
    object->abortValue = abortValue.val();
    object->terminatePtr = &terminate_init;
    cout << "Abort value = " << object->abortValue<<endl;
    bapcodInit().check(CPXsetinfocallbackfunc(cplexEnv, abortCallback, object),
                       "Failure to set abort callback");
}


#endif // CPLEX_FOUND
