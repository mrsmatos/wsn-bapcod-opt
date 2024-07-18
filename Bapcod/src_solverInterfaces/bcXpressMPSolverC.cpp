/**
 *
 * This file bcXpressMPSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _XPRESSMP_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMipSolverDef.hpp"
#include "bcXpressMPSolverC.hpp"
#include "bcPrintC.hpp"
#include "bcSwitches.h"
#include "bcTimeC.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 16;
const char endOfWord = '\0';
const std::string filling = "________________";


/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

using namespace std;
MipXpressMPInit::MipXpressMPInit(BapcodInit* bapcodInit) : MathProgSolverInit(bapcodInit)
{
  //printout("MipXpressMPInit::MipXpressMPInit()");
  /// Xpress
  std::cout << "XPRESS==" << getenv("XPRESS") << std::endl;
  check(XPRSinit(getenv("XPRESS")), "XPRESSMP Initialization failed");
  return;
}

/// New functionalities to add to mip inteface: begin
void MipXpressMPInit::setPresolve(const bool & flag)
{
  return;
}

void MipXpressMPInit::setTimeLimit(const double & seconds)
{
  return;
}

void MipXpressMPInit::setMultiThread(const int & maxNbThread)
{
  return;
}


void MipXpressMPInit::setSearchPriority(const int & flag)
{
  return;
}

void MipXpressMPInit::setRelativeMipGapLimit(const double & relativeGap)
{
  return;
}

void MipXpressMPInit::setWorkingMemorySpace(const double & sizeInMB)
{
  return;
}

/// New functionalities to add to mip inteface: end
MipXpressMPInit::~MipXpressMPInit()
{
  //printout("MipXpressMPInit::~MipXpressMPInit()");

  check(XPRSfree(),"could not free Xp prob memory", ProgStatus::terminate);
  return;
}

LpXpressMPInterface::LpXpressMPInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  MathProgSolverInterface(bapcodInit, ref, name),
  xpressProbPtr(NULL)
{
  //printout("LpXpressMPInterface::LpXpressMPInterface()", 6);

  /// Create problem pointer
  check(XPRScreateprob(&xpressProbPtr), "XPRScreateprob: Failed to create LP");
  check(xpressProbPtr == NULL, "XPRScreateprob: returns NULL pointer");

  /// Problem specific parameter setting
  char *logname = new char[name.size() + 5];

  sprintf(logname, "%s.log", name.c_str());
  logname[name.size()+ 4] = endOfWord;
  check(XPRSsetlogfile(xpressProbPtr, logname), "XPRESS log file Initialization failed");
  delete [] logname; logname = NULL;

  //   check(XPRSsetintcontrol(xpressProbPtr, 
  // 			  XPRS_SOLUTIONFILE , 
  // 			  1 ), 
  // 	"XPRESS setintcontrol failed");

  //   check(XPRSsetintcontrol(xpressProbPtr, 
  // 			  XPRS_PRESOLVE , 
  // 			  0 ), 
  // 	"XPRESS setintcontrol failed");

  //   check(XPRSsetintcontrol(xpressProbPtr, 
  // 			  XPRS_MIPPRESOLVE , 
  // 			  0 ), 
  // 	"XPRESS setintcontrol failed");

  //   check(XPRSsetintcontrol(xpressProbPtr, 
  // 			  XPRS_CUTSTRATEGY , 
  // 			  0 ), 
  // 	"XPRESS setintcontrol failed");

  //   check(XPRSsetintcontrol(xpressProbPtr, 
  // 			  XPRS_HEURSTRATEGY , 
  // 			  0 ),
  // 	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_EXTRAROWS,
			  param.MipSolverNbOfExtraRows ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_EXTRACOLS,
			  param.MipSolverNbOfExtraCols ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_EXTRAMIPENTS ,
			  param.MipSolverNbOfExtraGlobalEntities ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_EXTRAELEMS ,
			  param.MipSolverNbOfExtraMatrixElements ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_MAXNODE ,
			  param.MipSolverMaxBBNodes ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_MAXTIME ,
			  param.MipSolverMaxTime ),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_LPITERLIMIT,
			  param.MipSolverMaxNbLpIterations ),
	"XPRESS setintcontrol failed");

  /**
   *  = 0 use 2 phase LP solution
   *  = 1 use big M combined phase 1 and phase 2
   */
  // seticv( N_IFBIGM, 0);
  return;
}

LpXpressMPInterface::~LpXpressMPInterface()
{
  //printout("LpXpressMPInterface::~LpXpressMPInterface()");

  check(XPRSdestroyprob(xpressProbPtr),
	"xpress destroy problem failed",
	ProgStatus::terminate);
  return;
}

MipXpressMPInterface::MipXpressMPInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  LpXpressMPInterface(BapcodInit* bapcodInit, ref, name)
{
  //printout("MipXpressMPInterface::MipXpressMPInterface()");
  _pureLP = false;

  /// Intensive heuristic search
  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_HEURSTRATEGY ,
			  3),
	"XPRESS setintcontrol failed");

  check(XPRSsetintcontrol(xpressProbPtr,
			  XPRS_NODESELECTION,
			  3),
	"XPRESS setintcontrol failed");
  /**
   * 1= local frst
   * 2= best first 
   * 3 = local depth first search
   * 4 = best first, then local first
   * 5 = pure depth first
   */
  return;
}

void LpXpressMPInterface::loadFormulation(const std::string & name,
					  const int & minmaxStatus,
					  const int & ncol,
					  const int & nrow,
					  const std::map<int, std::string> & mapSeqnb2Cname,
					  const std::map<int, std::string> & mapSeqnb2Rname,
					  const std::set<ProbCoef> & objectRow,
					  const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
					  const std::set<ProbBound> & rhsv,
					  const std::set<ProbBound> & bounds,
					  const std::set<ProbType> & types,
					  const std::set<ProbIntC> & directs,
					  const std::set<ProbSetCoef> & sets)
{
  //printout("LpXpressMPInterface::loadFormulation");
  if (printL(3))
    std::cout << "LpXpressMPInterface::loadFormulation() ncol = " << ncol
	      << " , nrow = " << nrow
	      << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  char *probname = new char[name.size() + 1];
  sprintf(probname, "%s", name.c_str());
  probname[name.size()] = endOfWord;

  /// Assumes sorted colMatrix
  int nnze = colMatrix.size();
  char *rsense = new char[nrow];
  fill_n(rsense, nrow, ' ');
  double *rhs = new double[nrow];
  fill_n(rhs, nrow, 0);
  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      rsense[rvPtr->ref] = rvPtr->sense;
      rhs[rvPtr->ref] = zero(rvPtr->bound);
    }

  double *range = NULL;

  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);
  for (set<ProbCoef>::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  int *matbeg = new int[ncol + 1];
  fill_n(matbeg, ncol + 1, 0);
  int *matcnt = new int[ncol + 1];
  fill_n(matcnt, ncol + 1, 0);
  int *matind = new int[nnze];
  fill_n(matind, nnze, 0);
  double *matval = new double[nnze];
  fill_n(matval, nnze, 0);

  int cnt(0);
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();

  while (mPtr != colMatrix.end())
    {
      int curCol = mPtr->colRef;
      matbeg[curCol] = cnt;
      do
        {
          if (printL(7))
	    std::cout << "colMatrix = " << *mPtr;

	  matind[cnt] = mPtr->rowRef;
          matval[cnt] = zero(mPtr->coef);
          cnt++;
          mPtr++;
        } while ((mPtr != colMatrix.end()) && (mPtr->colRef == curCol));
      matcnt[curCol] = cnt - matbeg[curCol];
    }
  matbeg[ncol] = cnt;

  double *dlb = new double[ncol];
  fill_n(dlb, ncol, 0);
  double *dub = new double[ncol];
  fill_n(dub, ncol, BapcodInfinity);
  for (set<ProbBound>::const_iterator bPtr = bounds.begin();
       bPtr != bounds.end();
       bPtr++)
    {
      if (bPtr->sense == 'U')
	dub[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'L')
	dlb[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'F')
	{
	  dlb[bPtr->ref] = zero(bPtr->bound);
	  dub[bPtr->ref] = zero(bPtr->bound);
	}
    }

  std::string curName;

  char * place;
  char *rnames = new char[placeSize * nrow];
  fill_n(rnames, placeSize * nrow, endOfWord);

  for(int i = 0; i < nrow; i++)
    {
      double *dlb = new double[ncol];
      fill_n(dlb, ncol, 0);
      double *dub = new double[ncol];
      fill_n(dub, ncol, BapcodInfinity);
      for (set<ProbBound>::const_iterator bPtr = bounds.begin();
	   bPtr != bounds.end();
	   bPtr++)
	{
	  // zero(bPtr->bound);
	  if (bPtr->sense == 'U')
	    dub[bPtr->ref] = bPtr->bound;
	  //zero(bPtr->bound);
	  else if (bPtr->sense == 'L')
	    dlb[bPtr->ref] = bPtr->bound;
	  else if (bPtr->sense == 'F')
	    {
	      // zero(bPtr->bound);
	      dlb[bPtr->ref] = bPtr->bound;
	      // zero(bPtr->bound);
	      dub[bPtr->ref] = bPtr->bound;
	    }
	}

      std::string curName;
      char** rownames = new char*[nrow];
      fill_n(rownames, nrow, static_cast<char*>(NULL));
      for(int ir = 0; ir < nrow; ir++)
	{
	  rownames[ir] = new char[placeSize];
	  fill_n(rownames[ir], placeSize, endOfWord);
	  curName = (mapSeqnb2Rname^ir) + filling;
	  strncpy(rownames[ir], curName.c_str(), wordSize);
	}

      char** colnames = new char*[ncol];
      fill_n(colnames, ncol, static_cast<char*>(NULL));
      for(int ic = 0; ic < ncol; ic++)
	{
	  colnames[ic] = new char[placeSize];
	  fill_n(colnames[ic], placeSize, endOfWord);
	  curName = (mapSeqnb2Cname^ic) + filling;
	  strncpy(colnames[ic], curName.c_str(), wordSize);
	}
      place = & (rnames[placeSize * i]);
      // std::cout << "LpXpressMPInterface::loadFormulation(): mapSeqnb2Rname"
      //           << std::endl;
      curName = (mapSeqnb2Rname^i) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }

  if (printL(7))
    {
      for(int i = 0; i < nrow; i++)
	printf("rname %s\n", &rnames[placeSize * i]);
    }

  char *cnames = new char[placeSize * ncol];
  fill_n(cnames, placeSize * ncol, endOfWord);

  for(int i = 0; i < ncol; i++)
    {
      place = & (cnames[placeSize * i]);
      // std::cout << "LpXpressMPInterface::loadFormulation(): mapSeqnb2Cname"
      //           << std::endl;
      curName = (mapSeqnb2Cname^i) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }
  if (printL(7))
    {
      for(int i = 0; i < ncol; i++)
	printf("cname %s\n", &cnames[placeSize * i]);
    }

  check(XPRSloadlp(xpressProbPtr, probname,
		   ncol, nrow,
		   rsense, rhs, range,
		   obj,
		   matbeg, matcnt, matind, matval,
		   dlb, dub ),
        "LpXpressMPInterface::LpXpressMPInterface(): could not loadprob");

  //      int readNrow, readNcol;
  //      XPRSgetintattrib(xpressProbPtr, 
  // 		      XPRS_COLS, 
  // 		      &readNcol);

  //      XPRSgetintattrib(xpressProbPtr, 
  // 		      XPRS_ROWS, 
  // 		      &readNrow);

  //      std::cout << " readNcol = " << readNcol 
  // 	       << "    readNrow = " << readNrow 
  // 	       << std::endl;

  //      check(XPRSwriteprob(xpressProbPtr, "curprob.lp", "l"), 
  // 	   "could not write lp") ;

  //      exit(1);

  check(XPRSaddnames(xpressProbPtr,
		     1,
		     rnames,
		     0,
		     nrow-1),
        "LpXpressMPInterface::LpXpressMPInterface(): could not add row names",
	ProgStatus::terminate);

  check(XPRSaddnames(xpressProbPtr,
		     2,
		     cnames,
		     0,
		     ncol-1),
        "LpXpressMPInterface::LpXpressMPInterface(): could not add col names",
	ProgStatus::terminate);

  delete [] rnames;
  rnames = NULL;
  delete [] cnames;
  cnames = NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] rhs;
  rhs = NULL;
  delete [] range;
  range = NULL;
  delete [] obj;
  obj = NULL;
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
  delete [] probname;
  probname = NULL;
  return;
}

void LpXpressMPInterface::unLoadFormulation()
{
  //printout("LpXpressMPInterface::~unLoadFormulation()", 6);
  return;
  _ncol = 0;
  _nrow = 0;
  //printout("LpXpressMPInterface::unLoadFormulation() : END" );

  return;

}

void MipXpressMPInterface::loadFormulation(const std::string & name,
					   const int & minmaxStatus,
					   const int & ncol,
					   const int & nrow,
					   const std::map<int, std::string> & mapSeqnb2Cname,
					   const std::map<int, std::string> & mapSeqnb2Rname,
					   const std::set<ProbCoef> & objectRow,
					   const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
					   const std::set<ProbBound> & rhsv,
					   const std::set<ProbBound> & bounds,
					   const std::set<ProbType> & types,
					   const std::set<ProbIntC> & directs,
					   const std::set<ProbSetCoef> & sets)
{
  //printout("MipXpressMPInterface::loadFormulation");
  if (printL(5))
    std::cout << "MipXpressMPInterface::loadFormulation() ncol = " << ncol
	      << " , nrow = " << nrow
	      << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();

  const int EXTRA(0);
  _pureLP = false;
  char *probname = new char[name.size() + 1];
  sprintf(probname, "%s", name.c_str());
  probname[name.size()] = endOfWord;
  /// Assumes sorted colMatrix

  int nnze = colMatrix.size();
  char *rsense = new char[nrow + EXTRA];
  fill_n(rsense, nrow + EXTRA, ' ');

  double *rhs = new double[nrow + EXTRA];
  fill_n(rhs, nrow + EXTRA, 0);

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      rsense[rvPtr->ref] = rvPtr->sense;
      rhs[rvPtr->ref] = zero(rvPtr->bound);
    }

  double *range = NULL;
  double *obj = new double[ncol + EXTRA];
  fill_n(obj, ncol + EXTRA, 0);

  for (set<ProbCoef>::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  int *matbeg = new int[ncol + EXTRA];
  fill_n(matbeg, ncol + EXTRA, 0);

  int *matcnt = new int[ncol + EXTRA];
  fill_n(matcnt, ncol + EXTRA, 0);

  int *matind = new int[nnze + EXTRA];
  fill_n(matind, nnze + EXTRA, 0);

  double *matval = new double[nnze + EXTRA];
  fill_n(matval, nnze + EXTRA, 0);

  int cnt(0);
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();

  //   int curCol(0);
  //   for (curCol = 0; 
  //        curCol < ncol;  
  //        matcnt[curCol] = cnt - matbeg[curCol], curCol++)
  //     {
  //       matbeg[curCol] = cnt;
  //       if (printL(7)) 
  // 	std::cout << "ColMatrix = " << *mPtr 
  // 		  << "curCol = " << curCol
  // 		  << " mPtr->colRef = " << mPtr->colRef 
  // 		  << " matbeg[curCol] = " << matbeg[curCol] 
  // 		  << std::endl;

  //       if (mPtr == colMatrix.end()) 
  // 	continue;

  //       if (curCol < (mPtr->colRef)) 
  // 	continue;

  //       while (curCol == mPtr->colRef)
  // 	{
  // 	  matind[cnt] = mPtr->rowRef;
  // 	  matval[cnt] = zero(mPtr->coef);
  // 	  if (printL(7)) 
  // 	    std::cout << "colMatrix = " << *mPtr 
  // 		      << " cnt = " << cnt 
  // 		      << " matind[cnt] = " << matind[cnt] 
  // 		      << " matval[cnt] = " << matval[cnt] 
  // 		      << std::endl;
  // 	  cnt++;
  // 	  mPtr++;
  // 	  if (mPtr == colMatrix.end()) 
  // 	    break;
  // 	}
  //     }
  //   matbeg[curCol] = cnt;

  while (mPtr != colMatrix.end())
    {
      int curCol = mPtr->colRef;
      matbeg[curCol] = cnt;
      do
        {
          matind[cnt] = mPtr->rowRef;
          matval[cnt] = zero(mPtr->coef);
          cnt++;
          mPtr++;
        } while ((mPtr != colMatrix.end()) && (mPtr->colRef == curCol));
      matcnt[curCol] = cnt - matbeg[curCol];
    }

  double *dlb = new double[ncol];
  fill_n(dlb, ncol, 0);
  double *dub = new double[ncol];
  fill_n(dub, ncol, BapcodInfinity);
  for (set<ProbBound>::const_iterator bPtr = bounds.begin();
       bPtr != bounds.end();
       bPtr++)
    {
      if (bPtr->sense == 'U')
	dub[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'L')
	dlb[bPtr->ref] = zero(bPtr->bound);
      else if (bPtr->sense == 'F')
        {
          dlb[bPtr->ref] = zero(bPtr->bound);
          dub[bPtr->ref] = zero(bPtr->bound);
        }
    }

  int ngents = 0;

  char *qgtype = new char[types.size() + EXTRA];
  fill_n(qgtype, types.size() + EXTRA, ' ');

  int *mgcols = new int[types.size() + EXTRA];
  fill_n(mgcols, types.size() + EXTRA, -1);

  for (set<ProbType>::const_iterator typePt = types.begin();
       typePt != types.end();
       typePt++)
    {
      if (typePt->type == 'I')
        {
          qgtype[ngents] = 'I';
          mgcols[ngents] = typePt->ref;
          ngents++;
        }
      if (typePt->type == 'B')
        {
          qgtype[ngents] = 'B';
          mgcols[ngents] = typePt->ref;
          ngents++;
        }
    }

  int nsets = 0;
  double *mplim = NULL;
  char *qstype = NULL;
  int *msstart = NULL;
  int *mscols = NULL;
  double *dref = NULL;
  if (!sets.empty())
    {
      std::set< int > setNbs;
      {
        for(set<ProbSetCoef>::const_iterator sPtr = sets.begin();
	    sPtr != sets.end();
	    sPtr++)
          setNbs.insert(sPtr->setRef);
      }
      nsets = setNbs.size();
      qstype = new char[nsets];
      fill_n(qstype, nsets, ' ');
      msstart = new int[nsets + 1];
      fill_n(msstart, nsets + 1, 0);
      mscols = new int[sets.size()];
      fill_n(mscols, sets.size(), -1);
      dref = new double[sets.size()];
      fill_n(dref, sets.size(), 0);

      {
        int cnt(0);
        set<ProbSetCoef>::const_iterator sPtr = sets.begin();
        while (sPtr != sets.end())
          {
            int curSet = sPtr->setRef;
            qstype[curSet] = sPtr->setType;
            msstart[curSet] = cnt;
            do
              {
                mscols[cnt] = sPtr->colRef;
                dref[cnt] = sPtr->coef;
                cnt++;
                sPtr++;
              } while ((sPtr != sets.end()) && (sPtr->setRef == curSet));
          }
        msstart[nsets] = sets.size();
      }

    }

  char * place;
  std::string curName;
  char *rnames = new char[placeSize * (nrow + EXTRA)];
  fill_n(rnames, placeSize * (nrow + EXTRA), endOfWord);

  for(int ir = 0; ir < nrow; ir++)
    {
      place = & (rnames[placeSize * ir]);
      curName = (mapSeqnb2Rname^ir) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }

  if (printL(7))
    for(int i = 0; i < nrow; i++)
      printf("rname %s\n", &rnames[placeSize * i]);

  char *cnames = new char[placeSize * (ncol + EXTRA)];
  fill_n(cnames, placeSize * (ncol + EXTRA), endOfWord);

  for(int ic = 0; ic < ncol; ic++)
    {
      place = & (cnames[placeSize * ic]);
      curName = (mapSeqnb2Cname^ic) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }

  if (printL(7))
    for(int i = 0; i < ncol; i++)
      printf("cname %s\n", &cnames[placeSize * i]);

  check(XPRSloadglobal(xpressProbPtr, probname,
		       ncol, nrow,
		       rsense, rhs, range,
		       obj,
		       matbeg, matcnt, matind, matval,
		       dlb, dub,
		       ngents, nsets,
		       qgtype,
		       mgcols, mplim,
		       qstype,
		       msstart, mscols,
		       dref),
        "could not loadglobal");

  check(XPRSaddnames(xpressProbPtr,
		     1,
		     rnames,
		     0,
		     nrow-1),
        "could not add row names",
	ProgStatus::terminate);

  check(XPRSaddnames(xpressProbPtr,
		     2,
		     cnames,
		     0,
		     ncol-1),
        "could not add col names",
	ProgStatus::terminate);

  delete [] rnames;
  rnames = NULL;
  delete [] cnames;
  cnames = NULL;
  delete [] qgtype;
  qgtype = NULL;
  delete [] mgcols;
  mgcols = NULL;
  delete [] mplim;
  mplim = NULL;
  delete [] qstype;
  qstype = NULL;
  delete [] msstart;
  msstart = NULL;
  delete [] mscols;
  mscols = NULL;
  delete [] dref;
  dref = NULL;

  int ndir = 0;
  int *mcols = NULL;
  int *mpri = NULL;
  char *qbr = NULL;
  double *dupc = NULL;
  double *ddpc = NULL;

  if(!directs.empty())
    {
      mcols = new int[directs.size()];
      fill_n(mcols, directs.size(), -100);

      mpri = new int[directs.size()];
      fill_n(mpri, directs.size(), 0);

      qbr = new char[directs.size()];
      fill_n(qbr, directs.size(), ' ');

      for (set<ProbIntC>::const_iterator dPtr = directs.begin();
	   dPtr != directs.end();
	   dPtr++)
        {
          //cout << *dPtr;
          mcols[ndir] = dPtr->ref;
          mpri[ndir] = dPtr->val;



          qbr[ndir] = dPtr->type;
          //cout << "dir[" << ndir << "] = " 
	  //     << mcols[ndir] << " " 
	  //     << mpri[ndir] << " " 
	  //     << qbr[ndir] 
	  //     << std::endl;
          ndir++;
        }
    }
  if (ndir)
    check(XPRSloaddirs(xpressProbPtr,
		       ndir,
		       mcols,
		       mpri,
		       qbr,
		       dupc,
		       ddpc),
          "could not loaddir",
	  ProgStatus::terminate);

  delete [] mcols;
  mcols = NULL;
  delete [] mpri;
  mpri = NULL;
  delete [] qbr;
  qbr = NULL;
  delete [] dupc;
  dupc = NULL;
  delete [] ddpc;
  ddpc = NULL;
  delete [] probname;
  probname = NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] rhs;
  rhs = NULL;
  delete [] range;
  range = NULL;
  delete [] obj;
  obj = NULL;
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
  /// We work only on lp with lpsolve, nolt allow but could used with BaB of lpSolve
  return;
}


void LpXpressMPInterface::addCols(const std::set<ProbCoef> & objectiveRow,
				  const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
				  const std::set<ProbBound> & bounds,
				  const std::map<int, std::string> & mapSeqnb2Cname)
{
  //printout("LpXpressMPInterface::addCols()");
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  int readNcol;
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_COLS,
		   &readNcol);

  int newnz_xpressmp = colMatrix.size();
  // require(newnz_xpressmp > 0, "empty extra colMatrix for addCols");

  require(readNcol == _ncol,
	  "LpXpressMPInterface::addCols: readNcol != _ncol");

  double *eobjx = new double[newcol + 1];
  fill_n(eobjx, newcol, 0);
  for (set<ProbCoef>::const_iterator oPtr = objectiveRow.begin();
       oPtr != objectiveRow.end();
       oPtr++)
    eobjx[oPtr->colRef - _ncol] = zero(oPtr->coef);

  int *ematbeg = new int[newcol + 1];
  fill_n(ematbeg, newcol + 1, 0);

  int *ematind = new int[newnz_xpressmp];
  fill_n(ematind, newnz_xpressmp, 0);

  double *ematval = new double[newnz_xpressmp];
  fill_n(ematval, newnz_xpressmp, 0);

  int cnt(0);
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();
  int newColRef(0);
  for (newColRef = 0; newColRef < newcol; newColRef++)
    {
      ematbeg[newColRef] = cnt;
      if (printL(7)) std::cout << "ColMatrix = " << *mPtr
			       << " newColRef = " << newColRef
			       << " mPtr->colRef - _nrow = " << mPtr->colRef - _ncol
			       << " ematbeg[newColRef] = " << ematbeg[newColRef]
			       << std::endl;

      if (mPtr == colMatrix.end())
	continue;

      if (newColRef < (mPtr->colRef - _ncol))
	continue;

      while (newColRef == (mPtr->colRef - _ncol))
        {
          ematind[cnt] = mPtr->rowRef;
          ematval[cnt] = zero(mPtr->coef);
          if (printL(7))
	    std::cout << "colMatrix = " << *mPtr
		      << " cnt = " << cnt
		      << " ematind[cnt] = " << ematind[cnt]
		      << " ematval[cnt] = " << ematval[cnt]
		      << std::endl;

          cnt++;
          mPtr++;
          if (mPtr == colMatrix.end()) break;
        }
    }
  ematbeg[newColRef] = cnt;

  double * ebdl = new double[newcol];
  fill_n(ebdl, newcol, 0);

  double * ebdu = new double[newcol];
  fill_n(ebdu, newcol, BapcodInfinity);

  if (bounds.size() > 0)
    {
      for (set<ProbBound>::const_iterator bPtr = bounds.begin();
	   bPtr != bounds.end();
	   bPtr++)
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

  check(XPRSaddcols(xpressProbPtr,
		    newcol, newnz_xpressmp,
		    eobjx, ematbeg, ematind, ematval, ebdl, ebdu),
	"could not add new cols");

  char * place;

  char * ecnames = new char[placeSize * (newcol + 1)];
  fill_n(ecnames, placeSize * (newcol + 1), endOfWord);

  for(int ic = 0; ic < newcol; ic++)
    {
      place = & (ecnames[placeSize * ic]);
      // std::cout << "LpXpressMPInterface::addCol(): mapSeqnb2Cname" << std::endl;
      curName = (mapSeqnb2Cname^(ic+ _ncol)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
      if (printL(6))
	std::cout << " LpXpressMPInterface::addCols: cname " << place << std::endl;
    }

  check(XPRSaddnames(xpressProbPtr, 2, ecnames,
		     _ncol, _ncol + newcol - 1),
	"could not add new col names",
	ProgStatus::terminate);

  delete[] ecnames;
  ecnames = NULL;
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

void LpXpressMPInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  //printout("LpXpressMPInterface::delCols()");
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  int readNcol;
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_COLS,
		   &readNcol);

  require(readNcol <= _ncol,
	  "LpXpressMPInterface::delCols: readNcol > _ncol");

  require(nbCol2Delete <= readNcol,
	  "LpXpressMPInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);

  int cnt(0);

  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
    dindex[cnt++] = *iPtr;

  check(XPRSdelcols(xpressProbPtr,
		    nbCol2Delete,
		    dindex),
	"could not delete cols");
  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}

void LpXpressMPInterface::addRows( const std::set<ProbBound> & rhsv,
				   const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
				   const std::map<int, std::string> & mapSeqnb2Rname)
{
  //printout("LpXpressMPInterface::addRows()");
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;

  int newnz_xpressmp = rowMatrix.size();
  //require(newnz_xpressmp > 0, "empty extra rowMatrix for addRows");
  int readNrow;
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_ROWS,
		   &readNrow);

  if (printL(7))
    std::cout << "newrows = " << newrows
	      << "  newnz_xpressmp = " << newnz_xpressmp
	      << "  readNrow = " << readNrow
	      << "  _nrow = " << _nrow
	      << std::endl;

  require(readNrow == _nrow,
	  "LpXpressMPInterface::addRows: readNrow != _nrow");

  double *erhs = new double[newrows];
  fill_n(erhs, newrows, 0);

  char *esense = new char[newrows];
  fill_n(esense, newrows, ' ');

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      if (printL(7))
	std::cout << "rhsv = " << *rvPtr;

      esense[rvPtr->ref - _nrow] = rvPtr->sense;
      erhs[rvPtr->ref - _nrow] = zero(rvPtr->bound);
    }

  int *ematbeg = new int[newrows + 1];
  fill_n(ematbeg, newrows + 1, 0);

  int *ematind = new int[newnz_xpressmp];
  fill_n(ematind, newnz_xpressmp, 0);

  double *ematval = new double[newnz_xpressmp];
  fill_n(ematval, newnz_xpressmp, 0);

  int cnt(0);
  set<ProbCoef, ProbCoefRowSmallerThan>::const_iterator mPtr = rowMatrix.begin();
  int newRowRef(0);
  for (newRowRef = 0; newRowRef < newrows; newRowRef++)
    {
      ematbeg[newRowRef] = cnt;
      if (printL(7))
	std::cout << "rowMatrix = " << *mPtr
		  << " newRowRef = " << newRowRef
		  << " mPtr->rowRef - _nrow = " << mPtr->rowRef - _nrow
		  << " ematbeg[newRowRef] = " << ematbeg[newRowRef]
		  << std::endl;

      if (mPtr == rowMatrix.end())
	continue;

      if (newRowRef < (mPtr->rowRef - _nrow))
	continue;

      while (newRowRef == (mPtr->rowRef - _nrow))
        {
          ematind[cnt] = mPtr->colRef;
          ematval[cnt] = zero(mPtr->coef);
          if (printL(7))
	    std::cout << "rowMatrix = " << *mPtr
		      << " cnt = " << cnt
		      << " ematind[cnt] = " << ematind[cnt]
		      << " ematval[cnt] = " << ematval[cnt]
		      << std::endl;
          cnt++;
          mPtr++;
          if (mPtr == rowMatrix.end())
	    break;
        }
    }
  ematbeg[newRowRef] = cnt;

  MPSwrite();
  check(XPRSaddrows(xpressProbPtr,
		    newrows, newnz_xpressmp,
		    esense, erhs, NULL, ematbeg, ematind, ematval),
	"could not add new rows");

  std::string curName;
  char * place;

  char * ernames = new char[placeSize * (newrows + 1)];
  fill_n(ernames, placeSize * (newrows + 1), endOfWord);

  for(int ic = 0; ic < newrows; ic++)
    {
      place = & (ernames[placeSize * ic]);
      // std::cout << "LpXpressMPInterface::addRows(): mapSeqnb2Rname"
      //           << std::endl;

      curName = (mapSeqnb2Rname^(ic+ _nrow)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;

      if (printL(6))
	std::cout << " LpXpressMPInterface::addRows: cname "
		  << place
		  << std::endl;
    }

  check(XPRSaddnames(xpressProbPtr, 1,
		     ernames,
		     _nrow, _nrow + newrows - 1),
	"could not add new col names",
	ProgStatus::terminate);

  delete[] ernames;
  ernames = NULL;
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

void LpXpressMPInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  //printout("LpXpressMPInterface::delRows()");
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  int readNrow;
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_ROWS,
		   &readNrow);

  require(readNrow <= _nrow,
	  "LpXpressMPInterface::delRowss: readNrow > _nrow");
  require(nbRow2Delete <= readNrow,
	  "LpXpressMPInterface::delRows: nbRow2Delete > readNrow");

  int *dindex = new int[nbRow2Delete];
  fill_n(dindex, nbRow2Delete, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin();
       iPtr != indexSetOfRow2Delete.end();
       iPtr++)
    dindex[cnt++] = *iPtr;

  check(XPRSdelrows(xpressProbPtr,
		    nbRow2Delete,
		    dindex),
	"could not delete rows");

  delete[] dindex;
  dindex = NULL;
  _nrow -= nbRow2Delete;

  return;
}

void LpXpressMPInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("LpXpressMPInterface::addDirectives");
  check(1, "LP can not have directives");
  return;
}

/// Not yet for LPsolve

void MipXpressMPInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("MipXpressMPInterface::addDirectives");
  check(1,
	"MipXpressMPInterface::addDirectives cannot by used. Instead one should record directives at the outset");

  // #if (0)
  //   if(newDirects.empty()) return;
  //   int Lndir = 0;

  //   int *Lmcols = new int[newDirects.size()]; 
  //   fill_n(Lmcols, newDirects.size(), 0);

  //   int *Lmpri = new int[newDirects.size()]; 
  //   fill_n(Lmpri, newDirects.size(), 0);

  //   char *Lqbr = new char[newDirects.size()]; 
  //   fill_n(Lqbr, newDirects.size(), ' ');

  //    for (set<ProbIntC>::const_iterator dPtr = newDirects.begin(); 
  // 	dPtr != newDirects.end(); 
  // 	dPtr++)
  //      {
  //        if (printL(6)) 
  // 	 std::cout << *dPtr;

  //        Lmcols[Lndir] = dPtr->ref;
  //        Lmpri[Lndir] = dPtr->val;
  //        Lqbr[Lndir] = dPtr->type;
  //        Lndir++;
  //      }

  //    if (printL(6))
  //      {
  //        std::cout << "Directives\n";
  //        for (int iel = 0; iel < Lndir; iel++)
  // 	 cout <<  "col[" << Lmcols[iel] << "] " 
  // 	      << ", priority = " << Lmpri[iel] 
  // 	      << ", sense = " << Lqbr[iel] 
  // 	      << std::endl;
  //      }

  // #if _WITH_XPRESS
  //    check(iloaddir (Lndir, Lmcols, Lmpri, Lqbr, NULL, NULL),
  // 	 "could not loaddir",  
  // 	 ProgStatus::terminate);
  //  #endif

  //    // delete extra arrays
  //    delete[] Lmcols; 
  //    Lmcols = NULL;

  //    delete[] Lmpri; 
  //    Lmpri = NULL;

  //    delete[] Lqbr; 
  //    Lqbr = NULL;

  //  #endif

  return;
}


void LpXpressMPInterface::getObjCoef(ProbCoef & pc)
{
  //printout("LpXpressMPInterface::getObjCoef()", 6);
  double *coef = new double[1];
  check(XPRSgetobj(xpressProbPtr,
		   coef,
		   pc.colRef,
		   pc.colRef),
	"could not getobj");
  pc.coef = *coef;
  delete [] coef;
  return;
}

void LpXpressMPInterface::chgObjCoef(const ProbCoef & pc)
{
  //printout("LpXpressMPInterface::chgObjCoef()", 6);
  int *ref = new int[1];
  *ref = pc.colRef;
  double *coef = new double[1];
  *coef = pc.coef;
  check(XPRSchgobj(xpressProbPtr,
		   1,
		   ref,
		   coef),
	"could not chgobj");

  delete [] ref;
  delete [] coef;
  return;
}


void LpXpressMPInterface::chgMatCoef(const ProbCoef & pc)
{
  //printout("LpXpressMPInterface::chgMatCoef()", 6);

  check(XPRSchgcoef(xpressProbPtr,
		    pc.rowRef, pc.colRef,
		    double(pc.coef)),
	"could not chgcof");
  return;
}

void LpXpressMPInterface::chgRhs(const ProbBound & pb)
{
  //printout("LpXpressMPInterface::chgRhs()", 6);
  int *index = new int[1];
  index[0] = pb.ref;
  double *rhs = new double[1];
  rhs[0] = zero(pb.bound);
  check(XPRSchgrhs(xpressProbPtr, 1, index, rhs), "could not chgrhs");
  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}


void LpXpressMPInterface::chgRhs(const std::set<ProbBound> & newRhs)
{
  //printout("LpXpressMPInterface::chgRhs()", 6);
  if (newRhs.empty()) return;
  int nrhs(0);
  int *index = new int[newRhs.size()];
  fill_n(index, newRhs.size(), -1);
  double *rhs = new double[newRhs.size()];
  fill_n(rhs, newRhs.size(), 0);

  for (set<ProbBound>::const_iterator bPtr = newRhs.begin();
       bPtr != newRhs.end();
       bPtr++)
    {
      index[nrhs] = bPtr->ref;
      rhs[nrhs] = zero(bPtr->bound);
      nrhs++;
    }

  check(XPRSchgrhs(xpressProbPtr,
		   nrhs,
		   index,
		   rhs),
	"could not chgrhs");

  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}

void LpXpressMPInterface::chgBds(const std::set<ProbBound> & newBounds)
{
  //printout("LpXpressMPInterface::chgBds()");
  if (newBounds.empty())
    return;
  int Lnbnds(0);
  int * Lmindex = new int[newBounds.size()];
  fill_n(Lmindex, newBounds.size(), 0);

  char *Lqbtype = new char[newBounds.size()];
  fill_n(Lqbtype, newBounds.size(), ' ');

  double *Lbnd = new double[newBounds.size()];
  fill_n(Lbnd, newBounds.size(), 0);

  for (set<ProbBound>::const_iterator bPtr = newBounds.begin();
       bPtr != newBounds.end();
       bPtr++)
    {
      if (printL(7))
	std::cout << "bound = " << *bPtr;

      Lmindex[Lnbnds] = bPtr->ref;
      Lbnd[Lnbnds]= zero(bPtr->bound);
      if (bPtr->sense == 'L')
        {
          Lqbtype[Lnbnds] = 'L';
          // dlb[bPtr->ref] = Lbnd[Lnbnds];
        }
      else if ((bPtr->sense == 'U')
	       || (bPtr->sense == 'I')
	       || (bPtr->sense == 'B'))
        {
          Lqbtype[Lnbnds] = 'U';
          // dub[bPtr->ref] = Lbnd[Lnbnds];
        }
      else if (bPtr->sense == 'F')
        {
          Lqbtype[Lnbnds] = 'B';
          // dlb[bPtr->ref] = Lbnd[Lnbnds];
          // dub[bPtr->ref] = Lbnd[Lnbnds];
        }
      Lnbnds++;
    }
  check(XPRSchgbounds(xpressProbPtr,
		      Lnbnds,
		      Lmindex,
		      Lqbtype,
		      Lbnd),
	"could not chgbds");

  delete[] Lmindex;
  Lmindex = NULL;
  delete[] Lqbtype;
  Lqbtype = NULL;
  delete[] Lbnd;
  Lbnd = NULL;
  return;
}

void LpXpressMPInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("LpXpressMPInterface::chgColType()");
  check(1,"LP Form cannot have ColType");

  return;
}

void MipXpressMPInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("MipXpressMPInterface::chgColType()");
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
  check(XPRSchgcoltype(xpressProbPtr,
		       Lnels,
		       Lmindex,
		       Lqctype),
	"could not chgcoltype");

  /// Delete extra arrays
  delete[] Lmindex;
  Lmindex = NULL;
  delete[] Lqctype;
  Lqctype = NULL;
  /// Not yet for LPsolve
  return;
}

void LpXpressMPInterface::getObjVal(Double & objV)
{
  //printout("LpXpressMPInterface::getObjVal");
  double objval;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_LPOBJVAL,
			 &objval),
	"could not get objective value");
  objV = zero(objval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsolve
void MipXpressMPInterface::getObjVal(Double & objV)
{
  //printout("MipXpressMPInterface::getObjVal");
  double objval = 0.0;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_MIPOBJVAL,
			 &objval),
	"could not get objective value");
  objV = zero(objval, param.HIGHPRECISION);

  return;
}

void LpXpressMPInterface::getDualBound(Double & val)
{
  //printout("LpXpressMPInterface::getDualBound");
  double dualval = 0.0;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_LPOBJVAL,
			 &dualval),
	"could not get objective value");
  val = zero(dualval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsole
void MipXpressMPInterface::getDualBound(Double & val)
{
  //printout("MipXpressMPInterface::getDualBound");
  double dualval = 0.0;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_BESTBOUND,
			 &dualval),
	"could not get objective value");
  val = zero(dualval, param.HIGHPRECISION);

  return;
}


void LpXpressMPInterface::getPrimalBound(Double & val)
{
  //printout("LpXpressMPInterface::getPrimalBound");
  double primalval = 0.0;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_LPOBJVAL,
			 &primalval),
	"could not get primal bound");
  val = zero(primalval, param.HIGHPRECISION);
  return;
}


/// Not yet for LPsolve
void MipXpressMPInterface::getPrimalBound(Double & val)
{
  //printout("MipXpressMPInterface::getPrimalBound");
  double primalval = 0.0;
  check(XPRSgetdblattrib(xpressProbPtr,
			 XPRS_MIPOBJVAL,
			 &primalval),
	"could not get primal bound");
  val = zero(primalval, param.HIGHPRECISION);
  return;
}

bool LpXpressMPInterface::getOptimStatus(SolutionStatus & lpStatus,
					 SolutionStatus & mipStatus)
{
  //printout("LpXpressMPInterface::getOptimStatus");
  int status;
  check(XPRSgetintattrib(xpressProbPtr,
			 XPRS_LPSTATUS,
			 &status),
	"could not get LP status");

  if (status == XPRS_LP_OPTIMAL)
    {
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP optimal"
		  << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }
  if (printL(4)) MPSwrite();
  if (status == XPRS_LP_INFEAS)
    {
      lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP infeasible"
		  << std::endl;

      return(false);
    }
  if (status == XPRS_LP_CUTOFF)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP objective worse than cutoff"
		  << std::endl;

      return(false);
    }
  if (status == XPRS_LP_UNFINISHED)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP solving unfinished"
		  << std::endl;

      return(false);
    }
  if (status == XPRS_LP_UNBOUNDED)
    {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP unbounded"
		  << std::endl;
      return(false);
    }
  if (status == XPRS_LP_CUTOFF_IN_DUAL)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpXpressMPInterface::getOptimStatus: LP cutoff in dual"
		  << std::endl;

      return(false);
    }
  std::cout << "LpXpressMPInterface::getOptimStatus: undefined status"
	    << std::endl;

  return(false);
}

/// Not yet for LPsolve
bool MipXpressMPInterface::getOptimStatus( SolutionStatus & lpStatus,
					   SolutionStatus & mipStatus)
{
  //printout("MipXpressMPInterface::getOptimStatus");
  int status;
  check(XPRSgetintattrib(xpressProbPtr,
			 XPRS_MIPSTATUS,
			 &status),
	"could not get LP status");

  if (status == XPRS_MIP_OPTIMAL)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: global search complete and an integer solution has been found"
		  << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }

  /// An integer solution has been found
  if (status == XPRS_MIP_SOLUTION)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: an integer solution has been found"
		  << std::endl;

      /// An integer slution is to be retrieved
      return(true);
    }

  if (printL(4))
    MPSwrite();

  if (status == XPRS_MIP_NOT_LOADED)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: MIP not loaded"
		  << std::endl;

      return(false);
    }

  if (status == XPRS_MIP_LP_NOT_OPTIMAL)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: LP not optimized"
		  << std::endl;

      return(false);
    }

  if (status == XPRS_MIP_LP_OPTIMAL)
    {
      mipStatus = SolutionStatus::DualFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: LP has been optimized"
		  << std::endl;

      return(false);
    }

  if (status == XPRS_MIP_NO_SOL_FOUND)
    {
      mipStatus = SolutionStatus::DualFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: global search incomplete, no integer sol found"
		  << std::endl;

      return(false);
    }

  if (status == XPRS_MIP_INFEAS)
    {
      mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipXpressMPInterface::getOptimStatus: global search complete, no integer sol found"
		  << std::endl;

      return(false);
    }
  std::cout << "MipXpressMPInterface::getOptimStatus: undefined status"
	    << std::endl;

  return(false);
}

void LpXpressMPInterface::getSol(std::map<int, Double> & primSol,
				 const bool & ifPrint)
{
  //printout("LpXpressMPInterface::getSol(primSol)");
  primSol.clear();
  int readNcol = 0;
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_COLS,
		   &readNcol);

  require(readNcol <= _ncol,
	  "LpXpressMPInterface::getSol: readNcol > _ncol");

  double *x = new double[_ncol];
  fill_n(x, _ncol, 0);
  check(XPRSgetsol(xpressProbPtr,
		   x,
		   NULL,
		   NULL,
		   NULL),
	" could not get primal solution");

  if (printL(6))
    std::cout << "readNcol = " << readNcol
	      << "  _ncol = " << _ncol
	      << std::endl;

  for (int jcol = 0; jcol < _ncol; jcol++)
    if (zero(x[jcol], param.HIGHPRECISION) != 0)
      {
	if (printL(6))
	  std::cout << "primSol[" << jcol << "] = " << x[jcol]
		    << std::endl;

	primSol[jcol] = x[jcol];
      }

  delete[] x;
  x = NULL;
  return;
}

void MipXpressMPInterface::getSol(std::map<int, Double> & primSol,
				  const bool & ifPrint)
{
  //printout("MipXpressMPInterface::getSol(primSol)");
  primSol.clear();
  LpXpressMPInterface::getSol(primSol, ifPrint);

  return;
}

void LpXpressMPInterface::getSol(std::map<int, Double> & primSol,
				 std::map<int, Double> & dualSol,
				 const int & minmaxStatus,
				 const bool & ifPrint)
{
  //printout("LpXpressMPInterface::getSol(primSol,dualSol)");
  getSol(primSol, ifPrint);
  dualSol.clear();
  int readNrow(0);
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_ROWS,
		   &readNrow);

  require(readNrow <= _nrow,
	  "LpXpressMPInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[_nrow];
  fill_n(dual, _nrow, 0);

  char *rsense = new char[_nrow];
  fill_n(rsense, _nrow, ' ');

  flipsign = - 1;
  check(XPRSgetsol(xpressProbPtr,
		   NULL,
		   NULL,
		   dual,
		   NULL),
	" could not get dual solution");

  check(XPRSgetrowtype(xpressProbPtr,
		       rsense,
		       0,
		       _nrow-1),
	" could not get row type");

  for (int irow = 0; irow < _nrow; irow++)
    if (zero(dual[irow], param.HIGHPRECISION) != 0)
      {
	if (rsense[irow] == 'L')
	  dualSol[irow] = minmaxStatus * flipsign * dual[irow];
	else if (rsense[irow] == 'G')
	  dualSol[irow] = minmaxStatus * flipsign * dual[irow];
	else
	  ///@internal == 'E'
	  dualSol[irow] = flipsign * dual[irow];

	if (printL(6))
	  {
	    std::cout << "dual[" << irow << "] = " << dual[irow]
		      << std::endl;

	    std::cout << " dualSol[" << irow << "] = " << dualSol[irow]
		      << std::endl;
	  }
      }
  delete[] dual;
  dual = NULL;
  delete[] rsense;
  rsense = NULL;
  return;
}

void LpXpressMPInterface::getReducedCost(std::map<int, Double> & redCost, 
				  const bool & ifPrint)
{

}

//  const bool LpXpressMPInterface::preProcessorChangedForm(int * ptrProbRowCnt,
// 						       const bool & checkVarOnly)
//  {
//    //printout("LpXpressMPInterface::preProcessorFeedBack()");
//    require(*ptrProbRowCnt == _nrow, 
// 	   "LpXpressMPInterface::preProcessorChangedForm(): not (*ptrProbRowCnt == _nrow)");

//    // /// Since matrix is saved and restore in optimise()
//    // return(false);
//    int readNcol;
//  #if _WITH_XPRESS
//    getipv(N_NCOL, &readNcol);
//  #endif
//  #if _WITH_CPLEX
//    readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
//  #endif
//  #if _WITH_LPSOLVE
//    readNcol = get_ncolumns(lpStructPtr);
//  #endif
//    if (_ncol != readNcol)
//      {
//        if (printL(5)) 
// 	 std::cout << "LpXpressMPInterface::preProcessor has Changed Formulation (_ncol != readNcol): need to rebuild formulation"
// 		   << std::endl;

//        return(true);
//      }

//    if (checkVarOnly) 
//      return(false);

//    int readNrow;

//  #if _WITH_XPRESS
//    getipv(N_NROW, &readNrow);
//  #endif
//  #if _WITH_CPLEX
//    readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
//  #endif
//  #if _WITH_LPSOLVE
//    readNrow = get_nrows(lpStructPtr);
//  #endif
//    if (_nrow != readNrow)
//      {
//        if (printL(5)) 
// 	 std::cout << "LpXpressMPInterface::preProcessor has Changed Formulation (_nrow > readNrow): need to rebuild formulation"
// 		   << std::endl;
//        return(true);
//      }

// //    if (_nrow < readNrow)
// //       {
// //         if (printL(5)) 
// // 	  std::cout << "LpXpressMPInterface::preProcessor has Changed Formulation (_nrow < readNrow): need to rebuild formulation"
// // 		    << std::endl;
// //         *ptrProbRowCnt = _nrow = readNrow;
// //         return(false);
// //       }

//    if (printL(5)) 
//      std::cout << "LpXpressMPInterface::preProcessor has NOT changed Formulation"
// 	       << std::endl;

//    return(false);
//  }

void LpXpressMPInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  //printout("LpXpressMPInterface::LPwrite()");
  int irow, jcol, iel, nchar;
  int readNrow, readNcol, idone, readNglents, readNsets;
  char coltype;
  double readBdl, readBdu;
  char probname[64];
  check(XPRSgetprobname(xpressProbPtr, probname ),
	"could not get probname") ;

  XPRSgetintattrib(xpressProbPtr,
		   XPRS_COLS,
		   &readNcol);

  XPRSgetintattrib(xpressProbPtr,
		   XPRS_ROWS,
		   &readNrow);

  int *readMclind = new int[readNcol];
  fill_n(readMclind, readNcol, 0);
  char *readCnames = new char[placeSize * readNcol];
  fill_n(readCnames, placeSize * readNcol, ' ');

  XPRSgetnames(xpressProbPtr,
	       2,
	       readCnames,
	       0,
	       readNcol-1);
  printf("PROBLEM: %s whose formulation ref number is %d\n", probname, ref());

  /// Output Obj
  if ( minmaxStatus ==1 )
    printf("Minimize\n");
  else
    printf("Maximize\n");

  double *readDmatval = new double[readNcol];
  fill_n(readDmatval, readNcol, 0);

  XPRSgetobj(xpressProbPtr,
	     readDmatval,
	     0,
	     readNcol-1);
  nchar = printf(" Objectif:");

  for (iel = 0; iel < readNcol; iel++)
    if (Double(readDmatval[iel]) != 0)
      {
	///@todo Replace printf by cout
	nchar += printf(" %+g %s", readDmatval[iel], &readCnames[placeSize * iel]);
	if (nchar >= 80-9-13)
	  {
	    printf("\n");
	    nchar = 0;
	  }
      }

  /// Output Matrix
  printf("\nSubject To\n");

  for(irow = 0; irow < readNrow; irow++)
    printRow(irow, readNcol, readMclind, readDmatval, readCnames);
  printf("Bounds\n");

  for (jcol = 0; jcol < readNcol; jcol++)
    {
      XPRSgetlb(xpressProbPtr,
		&readBdl,
		jcol,
		jcol);

      XPRSgetub(xpressProbPtr,
		&readBdu,
		jcol,
		jcol);

      if (readBdu < XPRS_PLUSINFINITY/10)
	{
	  printf(" %g <= %s <= %g\n", readBdl, &readCnames[placeSize * jcol], readBdu);
	} else if (readBdl != 0.0)
	{
	  printf("      %s >= %g\n", &readCnames[placeSize * jcol], readBdl);
	}
    }

  /// Always Pure for LPSolve yet
  if (_pureLP)
    printf("End\n");
  else
    {
      XPRSgetglobal(xpressProbPtr,
		    &readNglents,
		    &readNsets,
		    NULL, NULL, NULL, NULL, NULL, NULL, NULL);

      if( readNglents )
        {
          printf("Integers\n");
          for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
            {
              XPRSgetcoltype(xpressProbPtr,
			     &coltype,
			     jcol,
			     jcol);

              if (coltype != 'C')
		{
		  ///@todo Replace printf by cout
		  printf(" %s ", &readCnames[placeSize * jcol]);
		  idone++;
		}
              if (idone >= wordSize)
		{
		  ///@todo Replace printf by cout
		  printf("\n");
		  idone = 0;
		}
            }
          if (idone != 0)
	    printf("\n");
        }
      printf("End\n");


      /// Add directives informations if any
      if( readNglents )
        {
          std::os << "Directives:" << std::endl;

          int Lndir = readNglents + readNsets;
          int *Lmcols = new int[Lndir];
	  fill_n(Lmcols, Lndir, 0);
          int *Lmpri = new int[Lndir];
	  fill_n(Lmpri, Lndir, 0);
          char *Lqbr = new char[Lndir];
	  fill_n(Lqbr, Lndir, ' ');
          char cname[placeSize];
          check(XPRSgetdirs(xpressProbPtr,
			    &Lndir,
			    Lmcols,
			    Lmpri,
			    Lqbr,
			    NULL,
			    NULL),
		"could not getdir");

          for (jcol = 0; jcol < Lndir; jcol++)
            {
              if (Lmcols[jcol] < 0)
		std::os << "set" << - Lmcols[jcol]
			  << " " << Lmpri[jcol]
			  << " " << Lqbr[jcol]
			  << std::endl;
              else
                {
                  XPRSgetnames(xpressProbPtr,
			       2,
			       cname,
			       Lmcols[jcol],
			       Lmcols[jcol]);

                  std::os << cname
			    << " " << Lmpri[jcol]
			    << " " << Lqbr[jcol]
			    << std::endl;
                }
            }
          delete[] Lqbr;
	  Lqbr = NULL;
          delete[] Lmcols;
	  Lmcols = NULL;
          delete[] Lmpri;
	  Lmpri = NULL;
        }

    }

  delete[] readMclind;
  readMclind = NULL;
  delete[] readDmatval;
  readDmatval = NULL;
  delete[] readCnames;
  readCnames = NULL;
  return;
}


void LpXpressMPInterface::MPSwrite()
{
  //printout("LpXpressMPInterface::MPSwrite()");
  check(XPRSwriteprob(xpressProbPtr,
		      "curprob.mps",
		      "p"),
	"could not get probname") ;

  check(XPRSwriteprob(xpressProbPtr,
		      "curprob.lp",
		      "l"),
	"could not get probname") ;
  return;
}


void LpXpressMPInterface::printRow(int irow,
				   const int & readNcol,
				   int *readMclind,
				   double *readDmatval,
				   char *readCnames)
{
  //printout("LpXpressMPInterface::printRow");
  int iel, readNels(0), nchar;
  double readDrhs;
  char readRtype, readSense[3];
  int readMstart[2];
  char rname[placeSize];
  fill_n(readMclind, readNcol, 0);
  fill_n(readDmatval, readNcol, 0);
  check(XPRSgetrows(xpressProbPtr,
		    readMstart,
		    readMclind,
		    readDmatval,
		    readNcol,
		    &readNels,
		    irow,
		    irow),
	"could not getrows");

  check(XPRSgetrhs(xpressProbPtr,
		   &readDrhs,
		   irow,
		   irow),
	"could not getrhs");

  check(XPRSgetrowtype(xpressProbPtr,
		       &readRtype,
		       irow,
		       irow),
	"could not getrowtype");

  check(XPRSgetnames(xpressProbPtr,
		     1,
		     rname,
		     irow,
		     irow),
	"could not getnames");

  strcpy(readSense,
         readRtype=='N'
	 ?"$" :(readRtype=='E'
		?"=" :(readRtype=='L'
		       ?"<=" :">=" ) ));

  nchar = printf(" %s:", rname);
  if (readNels > 0)
    for (iel = 0; iel < readNels; iel++)
      {
        nchar += printf(" %+g %s",
			readDmatval[iel],
			&readCnames[placeSize * readMclind[iel]]);
        if (nchar >= 80-9-13)
	  {
	    printf("\n");
	    nchar = 0;
	  }
      }

  // if (readRtype != 'N')
  printf(" %s %g", readSense, readDrhs);
  printf("\n");
  return;
}

void LpXpressMPInterface::optimise(const int & minmaxStatus,
				   const bool & preprocessorOn,
				   const bool & probingOn,
				   const bool & automaticCuttingPlanesOn,
				   const char & solverSelection)
{
  //printout("LpXpressMPInterface::optimise()", 6);
  require(_formCurrentlyLoaded, "Form not Currently Loaded",
	  ProgStatus::quit,
	  3);
  if (printL(8))
    MPSwrite();



  XPRSsetintcontrol(xpressProbPtr,
		    XPRS_PRESOLVE,
		    (preprocessorOn ? 1 : 0));

  XPRSsetintcontrol(xpressProbPtr,
		    XPRS_OUTPUTLOG,
		    1);

  if (minmaxStatus == 1)
    {
      /// "g" flag means call global
      if (solverSelection == 'p')
	check(XPRSminim(xpressProbPtr, "pl"), "could not solve LP");

      /// "g" flag means call global
      if (solverSelection == 'd')
	check(XPRSminim(xpressProbPtr, "dl"), "could not solve LP");

      /// "g" flag means call global
      if (solverSelection == 'b')
	check(XPRSminim(xpressProbPtr, "bl"), "could not solve LP");
    }
  else
    {
      if (solverSelection == 'p')
	check(XPRSmaxim(xpressProbPtr, "pl"), "could not solve LP");

      if (solverSelection == 'd')
	check(XPRSmaxim(xpressProbPtr, "dl"), "could not solve LP");

      if (solverSelection == 'b')
	check(XPRSmaxim(xpressProbPtr, "bl"), "could not solve LP");
    }
  int iterCounter(0);
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_SIMPLEXITER,
		   &iterCounter);

  if (printL(4))
    std::os << "LpXpressMPInterface::optimise(): iterCounter = " << iterCounter
	      << std::endl;

  if (iterCounter > bapcodInit().statistics().getCounter("bcMaxIterLp"))
    bapcodInit().statistics().setCounter("bcMaxIterLp",iterCounter);
  return;
}

void MipXpressMPInterface::optimise(const int & minmaxStatus,
				    const bool & preprocessorOn,
				    const bool & probingOn,
				    const bool & automaticCuttingPlanesOn,
				    const char & solverSelection)
{
  //printout("MipXpressMPInterface::optimise()", 6);
  require(_formCurrentlyLoaded,
	  "Form not Currently Loaded/2",
	  ProgStatus::quit,
	  3);
  if (printL(8))
    MPSwrite();





  XPRSsetintcontrol(xpressProbPtr,
		    XPRS_MIPPRESOLVE,
		    (preprocessorOn ? 3 : 0) + (probingOn ? 4 : 0));

  XPRSsetintcontrol(xpressProbPtr,
		    XPRS_CUTSTRATEGY,
		    (automaticCuttingPlanesOn ? 2 : 0));

  XPRSsetintcontrol(xpressProbPtr,
		    XPRS_OUTPUTLOG,
		    1);

  if (minmaxStatus == 1)
    {
      if (printL(5))
	std::cout << "MIP minim (solverSelection = " << solverSelection << ")"
		  << std::endl;

      if (solverSelection == 'p')
        {
          int retcod = XPRSminim(xpressProbPtr, "pg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/p " << retcod
		      << std::endl;

	  /// "g" flag means call global
          check(retcod, "could not solve MIP XPRESSMP p");
        }
      if (solverSelection == 'd')
        {
          int retcod = XPRSminim(xpressProbPtr, "dg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/d " << retcod
		      << std::endl;

	  // "g" flag means call global
          check(retcod, "could not solve MIP XPRESSMP d");
        }
      if (solverSelection == 'b')
        {
          int retcod = XPRSminim(xpressProbPtr, "bg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/b " << retcod
		      << std::endl;
	  /// "g" flag means call global
	  check(retcod, "could not solve MIP XPRESSMP b");
        }

      //       if (solverSelection == 'p') 
      // 	/// "g" flag means call global
      // 	check(XPRSminim(xpressProbPtr, "p"), 
      // 	      "could not solve LP relaxation of MIP ");
      //       if (solverSelection == 'd') 
      // 	/// "g" flag means call global
      // 	check(XPRSminim(xpressProbPtr, "d"), 
      // 	      "could not solve LP relaxation of MIP");

      //       if (solverSelection == 'b') 
      // 	/// "g" flag means call global
      // 	check(XPRSminim(xpressProbPtr, "b"), 
      // 	      "could not solve LP relaxation of MIP");

      //       int readNcol = 0;
      //       XPRSgetintattrib(xpressProbPtr, 
      // 		       XPRS_COLS, 
      // 		       &readNcol);

      //       require(readNcol <= _ncol, 
      // 	      "LpXpressMPInterface::getSol: readNcol > _ncol");

      //       double *x = new double[_ncol]; 
      //       fill_n(x, _ncol, 0);
      //       check(XPRSgetsol(xpressProbPtr, 
      // 		       x, 
      // 		       NULL, 
      // 		       NULL, 
      // 		       NULL), 
      // 	    " could not get primal solution");

      //       if (printL(5)) 
      // 	std::cout  << "readNcol = " << readNcol 
      // 		   << "  _ncol = " << _ncol 
      // 		   << std::endl;

      //       for (int jcol = 0; jcol < _ncol; jcol++)
      // 	if (printL(5)) 
      // 	  std::cout  << "MipXpressMPInterface::optimise(): xlp[" << jcol
      // 		     << "] = " << x[jcol] 
      // 		     << std::endl;

      //       check(XPRSglobal(xpressProbPtr), 
      // 	    "could not solve MIP XPRESSMP ?");
      //       check(XPRSgetsol(xpressProbPtr, 
      // 		       x, 
      // 		       NULL, NULL, NULL), 
      // 	    " could not get primal solution");

      //       if (printL(5)) 
      // 	std::cout  << "readNcol = " << readNcol 
      // 		   << "  _ncol = " << _ncol 
      // 		   << std::endl;

      //       for (int jcol = 0; jcol < _ncol; jcol++)
      // 	if (printL(5)) 
      // 	  std::cout  << "MipXpressMPInterface::optimise(): xmip[" << jcol
      // 		     << "] = " << x[jcol] 
      // 		     << std::endl;

      //       delete[] x; x = NULL;

    }
  else
    {
      if (printL(6))
	std::cout << "MIP maxim (solverSelection = " << solverSelection
		  << ")" << std::endl;

      if (solverSelection == 'p')
        {
          int retcod = XPRSmaxim(xpressProbPtr, "pg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/p2 " << retcod
		      << std::endl;

	  /// "g" flag means call global
          check(retcod, "could not solve MIP XPRESSMP p");
        }

      if (solverSelection == 'd')
        {
          int retcod = XPRSmaxim(xpressProbPtr, "dg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/d2 " << retcod << std::endl;

	  /// "g" flag means call global
          check(retcod, "could not solve MIP XPRESSMP d");
        }
      if (solverSelection == 'b')
        {
          int retcod = XPRSmaxim(xpressProbPtr, "bg");
          if (retcod)
	    std::cerr << " ERROR in XPRESSMP/b2 " << retcod << std::endl;

	  /// "g" flag means call global
	  check(retcod, "could not solve MIP XPRESSMP b");
        }
    }
  int bbNodeCounter(0);
  XPRSgetintattrib(xpressProbPtr,
		   XPRS_NODES,
		   &bbNodeCounter);

  if (printL(4))
    std::cout << " MipXpressMPInterface::optimise(): bbNodeCounter = " << bbNodeCounter
	      << std::endl;

  if (bbNodeCounter > bapcodInit().statistics().getCounter("bcMaxBBnodes"))
    bapcodInit().statistics().setCounter("bcMaxBBnodes", bbNodeCounter);
  return;
}

void LpXpressMPInterface::reset()
{
  //printout("LpXpressMPInterface::reset()");
  /// Do nothing

  return;
}



void MipXpressMPInterface::reset()
{
  //printout("MipXpressMPInterface::reset()", 6);
  XPRSinitglobal(xpressProbPtr);
  return;
}


void LpXpressMPInterface::makeSpaceForLoadingForm()
{
  /// Save loaded form if any
  //printout("LpXpressMPInterface::makeSpaceForLoadingForm()", 6);
  _formCurrentlyLoaded = true;
  return;
}

void LpXpressMPInterface::saveCopyOfCurForm()
{
  //printout("LpXpressMPInterface::saveCopyOfCurForm()", 6);
  /// Save loaded form if any
  return;
}

void LpXpressMPInterface::load()
{
  //printout("LpXpressMPInterface::load()", 6);
  return;
}

void LpXpressMPInterface::unload(const bool & deleteMat)
{
  //printout("LpXpressMPInterface::unload()");
  // if (deleteMat == false)

  return;
}

#endif // _XPRESSMP_FOUND
