/**
 *
 * This file bcXpressSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _WITH_XPRESS

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMipSolverDef.hpp"
#include "bcMathProgSolverIntenrfaceC.h"
#include "bcPrintC.hpp"
#include "bcSwitches.h"
#include "bcTimeC.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 12;
const char endOfWord = '\0';

/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;






using namespace std;
MipSolverInit::MipSolverInit()
{
  //printout("MipSolverInit::MipSolverInit()");
  /// Xpress
  check(initlz(getenv("XPRESS"), 0), "XPRESS Initialization failed");


  seticv( N_NRXTRA , param.MipSolverNbOfExtraRows );
  seticv( N_NCXTRA , param.MipSolverNbOfExtraCols );
  seticv( N_NGXTRA , param.MipSolverNbOfExtraGlobalEntities );
  seticv( N_NMXTRA , param.MipSolverNbOfExtraMatrixElements );
  seticv( N_NEXTRA , param.MipSolverNbOfExtraElementsInPresolve );

  /// Turn off printing solution to file
  seticv(N_IFMEM, 13);

  /// Turn on printing solution to file and post-processing of solution
  //seticv( N_IFMEM, 1);
  seticv(N_PRTMSG, 1);
  seticv(N_GPRINT,-100);
  // if(printL(3)) seticv(N_GPRINT, 1);
  // if(printL(4)) seticv(N_GPRINT, 2);
  // if(printL(5)) seticv(N_GPRINT, 3);
  seticv( N_MAXNOD , param.MipSolverMaxBBNodes );

  /**
   * 1= priority to descendent
   * 2= all nodes considered
   * 3 = pure depth first search
   */
  seticv( N_NDSEL1 , 3 );

  setdcv( N_ZTOLZE, param.MipSolverRightHAndSideZeroTol);
  setdcv( N_ZTCOST, param.MipSolverReducedCostTolerance);
  setdcv( N_ZTOLIS, param.MipSolverIntegerFeasibility);
  seticv( N_IFSCAL, 35 );

  /// Set the maxim number of interation in solving an LP
  seticv( N_ITRLIM, param.MipSolverMaxNbLpIterations);

  /**
   * = 0 use 2 phase LP solution
   * = 1 use big M combined phase 1 and phase 2
   */
  // seticv( N_IFBIGM, 0);
  return;
}

/// New functionalities to add to mip inteface: begin
void MipSolverInit::setPresolve(const bool & flag)
{
  return;
}

void MipSolverInit::setTimeLimit(const double & seconds)
{

  return;
}

void MipSolverInit::setMultiThread(const int & maxNbThread)
{

  return;
}


void MipSolverInit::setSearchPriority(const int & flag)
{
  return;
}


void MipSolverInit::setRelativeMipGapLimit(const double & relativeGap)
{






  return;
}

void MipSolverInit::setWorkingMemorySpace(const double & sizeInMB)
{






  return;
}

void MipSolverInit::setLPoptimalityTolerance(const double& tolerance)
{

}

void MipSolverInit::setLPfeasibilityTolerance(const double& tolerance)
{

}



/// New functionalities to add to mip inteface: end
MipSolverInit::~MipSolverInit()
{
  //printout("MipSolverInit::~MipSolverInit()");

  check(freexo(),"could not free Xp prob memory", ProgStatus::terminate);
  return;
}

LpSolverInterface::LpSolverInterface(const int & ref, const std::string & name):
  _ref(ref),
  _formCurrentlyLoaded(false),
  _formCurrentlySaved(false),
  _pureLP(true),
  _ncol(0),
  _nrow(0)

  , matRefNb(-1)
{
  //printout("LpSolverInterface::LpSolverInterface()", 6);
  return;
}

LpSolverInterface::~LpSolverInterface()
{
  //printout("LpSolverInterface::~LpSolverInterface()");
  return;
}

MipSolverInterface::MipSolverInterface(const int & ref, const std::string & name):
  LpSolverInterface(ref, name)
{
  //printout("MipSolverInterface::MipSolverInterface()");
  _pureLP = false;
return;
}

void LpSolverInterface::loadFormulation(const std::string & name,
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
  //printout("LpSolverInterface::loadFormulation");
  if (printL(3))
    std::cout << "LpSolverInterface::loadFormulation() ncol = " << ncol
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
          if (printL(7)) std::cout << "colMatrix = " << *mPtr;
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

  std::string filling = "____________";
  std::string curName;

  char * place;
  char *rnames = new char[placeSize * nrow];
  fill_n(rnames, placeSize * nrow, endOfWord);

  for(int i = 0; i < nrow; i++)
    {
      place = & (rnames[placeSize * i]);
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
      curName = (mapSeqnb2Cname^i) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }
  if (printL(7))
    {
      for(int i = 0; i < ncol; i++)
 printf("cname %s\n", &cnames[placeSize * i]);
    }

  check(loadprob(probname,
   ncol, nrow,
   rsense, rhs, range,
   obj,
   matbeg, matcnt, matind, matval,
   dlb, dub ),
        "LpSolverInterface::LpSolverInterface(): could not loadprob");

  check(addnames(1, rnames, 0, nrow-1),
        "LpSolverInterface::LpSolverInterface(): could not add row names",
 ProgStatus::terminate);

  check(addnames(2, cnames, 0, ncol-1),
        "LpSolverInterface::LpSolverInterface(): could not add col names",
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


void LpSolverInterface::unLoadFormulation()
{
  //printout("LpSolverInterface::~unLoadFormulation()", 6);

  /// If last loaded form delete it
  if (_formCurrentlyLoaded)
    {
      //printout("LpSolverInterface::unLoadFormulation(): _formCurrentlyLoaded at entry");
      require(ref() == refOfFormCurrentlyLoaded,
              "LpSolverInterface::unLoadFormulation(): ref currently laoded means ref() should be refOfFormCurrentlyLoaded");

      PtLastFormLoaded = NULL;
      refOfFormCurrentlyLoaded = -1;

      return;
    }
  if (_formCurrentlySaved)
    {
      /// I.e. if there is a current formulation
      if (PtLastFormLoaded != NULL)
        {
          saveCopyOfCurForm();
          PtLastFormLoaded->_formCurrentlyLoaded = false;
          PtLastFormLoaded = NULL;
          refOfFormCurrentlyLoaded = -1;
        }
      check(resmat(matRefNb),
     "LpSolverInterface::~LpSolverInterface: could not restore mat");

      matRefNb = -1;
      nbSavedFormMats--;
      _formCurrentlySaved = false;
      _formCurrentlyLoaded = false;
      if (printL(6))
 std::cout << "LpSolverInterface::nbSavedFormMats decreased to "
    << nbSavedFormMats
    << std::endl;

      return;
    }
  _ncol = 0;
  _nrow = 0;
  //printout("LpSolverInterface::unLoadFormulation() : END" );

  return;

}

void MipSolverInterface::loadFormulation(const std::string & name,
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
  //printout("MipSolverInterface::loadFormulation");
  if (printL(5))
    std::cout << "MipSolverInterface::loadFormulation() ncol = " << ncol
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

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
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
  //
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

  for (set<ProbType>::const_iterator typePt = types.begin(); typePt != types.end(); typePt++)
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
      std::set<int> setNbs;
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
  std::string filling = "____________";
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

  check(loadglobal (probname,
      ncol, nrow,
      rsense, rhs, range,
      obj,
      matbeg, matcnt, matind, matval,
      dlb, dub,
      ngents, nsets,
                    qgtype,
      mgcols, mplim, qstype,
      msstart, mscols,
      dref),
        "could not loadglobal");

  check(addnames(1, rnames, 0, nrow-1),
        "could not add row names", ProgStatus::terminate);

  check(addnames(2, cnames, 0, ncol-1),
        "could not add col names", ProgStatus::terminate);

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
          // cout << *dPtr;
          mcols[ndir] = dPtr->ref;
          mpri[ndir] = dPtr->val;
          qbr[ndir] = dPtr->type;
          // cout << "dir[" << ndir << "] = " 
   //      << mcols[ndir] << " " 
   //      << mpri[ndir] << " " 
   //      << qbr[ndir] 
   //      << std::endl;
          ndir++;
        }
    }
  /// Lower priority index means chosen first for branching
  if (ndir)
    check(loaddir (ndir,
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


void LpSolverInterface::addCols(const std::set<ProbCoef> & objectiveRow,
    const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
    const std::set<ProbBound> & bounds,
    const std::map<int, std::string> & mapSeqnb2Cname)
{
  //printout("LpSolverInterface::addCols()");
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;


  int readNcol;
  getipv(N_NCOL, &readNcol);
  int newnz_xpress = colMatrix.size();
  // require(newnz_xpress > 0, "empty extra colMatrix for addCols");
  require(readNcol == _ncol,
   "LpSolverInterface::addCols: readNcol != _ncol");

  double *eobjx = new double[newcol + 1];
  fill_n(eobjx, newcol, 0);
  for (set<ProbCoef>::const_iterator oPtr = objectiveRow.begin();
       oPtr != objectiveRow.end();
       oPtr++)
    eobjx[oPtr->colRef - _ncol] = zero(oPtr->coef);

  int *ematbeg = new int[newcol + 1];
  fill_n(ematbeg, newcol + 1, 0);

  int *ematind = new int[newnz_xpress];
  fill_n(ematind, newnz_xpress, 0);

  double *ematval = new double[newnz_xpress];
  fill_n(ematval, newnz_xpress, 0);

  int cnt(0);
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();
  int newColRef(0);
  for (newColRef = 0; newColRef < newcol; newColRef++)
    {
      ematbeg[newColRef] = cnt;
      if (printL(7))
 std::cout << "ColMatrix = " << *mPtr
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
  std::string filling = "____________";
  std::string curName;

  check(addcols(newcol, newnz_xpress,
  eobjx, ematbeg, ematind, ematval, ebdl, ebdu),
 "could not add new cols");

  char * place;

  char * ecnames = new char[placeSize * (newcol + 1)];
  fill_n(ecnames, placeSize * (newcol + 1), endOfWord);

  for(int ic = 0; ic < newcol; ic++)
    {
      place = & (ecnames[placeSize * ic]);
      curName = (mapSeqnb2Cname^(ic+ _ncol)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
      if (printL(6))
 std::cout << " LpSolverInterface::addCols: cname "
    << place
    << std::endl;
    }

  check(addnames(2, ecnames, _ncol, _ncol + newcol - 1),
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

void LpSolverInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  //printout("LpSolverInterface::delCols()");
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;


  int readNcol;
  getipv(N_NCOL, &readNcol);
  require(readNcol <= _ncol,
   "LpSolverInterface::delCols: readNcol > _ncol");

  require(nbCol2Delete <= readNcol,
   "LpSolverInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);
  int cnt(0);

  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
    dindex[cnt++] = *iPtr;

  check(delcols(nbCol2Delete, dindex),
 "could not delete cols");

  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}

void LpSolverInterface::addRows( const std::set<ProbBound> & rhsv,
     const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
     const std::map<int, std::string> & mapSeqnb2Rname)
{
  //printout("LpSolverInterface::addRows()");
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;

  int newnz_xpress = rowMatrix.size();
  // require(newnz_xpress > 0, "empty extra rowMatrix for addRows");







  int readNrow;
  getipv(N_NROW, &readNrow);
  if (printL(7))
    std::cout << "newrows = " << newrows
       << "  newnz_xpress = " << newnz_xpress
       << "  readNrow = " << readNrow
       << "  _nrow = " << _nrow
       << std::endl;

  require(readNrow == _nrow, "LpSolverInterface::addRows: readNrow != _nrow");

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

  int *ematind = new int[newnz_xpress];
  fill_n(ematind, newnz_xpress, 0);

  double *ematval = new double[newnz_xpress];
  fill_n(ematval, newnz_xpress, 0);

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

  std::string filling = "____________";
  std::string curName;
  check(addrows(newrows,
  newnz_xpress,
  esense,
  erhs,
  NULL,
  ematbeg,
  ematind,
  ematval),
 "could not add new rows");

  char * place;

  char * ernames = new char[placeSize * (newrows + 1)];
  fill_n(ernames, placeSize * (newrows + 1), endOfWord);

  for(int ic = 0; ic < newrows; ic++)
    {
      place = & (ernames[placeSize * ic]);
      curName = (mapSeqnb2Rname^(ic+ _nrow)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
      if (printL(6))
 std::cout << " LpSolverInterface::addRows: cname "
    << place
    << std::endl;
    }

  check(addnames(1, ernames, _nrow, _nrow + newrows - 1),
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

void LpSolverInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  //printout("LpSolverInterface::delRows()");
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;


  int readNrow;
  getipv(N_NROW, &readNrow);
  require(readNrow <= _nrow,
   "LpSolverInterface::delRowss: readNrow > _nrow");
  require(nbRow2Delete <= readNrow,
   "LpSolverInterface::delRows: nbRow2Delete > readNrow");

  int *dindex = new int[nbRow2Delete];
  fill_n(dindex, nbRow2Delete, -1);

  int cnt(0);

  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin();
       iPtr != indexSetOfRow2Delete.end();
       iPtr++)
    dindex[cnt++] = *iPtr;

  check(delrows(nbRow2Delete, dindex), "could not delete rows");
  delete[] dindex;
  dindex = NULL;
  _nrow -= nbRow2Delete;

  return;
}

void LpSolverInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("LpSolverInterface::addDirectives");
  check(1, "LP can not have directives");
  return;
}

/// Not yet for LPsolve

void MipSolverInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("MipSolverInterface::addDirectives");
  check(1,
 "MipSolverInterface::addDirectives cannot by used. Instead one should record directives at the outset");

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


void LpSolverInterface::getObjCoef(ProbCoef & pc)
{
  //printout("LpSolverInterface::getObjCoef()", 6);
  double *coef = new double[1];


  check(getobj(coef,
        pc.colRef,
        pc.colRef),
 "could not getobj");
  pc.coef = *coef;
  delete [] coef;
  return;
}

void LpSolverInterface::chgObjCoef(const ProbCoef & pc)
{
  //printout("LpSolverInterface::chgObjCoef()", 6);


  int *ref = new int[1];
  *ref = pc.colRef;
  double *coef = new double[1];
  *coef = pc.coef;
  check(chgobj(1,
        ref,
        coef),
 "could not chgobj");

  delete [] ref;
  delete [] coef;
  return;
}


void LpSolverInterface::chgMatCoef(const ProbCoef & pc)
{
  //printout("LpSolverInterface::chgMatCoef()", 6);

  check(chgcof(pc.rowRef,
        pc.colRef,
        double(pc.coef)),
 "could not chgcof");
  return;
}

void LpSolverInterface::chgRhs(const ProbBound & pb)
{
  //printout("LpSolverInterface::chgRhs()", 6);


  int *index = new int[1];
  index[0] = pb.ref;
  double *rhs = new double[1];
  rhs[0] = zero(pb.bound);
  check(chgrhs(1, index, rhs), "could not chgrhs");

  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}


void LpSolverInterface::chgRhs(const std::set<ProbBound> & newRhs)
{
  //printout("LpSolverInterface::chgRhs()", 6);
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

  check(chgrhs(nrhs,
        index,
        rhs),
 "could not chgrhs");

  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}

void LpSolverInterface::chgBds(const std::set<ProbBound> & newBounds)
{
  //printout("LpSolverInterface::chgBds()");
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
  check(chgbds(Lnbnds,
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

void LpSolverInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("LpSolverInterface::chgColType()");
  check(1,"LP Form cannot have ColType");

  return;
}

void MipSolverInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("MipSolverInterface::chgColType()");


  int Lnels = 0;
  int *Lmindex = new int[newTypes.size()];
  fill_n(Lmindex, newTypes.size(), 0);

  char *Lqctype = new char[newTypes.size()];
  fill_n(Lqctype, newTypes.size(), ' ');

  for (set<ProbType>::const_iterator bPtr = newTypes.begin();
       bPtr != newTypes.end();
       bPtr++)
    {
      /// Set data structure for chgcoltype
      Lmindex[Lnels] = bPtr->ref;
      Lqctype[Lnels] = bPtr->type;
      Lnels++;
    }
  check(chgcoltype(Lnels,
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

void LpSolverInterface::getObjVal(Double & objV)
{
  //printout("LpSolverInterface::getObjVal");
  double objval;


  check(getdpv(N_DOBJVL, &objval), "could not get objective value");
  objV = zero(objval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsolve
void MipSolverInterface::getObjVal(Double & objV)
{
  //printout("MipSolverInterface::getObjVal");
  double objval = 0.0;


  check(getdpv(N_DINTSL, &objval), "could not get objective value");
  objV = zero(objval, param.HIGHPRECISION);

  return;
}

void LpSolverInterface::getDualBound(Double & val)
{
  //printout("LpSolverInterface::getDualBound");
  double dualval = 0.0;


  check(getdpv(N_DBDOBJ, &dualval), "could not get dual bound");
  val = zero(dualval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsole
void MipSolverInterface::getDualBound(Double & val)
{
  //printout("MipSolverInterface::getDualBound");
  double dualval = 0.0;

  check(getdpv(N_DBESTB, &dualval), "could not get dual bound");
  val = zero(dualval, param.HIGHPRECISION);

  return;
}


void LpSolverInterface::getPrimalBound(Double & val)
{
  //printout("LpSolverInterface::getPrimalBound");
  double primalval = 0.0;


  check(getdpv(N_DBPOBJ, &primalval), "could not get primal bound");
  val = zero(primalval, param.HIGHPRECISION);
  return;
}


/// Not yet for LPsolve
void MipSolverInterface::getPrimalBound(Double & val)
{
  //printout("MipSolverInterface::getPrimalBound");
  double primalval = 0.0;


  check(getdpv(N_DINTSL, &primalval), "could not get primal bound");
  val = zero(primalval, param.HIGHPRECISION);
  return;
}

bool LpSolverInterface::getOptimStatus(SolutionStatus & lpStatus,
           SolutionStatus & mipStatus)
{
  //printout("LpSolverInterface::getOptimStatus");
  int status;


  ///@todo Add a switch for each "if" here.
  check(getipv(N_STATUS, &status), "could not get LP status");
  if (status == 1)
    {
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP optimal" << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }
  if (printL(4)) MPSwrite();
  if (status == 2)
    {
      lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP infeasible"
    << std::endl;

      return(false);
    }
  if (status == 3)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP objective worse than cutoff"
    << std::endl;

      return(false);
    }
  if (status == 4)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP solving unfinished"
    << std::endl;

      return(false);
    }
  if (status == 5)
    {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP unbounded"
    << std::endl;

      return(false);
    }
  if (status == 6)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP cutoff in dual"
    << std::endl;

      return(false);
    }
  std::cout << "LpSolverInterface::getOptimStatus: undefined status"
     << std::endl;

  return(false);
}

/// Not yet for LPsolve
bool MipSolverInterface::getOptimStatus( SolutionStatus & lpStatus,
      SolutionStatus & mipStatus)
{
  //printout("MipSolverInterface::getOptimStatus");
  int status;


  check(getipv(N_GLSTAT, &status), "could not get MIP status");

  if (status == 6)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: global search complete and an integer solution has been found"
    << std::endl;
      /// Optimal solution found => no problem encountered
      return(true);
    }

  /// An integer solution has been found
  if (status == 4)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3)) std::cout << "MipSolverInterface::getOptimStatus: an integer solution has been found"
          << std::endl;

      ///  An integer slution is to be retrieved
      return(true);
    }

  if (printL(4))
    MPSwrite();

  if (status == 0)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP not loaded"
    << std::endl;

      return(false);
    }

  if (status == 1)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: LP not optimized"
    << std::endl;

      return(false);
    }

  if (status == 2)
    {
      mipStatus = SolutionStatus::DualFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: LP has been optimized"
    << std::endl;

      return(false);
    }

  if (status == 3)
    {
      mipStatus = SolutionStatus::DualFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: global search incomplete, no integer sol found"
    << std::endl;

      return(false);
    }

  if (status == 5)
    {
      mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: global search complete, no integer sol found"
    << std::endl;

      return(false);
    }
  std::cout << "MipSolverInterface::getOptimStatus: undefined status"
     << std::endl;

  return(false);
}



void LpSolverInterface::getSol(std::map<int, Double> & primSol,
          const bool & ifPrint)
{
  //printout("LpSolverInterface::getSol(primSol)");
  primSol.clear();


  int readNcol = 0;
  getipv(N_NCOL, &readNcol);
  require(readNcol <= _ncol,
   "LpSolverInterface::getSol: readNcol > _ncol");

  double *x = new double[_ncol];
  fill_n(x, _ncol, 0);
  check(solution(x,
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

void MipSolverInterface::getSol(std::map<int, Double> & primSol,
    const bool & ifPrint)
{
  //printout("MipSolverInterface::getSol(primSol)");
  primSol.clear();


  LpSolverInterface::getSol(primSol, ifPrint);
  return;
}

void LpSolverInterface::getSol(std::map<int, Double> & primSol,
    std::map<int, Double> & dualSol,
    const int & minmaxStatus,
    const bool & ifPrint)
{
  //printout("LpSolverInterface::getSol(primSol,dualSol)");
  getSol(primSol, ifPrint);
  dualSol.clear();


  int readNrow(0);
  getipv(N_NROW, &readNrow);
  require(readNrow <= _nrow,
   "LpSolverInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[_nrow];
  fill_n(dual, _nrow, 0);

  char *rsense = new char[_nrow];
  fill_n(rsense, _nrow, ' ');

  flipsign = - 1;
  check(solution(NULL,
   NULL,
   dual,
   NULL),
 " could not get dual solution");

  check(getrowtype(rsense,
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

//  const bool LpSolverInterface::preProcessorChangedForm(int * ptrProbRowCnt, 
// 						       const bool & checkVarOnly)
//  {
//    //printout("LpSolverInterface::preProcessorFeedBack()");
//    require(*ptrProbRowCnt == _nrow, 
// 	   "LpSolverInterface::preProcessorChangedForm(): not (*ptrProbRowCnt == _nrow)");

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
// 	 std::cout << "LpSolverInterface::preProcessor has Changed Formulation (_ncol != readNcol): need to rebuild formulation" 
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
// 	 std::cout << "LpSolverInterface::preProcessor has Changed Formulation (_nrow > readNrow): need to rebuild formulation" 
// 		   << std::endl;
//        return(true);
//      }

// //    if (_nrow < readNrow)
// //       {
// //         if (printL(5)) 
// // 	  std::cout << "LpSolverInterface::preProcessor has Changed Formulation (_nrow < readNrow): need to rebuild formulation" 
// // 		    << std::endl;
// //         *ptrProbRowCnt = _nrow = readNrow;
// //         return(false);
// //       }

//    if (printL(5)) 
//      std::cout << "LpSolverInterface::preProcessor has NOT changed Formulation" 
// 	       << std::endl;

//    return(false);
//  }

void LpSolverInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  //printout("LpSolverInterface::LPwrite()");


  int irow, jcol, iel, nchar;
  int readNrow, readNcol, idone, readNglents, readNsets;
  char coltype;
  double readBdl, readBdu;
  char probname[64];

  /// Turn scaling off (less work)
  // seticv(N_IFSCAL, 0);
  check(getprob ( probname ), "could not get probname") ;
  getipv(N_NROW, &readNrow);
  getipv(N_NCOL, &readNcol);
  int *readMclind = new int[readNcol];
  fill_n(readMclind, readNcol, 0);
  char *readCnames = new char[placeSize * readNcol];
  fill_n(readCnames, placeSize * readNcol, ' ');
  getnames(2, readCnames, 0, readNcol - 1);
  printf("PROBLEM: %s whose formulation ref number is %d\n", probname, ref());

  /// Output Obj
  if ( minmaxStatus ==1 )
    printf("Minimize\n");
  else
    printf("Maximize\n");

  double *readDmatval = new double[readNcol];
  fill_n(readDmatval, readNcol, 0);
  getobj(readDmatval,
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
      getbdl(&readBdl, jcol, jcol);
      getbdu(&readBdu, jcol, jcol);
      if (readBdu < DPLINF/10)
 {
   ///@todo Replace printf by cout
   printf(" %g <= %s <= %g\n", readBdl, &readCnames[placeSize * jcol], readBdu);
 } else if (readBdl != 0.0)
 ///@todo Replace printf by cout
 printf("      %s >= %g\n", &readCnames[placeSize * jcol], readBdl);
    }

  /// Always pure LP for LPsolve yet
  if (_pureLP)
    printf("End\n");
  else
    {
      getglobal(&readNglents,
  &readNsets,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL);

      if( readNglents )
        {
          printf("Integers\n");
          for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
            {
              getcoltype(&coltype,
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
          check(getdir (&Lndir,
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
                  getnames(2,
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


void LpSolverInterface::MPSwrite()
{
  //printout("LpSolverInterface::MPSwrite()");


  char probname[64];
  check(getprob(probname), "could not get probname") ;





  check(output(probname, "p"),
 "could not output problem",
 ProgStatus::terminate);
  return;
}


void LpSolverInterface::printRow(int irow,
     const int & readNcol,
     int *readMclind,
     double *readDmatval,
     char *readCnames)
{
  //printout("LpSolverInterface::printRow");


  int iel, readNels(0), nchar;
  double readDrhs;
  char readRtype, readSense[3];
  int readMstart[2];
  char rname[placeSize];
  fill_n(readMclind, readNcol, 0);
  fill_n(readDmatval, readNcol, 0);
  readNels = 0;
  check(getrows(readMstart,
  readMclind,
  readDmatval,
  readNcol,
  &readNels,
  irow,
  irow),
 "could not getrows");

  check(getrhs(&readDrhs,
        irow,
        irow),
 "could not getrhs");

  check(getrowtype(&readRtype,
     irow,
     irow),
 "could not getrowtype");

  check(getnames(1,
   rname,
   irow,
   irow),
 "could not getnames");

  ///@todo Replace this type of expression by a more simple
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
        if (nchar >= 80-9-13) {
          printf("\n");
          nchar = 0;
        }
      }
  // if (readRtype != 'N')
  printf(" %s %g", readSense, readDrhs);
  printf("\n");
  return;
}

void LpSolverInterface::optimise(const int & minmaxStatus,
      const bool & preprocessorOn,
      const bool & probingOn,
      const bool & automaticCuttingPlanesOn,
                                  const char & solverSelection)
{
  //printout("LpSolverInterface::optimise()", 6);
  require(_formCurrentlyLoaded, "Form not Currently Loaded",
   ProgStatus::quit,
   3);
  if (printL(8))
    MPSwrite();


  /// preprocessorOn = 0 => Turn off preprocessor: lp/ip preprocessing
  seticv(N_IFPRES, (preprocessorOn ? 1 : 0));
  /// probingOn = 0 => Turn off preprocessor: probing at root node
  seticv(N_IFINTP, (probingOn ? 1 : 0));


  seticv(N_CUTSTRAT, (automaticCuttingPlanesOn ? 1 : 0));


  seticv(N_GPRINT,-100);



  if (minmaxStatus ==1)
    {
      /// "g" flag means call global
      if (solverSelection == 'p')
 check(minim("pl"), "could not solve LP");

      /// "g" flag means call global
      if (solverSelection == 'd')
 check(minim("dl"), "could not solve LP");

      /// "g" flag means call global
      if (solverSelection == 'b')
 check(minim("bl"), "could not solve LP");
    }
  else
    {
      if (solverSelection == 'p')
 check(maxim("pl"), "could not solve LP");

      if (solverSelection == 'd')
 check(maxim("dl"), "could not solve LP");

      if (solverSelection == 'b')
 check(maxim("bl"), "could not solve LP");
    }
  int iterCounter(0);
  getipv(N_ITCNT, &iterCounter);
  if (printL(4))
    std::cout << "LpSolverInterface::optimise(): iterCounter = "
       << iterCounter
       << std::endl;

  if (iterCounter > bapcodInit().statistics().getCounter("bcMaxIterLp"))
    bapcodInit().statistics().setCounter("bcMaxIterLp",iterCounter);
  return;
}

void MipSolverInterface::optimise(const int & minmaxStatus,
      const bool & preprocessorOn,
      const bool & probingOn,
      const bool & automaticCuttingPlanesOn,
                                  const char & solverSelection)
{
  //printout("MipSolverInterface::optimise()", 6);
  require(_formCurrentlyLoaded,
   "Form not Currently Loaded/2",
   ProgStatus::quit,
   3);
  if (printL(8))
    MPSwrite();


  // preprocessorOn = 0 => Turn off preprocessor: lp/ip preprocessing
  seticv(N_IFPRES, (preprocessorOn ? 1 : 0));
  // probingOn = 0 => Turn off preprocessor: probing at root node
  seticv(N_IFINTP, (probingOn ? 7 : 0));


  seticv(N_CUTSTRAT, (automaticCuttingPlanesOn ? - 1 : 0));


  seticv(N_GPRINT, -100);
  // if(printL(3)) 
  //   seticv(N_GPRINT, 1);
  // if(printL(4)) 
  //   seticv(N_GPRINT, 2);
  // if(printL(5)) 
  //   seticv(N_GPRINT, 3);



  if (minmaxStatus ==1)
    {
      if (solverSelection == 'p')
 /// "g" flag means call global
 check(minim("pg"), "could not solve MIP XPRESS p");

      if (solverSelection == 'd')
 /// "g" flag means call global
 check(minim("dg"), "could not solve MIP XPRESS d");

      if (solverSelection == 'b')
 // "g" flag means call global
 check(minim("bg"), "could not solve MIP XPRESS b");
    }
  else
    {
      if (solverSelection == 'p')
 check(maxim("pg"), "could not solve MIP XPRESS p2");

      if (solverSelection == 'd')
 check(maxim("dg"), "could not solve MIP XPRESS d2");

      if (solverSelection == 'b')
 check(maxim("bg"), "could not solve MIP XPRESS b2");
    }
  int bbNodeCounter(0);
  getipv(N_NODNUM, &bbNodeCounter);
  if (printL(4))
    std::cout << " MipSolverInterface::optimise(): bbNodeCounter = " << bbNodeCounter
       << std::endl;

  if (bbNodeCounter > bapcodInit().statistics().getCounter("bcMaxBBnodes"))
    bapcodInit().statistics().setCounter("bcMaxBBnodes",bbNodeCounter);
  return;
}

void LpSolverInterface::reset()
{
  //printout("LpSolverInterface::reset()");
  /// Do nothing

  return;
}



void MipSolverInterface::reset()
{
  //printout("MipSolverInterface::reset()", 6);


  require(_formCurrentlyLoaded,
   "MipSolverInterface::reset(): Form not Currently Loaded",
   ProgStatus::quit,
   3);
  iniglobal();
  return;
}


void LpSolverInterface::makeSpaceForLoadingForm()
{
  /// Save loaded form if any
  //printout("LpSolverInterface::makeSpaceForLoadingForm()", 6);


  /// I.e. if there is a current formulation
  if (PtLastFormLoaded != NULL)
    {
      saveCopyOfCurForm();
      PtLastFormLoaded->_formCurrentlyLoaded = false;
    }

  /// Set status of current formulation
  _formCurrentlyLoaded = true;
  refOfFormCurrentlyLoaded = ref();
  PtLastFormLoaded = this;
  return;
}

void LpSolverInterface::saveCopyOfCurForm()
{
  //printout("LpSolverInterface::saveCopyOfCurForm()", 6);
  /// Save loaded form if any

  require(refOfFormCurrentlyLoaded != -1,
   "LpSolverInterface::saveCopyOfCurForm(): cannot find LastFormLoaded");

  require(PtLastFormLoaded != NULL
   , "LpSolverInterface::saveCopyOfCurForm(): cannot find LastFormLoaded");
  require(PtLastFormLoaded->_formCurrentlyLoaded
   , "LpSolverInterface::saveCopyOfCurForm(): Last Form should be Currently Loaded");

  if (PtLastFormLoaded->_formCurrentlySaved)
    {
      /**
       * current form has already been saved
       * delete saved copy
       */
      check(1,
     "LpSolverInterface::makeSpaceForLoadingForm(): Last Form should not be Currently saved");
    }

  if (printL(6))
    std::cout << "LpSolverInterface::saveCopyOfCurForm() save problem " << PtLastFormLoaded->ref()
       << std::endl;

  check(savmat(&(PtLastFormLoaded->matRefNb)),
 "LpSolverInterface::saveCopyOfCurForm(): could not save mat");

  if (printL(7))
    std::cout << " after savmat PtLastFormLoaded->matRefNb = " << PtLastFormLoaded->matRefNb
       << std::endl;

  nbSavedFormMats++;

  if (printL(6))
    std::cout << "LpSolverInterface::saveCopyOfCurForm(): nbSavedFormMats inceased to "
       << nbSavedFormMats
       << std::endl;

  check(PtLastFormLoaded->matRefNb == -1,
 "LpSolverInterface::saveCopyOfCurForm(): matRefNb still = -1");

  PtLastFormLoaded->_formCurrentlySaved = true;
  std::cout << "LpSolverInterface::saveCopyOfCurForm() : END" << endl;
  return;
}

void LpSolverInterface::load()
{
  //printout("LpSolverInterface::load()", 6);


  if (printL(6))
    std::cout << "LpSolverInterface::reload() problem " << ref()
       << " _formCurrentlyLoaded = " << _formCurrentlyLoaded
       << " _formCurrentlySaved = " << _formCurrentlySaved
       << std::endl;

  if (_formCurrentlyLoaded)
    {
      //printout("LpSolverInterface::reload() _formCurrentlyLoaded at entry");
      require(ref() == refOfFormCurrentlyLoaded,
       "LpSolverInterface::reload(): ref currently laoded means ref() should be refOfFormCurrentlyLoaded");
      return;
    }
  require(_formCurrentlySaved,
   "LpSolverInterface::reload(): form should be either loaded or saved (since object are loaded at creation)");

  makeSpaceForLoadingForm();
  if (printL(7))
    std::cout << " before resmat matRefNb = " << matRefNb << std::endl;

  check(resmat(matRefNb),
 "LpSolverInterface::reload(): could not restore matrix");

  // check(delmat(matRefNb), "LpSolverInterface::reload(): could not delete saved matrix");
  matRefNb = -1;
  nbSavedFormMats--;

  if (printL(6))
    std::cout << "LpSolverInterface::load(): LpSolverInterface::nbSavedFormMats decreased to "
       << nbSavedFormMats
       << std::endl;

  _formCurrentlySaved = false;
  return;
}

void LpSolverInterface::unload(const bool & deleteMat)
{
  //printout("LpSolverInterface::unload()");
  // if (deleteMat == false)

  return;
}
#endif // _WITH_XPRESS
