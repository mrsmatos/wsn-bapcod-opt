/**
 *
 * This file bcLpSolveSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _WITH_LPSOLVE

int MathProgSolverInterface::nbSavedFormMats = 0;
int MathProgSolverInterface::refOfFormCurrentlyLoaded = -1;
MathProgSolverInterface * MathProgSolverInterface::PtLastFormLoaded = NULL;
long MathProgSolverInterface::WhereICanCopyFormMat = 0;

/**
 * Vector that contains the lprec, that we want to keep.
 * The first element is a temporary copie.
 */
std::vector<lprec * > LPform::lprecvect(3);

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMipSolverDef.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcPrintC.hpp"
#include "bcSwitches.h"
#include "bcTimeC.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 12;
const char endOfWord = '\0';

/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

CPXENVptr cplexEnv = NULL;
CPXFILEptr log_cplex = NULL;

using namespace std;
MipSolverInit::MipSolverInit()
{
  //printout("MipSolverInit::MipSolverInit()");
  /// Xpress
  int status(0);
  cplexEnv = CPXopenCPLEX(&status);
  if ( cplexEnv == NULL || status != 0 )
    {
      std::cout << "CPXopenCPLEX status = " << status <<endl;
      std::cout << "Could not open CPLEX environment.\n";
      char errmsg[1024];
      CPXgeterrorstring (cplexEnv,
    status,
    errmsg);
      std::cout << errmsg << std::endl;
      exit (1);
    }

  log_cplex = CPXfopen ("cplex.log", "w");

  status = CPXsetlogfile (cplexEnv, log_cplex);
  if ( status != 0 )
    {
      std::cout << "CPXsetlogfile status = " << status <<endl;
      std::cout << "Could not open CPLEX logfile.\n";
      char errmsg[1024];
      CPXgeterrorstring (cplexEnv,
    status,
    errmsg);
      std::cout << errmsg << std::endl;
      exit (1);
    }

  if (printL(1))
    {
      std::cout << "CPXopenCPLEXdevelop status = " << status <<endl;
      std::cout << "CPLEX version is " << CPXversion(cplexEnv) <<endl;
      check(CPXsetintparam (cplexEnv,
       CPX_PARAM_SCRIND,
       CPX_ON),
     "Failure to turn on screen indicator");
    }
  else
    {
      check(CPXsetintparam (cplexEnv,
       CPX_PARAM_SCRIND,
       CPX_OFF),
     "Failure to turn off screen indicator");
    }

  // /// Turn on/off preprocessor: lp/ip preprocessing
  //   check(CPXsetintparam (cplexEnv, 
  // 			CPX_PARAM_PREIND, 
  // 			CPX_ON), 
  // 	"Failure to turn on/off presolve ");

  //   check(CPXsetdblparam (cplexEnv, 
  // 			CPX_PARAM_TILIM, 
  // 			seconds), 
  // 	"Failure to set time limit ");

  //   check(CPXsetintparam (cplexEnv,
  // 			CPX_PARAM_MIPEMPHASIS, 
  // 			CPX_MIPEMPHASIS_FEASIBILITY), 
  // 	"Failure to set mipemphasis indicator ");

  //   check(CPXsetdblparam (cplexEnv, 
  // 			CPX_PARAM_WORKMEM, 
  // 			5000), 
  // 	"Failure to reset work memory ");

  //     check(CPXsetdblparam (cplexEnv, 
  // 			  CPX_PARAM_EPRHS, 
  // 			  1e-4), 
  // 	  "Failure to set LP feasibility tolerance ");

  //     check(CPXsetdblparam (cplexEnv, 
  // 			  CPX_PARAM_EPINT, 
  // 			  1e-4), 
  // 	  "Failure to set IP feasibility tolerance ");
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


  if ( log_cplex != NULL )
    CPXfclose (log_cplex);

  check(CPXcloseCPLEX(&cplexEnv),
 "could not close cplex and free memory",
 ProgStatus::terminate);







  return;
}

LpSolverInterface::LpSolverInterface(const int & ref, const std::string & name):
  _ref(ref),
  _formCurrentlyLoaded(false),
  _formCurrentlySaved(false),
  _pureLP(true),
  _ncol(0),
  _nrow(0)
  , cplexProbPtr(NULL)
{
  //printout("LpSolverInterface::LpSolverInterface()", 6);
  char *probname = new char[name.size() + 1];
  sprintf(probname, "%s", name.c_str());
  probname[name.size()] = endOfWord;
  int status(0);
  cplexProbPtr = CPXcreateprob(cplexEnv,
          &status,
          probname);
  check(cplexProbPtr == NULL,
 "CPXcreateprob: Failed to create LP");

  check(CPXchgprobtype(cplexEnv,
         cplexProbPtr,
         CPXPROB_LP),
 " CPXchgprobtype: Failed to change typ\e to LP.");

  delete [] probname; probname = NULL;
  return;
}

LpSolverInterface::~LpSolverInterface()
{
  //printout("LpSolverInterface::~LpSolverInterface()");

  check(CPXfreeprob(cplexEnv, &cplexProbPtr),
 "CPXfreeprob failed",
 ProgStatus::terminate);

  return;
}

MipSolverInterface::MipSolverInterface(const int & ref, const std::string & name):
  LpSolverInterface(ref, name)
{
  //printout("MipSolverInterface::MipSolverInterface()");
  _pureLP = false;
  check(cplexProbPtr == NULL,
 "MipSolverInterface: problem must have been  created");

  check(CPXchgprobtype(cplexEnv,
         cplexProbPtr,
         CPXPROB_MILP),
 " CPXchgprobtype: Failed to change typ\e to MILP.");
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
  /// Assumes sorted colMatrix
  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);

  for (set<ProbCoef>::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  double *rhs = new double[nrow];
  fill_n(rhs, nrow, 0);
  char *rsense = new char[nrow];
  fill_n(rsense, nrow, ' ');
  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
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
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();
  int curCol(0);
  for(;curCol<ncol;curCol++)
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
          dub[bPtr->ref] = dlb[bPtr->ref];
        }
    }

  double *range = NULL;
  std::string filling = "____________";
  std::string curName;
  char *colnamestore = new char[ncol * placeSize];
  fill_n(colnamestore, ncol * placeSize, endOfWord);
  char **colnames = new char*[ncol];
  fill_n(colnames, ncol, static_cast<char*>(NULL));
  for(int ic = 0; ic < ncol; ic++)
    {
      if (mapSeqnb2Cname.count(ic))
 curName = (mapSeqnb2Cname^ic) + filling;
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
 curName = (mapSeqnb2Rname^ir) + filling;
      else
 curName = filling;

      strncpy(&(rownamestore[ir * placeSize]), curName.c_str(), wordSize);
      rownamestore[(ir +1) * placeSize - 1] = endOfWord;
      rownames[ir]= &(rownamestore[ir*placeSize]);
    }

  int status;
  check(cplexProbPtr == NULL, "LpSolverInterface.load(): problem must have been  created");
  check(CPXcopylpwnames(cplexEnv, cplexProbPtr,
   ncol, nrow,
   minmaxStatus, obj,
   rhs, rsense,
   matbeg, matcnt, matind, matval,
                        dlb, dub,
   range, colnames, rownames),
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


void LpSolverInterface::unLoadFormulation()
{
  //printout("LpSolverInterface::~unLoadFormulation()", 6);
  //return;
  _ncol = 0;
  _nrow = 0;

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
  /// Assumes sorted colMatrix

  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);
  for (set<ProbCoef>::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  double *rhs = new double[nrow];
  fill_n(rhs, nrow, 0);

  char *rsense = new char[nrow];
  fill_n(rsense, nrow, ' ');

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
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
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();
  int curCol(0);
  for(;curCol<ncol;curCol++)
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
          dub[bPtr->ref] = dlb[bPtr->ref];
        }
    }

  double *range = NULL;
  std::string filling = "____________";
  std::string curName;
  char *colnamestore(NULL);
  char **colnames = NULL;
  char *rownamestore(NULL);
  char **rownames = NULL;

  if(param.MipSolverRecordNamesInFormulation)
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
     curName = (mapSeqnb2Cname^ic) + filling;
   else
     curName = filling;

   strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
   colnamestore[(ic +1) * placeSize - 1] = endOfWord;
   colnames[ic]= &(colnamestore[ic*placeSize]);
 }

      for(int ir = 0; ir < nrow; ir++)
 {
   if (mapSeqnb2Rname.count(ir))
     curName = (mapSeqnb2Rname^ir) + filling;
   else
     curName = filling;

   strncpy(&(rownamestore[ir * placeSize]), curName.c_str(), wordSize);
   rownamestore[(ir +1) * placeSize - 1] = endOfWord;
   rownames[ir]= &(rownamestore[ir*placeSize]);
 }

      /// Not required since it is done in CPXcopyctype
      check(cplexProbPtr == NULL,
     "MipSolverInterface.load(): problem must have been  created");

      check(CPXcopylpwnames(cplexEnv, cplexProbPtr,
       ncol, nrow,
       minmaxStatus, obj,
       rhs, rsense,
       matbeg, matcnt, matind, matval,
       dlb, dub, range,
       colnames, rownames),
     "CPXcopylpwnames: Failed to copy problem data.");
    }
  else
    {
      check(cplexProbPtr == NULL,
     "MipSolverInterface.load(): problem must have been  created");

      check(CPXcopylp(cplexEnv, cplexProbPtr,
        ncol, nrow,
        minmaxStatus, obj,
        rhs, rsense,
        matbeg, matcnt, matind, matval,
        dlb, dub,
        range),
     "CPXcopylp: Failed to copy problem data.");
    }

  char *ctype = new char[ncol];
  fill_n(ctype, ncol, 'C');
  for (set<ProbType>::const_iterator typePt = types.begin();
       typePt != types.end();
       typePt++)
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

  check(CPXcopyctype(cplexEnv,
       cplexProbPtr,
       ctype ),
 "CPXcopyctype: Failed to add integer type info");

  int ndir = 0;
  int *pricind = NULL;
  int *prival = NULL;
  int *pridir = NULL;

  if(!directs.empty())
    {
      pricind = new int[directs.size()];
      fill_n(pricind, directs.size(), -1);

      prival = new int[directs.size()];
      fill_n(prival, directs.size(), 0);

      pridir = new int[directs.size()];
      fill_n(pridir, directs.size(), 0);

      for (set<ProbIntC>::const_iterator dPtr = directs.begin();
    dPtr != directs.end();
    dPtr++)
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
      check(CPXcopyorder(cplexEnv, cplexProbPtr, ndir,
    pricind, prival, pridir ),
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
        for(set<ProbSetCoef>::const_iterator sPtr = sets.begin();
     sPtr != sets.end();
     sPtr++)
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
            } while ((sPtr != sets.end()) && (sPtr->setRef == curSet));
        }

      check(CPXcopysos(cplexEnv, cplexProbPtr,
         numsos, numsosnz,
                       sostype, sosbeg, sosind, sosref, NULL),
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




  /// We work only on lp with lpsolve, nolt allow but could used with BaB of lpSolve
  return;
}


void LpSolverInterface::addCols(const std::set<ProbCoef> & objectiveRow,
    const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
    const std::set<ProbBound> & bounds,
    const std::map<int, std::string> & mapSeqnb2Cname)
{
  //printout("LpSolverInterface::addCols");
  //printout("LpSolverInterface::addCols()");
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  int readNcol(0);
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  require(readNcol == _ncol,
   "LpSolverInterface::addCols: readNcol != _ncol");

  int newnzInCols = colMatrix.size();

  double *eobjx = new double[newcol];
  fill_n(eobjx, newcol, 0);

  for (set<ProbCoef>::const_iterator oPtr = objectiveRow.begin();
       oPtr != objectiveRow.end();
       oPtr++)
    eobjx[oPtr->colRef - _ncol] = zero(oPtr->coef);

  int *ematbeg = new int[newcol + 1];
  fill_n(ematbeg, newcol + 1, 0);

  int *ematind = new int[newnzInCols];
  fill_n(ematind, newnzInCols, 0);

  double *ematval = new double[newnzInCols];
  fill_n(ematval, newnzInCols, 0);

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
  char *colnamestore = new char[newcol * placeSize];
  fill_n(colnamestore, newcol * placeSize, endOfWord);

  char **colnames = new char*[newcol];
  fill_n(colnames, newcol, static_cast<char*>(NULL));

  for(int ic = 0; ic < newcol; ic++)
    {
      int colRefNb = ic + _ncol;
      if (mapSeqnb2Cname.count(colRefNb))
 curName = (mapSeqnb2Cname^colRefNb) + filling;
      else
 curName = filling;

      strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
      colnamestore[(ic +1) * placeSize - 1] = endOfWord;
      colnames[ic]= &(colnamestore[ic*placeSize]);
    }

  check(CPXaddcols(cplexEnv, cplexProbPtr,
     newcol, newnzInCols,
     eobjx, ematbeg, ematind, ematval, ebdl, ebdu,
     colnames),
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

void LpSolverInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  //printout("LpSolverInterface::delCols");
  //printout("LpSolverInterface::delCols()");
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  int readNcol;

  int begIndex, endIndex;

  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  require(readNcol <= _ncol, "LpSolverInterface::delCols: readNcol > _ncol");
  require(nbCol2Delete <= readNcol,
   "LpSolverInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
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
   check(CPXdelcols(cplexEnv,
      cplexProbPtr,
      begIndex,
      endIndex),
  "CPXdelcols: could not delete cols");

   /// Only available on Unix, because it doesn't compile on Windows.
   if ( log_cplex != NULL )
     fflush(log_cplex);

   begIndex = dindex[cpt];
   endIndex = dindex[cpt];
 }
    }
  check(CPXdelcols(cplexEnv,
     cplexProbPtr,
     begIndex,
     endIndex),
 "CPXdelcols: could not delete cols");


  /// Only available on Unix, because it doesn't compile on Windows.
  if ( log_cplex != NULL )
    fflush(log_cplex);


  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}

void LpSolverInterface::addRows( const std::set<ProbBound> & rhsv,
     const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
     const std::map<int, std::string> & mapSeqnb2Rname)
{
  //printout("LpSolverInterface::addRows");
  //printout("LpSolverInterface::addRows()");
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;
  int newnzInNewRows = rowMatrix.size();
  int readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  if (printL(7))
    std::cout << "newrows = " << newrows
       << "  newnzInNewRows = " << newnzInNewRows
       << "  readNrow = " << readNrow
       << "  _nrow = " << _nrow
       << std::endl;

  require(readNrow == _nrow,
   "LpSolverInterface::addRows: readNrow != _nrow");

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

  int *ematind = new int[newnzInNewRows];
  fill_n(ematind, newnzInNewRows, 0);

  double *ematval = new double[newnzInNewRows];
  fill_n(ematval, newnzInNewRows, 0);

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
          if (mPtr == rowMatrix.end()) break;
        }
    }

  /// For dummy element newrows + 1
  ematbeg[newRowRef] = cnt;

  std::string filling = "____________";
  std::string curName;

  char * emptyName = new char[newrows * placeSize];
  fill_n(emptyName, newrows * placeSize, endOfWord);

  char** rownames = new char*[newrows];
  fill_n(rownames, newrows, static_cast<char*>(NULL));

  for(int ir = 0; ir < newrows; ir++)
    {
      int newRowRef = ir + _nrow;
      // rownames[ir] =  new char[placeSize];
      if (mapSeqnb2Rname.count(newRowRef))
 {
   curName = (mapSeqnb2Rname^newRowRef) + filling;
 }
      else
 {
   curName = filling;
 }
      strncpy(&(emptyName[ir * placeSize]), curName.c_str(), wordSize);
      emptyName[(ir +1) * placeSize - 1] = endOfWord;
      rownames[ir]= &(emptyName[ir*placeSize]);
    }


  int newlyAddedColsAlongSideNewRows(0);
  char ** namesOfNewlyAddedCols = NULL;

  check(CPXaddrows(cplexEnv, cplexProbPtr,
     newlyAddedColsAlongSideNewRows, newrows, newnzInNewRows,
     erhs, esense, ematbeg, ematind, ematval,
     namesOfNewlyAddedCols, rownames),
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

void LpSolverInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  //printout("LpSolverInterface::delRows");
  //printout("LpSolverInterface::delRows()");
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  int readNrow;

  int begIndex, endIndex;

  readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  require(readNrow <= _nrow,
   "LpSolverInterface::delRowss: readNrow > _nrow");
  require(nbRow2Delete <= readNrow,
   "LpSolverInterface::delRows: nbRow2Delete > readNrow");

  int *dindex = new int[nbRow2Delete + 1];
  fill_n(dindex, nbRow2Delete, -1);
  int cnt(0);

  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin();
       iPtr != indexSetOfRow2Delete.end();
       iPtr++)
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
   check(CPXdelrows(cplexEnv,
      cplexProbPtr,
      begIndex,
      endIndex),
  "CPXdelrows: could not delete cols");

   /// Only available on Unix, because it doesn't compile on Windows.
   if ( log_cplex != NULL ) fflush(log_cplex);


   begIndex = dindex[cpt];
   endIndex = dindex[cpt];
 }
    }
  check(CPXdelrows(cplexEnv,
     cplexProbPtr,
     begIndex,
     endIndex),
 "CPXdelcols: could not delete cols");

  /// Only available on Unix, because it doesn't compile on Windows.
  if ( log_cplex != NULL )
    fflush(log_cplex);


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
  //printout("MipSolverInterface::addDirectives");
  check(1,
 "MipSolverInterface::addDirectives cannot by used. Instead one should record directives at the outset");

  // #ifdef (0)
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

  // #ifdef _WITH_XPRESS
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
  //printout("LpSolverInterface::getObjCoef");
  //printout("LpSolverInterface::getObjCoef()", 6);
  double *coef = new double[1];
  check(CPXgetobj(cplexEnv,
    cplexProbPtr,
    coef,
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
  check(CPXchgobj(cplexEnv,
    cplexProbPtr,
    1,
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
  check(CPXchgcoef(cplexEnv,
     cplexProbPtr,
     pc.rowRef,
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
  check(CPXchgrhs(cplexEnv,
    cplexProbPtr,
    1,
    index,
    rhs),
 "could not chgths");
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

  check(CPXchgrhs(cplexEnv,
    cplexProbPtr,
    nrhs,
    index,
    rhs),
 "could not chgths");

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
  check(CPXchgbds(cplexEnv,
		  cplexProbPtr,
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
  check(CPXchgctype(cplexEnv,
      cplexProbPtr,
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

void LpSolverInterface::getObjVal(Double & objV)
{
  //printout("LpSolverInterface::getObjVal");
  double objval;
  check(CPXgetobjval(cplexEnv,
       cplexProbPtr,
       &objval),
 "could not get objective value");
  objV = zero(objval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsolve
void MipSolverInterface::getObjVal(Double & objV)
{
  //printout("MipSolverInterface::getObjVal");
  double objval = 0.0;
  check(CPXgetobjval(cplexEnv,
       cplexProbPtr,
       &objval),
 "could not get objective value");





  objV = zero(objval, param.HIGHPRECISION);

  return;
}

void LpSolverInterface::getDualBound(Double & val)
{
  //printout("LpSolverInterface::getDualBound");
  double dualval = 0.0;
  // check(CPXgetdualval(cplexEnv, cplexProbPtr, &dualval), "could not get dual bound");
  int nrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  double * dualsol = new double[nrow];
  check(CPXgetpi(cplexEnv,
   cplexProbPtr,
   dualsol,
   0,
   nrow-1),
 "could not get dual bound");

  double * rhs = new double[nrow];
  int * rmatbeg = new int[nrow];
  int n;
  CPXgetrows(cplexEnv,
      cplexProbPtr,
      & n,
      rmatbeg,
      NULL,
      NULL,
      0,
      & n,
      0,
      nrow-1);

  for(int i = 0; i < nrow; i++)
    dualval = dualval + dualsol[i] * rhs[i];
  val = zero(dualval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsole
void MipSolverInterface::getDualBound(Double & val)
{
  //printout("MipSolverInterface::getDualBound");
  double dualval = 0.0;
  int status;
  status = CPXgetstat(cplexEnv, cplexProbPtr);
  if(status == CPXMIP_OPTIMAL)
    {
      require(1, "MIP optimal");
      check(CPXgetobjval(cplexEnv,
    cplexProbPtr,
    &dualval),
     "could not get dual bound");
    }
  else
    dualval = -BapcodInfinity;





  val = zero(dualval, param.HIGHPRECISION);

  return;
}


void LpSolverInterface::getPrimalBound(Double & val)
{
  //printout("LpSolverInterface::getPrimalBound");
  double primalval = 0.0;
  check(CPXgetobjval(cplexEnv,
       cplexProbPtr,
       &primalval),
 "could not get primal bound");
  val = zero(primalval, param.HIGHPRECISION);
  return;
}


/// Not yet for LPsolve
void MipSolverInterface::getPrimalBound(Double & val)
{
  //printout("MipSolverInterface::getPrimalBound");
  double primalval = 0.0;
  check(CPXgetobjval(cplexEnv,
       cplexProbPtr,
       &primalval),
 "could not get primal bound");





  val = zero(primalval, param.HIGHPRECISION);
  return;
}

bool LpSolverInterface::getOptimStatus(SolutionStatus & lpStatus,
           SolutionStatus & mipStatus)
{
  //printout("LpSolverInterface::getOptimStatus");
  int status;
  status = CPXgetstat(cplexEnv, cplexProbPtr);
  if (status == CPX_STAT_OPTIMAL)
    {
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP optimal"
    << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }

  if (printL(4))
    MPSwrite();

  if (status == CPX_STAT_INFEASIBLE)
    {
      lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP infeasible"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_UNBOUNDED)
    {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: LP unbounded"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_OBJ_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: objective limit exceeded in phase II"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_IT_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: iteration limit exceeded"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_TIME_LIM )
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: time limit exceeded"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_NUM_BEST)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: problem non-optimal"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_OPTIMAL_INFEAS)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: Optimal solution found, unscaled infeasibilities"
    << std::endl;

      return(false);
    }

  if (status == CPX_STAT_ABORT_USER)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "LpSolverInterface::getOptimStatus: Aborded"
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
  status = CPXgetstat(cplexEnv, cplexProbPtr);
  if (status == CPXMIP_OPTIMAL)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: global search complete and an integer solution has been found"
    << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }

  if (printL(4))
    MPSwrite();

  if (status == CPXMIP_OPTIMAL_TOL)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP optimal solution found within tolerance"
    << std::endl;

      // Optimal solution found
      return(true);
    }

  if (status == CPXMIP_INFEASIBLE)
    {
      mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP infeasible"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_SOL_LIM)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP solutions limit exceeded"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_NODE_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP node limit exceeded, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_NODE_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP node limit exceeded, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_TIME_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP time limit exceeded, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_TIME_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP time limit exceeded, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_FAIL_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP error termination, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_FAIL_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP error termination, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_MEM_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP Treememory limit exceeded, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_MEM_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP Treememory limit  exceeded, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_ABORT_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP aborded, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_ABORT_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP aborded, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_OPTIMAL_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP optimal with unscaled infeasibilities"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_FAIL_FEAS_NO_TREE)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP out of memory, no tree, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_FAIL_INFEAS_NO_TREE)
    {
      mipStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP out of memory, no tree, no integer solution exists"
    << std::endl;

      return(false);
    }

  if (status == CPXMIP_NODE_LIM_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP node file size limit exceeded, integer solution exists"
    << std::endl;

      return(true);
    }

  if (status == CPXMIP_NODE_LIM_INFEAS)
    {
      mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
 std::cout << "MipSolverInterface::getOptimStatus: MIP node file size limit  exceeded, no integer solution exists"
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
  readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
  require(readNcol <= _ncol,
   "LpSolverInterface::getSol: readNcol > _ncol");

  double *x = new double[_ncol];
  fill_n(x, _ncol, 0);
  check(CPXgetx(cplexEnv,
  cplexProbPtr,
  x,
  0,
  _ncol-1),
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
  readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  require(readNrow <= _nrow,
   "LpSolverInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[_nrow];
  fill_n(dual, _nrow, 0);

  char *rsense = new char[_nrow];
  fill_n(rsense, _nrow, ' ');

  check(CPXgetpi(cplexEnv,
   cplexProbPtr,
   dual,
   0,
   _nrow-1),
 " could not get dual solution");

  check(CPXgetsense(cplexEnv,
      cplexProbPtr,
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

//  const bool LpSolverInterface::preProcessorChangedForm(int * ptrProbRowCnt, 
// 						       const bool & checkVarOnly)
//  {
//    //printout("LpSolverInterface::preProcessorFeedBack()");
//    require(*ptrProbRowCnt == _nrow, 
// 	   "LpSolverInterface::preProcessorChangedForm(): not (*ptrProbRowCnt == _nrow)");

//    // /// Since matrix is saved and restore in optimise()
//    // return(false);
//    int readNcol;
//  #ifdef _WITH_XPRESS
//    getipv(N_NCOL, &readNcol);
//  #endif
//  #ifdef _WITH_CPLEX
//    readNcol = CPXgetnumcols(cplexEnv, cplexProbPtr);
//  #endif
//  #ifdef _WITH_LPSOLVE
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

//  #ifdef _WITH_XPRESS
//    getipv(N_NROW, &readNrow);
//  #endif
//  #ifdef _WITH_CPLEX
//    readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
//  #endif
//  #ifdef _WITH_LPSOLVE
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
  //printout("LpSolverInterface::LPwrite");
  //printout("LpSolverInterface::LPwrite()");
  //   check(CPXwriteprob(cplexEnv, 
  // 		     cplexProbPtr, 
  // 		     "lastprob.lp", 
  // 		     "LP"), 
  // 	" could not write prob",  
  // 	ProgStatus::terminate);

  int irow, jcol, iel, nchar;
  int readNrow, readNcol, idone, readNglents, readNsets;
  char coltype;
  double readBdl, readBdu;
  char probname[4 * wordSize];
  int pnsurplus(0);

  check(CPXgetprobname(cplexEnv,
         cplexProbPtr,
         probname,
         4 * wordSize,
         &pnsurplus),
 " could not getprobname");

  require(pnsurplus >= 0,
   "CPXgetprobname did not have enough space ");

  readNrow = CPXgetnumrows(cplexEnv,
      cplexProbPtr);

  readNcol = CPXgetnumcols(cplexEnv,
      cplexProbPtr);

  int cnsurplus(0);
  char *colnamestore = new char[readNcol * placeSize];
  fill_n(colnamestore, readNcol * placeSize, endOfWord);
  char **colnames = new char*[readNcol];
  fill_n(colnames, readNcol, static_cast<char*>(NULL));

  for(int ic = 0; ic < readNcol; ic++)
    colnames[ic] = &(colnamestore[ic*placeSize]);

  if(param.MipSolverRecordNamesInFormulation)
    check(CPXgetcolname(cplexEnv,
   cplexProbPtr,
   colnames,
   colnamestore,
   readNcol * placeSize,
   &cnsurplus,
   0,
   readNcol-1),
   " could not getcolname");

  require(cnsurplus >= 0,
   "CPXgetcolname did not have enough space ");
  printf("PROBLEM: %s whose formulation ref number is %d\n", probname, ref());

  /// Output Obj
  if ( minmaxStatus ==1 ) printf("Minimize\n");
  else printf("Maximize\n");

  int *readMclind = new int[readNcol];
  fill_n(readMclind, readNcol, 0);
  double *readDmatval = new double[readNcol];
  fill_n(readDmatval, readNcol, 0);
  check(CPXgetobj(cplexEnv,
    cplexProbPtr,
    readDmatval,
    0,
    readNcol-1),
 " could not getobj");

  nchar = printf(" Objectif:");
  for (iel = 0; iel < readNcol; iel++)
    if (Double(readDmatval[iel]) != 0)
      {
 nchar += printf(" %+.9g %s",
   readDmatval[iel],
   &(colnamestore[placeSize * iel]));
 if (nchar >= 80-9-13)
   {
     printf("\n");
     nchar = 0;
   }
      }

  /// Output Matrix
  printf("\nSubject To\n");
  for(irow = 0; irow < readNrow; irow++)
    printRow(irow,
      readNcol,
      readMclind,
      readDmatval,
      colnamestore);
  printf("Bounds\n");

  for (jcol = 0; jcol < readNcol; jcol++)
    {
      check(CPXgetlb(cplexEnv,
       cplexProbPtr,
       &readBdl,
       jcol,
       jcol),
     " could not getlb");

      check(CPXgetub(cplexEnv,
       cplexProbPtr,
       &readBdu,
       jcol,
       jcol),
     " could not getub");

      if (readBdu < CPX_INFBOUND / 10)
 printf(" %.9g <= %s <= %.9g\n", readBdl, &colnamestore[placeSize * jcol], readBdu);
      else if (readBdl != 0.0)
 printf("      %s >= %.9g\n", &colnamestore[placeSize * jcol], readBdl);
    }

  /// Always Pure for LPSolve yet
  if (_pureLP)
    printf("End\n");
  else
    {
      readNglents = CPXgetnumint(cplexEnv, cplexProbPtr)
 + CPXgetnumbin(cplexEnv, cplexProbPtr);
      readNsets = CPXgetnumsos(cplexEnv, cplexProbPtr);
      if( readNglents )
        {
   ///@todo Replace printf by cout
          printf("Integers\n");
          for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
            {
              check(CPXgetctype(cplexEnv,
    cplexProbPtr,
    &coltype,
    jcol,
    jcol),
      " could not getctype");

              if (coltype != 'C')
  {
    ///@todo Replace printf by cout
    printf(" %s ", &colnamestore[placeSize * jcol]);
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

      /// Add directives informations if any
      if( readNglents )
        {
          std::cout << "Directives:" << std::endl;

          int Lndir = readNglents + readNsets;
          int *Lmcols = new int[Lndir];
   fill_n(Lmcols, Lndir, 0);
          int *Lmpri = new int[Lndir];
   fill_n(Lmpri, Lndir, 0);

          int *LbrDir = new int[Lndir];
   fill_n(LbrDir, Lndir, ' ');
          int dsurplus(0);
          check(CPXgetorder(cplexEnv,
       cplexProbPtr,
       &Lndir,
       Lmcols,
       Lmpri,
       LbrDir,
       Lndir,
       &dsurplus),
  " could not getorder");

          for (jcol = 0; jcol < Lndir; jcol++)
            std::cout << &(colnamestore[placeSize * Lmcols[jcol]])
        << " " << Lmpri[jcol]
        << " " << LbrDir[jcol]
        << std::endl;

          delete[] LbrDir;
   LbrDir = NULL;
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
  delete[] colnamestore;
  colnamestore = NULL;
  delete[] colnames;
  colnames = NULL;
  return;
}


void LpSolverInterface::MPSwrite()
{
  //printout("LpSolverInterface::MPSwrite");
  //printout("LpSolverInterface::MPSwrite()");
  check(CPXwriteprob(cplexEnv,
       cplexProbPtr,
       "curprob.lp",
       "LP"),
 " could not write prob",
 ProgStatus::terminate);

  check(CPXwriteprob(cplexEnv,
       cplexProbPtr,
       "curprob.mps",
       "MPS"),
 " could not write prob",
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
  int iel, readNels(0), nchar, nbRows(1);
  double readDrhs;
  char readRtype, readSense[3];
  int readMstart[2];
  int rsurplus(0);
  check(CPXgetrows(cplexEnv,
     cplexProbPtr,
     &readNels,
     readMstart,
     readMclind,
     readDmatval,
     readNcol,
     &rsurplus,
     irow,
     irow),
        "could not getrows");

  check(CPXgetrhs(cplexEnv,
    cplexProbPtr,
    &readDrhs,
    irow,
    irow),
 "could not getrhs");

  check(CPXgetsense(cplexEnv,
      cplexProbPtr,
      &readRtype,
      irow,
      irow),
 "could not getrowtype");

  int rnsurplus(0);
  char *rownamestore = new char[placeSize * nbRows];
  fill_n(rownamestore, placeSize * nbRows, endOfWord);

  char **rownames = new char*[nbRows];
  fill_n(rownames, nbRows, static_cast<char*>(NULL));

  for(int ir = 0; ir < nbRows; ir++)
    rownames[ir] = &(rownamestore[ir*placeSize]);

  if(param.MipSolverRecordNamesInFormulation)
    check(CPXgetrowname(cplexEnv,
   cplexProbPtr,
   rownames,
   rownamestore,
   placeSize * nbRows,
   &rnsurplus,
   irow,
   irow),
   "could not getnames");

  require(rnsurplus >= 0,
   "CPXgetrowname did not get enough space");

  strcpy(readSense,
  readRtype=='O'
  ?"$" :(readRtype=='E'
  ?"=" :(readRtype=='L'
         ?"<=" :">=" ) ));

  nchar = printf(" %s:", rownamestore);
  if (readNels > 0)
    for (iel = 0; iel < readNels; iel++)
      {
        nchar += printf(" %+.9g %s",
   readDmatval[iel],
   &readCnames[placeSize * readMclind[iel]]);
        if (nchar >= 80-9-13) {
   printf("\n");
   nchar = 0;
        }
      }

  printf( " %s %.9g", readSense, readDrhs);
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
  check(CPXchgprobtype(cplexEnv,
         cplexProbPtr,
         CPXPROB_LP),
 "could not change type to LP");

  check(CPXlpopt(cplexEnv,
   cplexProbPtr),
 "could not solve LP by dualopt");
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

  const char * filename_str;
  const char * filetype_str;
  if(printL(1))
    CPXwriteprob (cplexEnv,
    cplexProbPtr,
    "myprob.mps",
    NULL);





  check(CPXmipopt(cplexEnv, cplexProbPtr),
 "could not solve MIP");
  return;
}

void LpSolverInterface::reset()
{
  //printout("LpSolverInterface::reset");
  //printout("LpSolverInterface::reset()");
  /// Do nothing

  return;
}

void MipSolverInterface::reset()
{
  //printout("MipSolverInterface::reset()", 6);
  // CPXsetdefaults(cplexEnv);
  return;
}

void LpSolverInterface::makeSpaceForLoadingForm()
{
  /// Save loaded form if any
  //printout("LpSolverInterface::makeSpaceForLoadingForm()", 6);
  _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << " _formCurrentlyLoaded = " << _formCurrentlyLoaded
       << std::endl;
  return;
}

void LpSolverInterface::saveCopyOfCurForm()
{
  //printout("LpSolverInterface::saveCopyOfCurForm()", 6);
  /// Save loaded form if any
  // _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << "LpSolverInterface::saveCopyOfCurForm(): formCurrentlyLoaded = "
       << _formCurrentlyLoaded
       << std::endl;
  return;
}

void LpSolverInterface::load()
{
  //printout("LpSolverInterface::load()", 6);
  if (printL(6))
    std::cout << "LpSolverInterface::load(): formCurrentlyLoaded = " << _formCurrentlyLoaded
       << std::endl;
  return;
}

void LpSolverInterface::unload(const bool & deleteMat)
{
  //printout("LpSolverInterface::unload");
  //printout("LpSolverInterface::unload()");
  // if (deleteMat == false)

  return;
}
#endif // _WITH_LPSOLVE
