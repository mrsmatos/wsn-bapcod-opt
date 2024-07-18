/**
 *
 * This file bcClpSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _CLP_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcClpSolverC.hpp"
#include "bcPrintC.hpp"

#include "bcModelC.hpp"
#include "bcMasterConfC.hpp"
#include "bcFormC.hpp"

#include "ClpSimplex.hpp"
#include "CoinWarmStartBasis.hpp"
#include "bcMathProgSolverException.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 16;
const char endOfWord = '\0';
const std::string filling = "________________";


/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

using namespace std;

//OK
LpClpInterface::LpClpInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  MathProgSolverInterface(bapcodInit, ref, name), clpPrecision(10e-8)
{
  char *probname = new char[name.size() + 1];
  snprintf(probname, 100, "%s", name.c_str());
  probname[name.size()] = endOfWord;

  clpModel.setStrParam(ClpProbName, probname);
  clpModel.setLogLevel(0);

  delete [] probname; probname = NULL;

  return;
}

//need debug
LpClpInterface::~LpClpInterface()
{
  return;
}

//OK
void LpClpInterface::loadFormulation(const std::string & name,
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
    std::cout << "LpClpInterface::loadFormulation() ncol = " << ncol
	      << " , nrow = " << nrow
	      << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  /// Assumes sorted colMatrix
  double *obj = new double[ncol];
  fill_n(obj, ncol, 0);

  for (ProbCoefContainer::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    obj[oPtr->colRef] = zero(oPtr->coef);

  double *rowUpper = new double[nrow];
  fill_n(rowUpper, nrow, COIN_DBL_MAX);
  double *rowLower = new double[nrow];
  fill_n(rowLower, nrow, -COIN_DBL_MAX);

  for (ProbRhsContainer::const_iterator rvPtr = rhsv.begin(); rvPtr != rhsv.end(); rvPtr++)
    {
	  if (rvPtr->sense != 'G')
		  rowUpper[rvPtr->ref] = zero(rvPtr->bound);
	  if (rvPtr->sense != 'L')
		  rowLower[rvPtr->ref] = zero(rvPtr->bound);
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

  clpModel.loadProblem(ncol, nrow, matbeg, matind, matval, dlb, dub, obj, rowLower, rowUpper);
  clpModel.setOptimizationDirection(minmaxStatus);

  for (int ir = 0; ir <nrow; ir++)
    {
      std::string rowName(rownames[ir]);
      clpModel.setRowName(ir, rowName);
    }

  for (int ic = 0; ic <ncol; ic++)
    {
      std::string colName(colnames[ic]);
      clpModel.setColumnName(ic, colName);
    }

  delete [] obj;
  obj = NULL;
  delete [] rowLower;
  rowLower = NULL;
  delete [] rowUpper;
  rowUpper = NULL;
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

//OK
void LpClpInterface::unLoadFormulation()
{
  return;
  _ncol = 0;
  _nrow = 0;

  return;
}



//OK
void LpClpInterface::addCols(const ProbCoefContainer & objectiveRow,
			     ProbColMatrixContainer & colMatrix,
			     const ProbBoundContainer & bounds,
			     const std::map < int, std::string > & mapSeqnb2Cname)
//                          (const std::set<ProbCoef> & objectiveRow,
//                           const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
//			     const std::set<ProbBound> & bounds,
//			     const std::map<int, std::string> & mapSeqnb2Cname)
{
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  int readNcol(0);
  readNcol = clpModel.getNumCols();
  bapcodInit().require(readNcol == _ncol,
	  "LpClpInterface::addCols: readNcol != _ncol");

  if (printL(7))
    std::cout << "_ncol = " <<  _ncol
	      << "newcol = " << newcol
	      << std::endl;

  int newnzInCols = colMatrix.size();

  double *eobjx = new double[newcol];
  fill_n(eobjx, newcol, 0);

  for (set<ProbCoef>::const_iterator oPtr = objectiveRow.begin();
       oPtr != objectiveRow.end();
       oPtr++)
    {
      if (printL(7))
	std::cout << "eobjx mPtr->colRef  = " << oPtr->colRef 
		  << ", oPtr->colRef - _ncol = " << oPtr->colRef - _ncol
		  << ", oPtr->coef= " << oPtr->coef
		  << std::endl;

      eobjx[oPtr->colRef - _ncol] = zero(oPtr->coef);
    }

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
	std::cout << " newColRef = " << newColRef
		  << " cnt = " << cnt
		  << " ematbeg[newColRef] = " << ematbeg[newColRef]
		  << std::endl;

      if (mPtr == colMatrix.end())
	continue;

      if (newColRef < (mPtr->colRef - _ncol))
	continue;

      if (printL(7))
	std::cout << "ColMatrix = " << *mPtr
		  << " newColRef = " << newColRef
		  << " mPtr->colRef = " << mPtr->colRef
		  << " _ncol = " <<  _ncol
		  << " newcol = " << newcol
		  << " mPtr->colRef - _ncol = " << mPtr->colRef - _ncol
		  << std::endl;

      while (newColRef == (mPtr->colRef - _ncol))
        {
          ematind[cnt] = mPtr->rowRef;
          ematval[cnt] = zero(mPtr->coef);

         if (printL(7))
	    std::cout << " cnt = " << cnt
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

  std::string curName;
  char *colnamestore = new char[newcol * placeSize];
  fill_n(colnamestore, newcol * placeSize, endOfWord);

  char **colnames = new char*[newcol];
  fill_n(colnames, newcol, static_cast<char*>(NULL));

  for(int ic = 0; ic < newcol; ic++)
    {
      int colRefNb = ic + _ncol;
      if (mapSeqnb2Cname.count(colRefNb))
	curName = mapSeqnb2Cname.at(colRefNb) + filling;
      else
	curName = filling;

      strncpy(&(colnamestore[ic * placeSize]), curName.c_str(), wordSize);
      colnamestore[(ic +1) * placeSize - 1] = endOfWord;
      colnames[ic]= &(colnamestore[ic*placeSize]);
    }

  clpModel.addColumns(newcol, ebdl, ebdu, eobjx, ematbeg, ematind, ematval);

  for (int ic = 0; ic<newcol; ic++)
    {
      string colName(colnames[ic]);
      clpModel.setColumnName(_ncol + ic, colName);
    }

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

//OK
void LpClpInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  int readNcol;

  int begIndex, endIndex;

  readNcol = clpModel.getNumCols();
  bapcodInit().require(readNcol <= _ncol, "LpClpInterface::delCols: readNcol > _ncol");
  bapcodInit().require(nbCol2Delete <= readNcol,
	  "LpClpInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
    {
      dindex[cnt++] = *iPtr;
    }

  clpModel.deleteColumns(cnt, dindex);

  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}

//OK
void LpClpInterface::addRows( const std::set<ProbBound> & rhsv,
				const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
				const std::map<int, std::string> & mapSeqnb2Rname)
{
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;

  int newnzInNewRows = rowMatrix.size();
  int readNrow = clpModel.getNumRows();

  if (printL(7))
    std::cout << " LpClpInterface::addRows() newrows = " << newrows
	      << "  newnzInNewRows = " << newnzInNewRows
	      << "  readNrow = " << readNrow
	      << "  _nrow = " << _nrow
	      << std::endl;


  bapcodInit().require(readNrow == _nrow,
	  "LpClpInterface::addRows: readNrow != _nrow");

  double *rowUpper = new double[newrows];
  fill_n(rowUpper, newrows, COIN_DBL_MAX);
  double *rowLower = new double[newrows];
  fill_n(rowLower, newrows, -COIN_DBL_MAX);

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      if (printL(7))
	std::cout << "rhsv = " << *rvPtr;

	  if (rvPtr->sense != 'G')
		  rowUpper[rvPtr->ref - _nrow] = zero(rvPtr->bound);
	  else if (rvPtr->sense != 'L')
		  rowLower[rvPtr->ref - _nrow] = zero(rvPtr->bound);
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
	  curName = mapSeqnb2Rname.at(newRowRef) + filling;
	}
      else
	{
	  curName = filling;
	}
      strncpy(&(emptyName[ir * placeSize]), curName.c_str(), wordSize);
      emptyName[(ir +1) * placeSize - 1] = endOfWord;
      rownames[ir]= &(emptyName[ir*placeSize]);
    }

  clpModel.addRows(newrows, rowLower, rowUpper, ematbeg, ematind, ematval);

  for (int ir = 0; ir<newrows; ir++)
    {
      string rowName(rownames[ir]);
      clpModel.setRowName(_nrow + ir, rowName);
    }

  delete[] rownames;
  rownames = NULL;
  delete[] emptyName;
  emptyName = NULL;
  delete[] rowLower;
  rowLower = NULL;
  delete[] rowUpper;
  rowUpper = NULL;
  delete[] ematbeg;
  ematbeg = NULL;
  delete[] ematind;
  ematind = NULL;
  delete[] ematval;
  ematval = NULL;
  _nrow += newrows;

  return;
}

//OK
void LpClpInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  int readNrow = clpModel.getNumRows();

  //int begIndex, endIndex;

  //readNrow = CPXgetnumrows(cplexEnv, cplexProbPtr);
  bapcodInit().require(readNrow <= _nrow,
	  "LpClpInterface::delRowss: readNrow > _nrow");
  bapcodInit().require(nbRow2Delete <= readNrow,
	  "LpClpInterface::delRows: nbRow2Delete > readNrow");

  //int *dindex = new int[nbRow2Delete + 1];
//  fill_n(dindex, nbRow2Delete, -1);
//  int cnt(0);
//
//  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin();
//       iPtr != indexSetOfRow2Delete.end();
//       iPtr++)
//    {
//      dindex[cnt++] = *iPtr;
//    }

//  clpModel.deleteRows(cnt, dindex);
  int* dindex = new int[indexSetOfRow2Delete.size()];
  std::copy(indexSetOfRow2Delete.begin(), indexSetOfRow2Delete.end(), dindex);

//  cerr << "EEEEEEEEEEEEEEEEE: indexSetOfRow2Delete.size() " << indexSetOfRow2Delete.size() <<  endl;
//
//  for (int i = 0; i < indexSetOfRow2Delete.size(); ++i) {
//      cerr << "EEE: dindex[" << i << "]: " << dindex[i]  << "; readNrow: "<< readNrow << " ; _nrow: " << _nrow << endl;
//  }

  clpModel.deleteRows(indexSetOfRow2Delete.size(), dindex);

//  cerr << "AAAAAAAAAAAAAAAAA" << endl;

  delete[] dindex;
  dindex = NULL;
  _nrow -= nbRow2Delete;

  return;
}

//OK
void LpClpInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  bapcodInit().check(1, "LP can not have directives");
  return;
}

//OK
void LpClpInterface::getObjCoef(ProbCoef & pc)
{
  const double *objCoef = clpModel.getObjCoefficients();
  pc.coef = objCoef[pc.colRef];

  return;
}

//OK
void LpClpInterface::chgObjCoef(const ProbCoef & pc)
{
  clpModel.setObjCoeff(pc.colRef, pc.coef);
 
  return;
}

//OK
void LpClpInterface::chgMatCoef(const ProbCoef & pc)
{
  clpModel.modifyCoefficient(pc.rowRef, pc.colRef, (double)pc.coef);
  
  return;
}

// OK
void LpClpInterface::chgRhs(const ProbBound & pb)
{
  double lower = -COIN_DBL_MAX;
  double upper = COIN_DBL_MAX;

  if (pb.sense != 'G')
	  upper = zero(pb.bound);
  else if (pb.sense != 'L')
	  lower = zero(pb.bound);

  clpModel.setRowBounds(pb.ref, lower, upper);

  return;
}

// OK
void LpClpInterface::chgRhs(const std::set<ProbBound> & newRhs)
{
  if (newRhs.empty()) return;
  
	double lower;
	double upper;

	for (set<ProbBound>::const_iterator bPtr = newRhs.begin();
		bPtr != newRhs.end();
		bPtr++)
	{
		lower = -COIN_DBL_MAX;
		upper = COIN_DBL_MAX;
		if (bPtr->sense != 'G')
			upper = zero(bPtr->bound);
		if (bPtr->sense != 'L')
			lower = zero(bPtr->bound);

		clpModel.setRowBounds(bPtr->ref, lower, upper);
	}

  return;
}

// need debug
void LpClpInterface::chgBds(const std::set<ProbBound> & newBounds)
{
  if (newBounds.empty())
    return;
 // int Lnbnds(0);
 // int * Lmindex = new int[newBounds.size()];
 // fill_n(Lmindex, newBounds.size(), 0);

 // char *Lqbtype = new char[newBounds.size()];
 // fill_n(Lqbtype, newBounds.size(), ' ');

 // double *Lbnd = new double[newBounds.size()];
 // fill_n(Lbnd, newBounds.size(), 0);

 // for (set<ProbBound>::const_iterator bPtr = newBounds.begin();
 //      bPtr != newBounds.end();
 //      bPtr++)
 //   {
 //     if (printL(7))
	//std::cout << "bound = " << *bPtr;

 //     Lmindex[Lnbnds] = bPtr->ref;
 //     Lbnd[Lnbnds]= zero(bPtr->bound);
 //     if (bPtr->sense == 'L')
 //       {
 //         Lqbtype[Lnbnds] = 'L';
 //         // dlb[bPtr->ref] = Lbnd[Lnbnds];
 //       }
 //     else if ((bPtr->sense == 'U')
	//       || (bPtr->sense == 'I')
	//       || (bPtr->sense == 'B'))
 //       {
 //         Lqbtype[Lnbnds] = 'U';
 //         // dub[bPtr->ref] = Lbnd[Lnbnds];
 //       }
 //     else if (bPtr->sense == 'F')
 //       {
 //         Lqbtype[Lnbnds] = 'B';
 //         // dlb[bPtr->ref] = Lbnd[Lnbnds];
 //         // dub[bPtr->ref] = Lbnd[Lnbnds];
 //       }
 //     Lnbnds++;
 //   }
 // bapcodInit().check(CPXchgbds(cplexEnv,
	//	  cplexProbPtr,
	//	  Lnbnds,
	//	  Lmindex,
	//	  Lqbtype,
	//	  Lbnd),
	//"could not chgbds");

    for (set<ProbBound>::const_iterator bPtr = newBounds.begin();
       bPtr != newBounds.end();
       bPtr++)
	{
		if (printL(7))
			std::cout << "bound = " << *bPtr;

		if (bPtr->sense == 'L')
			clpModel.setColLower(bPtr->ref, zero(bPtr->bound));
		else if (bPtr->sense == 'U' || bPtr->sense == 'I' || bPtr->sense == 'B')
			clpModel.setColumnUpper(bPtr->ref, zero(bPtr->bound));
		else if (bPtr->sense == 'F')
			clpModel.setColumnBounds(bPtr->ref, zero(bPtr->bound), zero(bPtr->bound));
	}

  return;
}

//OK
void LpClpInterface::chgColType(const std::set<ProbType> & newTypes)
{
  bapcodInit().check(1,"LP Form cannot have ColType");

  return;
}

//OK
void LpClpInterface::getObjVal(Double & objV)
{
  double objval = clpModel.getObjValue();
  objV = zero(objval, clpPrecision);

  return;
}

// need debug
void LpClpInterface::getDualBound(Double & val)
{
  double dualval = 0.0;
  int nrow; 
  const double * dualsol;
  const double * rowLower;
  const double * rowUpper;

  nrow = clpModel.getNumRows();
  dualsol = clpModel.getRowPrice();
  rowLower = clpModel.getRowLower();
  rowUpper = clpModel.getRowUpper();

  for (int ir = 0; ir < nrow; ir++)
  {
	// 'G'
	if (rowLower[ir] > - COIN_DBL_MAX && rowUpper[ir] == COIN_DBL_MAX)
		dualval += (dualsol[ir] * rowLower[ir]);
	// 'L'
	else if (rowLower[ir] == -COIN_DBL_MAX && rowUpper[ir] < COIN_DBL_MAX)
		dualval += (dualsol[ir] * rowUpper[ir]);
	// 'E'
	else
		dualval += (dualsol[ir] * rowUpper[ir]);
  }

  val = zero(dualval, clpPrecision);

  return;
}

//OK
void LpClpInterface::getPrimalBound(Double & val)
{
  double primalval = clpModel.getObjValue();
  val = zero(primalval, clpPrecision);

  return;
}

//OK
bool LpClpInterface::getOptimStatus(SolutionStatus & lpStatus,
				      SolutionStatus & mipStatus)
{
    /** Status of problem:
        -1 - unknown e.g. before solve or if postSolve says not optimal
        0 - optimal
        1 - primal infeasible
        2 - dual infeasible
        3 - stopped on iterations or time
        4 - stopped due to errors
        5 - stopped by event handler (virtual int ClpEventHandler::event())
    */

  int status = clpModel.status();

	if (status == 0)
	{
	lpStatus = SolutionStatus::Optimum;

	if (printL(3))
	std::cout << "LpClpInterface::getOptimStatus: LP optimal"
			<< std::endl;

	return true;
	}

    if (printL(4))
    MPSwrite();

	if (status == 1)
	{
		lpStatus = SolutionStatus::Infeasible;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: LP primal infeasible"
		  << std::endl;

		return false;
	}

	if (status == 2)
	{
		lpStatus = SolutionStatus::Infeasible;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: LP dual infeasible"
		  << std::endl;

		return false;
	}

	if (status == 3)
	{
		lpStatus = SolutionStatus::UnSolved;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: stopped on iterations or time"
		  << std::endl;

		return false;
	}

	if (status == 4)
	{
		lpStatus = SolutionStatus::UnSolved;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: stopped due to errors"
		  << std::endl;

		return false;
	}

	if (status == 5)
	{
		lpStatus = SolutionStatus::UnSolved;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: stopped by event handler (virtual int ClpEventHandler::event())"
		  << std::endl;

		return false;
	}

	if (status == -1)
	{
		lpStatus = SolutionStatus::Undefined;
		if (printL(3))
		std::cout << "LpClpInterface::getOptimStatus: unknown e.g. before solve or if postSolve says not optimal"
		  << std::endl;

		return false;
	}

	std::cout << "LpClpInterface::getOptimStatus: undefined status"
	<< std::endl;

  return false;
 
}

//OK
void LpClpInterface::getSol(std::map<int, Double> & primSol,
			      const bool & ifPrint)
{
  primSol.clear();
  int readNcol = clpModel.getNumCols();
  bapcodInit().require(readNcol <= _ncol,
	  "LpClpInterface::getSol: readNcol > _ncol");

  const double * x = clpModel.getColSolution();

  if (printL(6))
    std::cout << "readNcol = " << readNcol
	      << "  _ncol = " << _ncol
	      << std::endl;

  for (int jcol = 0; jcol < _ncol; jcol++)
    if (zero(x[jcol], clpPrecision) != 0)
      {
		if (printL(6))
		  std::cout << "primSol[" << jcol << "] = " << x[jcol]
				<< std::endl;

		primSol[jcol] = x[jcol];
      }

  return;
}

//OK
void LpClpInterface::getSol(std::map<int, Double> & primSol,
			      std::map<int, Double> & dualSol,
			      const int & minmaxStatus,
			      const bool & ifPrint)
{
  getSol(primSol, ifPrint);
  dualSol.clear();
  int readNrow(0);
  readNrow = clpModel.getNumRows();
  bapcodInit().require(readNrow <= _nrow,
	  "LpClpInterface::getSol: readNrow > _nrow");

  int flipsign = -1;

	const double *dual = clpModel.getRowPrice();
	const double *rowLower = clpModel.getRowLower();
	const double *rowUpper = clpModel.getRowUpper();

	for (int ir = 0; ir < _nrow; ir++)
	{
	  if (zero(dual[ir], clpPrecision) != 0)
		{
			// 'G'
			if (rowLower[ir] > - COIN_DBL_MAX && rowUpper[ir] == COIN_DBL_MAX)
				dualSol[ir] = minmaxStatus * flipsign * dual[ir];
			// 'L'
			else if (rowLower[ir] == -COIN_DBL_MAX && rowUpper[ir] < COIN_DBL_MAX)
				dualSol[ir] = minmaxStatus * flipsign * dual[ir];
			// 'E'
			else
				dualSol[ir] = flipsign * dual[ir];

			if (printL(6))
			{
			std::cout << "dual[" << ir << "] = " << dual[ir]
					<< std::endl;

			std::cout << " dualSol[" << ir << "] = " << dualSol[ir]
				  << std::endl;
			}
		}
	}

	return;
}

void LpClpInterface::getSolPool(std::vector< std::map<int, Double> > & primSolPool, const bool & ifPrint)
{
  //bapcodInit().check(1, "LpCplexInterface::getSolPool() should not be called since LP's do not have a solution pool");
  //return;
}

void LpClpInterface::getReducedCost(std::map<int, Double> & redCost, const bool & ifPrint)
{
    redCost.clear();
    int readNcol = clpModel.getNumCols();
    bapcodInit().require(readNcol <= _ncol, "LpClpInterface::getSol: readNcol > _ncol");

    const double * dj = clpModel.getReducedCost();

    for (int jcol = 0; jcol < _ncol; jcol++)
        if (zero(dj[jcol], clpPrecision) != 0)
        {
            if (printL(6))
                std::cout << "primSol[" << jcol << "] = " << dj[jcol] << std::endl;

            redCost[jcol] = dj[jcol];
        }

    return;
}

// by Laurent //to check later
void LpClpInterface::resetUpperCutOffValue()
{
    throw MathProgSolverException("LpClpInterface::resetUpperCutOffValue is not available", true);
}
void LpClpInterface::setUpperCutOffValue(const Double & cutOff)
{
    throw MathProgSolverException("LpClpInterface::setUpperCutOffValue is not available", true);
}
void LpClpInterface::setSolveFromScratch(const bool exactSolution)
{
    throw MathProgSolverException("LpClpInterface::setSolverFromScratch is not available", true);
}
void LpClpInterface::resetSolveFromScratch()
{
    throw MathProgSolverException("LpClpInterface::resetSolverFromScratch is not available", true);
}
void LpClpInterface::setAbortValue(const Double& abortValue)
{
 // DO NOTHING
}

void LpClpInterface::setTimeLimit(const double & seconds)
{
    clpModel.setMaximumSeconds(seconds);
}

void LpClpInterface::setMultiThread(const int & maxNbThread)
{
  //TODO
}
void LpClpInterface::setRelativeMipGapLimit(const double & relativeGap)
{
  //TODO
}
void LpClpInterface::setSearchPriority(const int & flag)
{
  //TODO
}
void LpClpInterface::setWorkingMemorySpace(const double & sizeInMB)
{
  //TODO
}
void LpClpInterface::setLPoptimalityTolerance(const double& tolerance)
{
  //TODO
  clpModel.setPrimalTolerance(tolerance);
  clpModel.setDualTolerance(tolerance);
  clpModel.setCurrentDualTolerance(tolerance);
  clpModel.setCurrentPrimalTolerance(tolerance);
}
void LpClpInterface::setLPfeasibilityTolerance(const double& tolerance)
{
  //TODO
}
void LpClpInterface::setScreenOutput(const bool & flag)
{
  //TODO
}
void LpClpInterface::setMaxBBNodes(const int & maxBBNodes)
{
  //TODO
}


// by Artur //need debug
void LpClpInterface::getBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
  int ncol = clpModel.getNumCols();
  int nrow = clpModel.getNumRows();
  CoinWarmStartBasis * basis = clpModel.getBasis();
  char * structural = basis->getStructuralStatus();
  char * artificial = basis->getArtificialStatus();	

  /* CLP
     isFree = 0x00,		///< Nonbasic free variable
     basic = 0x01,		///< Basic variable
     atUpperBound = 0x02,	///< Nonbasic at upper bound
     atLowerBound = 0x03		///< Nonbasic at lower bound
  */
  /* CPLEX 
     (Values of elements of cstat)
     CPX_AT_LOWER	0	variable at lower bound
     CPX_BASIC		1	variable is basic
     CPX_AT_UPPER	2	variable at upper bound
     CPX_FREE_SUPER	3	variable free and nonbasic

     (Values of elements of rstat in rows other than ranged rows)
     CPX_AT_LOWER	0	associated slack, surplus, or artificial variable is nonbasic at value 0.0 (zero)
     CPX_BASIC		1	associated slack, surplus, or artificial variable is basic
  */

	int n = 0;
	for (int ic = 0; ic<ncol; ic++)
	{
		const int st = (structural[ic>>2] >> ((ic&3)<<1)) & 3;

		if (st == 3)
			colStatus.push_back(0);
		else if (st == 1)
			colStatus.push_back(1); 
		else if (st == 2)
			colStatus.push_back(2);
		else if (st == 0)
			colStatus.push_back(3);
        else //else added by Issam.
        {
          exit(0);
        }
	}

	for (int ir = 0; ir<nrow; ir++)
	{
		const int st = (artificial[ir>>2] >> ((ir&3)<<1)) & 3;

        if (st == 3)
            rowStatus.push_back(0);
        else if (st == 1)
            rowStatus.push_back(1);
        else if (st == 2)
            rowStatus.push_back(2);
        else if (st == 0)
            rowStatus.push_back(3);
        else //else added by Issam.
        {
            exit(0);
        }

//        if (st == 1)
//			rowStatus.push_back(1);
//		//else if (st == 0) Commented by Issam (since it should be 3)
//        else if (st == 3)
//			rowStatus.push_back(0);
//		else
//		{
//			exit(0);
//		}
	}

  return;
}

// by Artur //need debug
void LpClpInterface::setBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
  int ncol = clpModel.getNumCols();
  int nrow = clpModel.getNumRows();

  for (int ic = 0; ic < ncol; ic++)
    {
      if (colStatus[ic] == 0)
	clpModel.setColumnStatus(ic, static_cast<ClpSimplex::Status>(3));
      else if (colStatus[ic] == 1)
	clpModel.setColumnStatus(ic, static_cast<ClpSimplex::Status>(1));
      else if (colStatus[ic] == 2)
	clpModel.setColumnStatus(ic, static_cast<ClpSimplex::Status>(2));
      else if (colStatus[ic] == 3)
	clpModel.setColumnStatus(ic, static_cast<ClpSimplex::Status>(0));
    }

  for (int ir = 0; ir < nrow; ir++)
    {
      if (rowStatus[ir] == 0)
	    clpModel.setRowStatus(ir, static_cast<ClpSimplex::Status>(3));
      else if (rowStatus[ir] == 1)
	    clpModel.setRowStatus(ir, static_cast<ClpSimplex::Status>(1));
      else if (rowStatus[ir] == 2)
        clpModel.setRowStatus(ir, static_cast<ClpSimplex::Status>(2));
      else if (rowStatus[ir] == 3)
          clpModel.setRowStatus(ir, static_cast<ClpSimplex::Status>(0));
    }

  return;
}

//OK
void LpClpInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  string probName, str;
  const vector<string> * colNames;
  const vector<string> * rowNames;
  int nRow, nCol, nchar;
  const double * obj;
  const double * colLower;
  const double * colUpper;
  const double * rowLower;
  const double * rowUpper;

  probName = clpModel.problemName();
  nRow = clpModel.getNumRows();
  nCol = clpModel.getNumCols();
  obj = clpModel.getObjCoefficients();
  colLower = clpModel.getColLower();
  colUpper = clpModel.getColUpper();
  rowLower = clpModel.getRowLower();
  rowUpper = clpModel.getRowUpper();

  if(param().MipSolverRecordNamesInFormulation)
    {
      colNames = clpModel.columnNames();
      rowNames = clpModel.rowNames();
    }

  os << "PROBLEM: "<< probName << " whose formulation ref number is " << ref() << endl;

  if ( minmaxStatus ==1 ) os << "Minimize" << endl;
  else os << "Maximize" << endl;

  nchar = printf(" Objectif:");
  for (int ic = 0; ic < nCol; ic++)
    {
      if (Double(obj[ic]) != 0)
	{
	  string s = (*colNames)[ic];
	  nchar += printf(" %+.9g %s", obj[ic], s.c_str());
	  if (nchar >= 80-9-13)
	    {
	      printf("\n");
	      nchar = 0;
	    }
	}
    }
  
  os << endl << "Subject To" << endl;

  CoinPackedMatrix * matrix = clpModel.matrix();
  
  for (int ir = 0; ir < nRow; ir++)
    {
      nchar = printf(" %s:", ((*rowNames)[ir]).c_str());
      for (int ic = 0; ic < nCol; ic++)
	{
	  double coef = matrix->getCoefficient(ir, ic);
	  if (Double(coef) != 0)
	    {
	      nchar += printf(" %+.9g %s", coef, ((*colNames)[ic]).c_str());
	      if (nchar >= 80-9-13) 
		{
		  printf("\n");
		  nchar = 0;
		}
	    }
	}
      if (rowUpper[ir] == rowLower[ir])
	printf( " %s %.9g", "=", rowUpper[ir]);
      else if (rowUpper[ir] < COIN_DBL_MAX)
	printf( " %s %.9g", "<=", rowUpper[ir]);
      else
	printf( " %s %.9g", ">=", rowLower[ir]);

      os << endl;
    }
  os << "Bounds" << endl;

  for (int ic = 0; ic < nCol; ic++)
    {
      os << " " << *colLower << " <= " << (*colNames)[ic] << " <= " << *colUpper << endl;
    }

  if (_pureLP)
    os << "END" << endl;
  else
    {
  
    }

  return;
}

//OK
void LpClpInterface::MPSwrite()
{
  clpModel.writeMps("curprob.mps");

  return;
}

//OK
void LpClpInterface::printRow(int irow,
			      const int & readNcol,
			      int *readMclind,
			      double *readDmatval,
			      char *readCnames)
{
  return;
}

//OK
void LpClpInterface::optimise(const int & minmaxStatus,
			      const double & BarrierConvergenceTolerance,
			      const double & rightHAndSideZeroTol,
			      const double & reducedCostTolerance,
			      const bool & preprocessorOn,
			      const bool & probingOn,
			      const bool & automaticCuttingPlanesOn,
			      const char & solverSelection)
{
  bapcodInit().require(_formCurrentlyLoaded, "Form not Currently Loaded",
		       ProgStatus::quit,
		       3);
  if (printL(8))
    MPSwrite();

  int status = clpModel.dual();

  return;
}


//OK
void LpClpInterface::makeSpaceForLoadingForm()
{
  _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << " _formCurrentlyLoaded = " << _formCurrentlyLoaded
	      << std::endl;
  return;
}

//OK
void LpClpInterface::saveCopyOfCurForm()
{
  if (printL(6))
    std::cout << "LpClpInterface::saveCopyOfCurForm(): formCurrentlyLoaded = "
	      << _formCurrentlyLoaded
	      << std::endl;
  return;
}

//OK
void LpClpInterface::load()
{
  if (printL(6))
    std::cout << "LpClpInterface::load(): formCurrentlyLoaded = " << _formCurrentlyLoaded
	      << std::endl;
  return;
}

//OK
void LpClpInterface::unload(const bool & deleteMat)
{
  return;
}

//OK
void LpClpInterface::reset()
{
  return;
}

void LpClpInterface::setPresolve(const bool& flag)
{
    throw MathProgSolverException("MipClpInit::setPresolve is not available", true);
}

// ************************************
// ***** Class MipClpInterface ******
// ************************************

MipClpInterface::MipClpInterface(BapcodInit* bapcodInit, const int& ref,
				 const std::string& name) : LpClpInterface(bapcodInit, ref, name)
{
#if _CBC_FOUND

  _pureLP = false;
  //bapcodInit->check(cplexProbPtr == NULL,"MipSolverInterface: problem must have been  created");
  
  //bapcodInit->check(CPXchgprobtype(cplexEnv, cplexProbPtr, CPXPROB_MILP),
  //                  " CPXchgprobtype: Failed to change typ\e to MILP.");
#else // _CBC_FOUND

#endif // _CBC_FOUND
}

void MipClpInterface::loadFormulation(const std::string& name,
				      const int& minmaxStatus, 
				      const int& ncol, 
				      const int& nrow,
				      const std::map<int, std::string>& mapSeqnb2Cname,
				      const std::map<int, std::string>& mapSeqnb2Rname,
				      const std::set<ProbCoef>& objectRow,
				      const std::set<ProbCoef, ProbCoefColSmallerThan>& colMatrix,
				      const std::set<ProbBound>& rhsv, 
				      const std::set<ProbBound>& bounds,
				      const std::set<ProbType>& types, 
				      const std::set<ProbIntC>& directs,
				      const std::set<ProbSetCoef>& sets)
{
#if _CBC_FOUND

#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);
#endif // _CBC_FOUND
}

void MipClpInterface::addDirectives(const std::set<ProbIntC>& varBrDirection)
{
  #if _CBC_FOUND

#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);
#endif // _CBC_FOUND
}

void MipClpInterface::chgColType(const std::set<ProbType>& newTypes)
{
#if _CBC_FOUND

#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);
#endif // _CBC_FOUND
}

void MipClpInterface::getObjVal(Double& objV)
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);
#endif // _CBC_FOUND
}

void MipClpInterface::getDualBound(Double& val)
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);  
#endif // _CBC_FOUND
}

void MipClpInterface::getPrimalBound(Double& val)
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);  
#endif // _CBC_FOUND
}

bool MipClpInterface::getOptimStatus(SolutionStatus& lpStatus,
				     SolutionStatus& mipStatus)
{
#if _CBC_FOUND

#else // _CBC_FOUND

#endif // _CBC_FOUND
  return true;
  //  throw MathProgSolverException("MipClpInterface is not available", true);
}

void MipClpInterface::optimise(const int& minmaxStatus,
			       const double & BarrierConvergenceTolerance,
			       const double & rightHAndSideZeroTol,
			       const double & reducedCostTolerance,			       
			       const bool& preprocessorOn,
			       const bool& probingOn,
			       const bool& automaticCuttingPlanesOn,
			       const char& solverSelection)
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);  
#endif // _CBC_FOUND
}

void MipClpInterface::getSol(std::map<int, Double>& primSol,
			     const bool& ifPrint)
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);  
#endif // _CBC_FOUND
}

void MipClpInterface::reset()
{
#if _CBC_FOUND
  
#else // _CBC_FOUND
  throw MathProgSolverException("MipClpInterface is not available", true);
#endif // _CBC_FOUND
}

#endif // CLP_FOUND

