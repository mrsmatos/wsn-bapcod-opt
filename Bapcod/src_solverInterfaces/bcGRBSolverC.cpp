/**
 *
 * This file bcGRBSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _GUROBI_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcGRBSolverC.hpp"
#include "bcMipSolverDef.hpp"
#include "bcPrintC.hpp"
#include "bcTimeC.hpp"

#include "bcMathProgSolverException.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 16;
const char endOfWord = '\0';
const std::string filling = "________________";


/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

using namespace std;


LpGRBInterface::LpGRBInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  MathProgSolverInterface(bapcodInit, ref, name)
{
  char *probname = new char[name.size() + 1];
  sprintf(probname, "%s", name.c_str());
  probname[name.size()] = endOfWord;

  if (GRBnewmodel(env, &model, probname, 0, NULL, NULL, NULL, NULL, NULL))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  delete [] probname;
  probname = NULL;
  
  return;
}


LpGRBInterface::~LpGRBInterface()
{
  if (GRBfreemodel(model))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  GRBfreeenv(env);

  return;
}


void LpGRBInterface::loadFormulation(const std::string & name,
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
  if (printL(3))
    std::cout << "LpGRBInterface::loadFormulation() ncol = " << ncol
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
      //rsense[rvPtr->ref] = rvPtr->sense;
	  if (rvPtr->sense == 'E')
		  rsense[rvPtr->ref] = '=';
	  else if (rvPtr->sense == 'L')
		  rsense[rvPtr->ref] = '<';
	  else if (rvPtr->sense == 'G')
		  rsense[rvPtr->ref] = '>';
	  else
	  {
		  printf("ERROR: %s\n", "the sense of constraint should be one of 'E','L' and 'G'");
		  exit(0);
	  }
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

  char *vtype = new char[ncol];
  fill_n(vtype, ncol, GRB_CONTINUOUS);

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


  if (GRBloadmodel(env, &model, "", ncol, nrow, minmaxStatus, 0.0, obj, rsense, 
							rhs, matbeg, matcnt, matind, matval, dlb, dub, vtype, colnames, rownames))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  return;
}


void LpGRBInterface::unLoadFormulation()
{
  return;
  _ncol = 0;
  _nrow = 0;

  return;
}




void LpGRBInterface::addCols(const std::set<ProbCoef> & objectiveRow,
			       const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
			       const std::set<ProbBound> & bounds,
			       const std::map<int, std::string> & mapSeqnb2Cname)
{
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  
  int readNcol(0);
  if (GRBgetintattr(model, "NumVars", &readNcol))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  bapcodInit().require(readNcol == _ncol,
	  "LpGRBInterface::addCols: readNcol != _ncol");

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

  char *vtype = new char[newcol];
  fill_n(vtype, newcol, GRB_CONTINUOUS);

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


  if (GRBaddvars(model, newcol, newnzInCols, ematbeg, ematind, ematval, eobjx, ebdl, ebdu, vtype, colnames))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
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


void LpGRBInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  
  int readNcol;
  if (GRBgetintattr(model, "NumVars", &readNcol))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  bapcodInit().require(readNcol <= _ncol, "LpClpInterface::delCols: readNcol > _ncol");
  bapcodInit().require(nbCol2Delete <= readNcol,
	  "LpGRBInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete];
  fill_n(dindex, nbCol2Delete, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
    {
      dindex[cnt++] = *iPtr;
    }

  if (GRBdelvars(model, cnt, dindex))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  delete[] dindex;
  dindex = NULL;
  _ncol -= nbCol2Delete;
  return;
}


void LpGRBInterface::addRows( const std::set<ProbBound> & rhsv,
				const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
				const std::map<int, std::string> & mapSeqnb2Rname)
{
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;

  int newnzInNewRows = rowMatrix.size();
  int readNrow;
  if (GRBgetintattr(model, "NumConstrs", &readNrow))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (printL(7))
    std::cout << "newrows = " << newrows
	      << "  newnzInNewRows = " << newnzInNewRows
	      << "  readNrow = " << readNrow
	      << "  _nrow = " << _nrow
	      << std::endl;

  bapcodInit().require(readNrow == _nrow,
	  "LpGRBInterface::addRows: readNrow != _nrow");

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
      //esense[rvPtr->ref - _nrow] = rvPtr->sense;
      erhs[rvPtr->ref - _nrow] = zero(rvPtr->bound);

	  if (rvPtr->sense == 'E')
		  esense[rvPtr->ref - _nrow] = '=';
	  else if (rvPtr->sense == 'L')
		  esense[rvPtr->ref - _nrow] = '<';
	  else if (rvPtr->sense == 'G')
		  esense[rvPtr->ref - _nrow] = '>';
	  else
	  {
		  printf("ERROR: %s\n", "the sense of constraint should be one of 'E','L' and 'G'");
		  exit(0);
	  }
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

  if (GRBaddconstrs(model, newrows, newnzInNewRows, 
					ematbeg, ematind, ematval, esense, erhs, rownames))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

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


void LpGRBInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  
  int readNrow;
  if (GRBgetintattr(model, "NumConstrs", &readNrow))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  bapcodInit().require(readNrow <= _nrow,
	  "LpGRBInterface::delRowss: readNrow > _nrow");
  bapcodInit().require(nbRow2Delete <= readNrow,
	  "LpGRBInterface::delRows: nbRow2Delete > readNrow");

  int* dindex = new int[nbRow2Delete];
  std::copy(indexSetOfRow2Delete.begin(), indexSetOfRow2Delete.end(), dindex);
 
  if (GRBdelconstrs(model, nbRow2Delete, dindex))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  delete[] dindex;
  dindex = NULL;
  _nrow -= nbRow2Delete;

  return;
}

//OK
void LpGRBInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  bapcodInit().check(1, "LP can not have directives");
  return;
}


void LpGRBInterface::getObjCoef(ProbCoef & pc)
{
  double *coef = new double[1];

  if (GRBgetdblattrelement(model, "Obj", pc.colRef, coef))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  pc.coef = *coef;
  delete [] coef;
  return;
}


void LpGRBInterface::chgObjCoef(const ProbCoef & pc)
{
  if (GRBsetdblattrelement(model, "Obj", pc.colRef, pc.coef))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  return;
}


void LpGRBInterface::chgMatCoef(const ProbCoef & pc)
{
  int *cind = new int[1];
  int *vind = new int[1];
  double *val = new double[1];
  cind[0] = pc.rowRef;
  vind[0] = pc.colRef;
  val[0] = double(pc.coef);

  if (GRBchgcoeffs(model, 1, cind, vind, val))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  return;
}


void LpGRBInterface::chgRhs(const ProbBound & pb)
{
  if (GRBsetdblattrelement(model, "RHS", pb.ref, zero(pb.bound)))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  return;
}


void LpGRBInterface::chgRhs(const std::set<ProbBound> & newRhs)
{
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

  if (GRBsetdblattrlist(model, "RHS", nrhs, index, rhs))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  delete [] index;
  index = NULL;
  delete [] rhs;
  rhs = NULL;
  return;
}


void LpGRBInterface::chgBds(const std::set<ProbBound> & newBounds)
{
	for (set<ProbBound>::const_iterator bPtr = newBounds.begin();
		bPtr != newBounds.end();
		bPtr++)
	{
		if (bPtr->sense == 'L')
			GRBsetdblattrelement(model, "LB", bPtr->ref, zero(bPtr->bound));
        else if ((bPtr->sense == 'U')
	       || (bPtr->sense == 'I')
	       || (bPtr->sense == 'B'))
			GRBsetdblattrelement(model, "UB", bPtr->ref, zero(bPtr->bound));
		else if (bPtr->sense == 'F')
		{
			GRBsetdblattrelement(model, "LB", bPtr->ref, zero(bPtr->bound));
			GRBsetdblattrelement(model, "UB", bPtr->ref, zero(bPtr->bound));
		}
	}

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  return;
}


void LpGRBInterface::chgColType(const std::set<ProbType> & newTypes)
{
  bapcodInit().check(1,"LP Form cannot have ColType");
  return;
}


void LpGRBInterface::getObjVal(Double & objV)
{
  double objVal;
  if (GRBgetdblattr(model, "ObjVal", &objVal))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  objV = zero(objVal, param().HIGHPRECISION);

  return;
}


void LpGRBInterface::getDualBound(Double & val)
{
  int nrow;
  if (GRBgetintattr(model, "NumConstrs", &nrow))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  double * dualsol = new double[nrow];
  if(GRBgetdblattrarray(model, "Pi", 0, nrow, dualsol))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  double * rhs = new double[nrow];
  if(GRBgetdblattrarray(model, "RHS", 0, nrow, rhs))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  double dualval = 0.0;
  for(int i = 0; i < nrow; i++)
  dualval = dualval + dualsol[i] * rhs[i];
  val = zero(dualval, param().HIGHPRECISION);

  return;
}


void LpGRBInterface::getPrimalBound(Double & val)
{
  double primalVal = 0.0;
  if (GRBgetdblattr(model, "ObjVal", &primalVal))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  val = zero(primalVal, param().HIGHPRECISION);

  return;
}


bool LpGRBInterface::getOptimStatus(SolutionStatus & lpStatus,
				      SolutionStatus & mipStatus)
{
  int status;
  if (GRBgetintattr(model, "Status", &status))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (status == GRB_OPTIMAL)
  {
	  lpStatus = SolutionStatus::Optimum;
      mipStatus = SolutionStatus::Optimum;
	  
	  if (printL(3))
		std::cout << "LpGRBInterface::getOptimStatus: LP optimal"
		  << std::endl;

	  return(true);
  }

  if (printL(4))
    MPSwrite();

  if (status == GRB_INFEASIBLE)
  {
	lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: LP infeasible"
		  << std::endl;

      return(false);
  }

  if (status == GRB_UNBOUNDED)
  {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: LP unbounded"
		  << std::endl;

      return(false);	
  }

  if (status == GRB_INF_OR_UNBD)
  {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: Model was proven to be either infeasible or unbounded"
		  << std::endl;

      return(false);	
  }

  if (status == GRB_LOADED)
  {
	 lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: Model is loaded, but no solution information is available"
		  << std::endl;

      return(false);
  }

  if (status == GRB_ITERATION_LIMIT)
  {
	 lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: the total number of simplex iterations performed exceeded"
		  << std::endl;

      return(false);
  }

  if (status == GRB_TIME_LIMIT)
  {
	 lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: the time expended exceeded"
		  << std::endl;

      return(false);
  }

  if (status == GRB_INTERRUPTED)
  {
	 lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: Optimization was terminated by the user"
		  << std::endl;

      return(false);
  }

  if (status == GRB_NUMERIC)
  {
	 lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGRBInterface::getOptimStatus: Optimization was terminated due to unrecoverable numerical difficulties"
		  << std::endl;

      return(false);
  }

  std::cout << "LpCplexInterface::getOptimStatus: undefined status"
	    << std::endl;

  return(false);
 
}


void LpGRBInterface::getSol(std::map<int, Double> & primSol,
			      const bool & ifPrint)
{
  primSol.clear();
  int readNcol = 0;
  if (GRBgetintattr(model, "NumVars", &readNcol))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  bapcodInit().require(readNcol <= _ncol,
	"LpGRBInterface::getSol: readNcol > _ncol");

  double *x = new double[_ncol];
  fill_n(x, _ncol, 0);

  if (GRBgetdblattrarray(model, "X", 0, _ncol, x))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

    if (printL(6))
    std::cout << "readNcol = " << readNcol
	      << "  _ncol = " << _ncol
	      << std::endl;

  for (int jcol = 0; jcol < _ncol; jcol++)
    if (zero(x[jcol], param().HIGHPRECISION) != 0)
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


void LpGRBInterface::getSol(std::map<int, Double> & primSol,
			      std::map<int, Double> & dualSol,
			      const int & minmaxStatus,
			      const bool & ifPrint)
{
  getSol(primSol, ifPrint);
  dualSol.clear();
  int readNrow(0);
  if (GRBgetintattr(model, "NumConstrs", &readNrow))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  bapcodInit().require(readNrow <= _nrow,
	  "LpGRBInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[_nrow];
  fill_n(dual, _nrow, 0);

  char *rsense = new char[_nrow];
  fill_n(rsense, _nrow, ' ');

  if(GRBgetdblattrarray(model, "Pi", 0, _nrow, dual))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if(GRBgetcharattrarray(model, "Sense", 0, _nrow, rsense))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

    for (int irow = 0; irow < _nrow; irow++)
    if (zero(dual[irow], param().HIGHPRECISION) != 0)
      {
	if (rsense[irow] == '<')
	  dualSol[irow] = minmaxStatus * flipsign * dual[irow];
	else if (rsense[irow] == '>')
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

void LpGRBInterface::getReducedCost(std::map<int, Double> & redCost, 
				  const bool & ifPrint)
{

}

void LpGRBInterface::getBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
	// GRB
	// variable status
	// 0 (basic), -1 (non-basic at lower bound), -2 (non-basic at upper bound), and -3 (super-basic)
	// constraint status
	// 0 (basic) or -1 (non-basic), A constraint is basic when its slack variable is in the simplex basis
  
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
	 
	
  int readNrow(0);
  int readNcol(0);
  
  if (GRBgetintattr(model, "NumConstrs", &readNrow))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
    
  if (GRBgetintattr(model, "NumVars", &readNcol))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  int *vbasis = new int[readNcol];
  int *cbasis = new int[readNrow];

  GRBgetintattrarray(model, "VBasis", 0, readNcol, vbasis);
  GRBgetintattrarray(model, "CBasis", 0, readNrow, cbasis);
	
  for (int i=0; i<readNcol; i++)
  {
	if (vbasis[i] == 0)
		colStatus.push_back(1);
	else if (vbasis[i] == -1)
		colStatus.push_back(0);
	else if (vbasis[i] == -2)
		colStatus.push_back(2);
	else if (vbasis[i] == -3)
		colStatus.push_back(3);
  }

  for (int i=0; i<readNrow; i++)
  {
  	if (cbasis[i] == 0)
		rowStatus.push_back(1);
	else if (cbasis[i] == -1)
		rowStatus.push_back(0);
  }

	return;
}


void LpGRBInterface::setBasis(std::vector<int> & colStatus, std::vector<int> & rowStatus)
{
  int readNrow(0);
  int readNcol(0);
  
  if (GRBgetintattr(model, "NumConstrs", &readNrow))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
    
  if (GRBgetintattr(model, "NumVars", &readNcol))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  	for (int ic = 0; ic < readNcol; ic++)
	{
		if (colStatus[ic] == 0)
			GRBsetintattrelement(model, "VBasis", ic, -1);
		else if (colStatus[ic] == 1)
			GRBsetintattrelement(model, "VBasis", ic, 0);
		else if (colStatus[ic] == 2)
			GRBsetintattrelement(model, "VBasis", ic, -2);
		else if (colStatus[ic] == 3)
			GRBsetintattrelement(model, "VBasis", ic, -3);
	}

	for (int ir = 0; ir < readNrow; ir++)
	{
		if (rowStatus[ir] == 0)
			GRBsetintattrelement(model, "CBasis", ir, -1);
		if (rowStatus[ir] == 1)
			GRBsetintattrelement(model, "CBasis", ir, 0);
	}

  if (GRBupdatemodel(model))
  {
  	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

	return;
}


void LpGRBInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  int nRow, nCol, nchar;
  GRBgetintattr(model, "NumVars", &nCol);
  GRBgetintattr(model, "NumConstrs", &nRow);

  char *probName;//, str; //Commented by Romain: variable not used.
  double * obj = new double[nCol];
  double *colLower = new double[nCol];
  double *colUpper = new double[nCol];
  double *rhs = new double[nRow];
  char *sense = new char[nRow];
  char **colNames = new char*[nCol];
  char **rowNames = new char*[nRow];


  GRBgetstrattr(model, "ModelName", &probName);
  GRBgetdblattrarray(model, "Obj", 0, nCol, obj);
  GRBgetdblattrarray(model, "LB", 0, nCol, colLower);
  GRBgetdblattrarray(model, "UB", 0, nCol, colUpper);
  GRBgetdblattrarray(model, "RHS", 0, nRow, rhs);
  GRBgetcharattrarray(model, "Sense", 0, nRow, sense);

  if(param().MipSolverRecordNamesInFormulation)
  {
	  GRBgetstrattrarray(model, "VarName", 0, nCol, colNames);
	  GRBgetstrattrarray(model, "ConstrName", 0, nRow, rowNames);
  }

  os << "PROBLEM: "<< probName << " whose formulation ref number is " << ref() << endl;

  if ( minmaxStatus ==1 ) os << "Minimize" << endl;
  else os << "Maximize" << endl;

  nchar = printf(" Objective:");
  for (int ic = 0; ic < nCol; ic++)
  {
	if (Double(obj[ic]) != 0)
	{
		string s = colNames[ic];
		nchar += printf(" %+.9g %s", obj[ic], s.c_str());
		if (nchar >= 80-9-13)
		{
			printf("\n");
			nchar = 0;
		}
	}
  }

  os << endl << "Subject To" << endl;

  for (int ir=0; ir<nRow; ir++)
  {
	nchar = printf(" %s:", rowNames[ir]);
	for (int ic=0; ic<nCol; ic++)
	{
		double coef;
		GRBgetcoeff(model, ir, ic, &coef);
		
		if (Double(coef) != 0)
		{
			nchar += printf(" %+.9g %s", coef, colNames[ic]);
			if (nchar >= 80-9-13) 
			{
			  printf("\n");
			  nchar = 0;
			}
		}
	}
	if (sense[ir] == '=')
		printf( " %s %.9g", "=", rhs[ir]);
	else if (sense[ir] == '>')
		printf( " %s %.9g", ">=", rhs[ir]);
	else
		printf( " %s %.9g", "<=", rhs[ir]);

	os << endl;
  }

  os << "Bounds" << endl;

  for (int ic=0; ic<nCol; ic++)
  {
    os << " " << colLower[ic] << " <= " << colNames[ic] << " <= " << colUpper[ic] << endl;
  }

  if (_pureLP)
	  os << "END" << endl;
  else
  {
  
  }

  return;
}


void LpGRBInterface::MPSwrite()
{
  if (GRBwrite(model, "curprob.lp"))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (GRBwrite(model, "curprob.mps"))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  return;
}


void LpGRBInterface::printRow(int irow,
				const int & readNcol,
				int *readMclind,
				double *readDmatval,
				char *readCnames)
{
  return;
}


void LpGRBInterface::optimise(const int & minmaxStatus,
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

  if (GRBoptimize(model))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }	

  return;
}



void LpGRBInterface::makeSpaceForLoadingForm()
{
  /// Save loaded form if any
  _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << " _formCurrentlyLoaded = " << _formCurrentlyLoaded
	      << std::endl;
  return;
}


void LpGRBInterface::saveCopyOfCurForm()
{
  /// Save loaded form if any
  // _formCurrentlyLoaded = true;
  if (printL(6))
    std::cout << "LpGRBInterface::saveCopyOfCurForm(): formCurrentlyLoaded = "
	      << _formCurrentlyLoaded
	      << std::endl;
  return;
}


void LpGRBInterface::load()
{
  if (printL(6))
    std::cout << "LpGRBInterface::load(): formCurrentlyLoaded = " << _formCurrentlyLoaded
	      << std::endl;
  return;
}


void LpGRBInterface::unload(const bool & deleteMat)
{

  return;
}


void LpGRBInterface::reset()
{
  return;
}

MipGRBInit::MipGRBInit(BapcodInit* bapcodInit) : MathProgSolverInit(bapcodInit)
{
  if (GRBloadenv(&env, "")) 
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if (printL(1))
  {
	int major, minor, technical;
	GRBversion(&major, &minor, &technical);
	printf("Gurobi library version %d.%d.%d\n", major, minor, technical);

	if (GRBsetintparam(env, "LogToConsole", 1))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}
  }
  else if (GRBsetintparam(env, "LogToConsole", 0))
  {
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  return;
}

MipGRBInit::~MipGRBInit()
{
  return;
}

void MipGRBInit::setPresolve(const bool& flag)
{
	// Controls the presolve level. 
	// A value of -1 corresponds to an automatic setting. 
	// Other options are off (0), conservative (1), or aggressive (2).
	if (flag)
	{
		if (GRBsetintparam(env, "Presolve", -1))
		{
			printf("ERROR: %s\n", GRBgeterrormsg(env));
			exit(0);
		}
	}
	else
	{
		if (GRBsetintparam(env, "Presolve", 0))
		{
			printf("ERROR: %s\n", GRBgeterrormsg(env));
			exit(0);
		}
	}

	return;
}

void MipGRBInit::setTimeLimit(const double& seconds)
{
	if (GRBsetdblparam(env, "TimeLimit", seconds))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}

	return;
}

void MipGRBInit::setMultiThread(const int& maxNbThread)
{
	if (GRBsetintparam(env, "Threads", maxNbThread))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}
	
	return;
}

void MipGRBInit::setRelativeMipGapLimit(const double& relativeGap)
{
	if (GRBsetdblparam(env, "MIPGap", relativeGap))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}
	
	return;
}

void MipGRBInit::setSearchPriority(const int& flag)
{
  switch (flag)
    {
    case 0 :
      {
		//default, balanced
		GRBsetintparam(env, "MIPFocus", 0);
		break;
      }
    case 1 :
      {
		//If you are more interested in finding feasible solution quickly
		GRBsetintparam(env, "MIPFocus", 1);
		break;
      }
    case 2 :
      {
		//if you wish to focus more attention on proving optimality
		GRBsetintparam(env, "MIPFocus", 2);
		break;
      }
    case 3 :
      {
		//If the best objective bound is moving very slowly, try to focus on the bound
		GRBsetintparam(env, "MIPFocus", 3);
		break;
      }
    case 4 :
      {
		throw MathProgSolverException("failure to set mipemphasis indicator 4", true);
		break;
      }
    }

  return;
}

void MipGRBInit::setWorkingMemorySpace(const double& sizeInMB)
{
	// measured in GBytes
	if (GRBsetdblparam(env, "NodefileStart", sizeInMB/1000))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}
	
	return;
}

void MipGRBInit::setLPoptimalityTolerance(const double& tolerance)
{
	if (GRBsetdblparam(env, "OptimalityTol", tolerance))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}
	
	return;
}

void MipGRBInit::setLPfeasibilityTolerance(const double& tolerance)
{
  throw MathProgSolverException("MipGRBInit::setLPfeasibilityTolerance is not available", true);
}

// ************************************
// ***** Class MipGRBInterface ******
// ************************************

MipGRBInterface::MipGRBInterface(BapcodInit* bapcodInit, const int& ref,
    const std::string& name) : LpGRBInterface(bapcodInit, ref, name)
{
  _pureLP = false;
  return;
}

void MipGRBInterface::loadFormulation(const std::string& name,
    const int& minmaxStatus, const int& ncol, const int& nrow,
    const std::map<int, std::string>& mapSeqnb2Cname,
    const std::map<int, std::string>& mapSeqnb2Rname,
    const std::set<ProbCoef>& objectRow,
    const std::set<ProbCoef, ProbCoefColSmallerThan>& colMatrix,
    const std::set<ProbBound>& rhsv, const std::set<ProbBound>& bounds,
    const std::set<ProbType>& types, const std::set<ProbIntC>& directs,
    const std::set<ProbSetCoef>& sets)
{
  if (printL(3))
    std::cout << "LpGRBInterface::loadFormulation() ncol = " << ncol
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
      //rsense[rvPtr->ref] = rvPtr->sense;
	  if (rvPtr->sense == 'E')
		  rsense[rvPtr->ref] = '=';
	  else if (rvPtr->sense == 'L')
		  rsense[rvPtr->ref] = '<';
	  else if (rvPtr->sense == 'G')
		  rsense[rvPtr->ref] = '>';
	  else
	  {
		  printf("ERROR: %s\n", "the sense of constraint should be one of 'E','L' and 'G'");
		  exit(0);
	  }
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

  char *vtype = new char[ncol];
  fill_n(vtype, ncol, GRB_CONTINUOUS);

   for (set<ProbType>::const_iterator typePt = types.begin();
       typePt != types.end();
       typePt++)
    {
      if (typePt->type == 'I')
		vtype[typePt->ref] = GRB_INTEGER;

      if (typePt->type == 'B')
		vtype[typePt->ref] = GRB_BINARY;

	  if (typePt->type == 'S')
		vtype[typePt->ref] = 'S';

      if (typePt->type == 'N')
		vtype[typePt->ref] = 'N';
    }

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
  
  if(param().MipSolverRecordNamesInFormulation)
    {

	  if (GRBloadmodel(env, &model, "", ncol, nrow, minmaxStatus, 0.0, obj, rsense, 
								rhs, matbeg, matcnt, matind, matval, dlb, dub, vtype, colnames, rownames))
	  {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	  }
  }
  else
  {
 	  if (GRBloadmodel(env, &model, "", ncol, nrow, minmaxStatus, 0.0, obj, rsense, 
								rhs, matbeg, matcnt, matind, matval, dlb, dub, vtype, NULL, NULL))
	  {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	  } 
  }

  //sos set
  if (!sets.empty())
    {
      int numsos = 0;
      int numsosnz = 0;
      int *sostype = NULL;
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
      sostype = new int[numsos];
      fill_n(sostype, numsos, 0);
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
		  if (sPtr->setType == '1')
			  sostype[curSet] = GRB_SOS_TYPE1;
		  else if (sPtr->setType == '2')
			  sostype[curSet] = GRB_SOS_TYPE2;
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

	  if(GRBaddsos(model, numsos, numsosnz, sostype, sosbeg, sosind, sosref))
	  {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	  } 

	  if (GRBupdatemodel(model))
	  {
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	  }

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

void MipGRBInterface::addDirectives(const std::set<ProbIntC>& varBrDirection)
{
	for (set<ProbIntC>::const_iterator dPtr = varBrDirection.begin(); 
   	dPtr != varBrDirection.end(); dPtr++)
    {
        if (printL(6)) 
   		std::cout << *dPtr;

		GRBsetintattrelement(model, "BranchPriority", dPtr->ref, (int) dPtr->val);

		throw MathProgSolverException("GRB does not allow to assign branch direction for each variable", true);
    }
}

void MipGRBInterface::chgColType(const std::set<ProbType>& newTypes)
{
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

	if(GRBgetcharattrlist(model, "VType", Lnels, Lmindex, Lqctype))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);	
	}

	if (GRBupdatemodel(model))
	{
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}

  /// Delete extra arrays
  delete[] Lmindex;
  Lmindex = NULL;
  delete[] Lqctype;
  Lqctype = NULL;

  return;
}

void MipGRBInterface::getObjVal(Double& objV)
{
  double objVal;
  if (GRBgetdblattr(model, "ObjVal", &objVal))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }
  
  objV = zero(objVal, param().HIGHPRECISION);

  return;
}

void MipGRBInterface::getDualBound(Double& val)
{
  double dualval = 0.0;
  int status;
  if (GRBgetintattr(model, "Status", &status))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

  if(status == 2) //optimal
    {
		bapcodInit().require(1, "MIP optimal");
	  if (GRBgetdblattr(model, "ObjVal", &dualval))
	  {  	
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	  }
    }
  else
    dualval = -BapcodInfinity;

  val = zero(dualval, param().HIGHPRECISION);

  return;
}

void MipGRBInterface::getPrimalBound(Double& val)
{
  double primalval = 0.0;
	if (GRBgetdblattr(model, "ObjVal", &primalval))
	{  	
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}

  val = zero(primalval, param().HIGHPRECISION);
  return;
}

bool MipGRBInterface::getOptimStatus(SolutionStatus& lpStatus,
    SolutionStatus& mipStatus)
{
  int status;
  if (GRBgetintattr(model, "Status", &status))
  {  	
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(0);
  }

 if (status == GRB_OPTIMAL)
 {
	mipStatus = SolutionStatus::Optimum;
    lpStatus = SolutionStatus::Optimum;
	if (printL(3))
	std::cout << "MipCplexInterface::getOptimStatus: global search complete and an integer solution has been found"
		  << std::endl;
	return true;
 }

   if (printL(4))
    MPSwrite();

   if (status == GRB_INFEASIBLE)
   {
      mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "MipCplexInterface::getOptimStatus: MIP infeasible"
		  << std::endl;

      return(false);   
   }

   if (status == GRB_INF_OR_UNBD)
   {
	  mipStatus = SolutionStatus::Infeasible;
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "MipCplexInterface::getOptimStatus: MIP infeasible"
		  << std::endl;

      return(false); 
   }

   if (status == GRB_CUTOFF)
   {
	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: cut off" << std::endl;
			 return(false); 
   }

   if (status == GRB_ITERATION_LIMIT)
   {
   	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: iteration limit";
	return(false); 
   }

   if (status == GRB_NODE_LIMIT)
   {
	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: node limit";
	return(false); 
   }

   if (status == GRB_TIME_LIMIT)
   {
	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: time limit";
	return(false); 
   }

   if (status == GRB_SOLUTION_LIMIT)
   {
	  mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: solution limit";
	return(false); 
   }

   if (status == GRB_INTERRUPTED)
   {
   	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: interrupted by the user";
	return(false); 
   }

    if (status == GRB_NUMERIC)
   {
   	  mipStatus = SolutionStatus::UnSolved;
      lpStatus = SolutionStatus::UnSolved;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: unrecoverable numerical difficulty";
	return(false); 
   }

	if (status == GRB_SUBOPTIMAL)
   {
   	  mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
	        if (printL(3))
	std::cout << 
	"MipCplexInterface::getOptimStatus: suboptimal solution is available";
	return(false); 
   }

	return LpGRBInterface::getOptimStatus(lpStatus, mipStatus);
}

void MipGRBInterface::optimise(const int& minmaxStatus,
    const bool& preprocessorOn, const bool& probingOn,
    const bool& automaticCuttingPlanesOn, const char& solverSelection)
{
  bapcodInit().require(_formCurrentlyLoaded,
	  "Form not Currently Loaded/2",
	  ProgStatus::quit,
	  3);
  if (printL(8))
    MPSwrite();

  if(printL(1))
	  GRBwrite(model, "myprob.mps");

  if (GRBoptimize(model))
	{  	
		printf("ERROR: %s\n", GRBgeterrormsg(env));
		exit(0);
	}

  return;
	
}

void MipGRBInterface::getSol(std::map<int, Double>& primSol,
    const bool& ifPrint)
{
  primSol.clear();
  LpGRBInterface::getSol(primSol, ifPrint);
  return;
}

void MipGRBInterface::reset()
{
  return;
}

#endif
