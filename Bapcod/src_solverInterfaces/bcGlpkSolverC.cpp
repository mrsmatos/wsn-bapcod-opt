/**
 *
 * This file bcGlpkSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if _GLPK_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcMipSolverDef.hpp"
#include "bcGlpkSolverC.hpp"
#include "bcPrintC.hpp"
#include "bcSwitches.h"
#include "bcTimeC.hpp"

/// Fixed length of names given to sub-solvers
const int wordSize = 16;
const char endOfWord = '\0';
const std::string filling = "________________";

//initialization of the counters
#if _GLPK_FOUND
int LpGlpkInterface::GLPKLpCounter = 1;
int LpGlpkInterface::GLPKMipCounter = 1;
#endif //_GLPK_FOUND

/// WordSize + a space for the end of word character 
const int placeSize = wordSize+ 1;

using namespace std;

static int hook(void *info = NULL, const char *s = NULL)
{
  //printout("hook");
  FILE *foo = fopen("GLPKexit","a");
  fputs(s, foo);
  fclose(foo);

  return(1);
}

MipGlpkInit::MipGlpkInit(BapcodInit* bapcodInit) : MathProgSolverInit(bapcodInit)
{
  //printout("MipGlpkInit::MipGlpkInit()");
  /// Xpress
  cout<<"GLPK solver : GLPK version "<< glp_version() <<endl;

  /// Clear GLPKexit
  if(printL(5)){
    FILE *foo = fopen("GLPKexit","w");
    fputs("////INIT SOLVER/////\n", foo);
    fclose(foo);
  }
  else{
    FILE *foo = fopen("GLPKexit","w");
    fclose(foo);
  }

  /// Init glpk hook
  if(printL(6)){
    FILE *file = NULL;
    glp_term_hook(hook, file);
  }
  else
    glp_term_hook(hook, NULL);
  // glp_term_hook(NULL, NULL);

  return;
}

/// New functionalities to add to mip inteface: begin
void MipGlpkInit::setPresolve(const bool & flag)
{
  return;
}

void MipGlpkInit::setTimeLimit(const double & seconds)
{
  return;
}

void MipGlpkInit::setMultiThread(const int & maxNbThread)
{
  return;
}


void MipGlpkInit::setSearchPriority(const int & flag)
{
  return;
}


void MipGlpkInit::setRelativeMipGapLimit(const double & relativeGap)
{
  return;
}

void MipGlpkInit::setWorkingMemorySpace(const double & sizeInMB)
{
  return;
}

/// New functionalities to add to mip inteface: end
MipGlpkInit::~MipGlpkInit()
{
  //printout("MipGlpkInit::~MipGlpkInit()");
  return;
}

LpGlpkInterface::LpGlpkInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  MathProgSolverInterface(bapcodInit, ref, name),
  glpkProbPtr(NULL)
{
  //printout("LpGlpkInterface::LpGlpkInterface()", 6);
  glpkProbPtr = glp_create_prob();
  glp_init_smcp(&glpkParamPtr);
  // glpkParamPtr.tm_lim = param.MipSolverMaxTime*1000;
  // glpkParamPtr.it_lim = param.MipSolverMaxNbLpIterations;

  return;
}

LpGlpkInterface::~LpGlpkInterface()
{
  //printout("LpGlpkInterface::~LpGlpkInterface");
  glp_delete_prob(glpkProbPtr);

  return;
}

MipGlpkInterface::MipGlpkInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name):
  LpGlpkInterface(bapcodInit, ref, name)
{
  //printout("MipGlpkInterface::MipGlpkInterface()");
  _pureLP = false;
  glp_init_iocp(&glpkParamMipPtr);
  // glpkParamMipPtr.tm_lim = param.MipSolverMaxBBNodes;
  glpkParamMipPtr.tm_lim = param.MipSolverMaxTime*1000;

  /// 0 = off
  // lpx_set_int_parm(glpkProbPtr,
  //LPX_K_PRESOL,
  //  0);
  lpx_set_int_parm(glpkProbPtr,
		   LPX_K_USECUTS,
		   0);
  lpx_set_int_parm(glpkProbPtr,
		   LPX_K_SCALE,
		   0);


  return;
}

void LpGlpkInterface::loadFormulation(const std::string & name,
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
  //printout("LpGlpkInterface::loadFormulation");
  if (printL(3))
    std::cout << "LpGlpkInterface::loadFormulation() ncol = " << ncol
	      << " , nrow = " << nrow
	      << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  char *probname= new char[name.size() + 1];
  sprintf(probname, "%s", name.c_str());
  probname[name.size()] = endOfWord;

  glp_set_prob_name(glpkProbPtr, probname);
  glp_set_obj_name(glpkProbPtr,"Objectif");
  glp_add_rows(glpkProbPtr, nrow);
  glp_add_cols(glpkProbPtr, ncol);

  std::string curName;
  char * place;

  char *rnames = new char[placeSize * nrow];
  fill_n(rnames, placeSize * nrow, endOfWord);

  for(int i = 0; i < nrow; i++){
    place = & (rnames[placeSize * i]);
    curName = (mapSeqnb2Rname^i) + filling;
    strncpy(place, curName.c_str(), wordSize);
    place[wordSize] = endOfWord;
  }

  if (printL(7))
    for(int i = 0; i < nrow; i++)
      printf("rname %s\n", &rnames[placeSize * i]);

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
    for(int i = 0; i < ncol; i++)
      printf("cname %s\n", &cnames[placeSize * i]);

  for(int i = 0; i < nrow; i++)
    glp_set_row_name(glpkProbPtr,
		     i + 1,
		     &rnames[placeSize * i]);

  for(int i = 0;i<ncol;i++)
    glp_set_col_name(glpkProbPtr,
		     i + 1,
		     &cnames[placeSize * i]);

  if(minmaxStatus==-1)
    glp_set_obj_dir(glpkProbPtr, GLP_MAX);
  else if(minmaxStatus==1)
    glp_set_obj_dir(glpkProbPtr, GLP_MIN);

  for (set<ProbCoef>::const_iterator oPtr = objectRow.begin();
       oPtr != objectRow.end();
       oPtr++)
    glp_set_obj_coef(glpkProbPtr,
		     oPtr->colRef + 1,
		     oPtr->coef);

  double *rhsub = new double[nrow];
  fill_n(rhsub, nrow, BapcodInfinity*10);
  double *rhslb = new double[nrow];
  fill_n(rhslb, nrow,-BapcodInfinity*10);
  char *rhss= new char[nrow];
  fill_n(rhss, nrow,' ');

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      if(rvPtr->sense=='G')
	{
	  rhslb[rvPtr->ref]=zero(rvPtr->bound);
	  //rhss[rvPtr->ref]='G';
	}
      else if(rvPtr->sense=='L')
	{
	  rhsub[rvPtr->ref]=zero(rvPtr->bound);
	  //rhss[rvPtr->ref]='L';
	}
      else if(rvPtr->sense=='E')
	{
	  rhsub[rvPtr->ref]=zero(rvPtr->bound);
	  rhslb[rvPtr->ref]=zero(rvPtr->bound);
	  rhss[rvPtr->ref]='E';
	}
    }

  for(int id = 0;id<nrow;id++)
    if(rhss[id]=='E')
      glp_set_row_bnds(glpkProbPtr,
		       id + 1,
		       GLP_FX,
		       rhslb[id],
		       rhsub[id]);
    else if(rhsub[id]>(BapcodInfinity))
      if(rhslb[id]<(-BapcodInfinity))
	glp_set_row_bnds(glpkProbPtr,
			 id + 1,
			 GLP_FR,
			 rhslb[id],
			 rhsub[id]);
      else
	glp_set_row_bnds(glpkProbPtr,
			 id + 1,
			 GLP_LO,
			 rhslb[id],
			 rhsub[id]);
    else
      if(rhslb[id]<(-BapcodInfinity))
	glp_set_row_bnds(glpkProbPtr,
			 id + 1,
			 GLP_UP,
			 rhslb[id],
			 rhsub[id]);
      else
	glp_set_row_bnds(glpkProbPtr,
			 id + 1,
			 GLP_DB,
			 rhslb[id],
			 rhsub[id]);

  int nnze = colMatrix.size();
  int *matrow = new int[nnze + 1];
  fill_n(matrow, nnze + 1, 0);
  int *matcol = new int[nnze + 1];
  fill_n(matcol, nnze + 1, 0);
  double *matval = new double[nnze + 1];
  fill_n(matval, nnze + 1, 0);
  int cnt(1);

  for(set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin() ;
      mPtr!=colMatrix.end() ;
      mPtr++)
    {
      matrow[cnt]=mPtr->rowRef + 1;;
      matcol[cnt]=mPtr->colRef + 1;
      matval[cnt]=zero(mPtr->coef);
      cnt++;
    }

  glp_load_matrix(glpkProbPtr,
		  nnze,
		  matrow,
		  matcol,
		  matval);

  double *dub = new double[ncol];
  fill_n(dub, ncol, BapcodInfinity*10);
  double *dlb = new double[ncol];
  fill_n(dlb, ncol, 0);
  char *rsense= new char[ncol];
  fill_n(rsense, ncol,' ');

  for (set<ProbBound>::const_iterator bPtr = bounds.begin();
       bPtr != bounds.end();
       bPtr++)
    {
      if(bPtr->sense=='U')
	{
	  dub[bPtr->ref]=zero(bPtr->bound);
	  //rsense[bPtr->ref]='U'; //@todo fv to verify 
	}
      else if(bPtr->sense=='L')
	{
	  dlb[bPtr->ref]=zero(bPtr->bound);
	  //rsense[bPtr->ref]='L'; //@todo  fv to verify 
	}
      else if(bPtr->sense=='F')
	{
	  dub[bPtr->ref]=zero(bPtr->bound);
	  dlb[bPtr->ref]=zero(bPtr->bound);
	  rsense[bPtr->ref]='F';
	  require(0,"LpGlpkInterface::loadFormulation() Variable free to verify");
	}
    }

  for(int id = 0;id<ncol;id++)
    if(zero(dub[id]-dlb[id])==0)
      glp_set_col_bnds(glpkProbPtr,
		       id + 1,
		       GLP_FX,
		       dlb[id],
		       dub[id]);
    else if(rsense[id]=='F')
      glp_set_col_bnds(glpkProbPtr,
		       id + 1,
		       GLP_FR,
		       dlb[id],
		       dub[id]);
    else if(dub[id]>(BapcodInfinity))
      glp_set_col_bnds(glpkProbPtr,
		       id + 1,
		       GLP_LO,
		       dlb[id],
		       dub[id]);
    else if(dlb[id]<-(BapcodInfinity))
      glp_set_col_bnds(glpkProbPtr,
		       id + 1,
		       GLP_UP,
		       dlb[id],
		       dub[id]);
    else
      glp_set_col_bnds(glpkProbPtr,
		       id + 1,
		       GLP_DB,
		       dlb[id],
		       dub[id]);

  delete [] rhsub;
  rhsub =NULL;
  delete [] rhslb;
  rhslb =NULL;
  delete [] rhss;
  rhss = NULL;
  delete [] dub;
  dub =NULL;
  delete [] dlb;
  dlb =NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] probname;
  probname = NULL;
  delete [] rnames;
  rnames = NULL;
  delete [] cnames;
  cnames = NULL;
  delete [] matrow;
  matrow = NULL;
  delete [] matcol;
  matcol = NULL;
  delete [] matval;
  matval = NULL;


  return;
}


void LpGlpkInterface::unLoadFormulation()
{
  //printout("LpGlpkInterface::~unLoadFormulation()", 6);
  // _formCurrentlySaved = false;
  _formCurrentlyLoaded = false;
  return;
}

void MipGlpkInterface::loadFormulation(const std::string & name,
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
  //printout("MipGlpkInterface::loadFormulation");
  if (printL(5))
    std::cout << "MipGlpkInterface::loadFormulation() ncol = " << ncol
	      << " , nrow = " << nrow
	      << std::endl;

  _ncol = ncol;
  _nrow = nrow;

  /// Free space for loading problem
  makeSpaceForLoadingForm();
  _pureLP = false;
  /// We work only on lp with lpsolve, nolt allow but could used with BaB of lpSolve


  LpGlpkInterface::loadFormulation(name,
				   minmaxStatus,
				   ncol,
				   nrow,
				   mapSeqnb2Cname,
				   mapSeqnb2Rname,
				   objectRow,
				   colMatrix,
				   rhsv,
				   bounds,
				   types,
				   directs,
				   sets);
  double lb, ub;
  for (set<ProbType>::const_iterator typePt = types.begin();
       typePt != types.end();
       typePt++)
    {
      if (typePt->type == 'I')
	{
	  glp_set_col_kind(glpkProbPtr,
			   typePt->ref + 1,
			   GLP_IV);
	}
      if (typePt->type == 'B')
	{
	  lb = glp_get_col_lb(glpkProbPtr,
			      typePt->ref + 1);
	  ub = glp_get_col_ub(glpkProbPtr,
			      typePt->ref + 1);
	  if(ub<1)
	    glp_set_col_kind(glpkProbPtr,
			     typePt->ref + 1,
			     GLP_IV);
	  else if(lb>0)
	    glp_set_col_kind(glpkProbPtr,
			     typePt->ref + 1,
			     GLP_IV);
	  else
	    glp_set_col_kind(glpkProbPtr,
			     typePt->ref + 1,
			     GLP_BV);
	}
    }

  return;
}


void LpGlpkInterface::addCols(const std::set<ProbCoef> & objectiveRow,
			      const std::set<ProbCoef, ProbCoefColSmallerThan> & colMatrix,
			      const std::set<ProbBound> & bounds,
			      const std::map<int, std::string> & mapSeqnb2Cname)
{
  //printout("LpGlpkInterface::addCols()");
  int newcol = objectiveRow.size();
  if (newcol <= 0)
    return;
  int readNcol = glp_get_num_cols(glpkProbPtr);
  int readNrow = glp_get_num_rows(glpkProbPtr);

  require(readNcol == _ncol,
	  "LpGlpkInterface::addCols: readNcol != _ncol");

  require(readNrow == _nrow,
	  "LpGlpkInterface::addCols: readNrow != _nrow");

  int rncol = glp_add_cols(glpkProbPtr, newcol);
  std::string curName;
  char* place;
  char* ecnames = new char[(wordSize + 1) * newcol];
  fill_n(ecnames, (wordSize + 1) * newcol, endOfWord);
  for(int i = 0; i < newcol; i++)
    {
      place = &(ecnames[(wordSize + 1)*i]);
      curName = (mapSeqnb2Cname^(i+_ncol)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
      if (printL(6))
	std::cout << " LpGlpkInterface::addCols: cname "
		  << place
		  << std::endl;
    }

  for(int i = 0;i < newcol; i++)
    glp_set_col_name(glpkProbPtr,
		     i + rncol,
		     &ecnames[(wordSize + 1)*i]);

  for (set<ProbCoef>::const_iterator oPtr = objectiveRow.begin();
       oPtr != objectiveRow.end();
       oPtr++)
    glp_set_obj_coef(glpkProbPtr,
		     oPtr->colRef + 1,
		     zero(oPtr->coef));

  double *ebdub = new double[newcol];
  fill_n(ebdub, newcol, BapcodInfinity*10);

  double *ebdlb = new double[newcol];
  fill_n(ebdlb, newcol, 0.);

  char *rsense= new char[newcol];
  fill_n(rsense, newcol,' ');

  for (set<ProbBound>::const_iterator bPtr = bounds.begin();
       bPtr != bounds.end();
       bPtr++)
    {
      if(bPtr->sense=='U')
	{
	  ebdub[bPtr->ref - readNcol] = zero(bPtr->bound);
	  //rsense[bPtr->ref - readNcol] = 'U';//@todo fv to verify 
	}
      else if(bPtr->sense == 'L')
	{
	  ebdlb[bPtr->ref - readNcol] = zero(bPtr->bound);
	  //rsense[bPtr->ref - readNcol] = 'L';//@todo fv to verify 
	}
      else if(bPtr->sense == 'F')
        {
          ebdub[bPtr->ref - readNcol] = zero(bPtr->bound);
          ebdlb[bPtr->ref - readNcol] = zero(bPtr->bound);
          rsense[bPtr->ref - readNcol] = 'F';
	  require(0,
		  "LpGlpkInterface::addCols(): Variable free to verify");
        }
    }

  for(int id = 0;id<newcol;id++)
    {
      if(zero(ebdub[id]-ebdlb[id]) == 0)
	glp_set_col_bnds(glpkProbPtr, id + rncol,
			 GLP_FX,
			 ebdlb[id],
			 ebdub[id]);
      else if(rsense[id]=='F')
	glp_set_col_bnds(glpkProbPtr,
			 id + rncol,
			 GLP_FR,
			 ebdlb[id],
			 ebdub[id]);
      else if(ebdub[id]>(BapcodInfinity))
        glp_set_col_bnds(glpkProbPtr,
			 id + rncol,
			 GLP_LO,
			 ebdlb[id],
			 ebdub[id]);
      else if(ebdlb[id]<-(BapcodInfinity))
	glp_set_col_bnds(glpkProbPtr,
			 id + rncol,
			 GLP_UP,
			 ebdlb[id],
			 ebdub[id]);
      else
        glp_set_col_bnds(glpkProbPtr,
			 id + rncol,
			 GLP_DB,
			 ebdlb[id],
			 ebdub[id]);
    }

  int *ematind = new int[readNrow + 1];
  double *ematval = new double[readNrow + 1];
  set<ProbCoef, ProbCoefColSmallerThan>::const_iterator mPtr = colMatrix.begin();
  int incol = rncol-1;
  int cnt;
  int newColRef;
  for (newColRef = 1; newColRef <= newcol; newColRef++)
    {
      incol++;
      if (mPtr == colMatrix.end())
	break;

      if (incol < mPtr->colRef + 1)
	continue;

      cnt = 0;
      fill_n(ematind, readNrow + 1, -1);
      fill_n(ematval, readNrow + 1, 0);
      while (incol == mPtr->colRef + 1)
	{
	  cnt++;
	  ematind[cnt] = mPtr->rowRef + 1;
	  ematval[cnt] = zero(mPtr->coef);
	  mPtr++;
	  if (mPtr == colMatrix.end())
	    break;
	}

      // require(incol==(newColRef + readNcol),"pb addCol");
      glp_set_mat_col(glpkProbPtr,
		      incol,
		      cnt,
		      ematind,
		      ematval);
    }

  delete [] ebdub;
  ebdub =NULL;
  delete [] ebdlb;
  ebdlb =NULL;
  delete [] rsense;
  rsense = NULL;
  delete [] ecnames;
  ecnames = NULL;
  delete [] ematind;
  ematind = NULL;
  delete [] ematval;
  ematval = NULL;

  _ncol += newcol;

  return;
}

void LpGlpkInterface::delCols(const std::set<int> & indexSetOfCol2Delete)
{
  //printout("LpGlpkInterface::delCols()");
  int nbCol2Delete = indexSetOfCol2Delete.size();
  if (nbCol2Delete <= 0)
    return;
  int readNcol;
  readNcol = glp_get_num_cols(glpkProbPtr);
  require(readNcol <= _ncol,
	  "LpGlpkInterface::delCols: readNcol > _ncol");

  require(nbCol2Delete <= readNcol,
	  "LpGlpkInterface::delCols: nbCol2Delete > readNcol");

  int *dindex = new int[nbCol2Delete + 1];
  fill_n(dindex, nbCol2Delete + 1, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfCol2Delete.begin();
       iPtr != indexSetOfCol2Delete.end();
       iPtr++)
    {
      cnt++;
      dindex[cnt] = *iPtr + 1;
    }

  require(cnt==nbCol2Delete,
	  "LpGlpkInterface::delCols: nbCol2Delete != cnt");

  glp_del_cols(glpkProbPtr,
	       nbCol2Delete,
	       dindex);

  delete[] dindex;
  dindex = NULL;


  _ncol -= nbCol2Delete;
  return;
}

void LpGlpkInterface::addRows( const std::set<ProbBound> & rhsv,
			       const std::set<ProbCoef, ProbCoefRowSmallerThan> & rowMatrix,
			       const std::map<int, std::string> & mapSeqnb2Rname)
{
  //printout("LpGlpkInterface::addRows()");
  int newrows = rhsv.size();
  if (newrows <= 0)
    return;
  int readNrow = glp_get_num_rows(glpkProbPtr);
  int readNcol = glp_get_num_cols(glpkProbPtr);
  require(readNrow == _nrow,
	  "LpGlpkInterface::addRows: readNrow != _nrow");
  require(readNcol == _ncol,
	  "LpGlpkInterface::addRows: readNcol != _ncol");

  // #if defined(DUMPING)
  // std::cout << "<>" <<  " add_row/1 " << newrows << std::endl;
  // #endif
  int rnrow = glp_add_rows(glpkProbPtr, newrows);
  std::string curName;
  char * place;

  char *ernames = new char[placeSize * newrows];
  fill_n(ernames, placeSize * newrows, endOfWord);

  for(int i = 0; i < newrows; i++)
    {
      place = & (ernames[placeSize * i]);
      curName = (mapSeqnb2Rname^(i+_nrow)) + filling;
      strncpy(place, curName.c_str(), wordSize);
      place[wordSize] = endOfWord;
    }
  if (printL(7))
    {
      for(int i = 0; i < newrows; i++)
	printf("ername %s\n", &ernames[placeSize * i]);
    }
  for(int i = 0;i<newrows;i++)
    {
      glp_set_row_name(glpkProbPtr,
		       i + rnrow,
		       &ernames[placeSize * i]);
    }

  double *rhsub = new double[newrows];
  fill_n(rhsub, newrows, BapcodInfinity*10);

  double *rhslb = new double[newrows];
  fill_n(rhslb, newrows,-BapcodInfinity*10);

  char *rhss= new char[newrows];
  fill_n(rhss, newrows,' ');

  for (set<ProbBound>::const_iterator rvPtr = rhsv.begin();
       rvPtr != rhsv.end();
       rvPtr++)
    {
      if(rvPtr->sense=='G')
	{
	  rhslb[rvPtr->ref - readNrow]=zero(rvPtr->bound);
	  //rhss[rvPtr->ref - readNrow]= 'G';
	}
      else if(rvPtr->sense=='L')
	{
	  rhsub[rvPtr->ref - readNrow]=zero(rvPtr->bound);
	  //rhss[rvPtr->ref - readNrow]= 'L';
	}
      else if(rvPtr->sense=='E')
	{
	  rhsub[rvPtr->ref - readNrow]=zero(rvPtr->bound);
	  rhslb[rvPtr->ref - readNrow]=zero(rvPtr->bound);
	  rhss[rvPtr->ref - readNrow]='E';
	}
    }

  for(int id = 0;id<newrows;id++)
    {
      if(rhss[id]=='E')
	{
	  glp_set_row_bnds(glpkProbPtr,
			   id + rnrow,
			   GLP_FX,
			   rhslb[id],
			   rhsub[id]);
	}
      else if(rhsub[id]>(BapcodInfinity))
	{
	  if(rhslb[id]<(-BapcodInfinity)) {
	    glp_set_row_bnds(glpkProbPtr,
			     id + rnrow,
			     GLP_FR,
			     rhslb[id],
			     rhsub[id]);
	  }
	  else {
	    glp_set_row_bnds(glpkProbPtr,
			     id + rnrow,
			     GLP_LO,
			     rhslb[id],
			     rhsub[id]);
	  }
	}
      else{
	if(rhslb[id]<(-BapcodInfinity)) {
	  glp_set_row_bnds(glpkProbPtr,
			   id + rnrow,
			   GLP_UP,
			   rhslb[id],
			   rhsub[id]);
	}
	else {
	  glp_set_row_bnds(glpkProbPtr,
			   id + rnrow,
			   GLP_DB,
			   rhslb[id],
			   rhsub[id]);
	}
      }
    }

  int *ematind = new int[readNcol + 1];
  double *ematval = new double[readNcol + 1];
  set<ProbCoef, ProbCoefRowSmallerThan>::const_iterator mPtr = rowMatrix.begin();

  int inrow = rnrow-1;
  int cnt;
  int newRowRef;
  for (newRowRef = 1; newRowRef <= newrows; newRowRef++)
    {
      inrow++;
      if (mPtr == rowMatrix.end())
	break;

      if (inrow < mPtr->rowRef + 1)
	continue;

      cnt = 0;//1;
      fill_n(ematind, readNcol + 1, -1);
      fill_n(ematval, readNcol + 1, 0);

      while (inrow == mPtr->rowRef + 1)
	{
	  cnt++;
	  ematind[cnt] = mPtr->colRef + 1;
	  ematval[cnt] = zero(mPtr->coef);
	  mPtr++;
	  if (mPtr == rowMatrix.end())
	    break;
	}

      // require(inrow==(newRowRef + readNrow),"pb addRow");
      glp_set_mat_row(glpkProbPtr,
		      inrow,
		      cnt,
		      ematind,
		      ematval);
    }

  delete[] rhsub;
  rhsub =NULL;
  delete[] rhslb;
  rhslb =NULL;
  delete[] rhss;
  rhss = NULL;
  delete[] ernames;
  ernames = NULL;
  delete[] ematind;
  ematind = NULL;
  delete[] ematval;
  ematval = NULL;

  _nrow += newrows;

  return;
}

void LpGlpkInterface::delRows(const std::set<int> & indexSetOfRow2Delete)
{
  //printout("LpGlpkInterface::delRows()");
  int nbRow2Delete = indexSetOfRow2Delete.size();
  if (nbRow2Delete <= 0)
    return;
  int readNrow;
  readNrow = glp_get_num_rows(glpkProbPtr);
  require(readNrow <= _nrow,
	  "LpGlpkInterface::delRows: readNrow > _nrow");
  require(nbRow2Delete <= readNrow,
	  "LpGlpkInterface::delRows: nbRow2Delete > readNrow");

  int *dindex = new int[nbRow2Delete + 1];
  fill_n(dindex, nbRow2Delete + 1, -1);

  int cnt(0);
  for (set<int>::const_iterator iPtr = indexSetOfRow2Delete.begin();
       iPtr != indexSetOfRow2Delete.end();
       iPtr++)
    {
      cnt++;
      dindex[cnt] = *iPtr +1;
    }

  require(cnt==nbRow2Delete,"LpGlpkInterface::delCols: nbCol2Delete != cnt");
  glp_del_rows(glpkProbPtr,
	       nbRow2Delete,
	       dindex);

  delete[] dindex; dindex = NULL;

  _nrow -= nbRow2Delete;

  return;
}

void LpGlpkInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("LpGlpkInterface::addDirectives");
  check(1, "LP can not have directives");
  return;
}

/// Not yet for LPsolve
void MipGlpkInterface::addDirectives(const std::set<ProbIntC> & newDirects)
{
  //printout("MipGlpkInterface::addDirectives");
  check(1,
	"MipGlpkInterface::addDirectives cannot by used. Instead one should record directives at the outset");

  /// It doesn't exist under GLPK
  require(0, "MipGlpkInterface::addDirectives(): don't exist in GLPK");

  return;
}


void LpGlpkInterface::getObjCoef(ProbCoef & pc)
{
  //printout("LpGlpkInterface::getObjCoef()", 6);
  double *coef = new double[1];
  *coef = glp_get_obj_coef(glpkProbPtr,
			   pc.colRef + 1);

  pc.coef = *coef;
  delete [] coef;
  return;
}

void LpGlpkInterface::chgObjCoef(const ProbCoef & pc)
{
  //printout("LpGlpkInterface::chgObjCoef()", 6);
  glp_set_obj_coef(glpkProbPtr,
		   pc.colRef + 1,
		   pc.coef);

  return;
}


void LpGlpkInterface::chgMatCoef(const ProbCoef & pc)
{
  //printout("LpGlpkInterface::chgMatCoef()", 6);
  int readNrow;
  readNrow = glp_get_num_rows(glpkProbPtr);
  require(readNrow <= _nrow,
	  "LpGlpkInterface::chgMatCoef: readNrow > _nrow");

  double * val= new double[readNrow + 1];
  fill_n(val, readNrow + 1, 0);

  int * ind = new int[readNrow + 1];
  fill_n(ind, readNrow + 1, -1);

  int len = glp_get_mat_col(glpkProbPtr,
			    pc.colRef + 1,
			    ind,
			    val);
  len = len + 1;
  ind[len] = pc.rowRef + 1;
  val[len] = pc.coef;

  glp_set_mat_col(glpkProbPtr,
		  pc.colRef + 1,
		  len,
		  ind,
		  val);

  delete [] val;
  val = NULL;
  delete [] ind;
  ind = NULL;


  return;
}

void LpGlpkInterface::chgRhs(const ProbBound & pb)
{
  //printout("LpGlpkInterface::chgRhs()", 6);
  //   int sense = glp_get_row_type(glpkProbPtr, pb.ref + 1);
  //   if (sense == GLP_FX)
  //   glp_set_row_bnds(glpkProbPtr, 
  // 		   pb.ref + 1, GLP_FX, 
  // 		   zero(pb.bound), 
  // 		   zero(pb.bound));
  //   else if (sense == GLP_UP)
  //     glp_set_row_bnds(glpkProbPtr,
  // 		     pb.ref + 1, GLP_UP,
  // 		     zero(pb.bound),
  // 		     zero(pb.bound));
  //   else if(sense == GLP_LO)
  //   glp_set_row_bnds(glpkProbPtr,
  // 		   pb.ref + 1,
  // 		   GLP_LO, zero(pb.bound),
  // 		   zero(pb.bound));

  //  if (pb.sense == 'E')
  //     glp_set_row_bnds(glpkProbPtr,
  // 		     pb.ref + 1, GLP_FX,
  // 		     zero(pb.bound),
  // 		     zero(pb.bound));
  //   else if (pb.sense == 'L')
  //     glp_set_row_bnds(glpkProbPtr,
  // 		     pb.ref + 1,
  // 		     GLP_UP, zero(pb.bound),
  // 		     zero(pb.bound));
  //   else if(pb.sense == 'G')
  //     glp_set_row_bnds(glpkProbPtr,
  // 		     pb.ref + 1,
  // 		     GLP_LO,
  // 		     zero(pb.bound),
  // 		     zero(pb.bound));

  require(0,"LpGlpkInterface::chgRhs(): chgRhs 1 not implemented");

  return;
}


void LpGlpkInterface::chgRhs(const std::set<ProbBound> & newRhs)
{
  //printout("LpGlpkInterface::chgRhs()", 6);
  if (newRhs.empty()) return;
  for (set<ProbBound>::const_iterator bPtr = newRhs.begin(); bPtr != newRhs.end(); bPtr++)
    if (bPtr->sense == 'E')
      glp_set_row_bnds(glpkProbPtr,
		       bPtr->ref + 1,
		       GLP_FX,
		       zero(bPtr->bound),
		       zero(bPtr->bound));
    else if (bPtr->sense == 'L')
      glp_set_row_bnds(glpkProbPtr,
		       bPtr->ref + 1,
		       GLP_UP,
		       zero(bPtr->bound),
		       zero(bPtr->bound));
    else if(bPtr->sense == 'G')
      glp_set_row_bnds(glpkProbPtr,
		       bPtr->ref + 1,
		       GLP_LO,
		       zero(bPtr->bound),
		       zero(bPtr->bound));

  return;
}

void LpGlpkInterface::chgBds(const std::set<ProbBound> & newBounds)
{
  //printout("LpGlpkInterface::chgBds()");
  if (newBounds.empty())
    return;
  double *dub = new double[newBounds.size()];
  fill_n(dub, newBounds.size(), BapcodInfinity * 10);

  double *dlb = new double[newBounds.size()];
  fill_n(dlb, newBounds.size(), 0);

  int * ind = new int[newBounds.size()];
  fill_n(ind, newBounds.size(), -1);

  char *rsense= new char[newBounds.size()];
  fill_n(rsense, newBounds.size(),' ');

  int nbCol = 0;
  for (set<ProbBound>::const_iterator bPtr = newBounds.begin();
       bPtr != newBounds.end();
       bPtr++)
    {
      if(nbCol>0 && ind[nbCol-1] == bPtr->ref + 1)
	{
	  if(bPtr->sense == 'U')
	    {
	    dub[nbCol-1] = zero(bPtr->bound);
	    // rsense[nbCol-1] = 'U';
	    }
	  else if(bPtr->sense == 'L')
	    {
	    dlb[nbCol-1] = zero(bPtr->bound);
	    // rsense[nbCol-1] = 'F';
	    }
	  else if(bPtr->sense == 'F')
	    {
	      dub[nbCol-1] = zero(bPtr->bound);
	      dlb[nbCol-1] = zero(bPtr->bound);
	      rsense[nbCol-1] = 'F';
	    }
	}
      else
	{
	  ind[nbCol] = bPtr->ref + 1;
	  if(bPtr->sense=='U')
	    {
	      dlb[nbCol] = glp_get_col_lb(glpkProbPtr,
					  bPtr->ref + 1);
	      dub[nbCol] = zero(bPtr->bound);
	    }
	  else if(bPtr->sense=='L')
	    {
	      dlb[nbCol] = zero(bPtr->bound);
	      dub[nbCol] = glp_get_col_ub(glpkProbPtr,
					  bPtr->ref + 1);
	    }
	  else if(bPtr->sense=='F')
	    {
	      dub[nbCol] = zero(bPtr->bound);
	      dlb[nbCol] = zero(bPtr->bound);
	      rsense[nbCol] = 'F';
	    }
	  nbCol++;
	}
    }

  for(int id = 0; id < nbCol; id++)
    {
      if(zero(dub[id]-dlb[id]) == 0)
	glp_set_col_bnds(glpkProbPtr,
			 ind[id],
			 GLP_FX,
			 dlb[id],
			 dub[id]);
      else if(rsense[id] == 'F')
	glp_set_col_bnds(glpkProbPtr,
			 ind[id],
			 GLP_FR,
			 dlb[id],
			 dub[id]);
      else if(dub[id] > (BapcodInfinity))
	glp_set_col_bnds(glpkProbPtr,
			 ind[id],
			 GLP_LO,
			 dlb[id],
			 dub[id]);
      else if(dlb[id] < -(BapcodInfinity))
	glp_set_col_bnds(glpkProbPtr,
			 ind[id],
			 GLP_UP,
			 dlb[id],
			 dub[id]);
      else
	glp_set_col_bnds(glpkProbPtr,
			 ind[id],
			 GLP_DB,
			 dlb[id],
			 dub[id]);
    }

  delete[] dub;
  dub = NULL;
  delete[] dlb;
  dlb = NULL;
  delete[] ind;
  ind = NULL;
  delete[] rsense;
  rsense = NULL;


  return;
}

void LpGlpkInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("LpGlpkInterface::chgColType()");
  check(1,"LP Form cannot have ColType");

  return;
}

void MipGlpkInterface::chgColType(const std::set<ProbType> & newTypes)
{
  //printout("MipGlpkInterface::chgColType()");
  /// Not yet for LPsolve


  double lb, ub;

  for (set<ProbType>::const_iterator bPtr = newTypes.begin();
       bPtr != newTypes.end();
       bPtr++)
    {
      if (bPtr->type == 'I')
	glp_set_col_kind(glpkProbPtr,
			 bPtr->ref + 1,
			 GLP_IV);
      if (bPtr->type == 'B')
	{
	  lb = glp_get_col_lb(glpkProbPtr,
			      bPtr->ref + 1);
	  ub = glp_get_col_ub(glpkProbPtr,
			      bPtr->ref + 1);
	  if(ub<1)
	    glp_set_col_kind(glpkProbPtr,
			     bPtr->ref + 1,
			     GLP_IV);
	  else if(lb>0)
	    glp_set_col_kind(glpkProbPtr,
			     bPtr->ref + 1,
			     GLP_IV);
	  else
	    glp_set_col_kind(glpkProbPtr,
			     bPtr->ref + 1,
			     GLP_BV);
	}
    }

  return;
}

void LpGlpkInterface::getObjVal(Double & objV)
{
  //printout("LpGlpkInterface::getObjVal");
  double objval;
  objval = glp_get_obj_val(glpkProbPtr);

  objV = zero(objval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsolve
void MipGlpkInterface::getObjVal(Double & objV)
{
  //printout("MipGlpkInterface::getObjVal");
  double objval = 0.0;
  objval = glp_mip_obj_val(glpkProbPtr);

  objV = zero(objval, param.HIGHPRECISION);

  return;
}

void LpGlpkInterface::getDualBound(Double & val)
{
  //printout("LpGlpkInterface::getDualBound");
  double dualval = 0.0;
  if(glp_get_status(glpkProbPtr) == GLP_OPT)
    dualval = glp_get_obj_val(glpkProbPtr);
  else
    {
      require(0, "LpGlpkInterface::getDualBound(): error getDualBound not implemented");
      if(glp_get_obj_dir(glpkProbPtr) == GLP_MAX)
	dualval = BapcodInfinity;
      else
	dualval = -BapcodInfinity;
    }


  val = zero(dualval, param.HIGHPRECISION);

  return;
}

/// Not yet for LPsole
void MipGlpkInterface::getDualBound(Double & val)
{
  //printout("MipGlpkInterface::getDualBound");
  double dualval = 0.0;
  dualval = glp_get_obj_val(glpkProbPtr);

  val = zero(dualval, param.HIGHPRECISION);

  return;
}


void LpGlpkInterface::getPrimalBound(Double & val)
{
  //printout("LpGlpkInterface::getPrimalBound");
  double primalval = 0.0;
  primalval = glp_get_obj_val(glpkProbPtr);

  val = zero(primalval, param.HIGHPRECISION);
  return;
}


/// Not yet for LPsolve
void MipGlpkInterface::getPrimalBound(Double & val)
{
  //printout("MipGlpkInterface::getPrimalBound");
  double primalval = 0.0;
  primalval = glp_mip_obj_val(glpkProbPtr);

  val = zero(primalval, param.HIGHPRECISION);
  return;
}

bool LpGlpkInterface::getOptimStatus(SolutionStatus & lpStatus,
				     SolutionStatus & mipStatus)
{
  //printout("LpGlpkInterface::getOptimStatus");
  int status;
  status = glp_get_status(glpkProbPtr);

  if (status == GLP_OPT) {
    lpStatus = SolutionStatus::Optimum;
    if (printL(3))
      std::cout << "LpGlpkInterface::getOptimStatus: LP solution is optimal"
		<< std::endl;

    /// Optimal solution found => no problem encountered
    return(true);
  }

  if (printL(4))
    MPSwrite();

  if (status == GLP_INFEAS)
    {
      lpStatus = SolutionStatus::Infeasible;
      if (printL(3))
	std::cout << "LpGlpkInterface::getOptimStatus: LP is infeasible"
		  << std::endl;

      return(false);
    }

  if (status == GLP_NOFEAS)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGlpkInterface::getOptimStatus: LP has no feasible solution"
		  << std::endl;

      return(false);
    }

  if (status == GLP_UNDEF)
    {
      lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "LpGlpkInterface::getOptimStatus: LP solution is undefined"
		  << std::endl;

      return(false);
    }

  if (status == GLP_UNBND)
    {
      lpStatus = SolutionStatus::Unbounded;
      if (printL(3))
	std::cout << "LpGlpkInterface::getOptimStatus: LP has unbounded solution"
		  << std::endl;

      return(false);
    }

  if (status == GLP_FEAS)
    {
      lpStatus = SolutionStatus::PrimalFeasSolFound;
      if (printL(3))
	std::cout << "LpGlpkInterface::getOptimStatus: LP solution is feasible"
		  << std::endl;

      return(false);
    }


  std::cout << "LpGlpkInterface::getOptimStatus: undefined status"
	    << std::endl;

  return(false);
}

/// Not yet for LPsolve
bool MipGlpkInterface::getOptimStatus( SolutionStatus & lpStatus,
				       SolutionStatus & mipStatus)
{
  //printout("MipGlpkInterface::getOptimStatus");
  int status;
  status = glp_mip_status(glpkProbPtr);
  int lpstat = glp_get_status(glpkProbPtr);

  if (status == GLP_OPT)
    {
      mipStatus = SolutionStatus::Optimum;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipGlpkInterface::getOptimStatus: MIP solution is integer optimal"
		  << std::endl;

      /// Optimal solution found => no problem encountered
      return(true);
    }

  if (status == GLP_FEAS)
    {
      mipStatus = SolutionStatus::PrimalFeasSolFound;
      lpStatus = SolutionStatus::Optimum;
      if (printL(3))
	std::cout << "MipGlpkInterface::getOptimStatus: MIP solution is integer feasible, however its optimality has not been proven, perhaps due to premature termination of the search"
		  << std::endl;
      return(true);
    }

  if (printL(4))
    MPSwrite();

  if (status == GLP_NOFEAS)
    {
      mipStatus = SolutionStatus::Infeasible;
      if(lpstat==GLP_OPT)
	lpStatus = SolutionStatus::Optimum;
      else
	lpStatus = SolutionStatus::UnSolved;
      if (printL(3))
	std::cout << "MipGlpkInterface::getOptimStatus: problem has no integer feasible solution (proven by the solver)"
		  << std::endl;

      return(false);
    }

  if (status == GLP_UNDEF)
    {
      mipStatus = SolutionStatus::UnSolved;
      if(lpstat==GLP_OPT)
	lpStatus = SolutionStatus::Optimum;
      else
	lpStatus = SolutionStatus::UnSolved;

      if (printL(3))
	std::cout << "MipGlpkInterface::getOptimStatus: MIP solution is undefined"
		  << std::endl;

      return(false);
    }

  std::cout << "MipGlpkInterface::getOptimStatus: undefined status"
	    << std::endl;

  return(false);
}



void LpGlpkInterface::getSol(std::map<int, Double> & primSol,
			     const bool & ifPrint)
{
  //printout("LpGlpkInterface::getSol(primSol)");
  primSol.clear();
  int readNcol = glp_get_num_cols(glpkProbPtr);
  require(readNcol <= _ncol,
	  "LpGlpkInterface::getSol: readNcol > _ncol");

  double *x = new double[readNcol];
  fill_n(x, readNcol, 0.);

  for(int j = 0; j < readNcol; j++)
    {
      x[j]=glp_get_col_prim(glpkProbPtr,
			    j + 1);
    }

  if (printL(6))
    std::cout << "readNcol = " << readNcol
	      << "  _ncol = " << _ncol
	      << std::endl;

  for (int jcol = 0; jcol < readNcol; jcol++)
    if (zero(x[jcol], param.HIGHPRECISION) != 0.0)
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

void MipGlpkInterface::getSol(std::map<int, Double> & primSol,
			      const bool & ifPrint)
{
  //printout("MipGlpkInterface::getSol(primSol)");
  primSol.clear();
  int readNcol = glp_get_num_cols(glpkProbPtr);
  require(readNcol <= _ncol,
	  "MipGlpkInterface::getSol: readNcol > _ncol");

  double *x = new double[readNcol];
  fill_n(x, readNcol, 0.);

  for(int j = 0; j < readNcol; j++)
    x[j] = glp_mip_col_val(glpkProbPtr, j + 1);

  if (printL(6))
    std::cout << "readNcol = " << readNcol
	      << "  _ncol = " << _ncol
	      << std::endl;

  for (int jcol = 0; jcol < readNcol; jcol++)
    if (zero(x[jcol], param.HIGHPRECISION) != 0.0)
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

void LpGlpkInterface::getSol(std::map<int, Double> & primSol,
			     std::map<int, Double> & dualSol,
			     const int & minmaxStatus,
			     const bool & ifPrint)
{
  //printout("LpGlpkInterface::getSol(primSol,dualSol)");
  getSol(primSol, ifPrint);
  dualSol.clear();
  int readNrow = glp_get_num_rows(glpkProbPtr);
  require(readNrow <= _nrow,
	  "LpGlpkInterface::getSol: readNrow > _nrow");

  int flipsign = -1;
  double *dual = new double[readNrow];
  fill_n(dual, readNrow, 0);
  int *rsense = new int[readNrow];
  fill_n(rsense, readNrow,-1);
  flipsign = - 1;
  for(int j = 0; j < readNrow; j++)
    {
      dual[j] = glp_get_row_dual(glpkProbPtr,
				 j + 1);

      rsense[j] = glp_get_row_type(glpkProbPtr,
				   j + 1);
    }

  for (int irow = 0; irow < readNrow; irow++)
    if (zero(dual[irow], param.HIGHPRECISION) != 0.0)
      {
	if (rsense[irow] == GLP_LO)
	  dualSol[irow] = minmaxStatus * flipsign * dual[irow];
	else if (rsense[irow] == GLP_UP)
	  dualSol[irow] = minmaxStatus * flipsign * dual[irow];
	else if (rsense[irow] == GLP_FX)
	  /// Was minmaxStatus*dual[irow];
	  dualSol[irow] = flipsign*dual[irow];
	else if (rsense[irow] == GLP_FR)
	  /// Was minmaxStatus*dual[irow];
	  dualSol[irow] = minmaxStatus*flipsign*dual[irow];
	else if (rsense[irow] == GLP_DB)
	  dualSol[irow] = minmaxStatus*flipsign*dual[irow];

	if (printL(6))
	  {
	    std::cout << "dual[" << irow << "] = " << dual[irow]
		      << std::endl;

	    std::cout << " dualSol[" << irow << "] = " << dualSol[irow]
		      << std::endl;
	  }
	/*#ifdef DUMPING    
	  #endif*/
      }

  delete[] dual;
  dual = NULL;
  delete[] rsense;
  rsense = NULL;

  return;
}

void LpGlpkInterface::getReducedCost(std::map<int, Double> & redCost, 
				  const bool & ifPrint)
{

}

//  const bool LpGlpkInterface::preProcessorChangedForm(int * ptrProbRowCnt,
// 						       const bool & checkVarOnly)
//  {
//    //printout("LpGlpkInterface::preProcessorFeedBack()");
//    require(*ptrProbRowCnt == _nrow, 
// 	   "LpGlpkInterface::preProcessorChangedForm(): not (*ptrProbRowCnt == _nrow)");

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
// 	 std::cout << "LpGlpkInterface::preProcessor has Changed Formulation (_ncol != readNcol): need to rebuild formulation"
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
// 	 std::cout << "LpGlpkInterface::preProcessor has Changed Formulation (_nrow > readNrow): need to rebuild formulation"
// 		   << std::endl;
//        return(true);
//      }

// //    if (_nrow < readNrow)
// //       {
// //         if (printL(5)) 
// // 	  std::cout << "LpGlpkInterface::preProcessor has Changed Formulation (_nrow < readNrow): need to rebuild formulation"
// // 		    << std::endl;
// //         *ptrProbRowCnt = _nrow = readNrow;
// //         return(false);
// //       }

//    if (printL(5)) 
//      std::cout << "LpGlpkInterface::preProcessor has NOT changed Formulation"
// 	       << std::endl;

//    return(false);
//  }

void LpGlpkInterface::LPwrite(const int & minmaxStatus, std::ostream& os)
{
  //printout("LpGlpkInterface::LPwrite()");
  int irow, jcol, iel, nchar;
  int readNrow, readNcol, idone, readNglents;
  double readBdl, readBdu;

  readNrow = glp_get_num_rows(glpkProbPtr);
  readNcol = glp_get_num_cols(glpkProbPtr);

  int *readMclind = new int[readNcol + 1];
  fill_n(readMclind, readNcol + 1, 0);
  char *readCnames = new char[placeSize * readNcol];
  fill_n(readCnames, placeSize * readNcol, ' ');
  double *readDmatval = new double[readNcol + 1];
  fill_n(readDmatval, readNcol + 1, 0);
  printf("PROBLEM: %s whose formulation ref number is %d\n",
	 glp_get_prob_name(glpkProbPtr),
	 ref());

  /// Output Obj
  if ( minmaxStatus ==1 )
    printf("Minimize\n");
  else
    printf("Maximize\n");

  for(int icol = 0; icol < readNcol; icol++)
    readDmatval[icol] = glp_get_obj_coef(glpkProbPtr,
					 icol + 1);

  nchar = printf(" Objectif:");

  for (iel = 0; iel < readNcol; iel++)
    if(Double(readDmatval[iel]) != 0)
      {
	///@todo Replace printf by cout
	nchar += printf(" %+g %s", readDmatval[iel], glp_get_col_name(glpkProbPtr, iel + 1));
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
  printf("\nBounds\n");
  for (jcol = 0; jcol < readNcol; jcol++)
    {
      readBdl = glp_get_col_lb(glpkProbPtr, jcol + 1);
      readBdu = glp_get_col_ub(glpkProbPtr, jcol + 1);
      if (glp_get_col_type(glpkProbPtr, jcol + 1) == GLP_DB)
	{
	  ///@todo Replace printf by cout
	  printf(" %g <= %s <= %g\n", readBdl, glp_get_col_name(glpkProbPtr, jcol + 1), readBdu);
	}
      else if(glp_get_col_type(glpkProbPtr, jcol + 1)==GLP_FX)
	{
	  ///@todo Replace printf by cout
	  printf(" %s = %g\n", glp_get_col_name(glpkProbPtr, jcol + 1), readBdl);
	}
      else if (glp_get_col_type(glpkProbPtr, jcol + 1)==GLP_LO) {
	///@todo Replace printf by cout
	printf(" %s <= %g\n", glp_get_col_name(glpkProbPtr, jcol + 1), readBdu);
      }
      else if(glp_get_col_type(glpkProbPtr, jcol + 1)==GLP_UP)
	{
	  ///@todo Replace printf by cout
	  printf(" %s >= %g\n", glp_get_col_name(glpkProbPtr, jcol + 1), readBdl);
	}
    }

  if (_pureLP) printf("End\n");
  else{
    int tmp;
    readNglents = glp_get_num_int(glpkProbPtr);
    if( readNglents!=0 )
      {
	printf("\nIntegers\n");
	for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
	  {
	    tmp = glp_get_col_kind(glpkProbPtr, jcol + 1);
	    if (tmp == GLP_IV)
	      {
		printf(" %s ", glp_get_col_name(glpkProbPtr, jcol + 1));
		// printf(" %s ", &readCnames[placeSize * jcol]);
		idone++;
	      }
	    if (idone >= wordSize)
	      {
		printf("\n");
		idone = 0;
	      }
	  }
	if (idone != 0)
	  printf("\n");
      }
    readNglents = glp_get_num_bin(glpkProbPtr);
    if( readNglents!=0 )
      {
	printf("\nBinary\n");
	for (idone = 0, jcol = 0; jcol < readNcol; jcol++)
	  {
	    tmp = glp_get_col_kind(glpkProbPtr, jcol + 1);
	    if (tmp == GLP_BV)
	      {
		printf(" %s ", glp_get_col_name(glpkProbPtr, jcol + 1));
		// printf(" %s ", &readCnames[placeSize * jcol]);
		idone++;
	      }
	    if (idone >= wordSize)
	      {
		printf("\n");
		idone = 0;
	      }
	  }
	if (idone != 0)
	  printf("\n");
      }

    printf("End\n");
  }

  /// Add directives informations if any

  delete[] readMclind;
  readMclind = NULL;
  delete[] readDmatval;
  readDmatval = NULL;
  delete[] readCnames;
  readCnames = NULL;


  return;
}

void LpGlpkInterface::MPSwrite()
{
  //printout("LpGlpkInterface::MPSwrite()");
  check(lpx_write_mps(glpkProbPtr,"curprob.mps") != 0,
	"could not write prob");
  check(lpx_write_cpxlp(glpkProbPtr,"curprob.lp") != 0,
	"could not write prob");
  check(lpx_print_prob(glpkProbPtr,"curprob.txt") != 0,
	"could not write prob");
  check(glp_write_mip(glpkProbPtr,"cursol.mip") != 0,
	"could not solve LP");
  check(glp_write_sol(glpkProbPtr,"cursol.sol") != 0,
	"could not solve LP");

  return;
}

void LpGlpkInterface::printRow(int irow,
			       const int & readNcol,
			       int *readMclind,
			       double *readDmatval,
			       char *readCnames)
{
  //printout("LpGlpkInterface::printRow");
  int iel, readNels(0), nchar(0);
  double readSrhs, readIrhs;
  int readType;
  fill_n(readMclind,
	 readNcol + 1,
	 -1);
  fill_n(readDmatval,
	 readNcol +1 ,
	 0);

  readNels = glp_get_mat_row(glpkProbPtr,
			     irow + 1,
			     readMclind,
			     readDmatval);

  readSrhs = glp_get_row_ub(glpkProbPtr,
			    irow + 1);

  readIrhs = glp_get_row_lb(glpkProbPtr,
			    irow + 1);

  readType = glp_get_row_type(glpkProbPtr,
			      irow + 1);

  nchar = printf(" %s:", glp_get_row_name(glpkProbPtr, irow + 1));

  if(readType==GLP_DB)
    {
      printf(" %g <= ", readIrhs);
    }

  if (readNels > 0)
    {
      for (iel = 1; iel <= readNels; iel++)
	{
	  nchar += printf(" %+g %s",
			  readDmatval[readNels-iel + 1],
			  glp_get_col_name(glpkProbPtr,
					   readMclind[readNels-iel + 1]));
	  if (nchar >= 80-9-13)
	    {
	      printf("\n");
	      nchar = 0;
	    }
	}
    }

  if(readType==GLP_FX)
    {
      printf(" = %g\n", readSrhs);
    }
  else if(readType==GLP_UP)
    {
      printf(" <= %g\n", readSrhs);
    }
  else if(readType==GLP_LO)
    {
      printf(" >= %g\n", readIrhs);
    }
  else if(readType==GLP_DB)
    {
      printf(" <= %g\n", readSrhs);
    }

  return;
}

void LpGlpkInterface::optimise(const int & minmaxStatus,
			       const bool & preprocessorOn,
			       const bool & probingOn,
			       const bool & automaticCuttingPlanesOn,
			       const char & solverSelection)
{
  //printout("LpGlpkInterface::optimise()", 6);
  require(_formCurrentlyLoaded, "Form not Currently Loaded",
	  ProgStatus::quit,
	  3);
  if (printL(8))
    MPSwrite();

  if(preprocessorOn==1)
    glpkParamPtr.presolve = GLP_ON;
  else
    glpkParamPtr.presolve = GLP_OFF;

  if (printL(5))
    {
      FILE *foo = fopen("GLPKexit","a");

      char * GLPKtmp= new char[100];
      sprintf(GLPKtmp,
	      "\n//////LP SOLVER num %d///////\n",
	      GLPKLpCounter);

      GLPKLpCounter++;
      fputs(GLPKtmp, foo);

      fclose(foo);
    }

  int verif = glp_bf_exists(glpkProbPtr);
  int ret = 0;

  if(verif == 0)
    ret = glp_factorize(glpkProbPtr);

  if (ret == GLP_EBADB)
    {
      if (printL(5))
	{
	  FILE *foo = fopen("GLPKexit","a");
	  char * GLPKtmp = new char[10];
	  sprintf(GLPKtmp, "WARM UP test\n");
	  fputs(GLPKtmp, foo);
	  fclose(foo);
	}
      lpx_std_basis(glpkProbPtr);
    }
  else if (ret == GLP_ESING)
    {
      /// The basis is invalid; build some valid basis
      if (printL(5))
	{
	  FILE *foo = fopen("GLPKexit","a");
	  char * GLPKtmp = new char[10];
	  sprintf(GLPKtmp, "WARM UP: initial basis matrix is singular\n");
	  fputs(GLPKtmp, foo);
	  fclose(foo);
	}
      lpx_std_basis(glpkProbPtr);
    }

  ret = glp_simplex(glpkProbPtr, &glpkParamPtr);
  check(ret != 0,"could not solve LP");

  //   if (ret == GLP_EBADB)
  //     {
  //       ret = lpx_warm_up(glpkProbPtr);
  //       /// The basis is invalid; build some valid basis
  //       if (ret == LPX_E_BADB)
  // 	{
  // 	  if (printL(5))
  // 	    {
  // 	      FILE *foo = fopen("GLPKexit","a");
  // 	      char * GLPKtmp= new char[10];
  // 	      sprintf(GLPKtmp,"WARM UP\n");
  // 	      fputs(GLPKtmp, foo);
  // 	      fclose(foo);
  // 	    }
  // 	  lpx_std_basis(glpkProbPtr);
  // 	  //lpx_adv_basis(glpkProbPtr);
  // 	  check(glp_simplex(glpkProbPtr, &glpkParamPtr) != 0,
  // 		"could not solve LP");
  // 	}
  //     }
  //   else if(ret == GLP_ESING)
  //     {
  //       ret = lpx_warm_up(glpkProbPtr);
  //       if(ret == LPX_E_SING)
  // 	{  
  // 	  /// The basis is invalid, build some valid basis
  // 	  if(printL(5))
  // 	    {
  // 	      FILE *foo = fopen("GLPKexit","a");
  // 	      char * GLPKtmp= new char[10];
  // 	      sprintf(GLPKtmp,"WARM UP: initial basis matrix is singular\n");
  // 	      fputs(GLPKtmp, foo);
  // 	      fclose(foo);
  // 	    }

  // 	  lpx_std_basis(glpkProbPtr);
  // 	  //lpx_adv_basis(glpkProbPtr);
  // 	}

  //       check(glp_simplex(glpkProbPtr,&glpkParamPtr) != 0,
  // 	    "could not solve LP");
  //   }
  //   else
  //   check(ret != 0,"could not solve LP");

  //   std::cout<<"obj_LP="<<glp_get_obj_val(glpkProbPtr)
  // 	   <<" obj_MIP="<<glp_mip_obj_val(glpkProbPtr)
  // 	   <<std::endl;

  //   std::cout<<"status_LP="<<glp_get_status(glpkProbPtr)
  // 	   <<" status_MIP="<<glp_mip_status(glpkProbPtr)
  // 	   <<std::endl;

  if(printL(5))
    {
      FILE *foo2 = fopen("GLPKexit", "a");
      fputs("//////END LP SOLVER///////\n\n", foo2);
      fclose(foo2);
    }
  return;
}

void MipGlpkInterface::optimise(const int & minmaxStatus,
				const bool & preprocessorOn,
				const bool & probingOn,
				const bool & automaticCuttingPlanesOn,
				const char & solverSelection)
{
  //printout("MipGlpkInterface::optimise()", 6);
  require(_formCurrentlyLoaded,
	  "Form not Currently Loaded/2",
	  ProgStatus::quit,
	  3);
  if (printL(8))
    MPSwrite();

  if(preprocessorOn == 1)
    glpkParamMipPtr.pp_tech = GLP_PP_ALL;
  else
    glpkParamMipPtr.pp_tech = GLP_PP_NONE;
  //   lpx_set_int_parm(glpkProbPtr,
  // 		   LPX_K_PRESOL,
  // 		   (preprocessorOn ? 1 : 0));

  if (printL(5))
    {
      FILE *foo = fopen("GLPKexit", "a");

      char * GLPKtmp= new char[100];
      sprintf(GLPKtmp,
	      "\n//////MIP SOLVER num %d///////\n",
	      GLPKMipCounter);
      GLPKMipCounter++;
      fputs(GLPKtmp, foo);

      fclose(foo);
    }

  //     check(lpx_write_cpxlp(glpkProbPtr, "test.cplx") != 0,"could not solve LP");
  //     check(lpx_print_prob(glpkProbPtr, "test.mod") != 0,"could not solve LP");
  //     check(glp_write_sol(glpkProbPtr,"test_sol.sol") != 0,"could not solve LP");
  //     check(glp_write_mip(glpkProbPtr,"test.sol") != 0,"could not solve LP");

  check(glp_simplex(glpkProbPtr,&glpkParamPtr) != 0
	,"could not solve MIP GLPK");

  int retcod = glp_intopt(glpkProbPtr,&glpkParamMipPtr);
  // lpx_intopt(glpkProbPtr);
  bool test = (retcod != 0 && retcod!=GLP_ETMLIM);
  // LPX_E_OK;
  if (test)
    {
      std::cerr << " ERROR in GLPK : " << retcod << std::endl;
      check(test,"could not solve MIP GLPK");
    }

  std::out << "obj_LP="<< glp_get_obj_val(glpkProbPtr)
	    << " obj_MIP=" << glp_mip_obj_val(glpkProbPtr)
	    <<std::endl;

  std::out << "status_LP="<< glp_get_status(glpkProbPtr)
	    << " status_MIP=" << glp_mip_status(glpkProbPtr)
	    << std::endl;

  if (printL(5))
    {
      FILE *foo2 = fopen("GLPKexit","a");
      fputs("//////END MIP SOLVER///////\n\n", foo2);
      fclose(foo2);
    }
  return;
}

void LpGlpkInterface::reset()
{
  //printout("LpGlpkInterface::reset()");
  /// Do nothing

  return;
}



void MipGlpkInterface::reset()
{
  //printout("MipGlpkInterface::reset()", 6);
  return;
}


void LpGlpkInterface::makeSpaceForLoadingForm()
{
  /// Save loaded form if any
  //printout("LpGlpkInterface::makeSpaceForLoadingForm()", 6);
  _formCurrentlyLoaded = true;

  return;
}

void LpGlpkInterface::saveCopyOfCurForm()
{
  //printout("LpGlpkInterface::saveCopyOfCurForm()", 6);
  /// Save loaded form if any
  return;
}

void LpGlpkInterface::load()
{
  //printout("LpGlpkInterface::load()", 6);
  return;
}

void LpGlpkInterface::unload(const bool & deleteMat)
{
  //printout("LpGlpkInterface::unload()");
  // if (deleteMat == false)

  return;
}

#endif // _WITH_GLPK
