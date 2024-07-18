/**
 *
 * This file bcFormC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcFormC.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcCplexSolverC.hpp"

#if _XPRESSMP_FOUND
#include "bcXpressMPSolverC.hpp"
#endif

#if _GLPK_FOUND
#include "bcGlpkSolverC.hpp"
#endif

#include "bcProblemC.hpp"
#include "bcPrintC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMathProgSolverBuilder.hpp"
#include "bcProbConfigC.hpp"

using namespace std;
using namespace VcIndexStatus;

BaseFormulation::BaseFormulation() :
    _problemPtr(NULL), _minmax(BcObjStatus::maxInt), _status(SolutionStatus::Undefined)
{

  return;
}

BaseFormulation::BaseFormulation(Problem * problemPtr) :
    _problemPtr(problemPtr), _minmax(((problemPtr->objStatus() == BcObjStatus::maxInt)
                                       || (problemPtr->objStatus() == BcObjStatus::maxFloat)) ? -1 : 1),
    _status(SolutionStatus::Undefined)
{
  return;
}

BapcodInit & BaseFormulation::bapcodInit() const
{
  return _problemPtr->probConfPtr()->bapcodInit();
}

const ControlParameters & BaseFormulation::param() const
{
  return _problemPtr->probConfPtr()->param();
}

ControlParameters & BaseFormulation::param()
{
  return _problemPtr->probConfPtr()->param();
}

LPform::LPform(Problem * probPtr, const bool & defineInterface) :
    BaseFormulation(probPtr), _interfacePtr(NULL), _probRowCnt(0), _probColCnt(0), _objScalFact(1)
{
  if (defineInterface)
  {
      MathProgSolverBuilder mathProgSolverBuilder;
      _interfacePtr = mathProgSolverBuilder.buildLpMathProgSolverInterface(probPtr->probConfPtr()->bapcodInitPtr(),
                                                                           probPtr->probConfPtr()->param().solverName(),
                                                                           probPtr->ref(), probPtr->name());
      if (_interfacePtr == nullptr)
      {
          std::string solverName = probPtr->probConfPtr()->param().solverName();
          std::string shortName = solverName.substr(0,solverName.length() - 7);
          std::cerr << "BaPCod error : solver " << shortName << " is not found!" << std::endl
                    << "Please define " << shortName << "_ROOT environment variable before running cmake" << std::endl;
          exit(1);
      }

    _interfacePtr->setLPoptimalityTolerance(param().MasterMipSolverReducedCostTolerance());
    _interfacePtr->setMultiThread(param().MipSolverMultiThread()); /// added by Ruslan
    _interfacePtr->setTimeLimit(param().MipSolverMaxTime()); /// added by Ruslan
    _interfacePtr->setMaxBBNodes(param().MipSolverMaxBBNodes());
    _interfacePtr->setRelativeMipGapLimit(param().relOptimalityGapTolerance());
  }
}

LPform::~LPform()
{
  _interfacePtr->unLoadFormulation();
  if (_interfacePtr != NULL)
  {
    delete _interfacePtr;
    _interfacePtr = NULL;
  }

  return;
}

void LPform::fillDataStruct(Constraint * constrPtr, int & sense, double & rhs, std::vector<int> & cutind,
                            std::vector<double> & cutval) const
{
    sense = (int)constrPtr->sense();
    rhs = constrPtr->curRhs();

    for (ConstVarConstrPtr2Double::iterator itm = constrPtr->member2coefMap().begin();
         itm != constrPtr->member2coefMap().end(); ++itm)
    {
        if (itm->first->inCurForm())
        {
            if (itm->first->isTypeOf(VcId::VariableMask) || itm->first->isTypeOf(VcId::MastColumnMask))
            {
                bapcodInit().require(itm->first->kind() == 'E',
                                     "LPform::fillDataStruct(constrPtr) contraint must be explicit if inCurForm");

                Double coeff = constrPtr->membCoef(itm->first);
                if (!zeroTest(coeff))
                {
                    cutind.push_back(itm->first->index());
                    cutval.push_back((double)coeff);
                }
            }
        }
    }
}

void LPform::retrievePrimalSol(const std::map<int, Double> & primalSolVect, VarPtrSet & inPrimalSol) const
{
    Variable * varPtr;
    for (std::map<int, Double>::const_iterator sPtr = primalSolVect.begin();
         sPtr != primalSolVect.end(); sPtr++)
    {
#if MAPCONSTR
        varPtr = _mapVarSeqNb2VarPtr.at(sPtr->first);
#else
        varPtr = _vectVarPtr.at(sPtr->first);
#endif
        const std::pair<VarPtrSet::iterator, bool> insertionResult = inPrimalSol.insert(varPtr);

        if (insertionResult.second == false)
        {
            bapcodInit().check( insertionResult.second == false,
                               "LPform::retrievePrimalSol(): Some variables are not inserted in inPrimalSol. "
                               "Check the output.",ProgStatus::run);
            std::cout
                    << "Var " << varPtr->name() << "\t ref:" << varPtr->ref()
                    << "\t is not inserted because of var \t"
                    << (*insertionResult.first)->name() << "\t ref:"  << (*insertionResult.first)->ref()
                    << std::endl;
        }

        varPtr->val(sPtr->second);
    }

}

/// Storing the solution from formulation

void LPform::retrieveSol(const char & flag, const bool & ifPrint,
                         VarPtrSet & inPrimalSol, ConstrPtrSet & inDualSol)
{
  std::map<int, Double> primalSolVect;
  std::map<int, Double> dualSolVect;
  Variable * varPtr;
  Constraint *constrPtr;

  // Flag 'r' added by Aurélien to retrieve the reduced cost in addition to the dual values
  if (flag == 'd' || flag == 'r')
    _interfacePtr->getSol(primalSolVect, dualSolVect, _minmax, ifPrint);

  else
    _interfacePtr->getSol(primalSolVect, ifPrint);

//    if (ifPrint)
//    {
//        std::vector<int> colStat, rowStat;
//        _interfacePtr->getBasis(colStat, rowStat);
//        for (int index = 0; index < colStat.size(); ++index)
//        {
//            Variable *varPtr = _mapVarSeqNb2VarPtr[index];
//            std::cout << "Basis status for var " << varPtr->name() << " is " << colStat[index] << std::endl;
//        }
//        for (int index = 0; index < rowStat.size(); ++index)
//        {
//            Constraint *constrPtr = _mapConstrSeqNb2ConstrPtr[index];
//            std::cout << "Basis status for constr " << constrPtr->name() << " is " << rowStat[index] << std::endl;
//        }
//    }

  for (std::map<int, Double>::const_iterator sPtr = primalSolVect.begin();
       sPtr != primalSolVect.end(); sPtr++)
  {
#if MAPCONSTR
    varPtr = _mapVarSeqNb2VarPtr.at(sPtr->first);
#else
    varPtr = _vectVarPtr.at(sPtr->first);
#endif
    const std::pair<VarPtrSet::iterator, bool> insertionResult = inPrimalSol.insert(varPtr);

    if (insertionResult.second == false)
    {
      bapcodInit().check(insertionResult.second == false,
                         "LPform::retrieveSol(): Some variables are not inserted in inPrimalSol. "
                         "Check the output.", ProgStatus::run);
      std::cout
        << "Var " << varPtr->name() << "\t ref:" << varPtr->ref()
        << "\t is not inserted because of var \t"
        << (*insertionResult.first)->name() << "\t ref:"  << (*insertionResult.first)->ref()
        << std::endl;
    }

    varPtr->val(sPtr->second);
    //std::cout << "$$$ primalSolVect[" << varPtr->name() << "] = " << varPtr->val() << std::endl;
    if (ifPrint)
    {
      std::cout << "primalSolVect[" << varPtr->name() << "] = " << varPtr->val()
      << std::endl;

      if(printL(5))
      {
        if (varPtr->isTypeOf(VcId::MastColumnMask))
        {
          MastColumn* col = static_cast<MastColumn*> (varPtr);
          col->spSol()->print();
        }
      }
    }
  }

  for (std::map<int, Double>::const_iterator dsPtr = dualSolVect.begin();
       dsPtr != dualSolVect.end(); dsPtr++)
  {
#if MAPCONSTR
    constrPtr = _mapConstrSeqNb2ConstrPtr.at(dsPtr->first);
#else
    constrPtr = _vectConstrPtr.at(dsPtr->first);
#endif
    inDualSol.insert(constrPtr);
    constrPtr->val(dsPtr->second * _objScalFact);
    if (ifPrint)
      std::cout << "dualSol[" << constrPtr->name() << "] = " << constrPtr->val() << " ref=" << constrPtr->VCref()
      << std::endl;
  }
  // Added by Aurélien to retrieve the reduced cost of the variables
    if (flag == 'r')
    {
        std::map<int, Double> redCostVect;
        _interfacePtr->getReducedCost(redCostVect);
        Variable * varPtr;

        for (std::map<int, Double>::const_iterator sPtr = redCostVect.begin(); sPtr != redCostVect.end(); sPtr++)
        {
#if MAPCONSTR
            varPtr = _mapVarSeqNb2VarPtr.at(sPtr->first);
#else
            varPtr = _vectVarPtr.at(sPtr->first);
#endif
            varPtr->reducedCost(sPtr->second);
            if (ifPrint)
                std::cout << "redCostVect[" << varPtr->name() << "] = " << varPtr->reducedCost() << std::endl;
        }
    }

  return;
}

void LPform::retrieveRedCosts(const bool & ifPrint, VarPtrSet & posRedCostVars)
{
  std::map<int, Double> redCostVect;
  _interfacePtr->getReducedCost(redCostVect);
  Variable * varPtr;

  for (std::map<int, Double>::const_iterator sPtr = redCostVect.begin();
       sPtr != redCostVect.end(); sPtr++)
  {
#if MAPCONSTR
    varPtr = _mapVarSeqNb2VarPtr.at(sPtr->first);
#else
    varPtr = _vectVarPtr.at(sPtr->first);
#endif
    posRedCostVars.insert(varPtr);
    varPtr->reducedCost(sPtr->second);
    if (ifPrint)
      std::cout << "redCostVect[" << varPtr->name() << "] = "
      << varPtr->reducedCost() << std::endl;
  }
}

bool LPform::solve(const double & BarrierConvergenceTolerance,
                   const double & rightHAndSideZeroTol,
                   const double & reducedCostTolerance,
                   const char & flag, const bool & ifPrint,
                   const SolutionStatus & requiredStatus, Double & objVal,
                   Double & primalBound, Double & dualBound,
                   VarPtrSet & inPrimalSol, ConstrPtrSet & inDualSol,
                   const bool & preprocessorOn, const bool & probingOn,
                   const bool & MIPautomaticCuttingPlanesOn,
                   const char & solverSelection)
{
  bool foundSol(false);

  _interfacePtr->load();
  if (ifPrint)
    _interfacePtr->LPwrite(_minmax);

  if (printL(7))
    _interfacePtr->MPSwrite();

  if (printL(2))
    std::cout << "We are just before _interfacePtr->optimiseLp " << std::endl;

  _interfacePtr->optimiseLp(_minmax, BarrierConvergenceTolerance, rightHAndSideZeroTol, reducedCostTolerance,
                            preprocessorOn, probingOn, MIPautomaticCuttingPlanesOn, solverSelection);

  SolutionStatus fakeMipStatus;

  _interfacePtr->getOptimStatus(_status, fakeMipStatus);

  if (printL(6))
  {
    std::cout << "status() = " << status() << std::endl;
    std::cout << "requiredStatus = " << requiredStatus << std::endl;
  }

  bool statusIntersects = _status.intersects(requiredStatus);
  if (!statusIntersects)
  {
      std::cout << requiredStatus << "Current status is : " << BaseFormulation::status();
      _interfacePtr->MPSwrite();
      printForm();
  }

  if (bapcodInit().require(statusIntersects,
                           "LPform::solve(): Formulation could not be solved according to prescribed status"))
    //, ProgStatus::run))
  {
    //foundSol = true;
    LPform::setBounds(objVal, primalBound, dualBound);

    foundSol = _status.count(SolutionStatus::Optimum) || _status.count(SolutionStatus::OptimumUnscalInfeas);
    if (foundSol)
    {
      retrieveSol(flag, ifPrint, inPrimalSol, inDualSol);
    }
  }

  _interfacePtr->unload();

  return (foundSol);
}

void LPform::setBounds(Double & objVal, Double & primalBound,
                       Double & dualBound)
{
  _interfacePtr->getObjValLp(objVal);

  /// Scaling
  objVal *= _objScalFact;

  if (BaseFormulation::status().count(SolutionStatus::Optimum)
      || BaseFormulation::status().count(SolutionStatus::OptimumUnscalInfeas))
  {
    primalBound = dualBound = objVal;
    return;
  }

  if (BaseFormulation::status().count(SolutionStatus::Infeasible) || status()
      .count(SolutionStatus::UnSolved))
  {
    primalBound = dualBound = _minmax * BapcodInfinity;
    return;
  }

  if (BaseFormulation::status().count(SolutionStatus::Unbounded))
  {
    primalBound = dualBound = -_minmax * BapcodInfinity;
    return;
  }

  return;
}

void LPform::fillDataStruct(Constraint * constrPtr,
                            const bool & fillConstrMatrix)
{
  if (printL(7))
    std::cout << "LPform::fillDataStruct(Constraint * constrPtr) name = "
              << constrPtr->name() << " index = " << constrPtr->index()
              << " sense = " << constrPtr->sense() << " rhs = "
              << constrPtr->curRhs() << std::endl;

  _mapSeqnb2Rname[constrPtr->index()] = constrPtr->name();

#ifdef PROBCONTAINERSET
  _rhsv.insert(ProbBound(constrPtr->index(), constrPtr->sense(), constrPtr->curRhs()));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _rhsv.push_back(ProbBound(constrPtr->index(), constrPtr->sense(), constrPtr->curRhs()));
#endif

  if (!fillConstrMatrix)
    return;

  for (ConstVarConstrPtr2Double::iterator itm = constrPtr->member2coefMap()
       .begin(); itm != constrPtr->member2coefMap().end(); ++itm)
  {
    if (printL(7))
      std::cout << "member var name = " << itm->first->name() << " in form = "
      << itm->first->inCurForm() << std::endl;

    if (itm->first->inCurForm())
    {
      //Begin changed by Romain
      if (itm->first->isTypeOf(VcId::VariableMask) || itm->first->isTypeOf(VcId::MastColumnMask))
      {
        bapcodInit().require(itm->first->kind() == 'E',
                             "LPform::fillDataStruct(constrPtr)  contraint must be explicit if inCurForm");

        Double coef = constrPtr->membCoef(itm->first);

        if (printL(7))
          std::cout << "includes var name = " << itm->first->name()
                     << " index = " << itm->first->index() << " coef = " << coef << std::endl;

        if (!zeroTest(coef))
        {
#ifdef PROBMATRIXCONTAINERSET
          _rowMatrix.insert(ProbCoef(constrPtr->index(), itm->first->index(), coef));
#endif

#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
          _rowMatrix.push_back(ProbCoef(constrPtr->index(), itm->first->index(), coef));
#endif
        }
      } else
      {
        bapcodInit().require(false,
                             "LPform::fillDataStruct(constrPtr) constr membership is not a variable");
      }
      //End changed by Romain
    }
  }

  return;
}

void LPform::addConstr2Formulation()
{
  if (_rhsv.empty() && _rowMatrix.empty() && _mapSeqnb2Rname.empty())
  {
    if (printL(5))
      std::cout << "LPform::addConstr2Formulation(): empty constraint: nothing to add in formulation to update" << endl;

    /// Nothing to add in formulation to update
    return;
  }

  if (printL(5))
    std::cout << "LPform::addConstr2Formulation(): add  " << _rhsv.size()
              << " rows with a total number of coef " << _rowMatrix.size() << endl;

  _interfacePtr->load();
  _interfacePtr->addRows(_rhsv, _rowMatrix, _mapSeqnb2Rname);
  if (printL(6))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  clearRowFormulationDataStruct();

  return;
}

void LPform::addVar2Formulation()
{
  if (printL(5))
    std::cout << "LPform::addVar2Formulation(): add  " << _mapSeqnb2Cname.size() << " cols " << endl;

  _interfacePtr->load();
  _interfacePtr->addCols(_objectRow, _colMatrix, _bounds, _mapSeqnb2Cname);

  if (printL(6))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  clearColFormulationDataStruct();

  return;
}

void LPform::chgObjCoef(const ProbCoefContainer & newObjCoef)
{
  for (ProbCoefContainer::const_iterator cPtr = newObjCoef.begin(); cPtr != newObjCoef.end(); cPtr++)
  {
    if (printL(6))
      std::cout << *cPtr;

    _interfacePtr->chgObjCoef(*cPtr);
  }

  return;
}

void LPform::setConstr2Form(Constraint * constrPtr,
                            const bool & fillConstrMatrix)
{
  /// Defined index and update two-way reference map
  bapcodInit().require(constrPtr->inCurForm(),
                       "LPform::setConstr2Form(): constr shod be marqued as to be inclued in the explicit formulation");

  if (printL(6))
    std::cout << "LPform::setConstr2Form(): constr " << constrPtr->name()
              << " fillConstrMatrix " << fillConstrMatrix << endl;

  constrPtr->index(_probRowCnt);
  constrPtr->val(0);
#if MAPCONSTR
  _mapConstrSeqNb2ConstrPtr[_probRowCnt] = constrPtr; //Commented by Romain
#else
  _vectConstrPtr.push_back(constrPtr);
#endif

  _probRowCnt++;
  fillDataStruct(constrPtr, fillConstrMatrix);

  if (printL(7))
    printMatrix();

  return;
}

void LPform::unsetConstr2Form(Constraint * constrPtr)
{
  _indexSetOfRow2Delete.insert(constrPtr->index());
  constrPtr->index(-1);
  _probRowCnt--;

  return;
}

void LPform::setVar2Form(Variable * varPtr)
{
  bapcodInit().require(varPtr->inCurForm(),
                       "LPform::setVar2Form(): var should be marqued as to be included in the explicit formulation");

  /// Define index and update two-way reference map
  varPtr->index(_probColCnt);
  varPtr->val(0);
#if MAPCONSTR
  _mapVarSeqNb2VarPtr[_probColCnt] = varPtr; //Commented by Romain
#else
  _vectVarPtr.push_back(varPtr);
#endif
  _probColCnt++;
  fillDataStruct(varPtr);

  return;
}

void LPform::unsetVar2Form(Variable * varPtr)
{
  if (printL(6))
    std::cout << "LPform::unsetVar2Form(" << varPtr->name() << ")" << std::endl;

  _indexSetOfCol2Delete.insert(varPtr->index());
  varPtr->index(-1);
  _probColCnt--;

  return;
}

void LPform::fillDataStruct(Variable * varPtr)
{
  if (printL(6))
    std::cout << "LPform::fillDataStruct(Variable * varPtr) name = "
              << varPtr->name() << " cost = " << varPtr->curCost() << std::endl;

  _mapSeqnb2Cname[varPtr->index()] = varPtr->name();

#ifdef PROBCONTAINERSET
  _objectRow.insert(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _objectRow.push_back(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif
  /**
   * Variables:
   * 'P' = positive
   * 'N' = negative
   * 'F' = free
   */

  if (varPtr->sense() == 'P')
    bapcodInit().require(varPtr->lb() >= 0,
                         "LPform::fillDataStruct(): ERROR sense() == 'P' && lb() < 0",
                         ProgStatus::run);

  if (varPtr->sense() == 'N')
    bapcodInit().require(varPtr->ub() <= 0,
                         "LPform::fillDataStruct(): ERROR sense() == 'N' && ub() > 0",
                         ProgStatus::run);
#ifdef PROBCONTAINERSET
  _bounds.insert(ProbBound(varPtr->index(), 'U', varPtr->curUb()));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _bounds.push_back(ProbBound(varPtr->index(), 'U', varPtr->curUb()));
#endif
  if (printL(6))
    std::cout << " LPform::fillDataStruct() var " << varPtr->name() << " ub = "
    << varPtr->curUb() << std::endl;

#ifdef PROBCONTAINERSET
  _bounds.insert(ProbBound(varPtr->index(), 'L', varPtr->curLb()));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _bounds.push_back(ProbBound(varPtr->index(), 'L', varPtr->curLb()));
#endif
  if (printL(6))
    std::cout << " LPform::fillDataStruct() var " << varPtr->name() << " lb = "
    << varPtr->curLb() << std::endl;

  for (ConstVarConstrPtr2Double::iterator itm =
       varPtr->member2coefMap().begin(); itm != varPtr->member2coefMap().end();
       ++itm)
    if (itm->first->inCurForm())
    {
      //Begin changed by Romain
      if (itm->first->isTypeOf(VcId::ConstraintMask))
	{
	  bapcodInit().require(itm->first->kind() == 'E',
                           "LPform::fillDataStruct()  contraint must be explicit if inCurForm");
    Double coef = itm->first->membCoef(varPtr);

	  if (!zeroTest(coef))
	    {
#ifdef PROBMATRIXCONTAINERSET
	      _colMatrix.insert(ProbCoef(itm->first->index(), varPtr->index(), coef));
#endif

#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
        _colMatrix.push_back(ProbCoef(itm->first->index(), varPtr->index(), coef));
#endif
	    }
	}
      else
	{
	  bapcodInit().require(false,
                           "LPform::fillDataStruct() var membership is not a contraint");
	}
      //End changed by Romain
    }

  return;
}

void LPform::delConstrFromFormulation()
{
  /// No constr to delete
  if (_indexSetOfRow2Delete.empty())
    return;

  /// Reset index of existing constraints (shift them by nb of deleted columns before them)
  int shift(0);

  // Commented by Romain.
#if MAPCONSTR
 for (std::map<int, Constraint *>::const_iterator it = _mapConstrSeqNb2ConstrPtr.begin();
      it != _mapConstrSeqNb2ConstrPtr.end(); ++it)
 {
   if (_indexSetOfRow2Delete.count(it->first))
   {
     shift++;
   }
   else
   {
     if(printL(7))
        cout << "LPform::delConstrFromFormulation constr stays in prob "
             << hex << (long) (it->second) << dec << endl;
     it->second->index(it->second->index() - shift);
   }

 }
#else
  for (std::deque<Constraint*>::const_iterator it = _vectConstrPtr.begin();
       it != _vectConstrPtr.end(); ++it)
  {
    if (_indexSetOfRow2Delete.count((*it)->index()))
    {
      ++shift;
    } else
    {
      (*it)->index((*it)->index() - shift);
    }
  }
#endif
  /// Rebuild mapConstrSeqNb2VarPtr

  //Commented by Romain
#if MAPCONSTR
  _mapConstrSeqNb2ConstrPtr.clear();
  for (ConstrIndexManager::iterator constrPt = problemPtr()->probConstrSet()
       .begin(Active, 's');
       constrPt != problemPtr()->probConstrSet().end(Active, 's'); ++constrPt)
  {
    _mapConstrSeqNb2ConstrPtr[(*constrPt)->index()] = *constrPt;
  }

  //Added by Romain: Comment this part if you don't want to use dynamic constraint.
  for (ConstrIndexManager::iterator constrPt = problemPtr()->probConstrSet()
       .begin(Active, 'd');
       constrPt != problemPtr()->probConstrSet().end(Active, 'd'); ++constrPt)
  {
    _mapConstrSeqNb2ConstrPtr[(*constrPt)->index()] = *constrPt;
  }
#else
  _vectConstrPtr.clear();
  for (ConstrIndexManager::iterator constrPt =
       problemPtr()->probConstrSet().begin(Active, 's');
       constrPt != problemPtr()->probConstrSet().end(Active, 's');
       ++constrPt)
  {
    _vectConstrPtr.push_back(*constrPt);
  }
#endif

  /// Remove constr from formulation
  _interfacePtr->load();
  _interfacePtr->delRows(_indexSetOfRow2Delete);
  if (printL(6))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  _indexSetOfRow2Delete.clear();

  return;
}

void LPform::delVarFromFormulation()
{
  if (_indexSetOfCol2Delete.empty())
    /// No var to delete
    return;

  /// Reset index of existing variables (shift them by nb of deleted columns before them)
  int shift(0);
  //Commented by Romain
#if MAPCONSTR
  for (std::map<int, Variable *>::const_iterator it =
       _mapVarSeqNb2VarPtr.begin(); it != _mapVarSeqNb2VarPtr.end(); ++it)
  {
    if (_indexSetOfCol2Delete.count(it->first))
      shift++;
    else
      it->second->index(it->second->index() - shift);
  }
#else
  for (std::deque<Variable*>::const_iterator it = _vectVarPtr.begin();
       it != _vectVarPtr.end(); ++it)
  {
    if (_indexSetOfCol2Delete.count((*it)->index()))
    {
      ++shift;
    } else
    {
      (*it)->index((*it)->index() - shift);
    }
  }
#endif
#if MAPCONSTR
  _mapVarSeqNb2VarPtr.clear();
  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 's');
       varPt != problemPtr()->probVarSet().end(Active, 's'); varPt++)
  {
    _mapVarSeqNb2VarPtr[(*varPt)->index()] = *varPt;
  }

  //Added by Romain: Comment this part if you don't want to use dynamic or artificial var.
  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 'd');
       varPt != problemPtr()->probVarSet().end(Active, 'd'); varPt++)
  {
    _mapVarSeqNb2VarPtr[(*varPt)->index()] = *varPt;
  }

  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 'a');
       varPt != problemPtr()->probVarSet().end(Active, 'a'); varPt++)
  {
    _mapVarSeqNb2VarPtr[(*varPt)->index()] = *varPt;
  }

#else
  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 's');
       varPt != problemPtr()->probVarSet().end(Active, 's'); varPt++)
  {
    _vectVarPtr.push_back(*varPt);
  }

  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 'd');
       varPt != problemPtr()->probVarSet().end(Active, 'd'); varPt++)
  {
    _vectVarPtr.push_back(*varPt);
  }

  for (VarIndexManager::iterator varPt = problemPtr()->probVarSet().begin(Active, 'a');
       varPt != problemPtr()->probVarSet().end(Active, 'a'); varPt++)
  {
    _vectVarPtr.push_back(*varPt);
  }
#endif

  /// Remove var from formulation
  _interfacePtr->load();
  _interfacePtr->delCols(_indexSetOfCol2Delete);
  if (printL(6))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  _indexSetOfCol2Delete.clear();

  return;
}

void LPform::resetRhs(Constraint * constrPtr)
{
#ifdef PROBCONTAINERSET
  _rhsv.insert(ProbBound(constrPtr->index(), constrPtr->sense(), constrPtr->curRhs()));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _rhsv.push_back(ProbBound(constrPtr->index(), constrPtr->sense(), constrPtr->curRhs()));
#endif
  return;
}

void LPform::updateConstrRhsInFormulation()
{
  // No rhs to change
  if (_rhsv.empty())
    return;

  _interfacePtr->load();
  _interfacePtr->chgRhs(_rhsv);
  if (printL(7))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  _rhsv.clear();

  return;
}

void LPform::resetObjCoef(Variable * varPtr)
{
  bapcodInit().require(varPtr->index() >= 0,  "LPform::resetObjCoef(): colRef < 0");
  bapcodInit().require(varPtr->index() < _probColCnt, "LPform::resetObjCoef(): colRef >= _probColCnt");

#ifdef PROBCONTAINERSET
  _objectRow.insert(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _objectRow.push_back(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif

  return;
}

void LPform::replaceObjCoef(Variable * varPtr)
{
  bapcodInit().require(varPtr->index() >= 0,  "LPform::resetObjCoef(): colRef < 0");
  bapcodInit().require(varPtr->index() < _probColCnt, "LPform::resetObjCoef(): colRef >= _probColCnt");

#ifdef PROBCONTAINERSET
  _objectRow.erase(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
  _objectRow.insert(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _objectRow.erase(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
  _objectRow.push_back(ProbCoef(-1, varPtr->index(), varPtr->curCost() / _objScalFact));
#endif

  return;

}

void LPform::updateObjectiveInFormulation()
{
  if (_objectRow.empty())
    /// No rhs to change
    return;

  _interfacePtr->load();
  chgObjCoef(_objectRow);
  if (printL(7))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  _objectRow.clear();

  return;
}

void LPform::resetBounds(Variable * varPtr)
{
  if (printL(6))
    std::cout << " LPform::resetBounds var " << varPtr->name() << " lb = "
              << varPtr->curLb() << " ub = " << varPtr->curUb() << std::endl;

#ifdef PROBCONTAINERSET
  _bounds.insert(ProbBound(varPtr->index(), 'L', varPtr->curLb()));

  _bounds.insert(ProbBound(varPtr->index(), 'U', varPtr->curUb()));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
  _bounds.push_back(ProbBound(varPtr->index(), 'L', varPtr->curLb()));

  _bounds.push_back(ProbBound(varPtr->index(), 'U', varPtr->curUb()));
#endif
  return;
}

void LPform::updateBoundsInFormulation()
{
  if (_bounds.empty())
    /// No bds to change
    return;

  _interfacePtr->load();
  _interfacePtr->chgBds(_bounds);

  if (printL(7))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  _bounds.clear();

  return;
}

void LPform::scaleObjectif()
{
  ///@todo: implement scaling done a the formulation level

  /// Scaled Objectif coefficient
  if (_interfacePtr->param().UseObjScalingFact)
  {
    Double minC = BapcodInfinity;
    Double maxC = 0;
    Double averC = 0;
    Double absCoef = 0;
    {
      for (ProbCoefContainer::const_iterator obPtr = _objectRow.begin(); obPtr != _objectRow.end(); obPtr++)
      {
        absCoef = Fabs(obPtr->coef);
        averC += absCoef;
        if (absCoef < minC)
          minC = absCoef;

        if (absCoef > maxC)
          maxC = absCoef;
      }
      averC /= double(_objectRow.size());
    }

    const double threshold = 1e-4; /// used to be LOWPRECISION parameter

    /// Bring averC to 1
    _objScalFact = averC;
    if (minC / _objScalFact < threshold)
    {
      /// Bring averC to 10
      _objScalFact /= 10;
      if (minC / _objScalFact < threshold)
      {
        /// Bring averC to 100
        _objScalFact /= 10;
        if (minC / _objScalFact < threshold)
        {
          /// Bring averC to 1000
          _objScalFact /= 10;
          if (minC / _objScalFact < threshold)
          {
            /// Bring averC to 10000
            _objScalFact /= 10;
            if (minC / _objScalFact < threshold)
              /// Bring averC to 100000
              _objScalFact /= 10;
          }
        }
      }
    }

    if (printL(6))
      std::cout << "_objScalFact = " << _objScalFact << std::endl;

    ProbCoefContainer newObjRow;
    for (ProbCoefContainer::const_iterator obPtr = _objectRow.begin();
         obPtr != _objectRow.end(); ++obPtr)
    {
#ifdef PROBCONTAINERSET
      newObjRow.insert(ProbCoef(-1, obPtr->colRef, obPtr->coef / _objScalFact));
#endif

#if defined(PROBCONTAINERVECTOR) || defined(PROBCONTAINERDEQUE)
      newObjRow.push_back(ProbCoef(-1, obPtr->colRef, obPtr->coef / _objScalFact));
#endif
    }

    _objectRow = newObjRow;
  } else
    _objScalFact = 1;

  return;
}

void LPform::checkFormulation()
{
  if (printL(5))
    std::cout << "_probColCnt = " << _probColCnt << "  _probRowCnt = " << _probRowCnt << std::endl;

  if (printL(7))
    printMatrix();

  for (ProbCoefContainer::const_iterator obPtr = _objectRow.begin(); obPtr != _objectRow.end(); obPtr++)
  {
    if (printL(6))
      std::cout << "obPtr->colRef = " << obPtr->colRef << ", _probColCnt = "
                << _probColCnt << std::endl;

    bapcodInit().require(obPtr->colRef >= 0,
                         "Problem::checkProblem: negative obj colRef");

    bapcodInit().require(obPtr->colRef < _probColCnt,
                         "Problem::checkProblem: obj colRef out of range");
  }
  for (ProbRowMatrixContainer::const_iterator maPtr = _rowMatrix.begin(); maPtr != _rowMatrix.end(); maPtr++)
  {
    bapcodInit().require(maPtr->colRef >= 0,
                         "Problem::checkProblem: negative mat colRef");
    bapcodInit().require(maPtr->colRef < _probColCnt,
                         "Problem::checkProblem: mat colRef out of range");
    bapcodInit().require(maPtr->rowRef >= 0,
                         "Problem::checkProblem: negative mat rowRef");
    bapcodInit().require(maPtr->rowRef < _probRowCnt,
                         "Problem::checkProblem: mat rowRef out of range");
  }
  for (ProbColMatrixContainer::const_iterator maPtr = _colMatrix.begin(); maPtr != _colMatrix.end(); maPtr++)
  {
    bapcodInit().require(maPtr->colRef >= 0,
                         "Problem::checkProblem: negative mat colRef");
    bapcodInit().require(maPtr->colRef < _probColCnt,
                         "Problem::checkProblem: mat colRef out of range");
    bapcodInit().require(maPtr->rowRef >= 0,
                         "Problem::checkProblem: negative mat rowRef");
    bapcodInit().require(maPtr->rowRef < _probRowCnt,
                         "Problem::checkProblem: mat rowRef out of range");
  }
  for (ProbRhsContainer::const_iterator rhPtr = _rhsv.begin(); rhPtr != _rhsv.end(); rhPtr++)
  {
    bapcodInit().require(rhPtr->ref >= 0,
                         "Problem::checkProblem: negative rhs rowRef");
    bapcodInit().require(rhPtr->ref < _probRowCnt,
                         "Problem::checkProblem: rhs rowRef out of range");
    bapcodInit().require((rhPtr->sense == 'O') || (rhPtr->sense == 'G') || (rhPtr->sense == 'L')
                         || (rhPtr->sense == 'E'),
                         "Problem::checkProblem: Problem::checkProblem: constraint sense should be"
                         " 'O', 'G', 'L', or 'E' ");
  }
  for (ProbBoundContainer::const_iterator bdPtr = _bounds.begin(); bdPtr != _bounds.end(); bdPtr++)
  {
    bapcodInit().require(bdPtr->ref >= 0,
                         "Problem::checkProblem: negative bound colRef");
    bapcodInit().require(bdPtr->ref < _probColCnt,
                         "Problem::checkProblem: bound colRef out of range");
  }

  return;
}

void LPform::buildFormulation()
{
  scaleObjectif();

  /// Check consistency of indexing
  if (bapcodInit().testLevel() >= 6)
    checkFormulation();

  for (ProbRowMatrixContainer::const_iterator mPtr = _rowMatrix.begin(); mPtr != _rowMatrix.end(); mPtr++)
  {
#ifdef PROBMATRIXCONTAINERSET
    _colMatrix.insert(*mPtr);
#endif

#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
    _colMatrix.push_back(*mPtr);
#endif
  }
  _interfacePtr->loadFormulation(problemPtr()->name(), _minmax, _probColCnt, _probRowCnt,
                                 _mapSeqnb2Cname, _mapSeqnb2Rname, _objectRow,
                                 _colMatrix, _rhsv, _bounds, _varTypes,
                                 _varBrDirective, _soSets);

  /// Erase problem matrix representation
  clearColFormulationDataStruct();
  clearRowFormulationDataStruct();

  return;
}

void LPform::clearRowFormulationDataStruct()
{
  /// Erase current modifications to matrix
  _indexSetOfCol2Delete.clear();
  _rowMatrix.clear();
  _rhsv.clear();
  _mapSeqnb2Rname.clear();

  return;
}

void LPform::clearColFormulationDataStruct()
{
  /// Erase current modifications to matrix
  _indexSetOfCol2Delete.clear();
  _objectRow.clear();
  _colMatrix.clear();
  _bounds.clear();
  _mapSeqnb2Cname.clear();

  return;
}

void LPform::printForm(std::ostream& os)
{
  _interfacePtr->load();
#ifndef ADLIB
    _interfacePtr->LPwrite(_minmax, os);
#endif

  _interfacePtr->unload();
}

std::ostream& LPform::printMatrix(std::ostream& os) const
{
  os << " LPform::printMatrix(Problem name= " << problemPtr()->name() << "), objStatus= "
     << (int) _minmax << std::endl;
  os << "    _probColCnt = " << _probColCnt << "    _probRowCnt = " << _probRowCnt << std::endl;
  for (std::map<int, std::string>::const_iterator itr = _mapSeqnb2Rname.begin(); itr != _mapSeqnb2Rname.end(); itr++)
    std::cout << " row name[" << itr->first << "] = " << itr->second << std::endl;

  for (std::map<int, std::string>::const_iterator itc = _mapSeqnb2Cname.begin(); itc != _mapSeqnb2Cname.end(); itc++)
    std::cout << " col name[" << itc->first << "] = " << itc->second << std::endl;

  printSeq(_objectRow, "objective: ", 1, os);
  printSeq(_rowMatrix, "_rowMatrix: ", 1, os);
  printSeq(_colMatrix, "_colMatrix: ", 1, os);
  printSeq(_rhsv, "_rhsv: ", 1, os);
  printSeq(_bounds, "_bounds: ", 1, os);

  return (os);
}

void LPform::retrieveBasis(LpBasisRecord & basis, const bool markInVars, const bool markInConstrs)
{
  // get the basis from the solver
  std::vector<int> colStat, rowStat;
  _interfacePtr->getBasis(colStat, rowStat);

  // insert into the basis._varInBasis all variables not at the lower bound
  /// TO DO : take into account variables which are at their upper bounds
  int nbVarInBasicSol(colStat.size());
  for (int index = 0; index < nbVarInBasicSol; ++index)
  {
    if (colStat[index] != MathProgSolverInterface::AtLowerVarConstr)
    {
      Variable * varPtr;
#if MAPCONSTR
        varPtr = _mapVarSeqNb2VarPtr[index];
#else
        varPtr = _vectVarPtr[index];
#endif
        basis._varInBasis.push_back(VarPtr_MpFormIndexStatus(varPtr, colStat[index]));
        if (markInVars)
            varPtr->isInBasis(true);
      if (printL(5))
        std::cout << "LPform::retrieveBasis() ############## " << std::endl
                  << " index = " << index << " var = " << varPtr->name()
                  << " status = " << colStat[index] << std::endl;
    }
  }

  // insert into the basis._constrInBasis all constraints not at the lower bound
  int nbConstrInBasicSol(rowStat.size());
  for (int index = 0; index < nbConstrInBasicSol; ++index)
  {
    if (rowStat[index] != MathProgSolverInterface::AtLowerVarConstr)
    {
      Constraint * constrPtr;
#if MAPCONSTR
      constrPtr = _mapConstrSeqNb2ConstrPtr[index];
#else
      constrPtr = _vectConstrPtr[index];
#endif
      basis._constrInBasis.push_back(ConstrPtr_MpFormIndexStatus(constrPtr, rowStat[index]));
      if (markInConstrs)
        constrPtr->isInBasis(true);
      if (printL(5))
        std::cout << "LPform::retrieveBasis() ############## " << std::endl
                  << " index = " << index << " constr = " << constrPtr->name()
                  << " status = " << rowStat[index] << std::endl;
    }
  }

}

void LPform::loadBasis(const LpBasisRecord & basis)
{
  // allocate the status arrays
  std::vector<int> colStat(_interfacePtr->ncol(),
                           MathProgSolverInterface::AtLowerVarConstr);
  std::vector<int> rowStat(_interfacePtr->nrow(),
                           MathProgSolverInterface::AtLowerVarConstr);
  bool isValid = true;

  // set the status for all variables not at the lower bound
  int nbVarInBasis = basis._varInBasis.size();
  int index(-1);

  for (int cnt = 0; cnt < nbVarInBasis; ++cnt)
  {
    if (basis._varInBasis[cnt]._varPtr->vcIndexStatus() == Active)
    {
      index = basis._varInBasis[cnt]._varPtr->index();
      if ((index < 0) || (index >= _interfacePtr->ncol()))
      {
        std::cerr << "LPform::loadBasis() "
                  << basis._varInBasis[cnt]._varPtr->name()
                  << " has index out of bounds "
                  << ": varindex = "
                  << basis._varInBasis[cnt]._varPtr->index()
                  << " while ncol = "
                  << _interfacePtr->ncol()
#if MAPCONSTR
                  << " and mapVarSeqNb2VarPtr.size() = "
                  << _mapVarSeqNb2VarPtr.size()
                  << std::endl; //Commented by Romain
#else
                  << " and mapVarSeqNb2VarPtr.size() = " << _vectVarPtr.size() << std::endl;
#endif
        exit(1);
      } else
      {
        if (printL(5))
          std::cout << "LPform::loadBasis() ############## " << std::endl
                << " index = " << index << " var = "
                << basis._varInBasis[cnt]._varPtr->name() << " status = "
          << basis._varInBasis[cnt]._statusInBasicSol << std::endl;
      }

      colStat[index] = basis._varInBasis[cnt]._statusInBasicSol;
    } else
    {
      if (printL(5))
        std::cout << "LPform::loadBasis() ############## NOT ACTIVE "
              << std::endl << " index = " << index << " var = "
              << basis._varInBasis[cnt]._varPtr->name() << " status = "
        << basis._varInBasis[cnt]._statusInBasicSol << std::endl;

      isValid = false;
      break;
    }
  }

  if (isValid)
  {
    // set the status for all constraint not at the lower bound
    int nbConstrInBasis = basis._constrInBasis.size();
    for (int cnt = 0; cnt < nbConstrInBasis; ++cnt)
    {
      if (basis._constrInBasis[cnt]._constrPtr->vcIndexStatus() == Active)
      {
        index = basis._constrInBasis[cnt]._constrPtr->index();
        if ((index < 0) || (index >= _interfacePtr->nrow()))
        {
          std::cerr << "LPform::loadBasis() " << basis._constrInBasis[cnt]._constrPtr->name()
                    << " has index out of bounds " << std::endl;
          exit(1);
        }
        if (printL(5))
          std::cout << "LPform::loadBasis() ##############  " << std::endl
                    << " index = " << index << " constr = "
                    << basis._constrInBasis[cnt]._constrPtr->name()
                    << " status = " << basis._constrInBasis[cnt]._statusInBasicSol << std::endl;
        rowStat[index] = basis._constrInBasis[cnt]._statusInBasicSol;
      } else
      {
        if (printL(5))
          std::cout << "LPform::loadBasis() ############## NOT ACTIVE "
                    << std::endl << " index = " << index << " constr = "
                    << basis._constrInBasis[cnt]._constrPtr->name()
                    << " status = " << basis._constrInBasis[cnt]._statusInBasicSol << std::endl;
        isValid = false;
        break;
      }
    }
  }

  // if the basis is not valid, then it is better setting the default basis
  if (!isValid)
  {
    return;
    /// the code here causes Cplex to crash,
    /// TO DO : repair it
    for (int j = 0; j < (int) colStat.size(); j++)
      colStat[j] = MathProgSolverInterface::AtLowerVarConstr;
    for (int i = 0; i < (int) rowStat.size(); i++)
      rowStat[i] = MathProgSolverInterface::BasicVarConstr;
    if (printL(1))
        std::cout << "LPform::loadBasis: current basis is not valid" << std::endl;
  }

  // load the basis
  _interfacePtr->setBasis(colStat, rowStat);
}

// *************************
// ***** class MIPform *****
// *************************

MIPform::MIPform(Problem * problemPtr) :
    LPform(problemPtr, false)
{
    MathProgSolverBuilder mathProgSolverBuilder;
    _interfacePtr = mathProgSolverBuilder.buildMipMathProgSolverInterface(problemPtr->probConfPtr()->bapcodInitPtr(),
                                                                          problemPtr->probConfPtr()->param().solverName,
                                                                          problemPtr->ref(), problemPtr->name());

    if (_interfacePtr == nullptr)
    {
        std::string solverName = problemPtr->probConfPtr()->param().solverName();
        std::string shortName = solverName.substr(0,solverName.length() - 7);
        std::cerr << "BaPCod error : solver " << shortName << " is not found!" << std::endl
                  << "Please define " << shortName << "_ROOT environment variable before running cmake" << std::endl;
        exit(1);
    }

    _interfacePtr->setLPoptimalityTolerance(param().MasterMipSolverReducedCostTolerance());
    _interfacePtr->setMultiThread(param().MipSolverMultiThread()); /// added by Ruslan
    _interfacePtr->setTimeLimit(param().MipSolverMaxTime()); /// added by Ruslan
    _interfacePtr->setMaxBBNodes(param().MipSolverMaxBBNodes());
    _interfacePtr->setRelativeMipGapLimit(param().relOptimalityGapTolerance());
}

bool MIPform::solve(const double & BarrierConvergenceTolerance,
                    const double & rightHAndSideZeroTol,
                    const double & reducedCostTolerance,
                    const char & flag, const bool & ifPrint,
                    const SolutionStatus & requiredStatus, Double & objVal,
                    Double & primalBound, Double & dualBound,
                    VarPtrSet & inPrimalSol, ConstrPtrSet & inDualSol,
                    const bool & preprocessorOn, const bool & probingOn,
                    const bool & automaticCuttingPlanesOn,
                    const char & solverSelection)
{
  bool foundSol(false);
  _interfacePtr->load();

  if (ifPrint)
    _interfacePtr->LPwrite(_minmax);

  if (printL(7))
    _interfacePtr->MPSwrite();

  if (flag == 'r')
    _interfacePtr->reset();

   _interfacePtr->setSolveFromScratch();

  _interfacePtr->optimise(_minmax, BarrierConvergenceTolerance, rightHAndSideZeroTol, reducedCostTolerance,
                          preprocessorOn, probingOn, automaticCuttingPlanesOn, solverSelection);

  _interfacePtr->getOptimStatus(_status, _mipStatus);
  if (printL(6))
  {
    std::cout << "status() = " << status() << std::endl;
    std::cout << "requiredStatus = " << requiredStatus << std::endl;
  }

  bool statusIsOk = bapcodInit().require(status().intersects(requiredStatus),
                                         "MIPform::solve(): Formulation could not be solved according "
                                         "to prescribed status", ProgStatus::run);

  if (statusIsOk)
    {
      foundSol = status().intersects(SolutionStatus(SolutionStatus::Optimum,
                                                          SolutionStatus::PrimalFeasSolFound));
      if (foundSol)
      {
        LPform::retrieveSol('p', ifPrint, inPrimalSol, inDualSol);
      }

      setBounds(objVal, primalBound, dualBound);
    }
  else
    {
      _interfacePtr->MPSwrite();
      std::cout << "MIPform::solve() status = " << status() << std::endl;
      return false;
    }

  _interfacePtr->unload();

  return (foundSol);
}

void MIPform::setBounds(Double & objVal, Double & primalBound, Double & dualBound)
{
  if (status().count(SolutionStatus::Optimum))
  {
    _interfacePtr->getObjVal(objVal);
    _interfacePtr->getPrimalBound(primalBound);
    _interfacePtr->getDualBound(dualBound);

    /// Scaling
    objVal *= _objScalFact;
    primalBound *= _objScalFact;
    dualBound *= _objScalFact;
    if (printL(4))
      std::cout << "Solution MIP status " << status() << " objVal = " << objVal
                << " primalBound = " << primalBound << " dualBound = " << dualBound << std::endl;

    primalBound = dualBound = objVal;
    if (printL(4))
      std::cout << "Solution MIP status " << status() << " objVal = " << objVal
                << " primalBound = " << primalBound << " dualBound = " << dualBound << std::endl;

    return;
  }

  if (status().count(SolutionStatus::PrimalFeasSolFound))
  {
    _interfacePtr->getObjVal(objVal);
    _interfacePtr->getPrimalBound(primalBound);

    /// Scaling
    objVal *= _objScalFact;
    primalBound *= _objScalFact;
    if (printL(4))
      std::cout << "Solution MIP status " << status() << " objVal = " << objVal
                << " primalBound = " << primalBound << " dualBound = " << dualBound << std::endl;
    return;
  }

  if (status().count(SolutionStatus::DualFeasSolFound))
  {
    _interfacePtr->getObjVal(objVal);
    _interfacePtr->getDualBound(dualBound);

    /// Scaling
    objVal *= _objScalFact;
    dualBound *= _objScalFact;
    if (printL(4))
      std::cout << "Solution MIP status " << status() << " objVal = " << objVal
                << " primalBound = " << primalBound << " dualBound = " << dualBound << std::endl;
    return;
  }
  if (status().count(SolutionStatus::Infeasible) || status().count(SolutionStatus::UnSolved))
  {
    primalBound = dualBound = _minmax * BapcodInfinity;
    return;
  }
  if (status().count(SolutionStatus::Unbounded))
  {
    primalBound = dualBound = -_minmax * BapcodInfinity;
    return;
  }

  return;
}

void MIPform::fillDataStruct(Variable * varPtr)
{
  /// Set the lp part of the variable
  LPform::fillDataStruct(varPtr);
  _varTypes.insert(ProbType(varPtr->index(), varPtr->type()));
  if (varPtr->type() != 'C')
  {
    _varBrDirective.insert(ProbIntC(varPtr->index(), varPtr->directive(), varPtr->priority()));
  }

  return;
}

void MIPform::setLazyConstraintsCallback(MasterConf * masterConfPtr)
{
    _interfacePtr->setLazyConstraintsCallback(masterConfPtr);
}

//added by Issam (Only tested for cplex).
void MIPform::resetMIPpartOfFormulation(const Double & cutOffValue, const bool exactSolution)
{
#if MAPCONSTR
  for (int varSeqNb = 0; varSeqNb < _mapVarSeqNb2VarPtr.size(); varSeqNb++)
  {
    Variable * varPtr = _mapVarSeqNb2VarPtr[varSeqNb];
#else
  for (int varSeqNb = 0; varSeqNb < _vectVarPtr.size(); varSeqNb++)
  {
    Variable * varPtr = _vectVarPtr[varSeqNb];
#endif
    _varTypes.insert(ProbType(varPtr->index(), varPtr->type()));
    if (varPtr->type() != 'C')
    {
      _varBrDirective.insert(
                             ProbIntC(varPtr->index(), varPtr->directive(), varPtr->priority()));
    }
  }

  // This instruction automatically transforms an LP object of cplex to and MIP object.
  _interfacePtr->chgColType(_varTypes);

  /// commented by Ruslan, as this may cause miss of optimal solution in enumerated MIP, we need to see why
  //_interfacePtr->addDirectives(_varBrDirective);

  _interfacePtr->setUpperCutOffValue(cutOffValue);
  _interfacePtr->setSolveFromScratch(exactSolution);

  _varTypes.clear();
  _varBrDirective.clear();
}

void MIPform::resetAfterMIP()
{
  _interfacePtr->resetUpperCutOffValue();
  _interfacePtr->resetSolveFromScratch();
  _interfacePtr->removeLazyConstraintsCallback();
}

void MIPform::setTimeLimit(const double timeLimit)
{
  _interfacePtr->setTimeLimit(timeLimit);
}

void MIPform::addVar2Formulation()
{
  if (_objectRow.empty())
    /// Nothing to add in formulation to update
    return;

  _interfacePtr->load();
  _interfacePtr->addCols(_objectRow, _colMatrix, _bounds, _mapSeqnb2Cname);
  if (!_varTypes.empty())
    _interfacePtr->chgColType(_varTypes);

  /// commented by Ruslan, as this may cause miss of optimal solution in enumerated MIP, we need to see why
//  if (!_varBrDirective.empty())
//    _interfacePtr->addDirectives(_varBrDirective);

  if (printL(6))
    _interfacePtr->LPwrite(_minmax);

  _interfacePtr->unload();
  clearColFormulationDataStruct();

  return;
}

void MIPform::checkFormulation()
{
  LPform::checkFormulation();

  for (set<ProbType>::const_iterator tPtr = _varTypes.begin(); tPtr != _varTypes.end(); tPtr++)
  {
    bapcodInit().require(tPtr->ref >= 0,
                         "MipProblem::checkProblem: negative type colRef");
    bapcodInit().require(tPtr->ref < _probColCnt,
                         "MipProblem::checkProblem: type colRef out of range");
  }
  for (set<ProbIntC>::const_iterator dPtr = _varBrDirective.begin(); dPtr != _varBrDirective.end(); dPtr++)
  {
    bapcodInit().require(dPtr->ref >= 0,
                         "MipProblem::checkProblem: negative _varBrDirective colRef");
    bapcodInit().require(dPtr->ref < _probColCnt,
                         "MipProblem::checkProblem: _varBrDirective colRef out of range");
    bapcodInit().require((dPtr->type == 'U') || (dPtr->type == 'D'),
                         "MipProblem::checkProblem: directive should be 'U' or 'D'");
  }

  for (set<ProbSetCoef>::const_iterator sPtr = _soSets.begin(); sPtr != _soSets.end(); sPtr++)
  {
    if (printL(6))
      cout << "MipProblem::checkProblem:  set[" << sPtr->setRef
           << "] includes col " << sPtr->colRef << " with coef " << sPtr->coef << std::endl;

    bapcodInit().require(sPtr->setRef >= 0,
                         "MipProblem::checkProblem: negative set index");
    bapcodInit().require(sPtr->colRef >= 0,
                         "MipProblem::checkProblem: negative col index in set description");
    bapcodInit().require(sPtr->colRef < _probColCnt,
                         "MipProblem::checkProblem: set col index out of range");
  }

  return;
}

void MIPform::buildFormulation()
{
  scaleObjectif();

  /// Check consistency of indexing
  if (bapcodInit().testLevel() >= 6)
    checkFormulation();

  for (ProbRowMatrixContainer::const_iterator mPtr = _rowMatrix.begin(); mPtr != _rowMatrix.end(); mPtr++)
  {
#ifdef PROBMATRIXCONTAINERSET
    _colMatrix.insert(*mPtr);
#endif

#if defined(PROBMATRIXCONTAINERVECTOR) || defined(PROBMATRIXCONTAINERDEQUE)
    _colMatrix.push_back(*mPtr);
#endif
  }

  _interfacePtr->loadFormulation(problemPtr()->name(), _minmax, _probColCnt, _probRowCnt,
                                 _mapSeqnb2Cname, _mapSeqnb2Rname, _objectRow,
                                 _colMatrix, _rhsv, _bounds, _varTypes,
                                 _varBrDirective, _soSets);

  /// Erase problem matrix representation
  clearColFormulationDataStruct();
  clearRowFormulationDataStruct();

  return;
}

void MIPform::clearColFormulationDataStruct()
{
  LPform::clearColFormulationDataStruct();

  /// Erase current modifications to matrix
  _varTypes.clear();
  _varBrDirective.clear();
  _soSets.clear();

  return;
}

std::ostream& MIPform::printMatrix(std::ostream& os) const
{
  LPform::printMatrix(os);
  printSeq(_varBrDirective, "_varBrDirective: ", 1, os);

  return (os);
}
