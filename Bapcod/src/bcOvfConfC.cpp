/**
 *
 * This file bcOvfConfC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcOvfConfC.hpp"

using namespace std;

/**
 * Generic code
 */

OvfConf::OvfConf(Model * modelPtr, Problem * problemClassPtr) :
    ProbConfig(ProbConfig::ovf, modelPtr, "ovf", IndexCell(), modelPtr->bestPrimalBound(),
               modelPtr->bestDualBound(), problemClassPtr)
{
  return;
}

OvfConf::~OvfConf()
{
}

Constraint * OvfConf::castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddConstraint() should not be called");
  return NULL;
}

InstanciatedConstr * OvfConf::castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddConstraint() should not be called");

  return NULL;

}

Variable * OvfConf::castAndAddVariable(Variable * varPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddVariable() should not be called");
  return NULL;
}

InstanciatedVar * OvfConf::castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddVariable() should not be called");
  return NULL;
}

/// Set initial variables and constraints
  
void OvfConf::prepareProbConfig()
{
  map<string, int> varCountMap; //gen name to count
  map<string, int> constrCountMap; //gen name to count

  if (!progStatus().doRun()) return;

  if (param().ovfSolMode().status() == SolutionMethod::none)
    return;

  if (_isPrepared) return;
  _isPrepared = true;


  if (printL(5))
    std::cout << "OvfConf::prepareProbConfig() " << name() << " _name2GenericVarPtrMap.size() "
              << _name2GenericVarPtrMap.size() << std::endl;

  for (std::map < std::string, GenericVar * >::const_iterator it = _name2GenericVarPtrMap.begin();
          it != _name2GenericVarPtrMap.end();
          ++it)
  {
    if (printL(5))
      std::cout << " setup BranchingConstr for "
            << it->first
            << std::endl;
    it->second->setupGenericBranchingConstr();
  }

  int lowPrintLevel = 6;

  if (_modelPtr->masterConfPtr() != NULL)
  {
    /// Set pure master variable
    for (InstMasterVarPtrSet::const_iterator pmvPt = _modelPtr->masterConfPtr()->setOfPureMasterVar().begin();
            pmvPt != _modelPtr->masterConfPtr()->setOfPureMasterVar().end();
            pmvPt++)
    {
      if (printL(lowPrintLevel))
        std::cout << "OvfConf::prepareProbConfig() add mastVar " << (*pmvPt)->name() << std::endl;
      
      varCountMap[(*pmvPt)->genericName()] += 1;
      _pcVarPtrList.push_back(new OvfVar(_modelPtr->masterConfPtr(), *pmvPt));
    }

    /// Set pure master constraint
    for (InstMasterConstrPtrSet::const_iterator pmcPt = _modelPtr->masterConfPtr()->setOfPureMasterConstr().begin();
            pmcPt != _modelPtr->masterConfPtr()->setOfPureMasterConstr().end();
            pmcPt++)
    {
      if (printL(lowPrintLevel))
        std::cout << "OvfConf::prepareProbConfig() add mastConstr " << (*pmcPt)->name() << std::endl;

      constrCountMap[(*pmcPt)->genericName()] += 1;
      _pcConstrPtrList.push_back(new OvfConstr(_modelPtr->masterConfPtr(), *pmcPt));
    }

    /// Set colGenSp variable and constr
    for (std::vector<ColGenSpConf *>::const_iterator spconfPt = _modelPtr->masterConfPtr()->colGenSubProbConfPts().begin();
            spconfPt != _modelPtr->masterConfPtr()->colGenSubProbConfPts().end();
            ++spconfPt)
    {

      bool setupVarNeeded(!zeroTest((*spconfPt)->fixedCost())
              || (((*spconfPt)->lowerBoundPtr() != NULL) && (*((*spconfPt)->lowerBoundPtr()) > 0)));

      int nbCopy(0);
      if ((*spconfPt)->upperBoundPtr() != NULL)
        nbCopy = *((*spconfPt)->upperBoundPtr());
      else
      {
        bapcodInit().require((*spconfPt)->lowerBoundPtr() != NULL
                , "OvfConf::prepareProbConfig() spConf should be bounded");
        nbCopy = *((*spconfPt)->lowerBoundPtr());
      }

      if (printL(lowPrintLevel))  
      std::cout << "OvfConf::prepareProbConfig() consider sp " << (*spconfPt)->name()
              << " nbCopy = " << nbCopy << std::endl;


      for (int ct = 0; ct < nbCopy; ct++)
      {
        if (setupVarNeeded) 
        {
          varCountMap["_setupVar"] += 1;
          _pcVarPtrList.push_back(new SpSetupOvfVar(*spconfPt, ct));
        }

        for (std::list< InstanciatedVar *>::const_iterator spvPt = (*spconfPt)->iVarPtrList().begin();
                spvPt != (*spconfPt)->iVarPtrList().end(); spvPt++)
        {
          if (printL(lowPrintLevel))  
          std::cout << "OvfConf::prepareProbConfig() copy " << ct << " add spVar " << (*spvPt)->name() << std::endl;

          varCountMap[(*spvPt)->genericName()] += 1;
          _pcVarPtrList.push_back(new OvfVar(*spconfPt, *spvPt, ct));

          if (setupVarNeeded)
          {
            /// Add component bound constraint
            if ((*spvPt)->curLb() != 0)
            {
              constrCountMap["_spVarLb"] += 1;
              _pcConstrPtrList.push_back(new SpVarLbOvfConstr(*spconfPt, *spvPt, ct));
            }

            if ((*spvPt)->curUb() != 0)
            {
              constrCountMap["_spVarUb"] += 1;
              _pcConstrPtrList.push_back(new SpVarUbOvfConstr(*spconfPt, *spvPt, ct));
            }
          }
        }

        for (std::list< InstanciatedConstr *>::const_iterator spcPt = (*spconfPt)->iConstrPtrList().begin();
                spcPt != (*spconfPt)->iConstrPtrList().end();
                spcPt++)
        {
          if (printL(lowPrintLevel))  
          std::cout << "OvfConf::prepareProbConfig() copy " << ct << " add spConstr " << (*spcPt)->name() << std::endl;

          constrCountMap[(*spcPt)->genericName()] += 1;
          _pcConstrPtrList.push_back(new OvfConstr(*spconfPt, *spcPt, ct));
        }

      }
      /// Add sp lowerBound constraint
      if (((*spconfPt)->lowerBoundPtr() != NULL) && (*((*spconfPt)->lowerBoundPtr()) > 0))
      {
        constrCountMap["_spLb"] += 1;
        _pcConstrPtrList.push_back(new SpLbOvfConstr(*spconfPt));
      }
    }
  }    /// Set pure ovf variable and constraints
  else
  {
    /// Cast master instantiated variable in type InstOvfVar
    OvfVar * instOvfVarPtr(NULL);
    for (std::list<InstanciatedVar *>::const_iterator iVarIterator = _iVarPts.begin();
            iVarIterator != _iVarPts.end();
            ++iVarIterator)
    {
      {
        if (printL(lowPrintLevel))
          std::cout << "OvfConf::prepareProbConfig(): Add Pure Var "
                << (*iVarIterator)->name() << std::endl;

        instOvfVarPtr = new OvfVar(this, *iVarIterator);
        varCountMap[instOvfVarPtr->genericName()] += 1;
        _pcVarPtrList.push_back(instOvfVarPtr);
        delete (*iVarIterator);
      }
    }
    _iVarPts.clear();

    /// Cast master instantiated constraint in type InstOvfConstr
    OvfConstr * instOvfConstrPtr(NULL);
    for (std::list<InstanciatedConstr *>::const_iterator iConstrIterator = _iConstrPts.begin();
            iConstrIterator != _iConstrPts.end();
            ++iConstrIterator)
    {
        if (printL(lowPrintLevel))
          std::cout << "OvfConf::prepareProbConfig(): Add Pure Constr " << (*iConstrIterator)->name() << std::endl;

        instOvfConstrPtr = new OvfConstr(this, *iConstrIterator);
        constrCountMap[instOvfConstrPtr->genericName()] += 1;
        _pcConstrPtrList.push_back(instOvfConstrPtr);
        delete (*iConstrIterator);
    }
    _iConstrPts.clear();
  }

  probPtr()->defineFormulation();

  _pcVarPtrList.insert(_pcVarPtrList.end(), _iVarPts.begin(), _iVarPts.end());
  _pcConstrPtrList.insert(_pcConstrPtrList.end(), _iConstrPts.begin(), _iConstrPts.end());
  probPtr()->addVarSet(_pcVarPtrList, 1, 0);
  //std::cout << "_pcConstrPtrList.size() " << _pcConstrPtrList.size() << std::endl;
  probPtr()->addConstrSet(_pcConstrPtrList, 1, 0);

  probPtr()->buildProblem();

  if (printL(3)) probPtr()->print();

  if(printL(1))
  {
    
    cout << "OvfConf::prepareProbConfig() Number of variables and constraints: " << endl;
    cout << "  Total Number of variables : " << _pcVarPtrList.size() << endl;
    cout << "  Total Number of constraints : " << _pcConstrPtrList.size() << endl;
    cout << "  --Details--" << endl;
    for(map<string, int>::iterator it = varCountMap.begin(); it != varCountMap.end(); it++)
    {
      cout << "  Number of variables " << it->first << " : " << it->second << endl;
    }
    for(map<string, int>::iterator it = constrCountMap.begin(); it != constrCountMap.end(); it++)
    {
      cout << "  Number of constraints " << it->first << " : " << it->second << endl;
    }
  }
  

  return;
}

Solution * OvfConf::solvePC()
{
  if (param().ovfSolMode().status() == SolutionMethod::none)
    return (NULL);

  if (_probPtr == NULL)
    return (NULL);

  Time start;

  _primalSolPtr = ProbConfig::solvePC(printL(-1));


  if (_primalSolPtr != NULL)
  {
    if (printL(1))
      std::cout << "bcRecBestDb " << _primalSolPtr->cost() << std::endl;

    updatePrimalIncBound(Bound(_primalSolPtr->cost(), _modelPtr->objectiveSense()));
    /// The following two lines are not necessarily correct! (hsen) 20/04/2017
    /// If Ovf::solvePC stops due to a limit/tolerance/etc. dualbound might be different.
    /// I didn't modified it since it was like this. But the dualbound (bestbound)
    /// should be retrieved directly from the solver.
    bapcodInit().statistics().incrValue("bcRecBestDb", _primalSolPtr->cost());
    updateDualIncBound(Bound(_primalSolPtr->cost(), _modelPtr->objectiveSense())); // Compelementary slackness
  
    _primalSolPtr = getDissagregatedSolution(_primalSolPtr);
  }

  double time = start.getElapsedTime_dbl();
  if (printL(3))
    std::cout << "bcTimeOvfMPsol" << time << std::endl;

  bapcodInit().statistics().incrTimer("bcTimeOvfMPsol", time);

  return _primalSolPtr;
}

Solution * OvfConf::getDissagregatedSolution(Solution * solPtr)
{
  VarPtr2DoubleMap mastSolutionMap;
  
  //In the following map,
  //the keys are pairs containing the sp probConf pointer and ovf cnt.
  //the values are solVarVal maps with the vars being original instanciated vars.
  map< pair<ProbConfig*, int>, VarPtr2DoubleMap > spSolutionsMap;
  
  for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin();
                it != solPtr->solVarValMap().end(); it++)
  {
    if (printL(5))
            std::cout << "curDisaggregateSol[" << it->first->name() << "] = " << it->second << std::endl;
    
    OvfVar * ovfVarPtr = NULL;
    if (printL(5))
        cout << "OVF variable has name " << it->first->name() << " and type " << it->first->name() << endl;
    if(it->first->isTypeOf(VcId::OvfVarMask))
    {      
      ovfVarPtr = static_cast<OvfVar*>(it->first);    
    }
    else
    {
      std::cerr << "OvfConf::getDissagregatedSolution Error : all variables should have" << "ovfVarPtr type " << endl;
      exit(1);
    }            
    Double val = it->second;
    
    Variable* orignalVarPtr = ovfVarPtr->originatingVarPtr();    
    
    if(orignalVarPtr->isTypeOf(VcId::InstMasterVarMask))
    {
      mastSolutionMap[orignalVarPtr] = val;
    }
    else
    {
      ProbConfig* probConfPtr = orignalVarPtr->probConfPtr();
      int counter = ovfVarPtr->cnt();    
      
      spSolutionsMap[make_pair(probConfPtr, counter)][orignalVarPtr] = val;
    }
  }
  
  Solution * pureMasterSolPtr = new Solution(this->mastConfPtr(), NULL);
  for (VarPtr2DoubleMap::const_iterator it = mastSolutionMap.begin(); it != mastSolutionMap.end(); it++)
  {
    pureMasterSolPtr->includeVar(it->first, it->second, false);    
  }
  pureMasterSolPtr->resetCost();

  Solution * previousSolptr = pureMasterSolPtr;
  
  for (map< pair<ProbConfig*, int>, VarPtr2DoubleMap >::const_iterator solIt = 
          spSolutionsMap.begin() ; solIt != spSolutionsMap.end() ; solIt++)
  {
    ProbConfig* spConfPtr = solIt->first.first;
    
    Solution * projectedSPSolPtr = new Solution(spConfPtr, previousSolptr);
    previousSolptr = projectedSPSolPtr;
    
    for (VarPtr2DoubleMap::const_iterator it = solIt->second.begin();
                it != solIt->second.end(); it++)
    {
      projectedSPSolPtr->includeVar(it->first, it->second, false);    
    }
    
    projectedSPSolPtr->resetCost();

    projectedSPSolPtr->multiplicity(1);
  }
  
  return pureMasterSolPtr;
}

std::ostream& OvfConf::print(std::ostream& os) const
{
  os << " OvfConf name = " << name() << std::endl;

  if (probPtr() != NULL)
    probPtr()->print(os);

  return (os);
}
