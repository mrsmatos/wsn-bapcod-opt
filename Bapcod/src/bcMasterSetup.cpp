/**
 *
 * This file bcMasterSetup.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMasterConfC.hpp"
#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcGenMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcModelC.hpp"

#include "bcSpVarConstrC.hpp"
#include "bcOvfVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcModelFormulationC.hpp"

#include "bcNetworkFlowC.hpp"

using namespace std;
using namespace VcIndexStatus;

Constraint * MasterConf::castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately)
{
  bapcodInit().check(1, "ProbConfig::castAndAddConstraint() should not be called");
  return constrPtr;
}

InstanciatedConstr * MasterConf::castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately)
{
  if (printL(6))
  {
    std::cout << "cast instantiated constraint of the master into type InstMasterConstr for constr : "
              << iconstrPtr->name() << std::endl;

    if (printL(7))
      std::cout << iconstrPtr << std::endl;
  }

  InstMasterConstr * instMasterConstrPtr =
    dynamic_cast<InstMasterConstr *> (iconstrPtr);

  if (instMasterConstrPtr == NULL) /// cast the constraint
  {
    if (printL(6))
    {
      std::cout << "castint required for : " << iconstrPtr->name() << std::endl;
    }
    NonLinearInstConstr * nonLinearInstConstrPtr =
      dynamic_cast<NonLinearInstConstr *> (iconstrPtr);

    if (nonLinearInstConstrPtr != NULL)
      instMasterConstrPtr = new NonLinearInstMastConstr(nonLinearInstConstrPtr);
    else
      instMasterConstrPtr = new InstMasterConstr(iconstrPtr);

    delete (iconstrPtr);
  }
  else
  {
    if (printL(6))
    {
      std::cout << "castint NOT required for : " << iconstrPtr->name()
                << std::endl;
    }
  }

  if (instMasterConstrPtr->flag() == 's')
    _setOfPureMasterConstr.insert(instMasterConstrPtr);

  if (insertImmediately)
  {
    std::list<Constraint *> constrPts;
    constrPts.push_back(instMasterConstrPtr);
    probPtr()->addConstrSet(constrPts, 1, 2);
  }
  else
  {
    if (instMasterConstrPtr->kind() == 'R')
      _pcDelayedConstrPtrList.push_back(instMasterConstrPtr);
    else if (instMasterConstrPtr->flag() == 's')
      _pcConstrPtrList.push_back(instMasterConstrPtr);
  }
  return instMasterConstrPtr;
}

void MasterConf::insertPureVar(InstanciatedVar * ivarPtr)
{
  _setOfPureMasterVar.insert(dynamic_cast<InstMasterVar *> (ivarPtr));
  return;
}

void MasterConf::insertPureConstr(InstanciatedConstr * iconstrPtr)
{
  _setOfPureMasterConstr.insert(dynamic_cast<InstMasterConstr *> (iconstrPtr));
  return;
}

Variable * MasterConf::castAndAddVariable(Variable * varPtr,
    const bool & insertImmediately)
{
  bapcodInit().check(1,"MasterConf::castAndAddVariable() should not be called");
  return varPtr;
}

InstanciatedVar * MasterConf::castAndAddVariable(InstanciatedVar * ivarPtr,
    const bool & insertImmediately)
{
  InstMasterVar * instMasterVarPtr = dynamic_cast<InstMasterVar *> (ivarPtr);

  if (instMasterVarPtr == NULL) /// cast the variable
  {
    instMasterVarPtr = new InstMasterVar(ivarPtr);
    delete (ivarPtr);
  }

  if (insertImmediately)
  {
    probPtr()->addVar(instMasterVarPtr);
  }
  else
  {
    _pcVarPtrList.push_back(instMasterVarPtr);
  }

  return instMasterVarPtr;
}

void MasterConf::insertColGenSpConf(ColGenSpConf * cgSpConfPtr)
{
  _colGenSubProbConfPts.push_back(cgSpConfPtr);
  _colGenSubProbFormList.push_back(BcFormulation(cgSpConfPtr));

  return;
}

void MasterConf::addArtVar(MissingColumn * artVarPtr)
{
  _nonStabStaticArtVarPtrList.push_back(artVarPtr);

  return;
}

void MasterConf::addGlobalArtVar(GlobalArtificialVar * artVarPtr)
{
  _nonStabStaticArtVarPtrList.push_back(artVarPtr);

  return;
}

void MasterConf::addNonStabilizationLocalArtVar(LocalArtificialVar * artVarPtr)
{
  _nonStabStaticArtVarPtrList.push_back(artVarPtr);
  return;
}

void MasterConf::prepareProbConfig()
{
  Time bcTimeMastPrepareProbConfig;

  if (!progStatus().doRun())
    return;

  if (param().masterSolMode().status() == SolutionMethod::none)
  {

    cout << "_isPrepared " << _isPrepared << endl;

    if (!_isPrepared) {
      _pcConstrPtrList.insert(_pcConstrPtrList.end(), _iConstrPts.begin(), _iConstrPts.end());
      _pcVarPtrList.insert(_pcVarPtrList.end(), _iVarPts.begin(), _iVarPts.end());
    }

    _isPrepared = true;
    return;
  }

  if (_isPrepared)
    return;

  _isPrepared = true;

  if (printL(5))
    std::cout << " MasterConf::prepareProbConfig() " << name()
              << " _name2GenericVarPtrMap.size() " << _name2GenericVarPtrMap.size() << std::endl;

  for (std::map<std::string, GenericVar *>::const_iterator it = _name2GenericVarPtrMap.begin();
       it != _name2GenericVarPtrMap.end(); ++it)
  {
    if (printL(5))
      std::cout << " setup master BranchingConstr for " << it->first << std::endl;
    it->second->setupGenericBranchingConstr();
  }

  if ((param().ovfSolMode().status() == SolutionMethod::none)
       && (param().masterSolMode().status() == SolutionMethod::none))
    return;

  probPtr()->defineFormulation();

    /// Cast master instantiated variable in type InstMasterVar
  for (std::list<InstanciatedVar *>::const_iterator it = _iVarPts.begin(); it != _iVarPts.end(); ++it)
    castAndAddVariable(*it, false);

  _iVarPts.clear();

  /// Cast instantiated constraint of the master into type InstMasterConstr

  std::list<InstanciatedConstr *> tempListOfMasterConstr2Upcast(_iConstrPts.begin(), _iConstrPts.end());
  _iConstrPts.clear();
  std::list<InstanciatedConstr *>::iterator iCnstrPtrIt;

  /// Add global artificial columns to master problem
  GlobalArtificialVar *posGlobalArtVarPtr(NULL);
  GlobalArtificialVar *negGlobalArtVarPtr(NULL);

  if ((param().mastInitMode().status() == MasterInitMode::globalArtCol)
      || (param().mastInitMode().status() == MasterInitMode::incSolColAndGac)
      || (param().mastInitMode().status() == MasterInitMode::localAndGlobAc)
      || (param().mastInitMode().status() == MasterInitMode::defaultInit)) {
    posGlobalArtVarPtr = new GlobalArtificialVar(_modelPtr, _modelPtr->artVarCost(), 'G',
                                                 _modelPtr->objectiveSense());
    if (printL(6))
      std::cout << "add globalArtVar to master  " << posGlobalArtVarPtr->name() << std::endl;
    addGlobalArtVar(posGlobalArtVarPtr);

    probPtr()->recordPosGlobalArtVar(posGlobalArtVarPtr);

    negGlobalArtVarPtr = new GlobalArtificialVar(_modelPtr, _modelPtr->artVarCost(), 'L',
                                                 _modelPtr->objectiveSense());
    if (printL(6))
      std::cout << "add globalArtVar to master  " << negGlobalArtVarPtr->name() << std::endl;
    addGlobalArtVar(negGlobalArtVarPtr);

    probPtr()->recordNegGlobalArtVar(negGlobalArtVarPtr);
  }

  ///  Add master convexity constraints
  _convexityGenConstrPtrList.push_back(new ConvexityGenConstr(_modelPtr, this,
                                                                 colGenSubProbConfPts()));

  /// Add convexity constraint UB into _iConstrPt
  for (iCnstrPtrIt = _iConstrPts.begin(); iCnstrPtrIt != _iConstrPts.end(); iCnstrPtrIt++) {
    if (printL(6))
      std::cout << "MasterConf::prepareProbConfig(): add convexity constraint and the like " << (*iCnstrPtrIt)->name()
                << std::endl;

    _pcConstrPtrList.push_back(*iCnstrPtrIt);
  }

  for (iCnstrPtrIt = tempListOfMasterConstr2Upcast.begin(); iCnstrPtrIt != tempListOfMasterConstr2Upcast.end();
       iCnstrPtrIt++)
    castAndAddConstraint(*iCnstrPtrIt);

  bapcodInit().statistics().setCounter("bcMaxPrimalSpaceSize", _pcConstrPtrList.size());

  if (printL(6))
    std::cout << "add master Constr: nb = " << _pcConstrPtrList.size() << std::endl;


  /// Prepare subproblems which include upcasting into subProblem variables and constraints
  /// (and hence should be done before the buildMemebership done in addConstr and addVar)

  _problemList.push_back(probPtr()); /// added by Ruslan
  std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt;
  for (cgSpConfPtrIt = colGenSubProbConfPts().begin(); cgSpConfPtrIt != colGenSubProbConfPts().end(); ++cgSpConfPtrIt) {
    (*cgSpConfPtrIt)->prepareProbConfig();
    _problemList.push_back((*cgSpConfPtrIt)->probPtr()); /// added by Ruslan
  }

  /// prepare cut and branching separation procedures
  /// should be called after preparation of col. gen. subproblems, as some cut separation procedures
  /// use network flow information of these subproblems
  std::map<std::string, GenericCutConstr *>::iterator cutMapIt;
  for (cutMapIt = _name2GenericCutConstrPtrMap.begin(); cutMapIt != _name2GenericCutConstrPtrMap.end(); ++cutMapIt)
  {
      if (!cutMapIt->second->prepareSeparation())
          progStatus().setStat(ProgStatus::terminate);
  }

  std::map<std::string, GenericBranchingConstr *>::iterator brMapIt;
  for (brMapIt = _name2GenericBranchingConstrPtrMap.begin(); brMapIt != _name2GenericBranchingConstrPtrMap.end();
       ++brMapIt)
  {
      if (!brMapIt->second->prepareSeparation())
          progStatus().setStat(ProgStatus::terminate);
  }

  /// Add local artificial columns to master problem
  InstMasterConstr *instMasterConstrPtr(NULL);

  if ((param().colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
      || (param().colGenDualPriceSmoothingAlphaFactor() > 0))
  {
    for (std::list<Constraint *>::const_iterator it = _pcConstrPtrList.begin(); it != _pcConstrPtrList.end(); ++it) {
      if (!(*it)->isTypeOf(VcId::InstMasterConstrMask))
        continue;

      instMasterConstrPtr = static_cast<InstMasterConstr *>(*it);
      instMasterConstrPtr->createStabInfo(modelPtr()->objectiveSense());
    }
  }

  int countDualizedConstr(0);
  if ((param().mastInitMode().status() == MasterInitMode::localArtCol)
      || (param().mastInitMode().status() == MasterInitMode::defaultInit)
      || (param().mastInitMode().status() == MasterInitMode::incSolColAndLac))
  {
    for (std::list<Constraint *>::const_iterator it = _pcConstrPtrList.begin(); it != _pcConstrPtrList.end(); ++it)
    {
      if (!(*it)->isTypeOf(VcId::InstMasterConstrMask))
        continue;

      instMasterConstrPtr = static_cast<InstMasterConstr *> (*it);

      if (printL(6))
        std::cout << " add localArtVar in instMasterConstr name  " << instMasterConstrPtr->name()
                  << "  subProbVarMember2coefMap().size()  "
                  << instMasterConstrPtr->subProbVarMember2coefMap().size() << std::endl;

      if (instMasterConstrPtr->addLocalArtVar(modelPtr()->objectiveSense()))
        countDualizedConstr++;
    }
  }
  if ((param().mastInitMode().status() == MasterInitMode::localAndGlobAc)) {
    for (std::list<Constraint *>::const_iterator it = _pcConstrPtrList.begin(); it != _pcConstrPtrList.end(); ++it)
    {
      if ((*it)->sense() != 'E')
        continue;

      if (!(*it)->isTypeOf(VcId::InstMasterConstrMask))
        continue;

      instMasterConstrPtr = static_cast<InstMasterConstr *> (*it);

      if (printL(6))
        std::cout << " add localArtVar in instMasterConstr name  " << instMasterConstrPtr->name()
                  << "  subProbVarMember2coefMap().size()  "
                  << instMasterConstrPtr->subProbVarMember2coefMap().size() << std::endl;

      if (instMasterConstrPtr->addLocalArtVar(modelPtr()->objectiveSense()))
        countDualizedConstr++;
    }
  }

  /// Record number of dualised master constraints
  bapcodInit().statistics().setCounter("bcMaxDualSpaceSize", countDualizedConstr);

  /// Add missing columns to master problem
  MissingColumn *misColPtr(NULL);
  for (cgSpConfPtrIt = colGenSubProbConfPts().begin(); cgSpConfPtrIt != colGenSubProbConfPts().end(); ++cgSpConfPtrIt) {
    if (param().mastInitMode().status() == MasterInitMode::subProbArtCol)
    {
      misColPtr = new MissingColumn(this, *cgSpConfPtrIt, _modelPtr->artVarCost());
      misColPtr->setAggregateVariable(misColPtr); /// moved here from the constructor of MastColumn by Ruslan

      if (printL(6))
        std::cout << "add MissingColumn to master  " << misColPtr->name() << std::endl;

      addArtVar(misColPtr);
    }
  }

  probPtr()->addVarSet(_pcVarPtrList, 1, 0);
  probPtr()->addVarSet(_nonStabStaticArtVarPtrList, 1, 0);

  /// Add master Constr
  probPtr()->addConstrSet(_pcConstrPtrList, 1, 0);
  probPtr()->addConstrSet(_pcDelayedConstrPtrList, 3, 0);

  probPtr()->buildProblem();

  if (printL(5))
  {
    probPtr()->print(std::cout);
  }
  else if (printL(2))
  {
    nicePrintAllConstraints(std::cout);
  }

  if (printL(-1) && progStatus().doRun())
  {
    std::cout << "Model is built ";
    printTime(bapcodInit().startTime().getElapsedTime(), std::cout);
  }

  bapcodInit().statistics().incrTimer("bcTimeMastPrepareProbConfig",
                                      bcTimeMastPrepareProbConfig.getElapsedTime_dbl());

  readDebugSolutionFromFile(param().debugsolution_file());

  return;
}

void MasterConf::addVariablesToForm()
{
    if (!_iVarPts.empty())
    {
        for (std::list<InstanciatedVar *>::const_iterator it = _iVarPts.begin(); it != _iVarPts.end(); ++it)
        {
            castAndAddVariable(*it, false);
        }
        _probPtr->addVarsSimplyInForm(_iVarPts);
        _iVarPts.clear();
    }
}

template <typename Out>
void split(const std::string &s, char delim, Out result)
{
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim))
    {
        *result++ = item;
    }
}

void MasterConf::readDebugSolutionFromFile(const std::string & fileName)
{
    if (fileName == "")
        return;

    if (_debugSolutionPtr != NULL)
    {
        _debugSolutionPtr->deleteSolutionsChain();
        delete _debugSolutionPtr;
    }
    _debugSolutionPtr = new Solution(this);

    std::map<std::string, double> pureMastVarNamesMap;
    std::map<ColGenSpConf *, std::vector<std::vector<int> > > spConfOrderedSolutions;
    for (auto * cgSpConfPtr : _colGenSubProbConfPts)
        spConfOrderedSolutions[cgSpConfPtr] = std::vector<std::vector<int> >();

    bool arcIdBased = false;
    bool vertexIdBased = false;

    std::map<int, std::map<std::string, int> > arcNameToIdMaps;

    std::vector<std::string> stringElems;
    std::ifstream file(fileName);
    std::string stringFromFile;
    while (getline(file, stringFromFile))
    {
        stringElems.clear();
        split(stringFromFile, ' ', std::back_inserter(stringElems));

        if ((stringElems[0] == "M") && (stringElems.size() >= 3))
        {
            try {
                double value = std::stod(stringElems[2]);
                pureMastVarNamesMap[stringElems[1]] = value;
            }
            catch (const std::invalid_argument &ia) {
                if (printL(0))
                    std::cout << "BaPCod WARNING : cannot read debug value for pure master variable : "
                              << stringElems[1] << std::endl;
            }
        }

        if ((stringElems[0] == "AN") && (stringElems.size() >= 3))
        {
            int subProbId = std::stoi(stringElems[1]);

            ColGenSpConf * subProbPtr = nullptr;
            for (ColGenSpConf * cgSpConfPtr : _colGenSubProbConfPts)
            {
                if (cgSpConfPtr->id().first() == subProbId)
                    subProbPtr = cgSpConfPtr;
            }

            if ((subProbPtr == nullptr) || (subProbPtr->networkFlowPtr() == nullptr))
            {
                if (printL(0))
                    std::cout << " BaPCod WARNING : cannot find RCSP subproblem with id " << subProbId
                              << " while reading a suproblem debug solution" << std::endl;
                continue;
            }

            if (arcNameToIdMaps.find(subProbId) == arcNameToIdMaps.end())
                arcNameToIdMaps[subProbId] = std::map<std::string, int>();
            auto & arcNameToIdMap = arcNameToIdMaps[subProbId];

            bool success = true;

            if (arcNameToIdMap.empty())
            {
                NetworkFlow *networkFlowPtr = subProbPtr->networkFlowPtr();
                for (lemon::ListDigraph::ArcIt lemonArc(networkFlowPtr->digraph());
                     lemonArc != lemon::INVALID; ++lemonArc)
                {
                    NetworkArc *netArcPtr = networkFlowPtr->netArcPtr(lemonArc);
                    auto pair = arcNameToIdMap.insert(std::make_pair(netArcPtr->name(), netArcPtr->id()));
                    if (!pair.second)
                    {
                        if (printL(0))
                            std::cout << " BaPCod WARNING : two arcs with the same name " << netArcPtr->name()
                                      << " in RCSP subproblem with id " << subProbId
                                      << " while reading a suproblem debug solution" << std::endl;
                        success = false;
                        break;
                    }
                }
                if (!success)
                    continue;
            }

            std::vector<int> arcIds;
            for (int elemPos = 2; elemPos < (int)stringElems.size(); ++elemPos)
            {
                auto findIt = arcNameToIdMap.find(stringElems[elemPos]);
                if (findIt == arcNameToIdMap.end())
                {
                    if (printL(0))
                        std::cout << " BaPCod WARNING : cannot find arc with name " << stringElems[elemPos]
                                  << " in RCSP subproblem with id " << subProbId
                                  << " while reading a suproblem debug solution" << std::endl;
                    success = false;
                    break;
                }
                arcIds.push_back(findIt->second);
            }
            if (!success)
                continue;

            arcIdBased = true;

            spConfOrderedSolutions[subProbPtr].push_back(arcIds);
        }

        if ((stringElems[0] == "V") && (stringElems.size() >= 3))
        {
            int subProbId = -1;
            std::vector<int> spOrderedSolution;
            bool success = true;
            for (int elemPos = 1; elemPos < (int)stringElems.size(); ++elemPos)
            {
                try {
                    int value = std::stoi(stringElems[elemPos]);
                    if (elemPos == 1)
                        subProbId = value;
                    else
                        spOrderedSolution.push_back(value);
                }
                catch (const std::invalid_argument &ia) {
                    if (printL(0))
                        std::cout << "BaPCod WARNING : cannot convert "<< stringElems[elemPos]
                                  << " to an integer while reading a suproblem debug solution " << std::endl;
                    success = false;
                    break;
                }
            }

            if (!success)
                continue;

            vertexIdBased = true;

            bool foundCgSp = false;
            for (auto * cgSpConfPtr : _colGenSubProbConfPts)
            {
                if (cgSpConfPtr->id().first() == subProbId)
                {
                    spConfOrderedSolutions[cgSpConfPtr].push_back(spOrderedSolution);
                    foundCgSp = true;
                }
            }
            if (!foundCgSp && printL(0))
                std::cout << " BaPCod WARNING : cannot find subproblem with id " << subProbId
                          << " while reading a suproblem debug solution" << std::endl;
        }
    }

    if (arcIdBased && vertexIdBased)
    {
        if (printL(0))
            std::cout << " BaPCod WARNING : cannot mix arcs and vertices in the debug solution" << std::endl;
        _debugSolutionPtr->deleteSolutionsChain();
        delete _debugSolutionPtr;
        _debugSolutionPtr = nullptr;
        return;
    }

    for (auto * cgSpConfPtr : _colGenSubProbConfPts)
    {
        if (cgSpConfPtr->probPtr()->setDebugSolution(spConfOrderedSolutions[cgSpConfPtr], vertexIdBased))
        {
            Solution * solutionPtr = new Solution(cgSpConfPtr);
            cgSpConfPtr->probPtr()->getDebugSolution(solutionPtr);
            Solution * solPtr = solutionPtr;
            while (solPtr != NULL)
            {
                if (!solPtr->solVarValMap().empty())
                    cgSpConfPtr->recordSubproblemSolution(solPtr, false, 3, _debugSolutionPtr);
                solPtr = solPtr->nextSolPtr();
            }
            cgSpConfPtr->clearColPtrList4Insertion(); /// otherwise columns will be added to the master on the first
                                                      /// col.gen. iteration (and we do not want that in order not to
                                                      /// modify the master in comparison with the run without
                                                      /// giving the debug solution)
            solutionPtr->deleteSolutionsChain();
            delete solutionPtr;
        }
    }

    for (VarIndexManager::iterator varPtrIt = _probPtr->probVarSet().begin(VcIndexStatus::Active, 's');
         varPtrIt != _probPtr->probVarSet().end(VcIndexStatus::Active, 's'); ++varPtrIt)
    {
        auto mapIt = pureMastVarNamesMap.find((*varPtrIt)->name());
        if (mapIt != pureMastVarNamesMap.end())
            _debugSolutionPtr->includeVar(*varPtrIt, mapIt->second, true);
    }
    _debugSolutionPtr->resetCost();
    if (printL(-1))
        std::cout << "Debug solution with cost " << _debugSolutionPtr->cost() << " is defined." << std::endl;
}
