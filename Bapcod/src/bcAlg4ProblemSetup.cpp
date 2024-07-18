/**
 *
 * This file bcAlg4ProblemSetup.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4ProblemSetup.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcStabilizationInfo.hpp"
#include "bcGlobalException.hpp"
#include "bcApplicationException.hpp"
#include "bcModelFormulationC.hpp"

#define ProblemSetupPrintLevel 3
#define ProblemSetupPrintLevelPlus 6

using namespace std;

SubProblemInfo::SubProblemInfo(ColGenSpConf * spConfPtrV) :
    spConfPtr(spConfPtrV), lb(InstMastConvexityConstr::lowerBoundWhenInactive),
        ub(InstMastConvexityConstr::upperBoundWhenInactive)
{
  if (spConfPtr->lowerBoundMastConstrPtr() != NULL)
    lb = spConfPtr->lowerBoundMastConstrPtr()->curRhs();
  if (spConfPtr->upperBoundMastConstrPtr() != NULL)
    ub = spConfPtr->upperBoundMastConstrPtr()->curRhs();
}

ProblemSetupInfo::ProblemSetupInfo(const int treatOrder) :
  treatOrderId(treatOrder), numberOfNodes(0), fullSetupIsObligatory(false), suitableMasterColumnsInfo(),
  suitableMasterCutsInfo(), activeBranchingConstraintsInfo(), subProblemsInfo(),
  masterPartialSolutionInfo(), modifiedStaticVarsInfo(), modifiedStaticConstrsInfo(),
  solverOracleInfo()
{
}

ProblemSetupInfo::ProblemSetupInfo(const int treatOrder,
                                   const std::vector<Variable *> & initialSetOfActiveColumns,
                                   const std::vector<Variable *> & initialSetOfInactiveColumns) :
  treatOrderId(treatOrder), numberOfNodes(0), fullSetupIsObligatory(false), suitableMasterColumnsInfo(),
  suitableMasterCutsInfo(), activeBranchingConstraintsInfo(), subProblemsInfo(),
  masterPartialSolutionInfo(), modifiedStaticVarsInfo(), modifiedStaticConstrsInfo(),
  solverOracleInfo()
{
  std::vector<Variable *>::const_iterator varPtrIt;
  for (varPtrIt = initialSetOfActiveColumns.begin();
      varPtrIt != initialSetOfActiveColumns.end(); ++varPtrIt)
  {
    (*varPtrIt)->incrParticipation(2);
        if(printL(7))
            cout <<"participation of "<<hex<<(long)(*varPtrIt)<<dec<<" incremeneted to : "
                    << (*varPtrIt)->participation() << endl;
    suitableMasterColumnsInfo.push_back(VariableSmallInfo(*varPtrIt, Active));
  }

  for (varPtrIt = initialSetOfInactiveColumns.begin();
      varPtrIt != initialSetOfInactiveColumns.end(); ++varPtrIt)
  {
    (*varPtrIt)->incrParticipation(2);
        if(printL(7))
            cout <<"participation of "<<hex<<(long)(*varPtrIt)<<dec<<" incremeneted to : "
                    << (*varPtrIt)->participation() << endl;

    suitableMasterColumnsInfo.push_back(VariableSmallInfo(*varPtrIt, Inactive));
  }

}

void ProblemSetupInfo::setFullSetupIsObligatory(const bool & value)
{
  fullSetupIsObligatory = value;
}

ProblemSetupInfo::~ProblemSetupInfo()
{
  int printLevel = 7;

  while (!modifiedStaticVarsInfo.empty())
    {
      delete modifiedStaticVarsInfo.front();
      modifiedStaticVarsInfo.pop_front();
    }

  for (std::list<ConstraintInfo>::iterator constrInfoIt= activeBranchingConstraintsInfo.begin() ;
          constrInfoIt != activeBranchingConstraintsInfo.end() ; constrInfoIt++)
  {
      (*constrInfoIt).constrPtr->decrParticipation(1);
      if (printL(printLevel))
       std::cout << "ProblemSetupInfo::~ProblemSetupInfo() participation of brConstr "
                 << (*constrInfoIt).constrPtr->name()
                 <<  "at " << static_cast<void*>((*constrInfoIt).constrPtr)
                 << " was decremented to " << (*constrInfoIt).constrPtr->participation() << std::endl;
      if ((*constrInfoIt).constrPtr->participation() < 0)
      {
          if(printL(printLevel))
            {
              std::cout << "A dynamic var/constr (and a branching constr in particular) "
                        << "can only participate in a positive number of nodes" << std::endl;
              std::cout << (*constrInfoIt).constrPtr->name() << " has a participation = "
                        << (*constrInfoIt).constrPtr->participation() << std::endl;
            }
          exit(1);
      }
  }

  for(std::list<ConstraintInfo>::iterator constrInfoIt = suitableMasterCutsInfo.begin() ;
          constrInfoIt != suitableMasterCutsInfo.end() ; constrInfoIt++)
  {
      (*constrInfoIt).constrPtr->decrParticipation(1);
      if (printL(printLevel))
        std::cout << "ProblemSetupInfo::~ProblemSetupInfo() participation of cut "
                  << (*constrInfoIt).constrPtr->name()
                  <<  "at " << static_cast<void*>((*constrInfoIt).constrPtr)
                  << " was decremented to " << (*constrInfoIt).constrPtr->participation() << std::endl;
      if ((*constrInfoIt).constrPtr->participation() < 0)
        {
          if (printL(printLevel))
            {
              std::cout << "A dynamic var/constr (and a cut in particular) "
                        << "can only participate in a positive number of nodes" << std::endl;
              std::cout << (*constrInfoIt).constrPtr->name() << " has a participation = "
                        << (*constrInfoIt).constrPtr->participation() << std::endl;
            }
          exit(1);
        }
  }

  for(std::list<VariableSmallInfo>::iterator varInfoIt = suitableMasterColumnsInfo.begin() ;
          varInfoIt != suitableMasterColumnsInfo.end() ; varInfoIt++)
  {
      (*varInfoIt).varPtr->decrParticipation(3);
      if (printL(printLevel))
            cout <<"participation of "<<hex<<(long)((*varInfoIt).varPtr)<<dec<<" decremeneted to : "
                    << ((*varInfoIt).varPtr)->participation() << endl;
      if ((*varInfoIt).varPtr->participation() < 0)
      {
          if (printL(0))
            std::cout << "BaPCod WARNING : "<< (*varInfoIt).varPtr->name() << " participation is negative" << std::endl;
          if(printL(printLevel))
          {
              std::cout << "A dynamic var/constr (and a branching constr in particular) "
                        << "can only participate in a positive number of nodes" << std::endl;
              std::cout << (*varInfoIt).varPtr->name() << " has a participation = "
                        << (*varInfoIt).varPtr->participation() << std::endl;
          }
          exit(1);
      }
  }
  for (std::map<const Problem *, BcSolverOracleInfo *>::iterator mapIt = solverOracleInfo.begin();
       mapIt != solverOracleInfo.end(); ++mapIt)
    delete mapIt->second;
}

Alg4ProblemSetDownOfNode::Alg4ProblemSetDownOfNode(
    MasterCommons4ProblemSetup & masterCommons) :
    _masterCommons(masterCommons), _masterProbPtr(_masterCommons.problemList().front())
{
}

Alg4ProblemSetDownOfNode::~Alg4ProblemSetDownOfNode()
{
}

void Alg4ProblemSetDownOfNode::run()
{
  if (printL(ProblemSetupPrintLevel))
    std::cout << "ProblemSetDownAlgorithm::run()" << std::endl;

  /// removal of unused dynamic vars and constrs from memory
  /// is moved to ProblemFullSetDownAlgorithm::run()
  /// (this is to make lighter the phase 1 of strong branching)

  for (std::list<Problem *>::const_iterator probPtrIt =
          _masterCommons.problemList().begin();
          probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    (*probPtrIt)->curNodePtr(NULL);
}

ProblemSetupInfo * Alg4ProblemSetDownOfNode::recordProblemInfo(int globalTreatOrder)
{
  Problem * masterProbPtr = _masterCommons.problemList().front();
  return new ProblemSetupInfo(masterProbPtr->curNodePtr()->treatOrder());
}

ProblemFullSetDownAlgorithm::ProblemFullSetDownAlgorithm(
    MasterCommons4ProblemSetup & masterCommons) :
    Alg4ProblemSetDownOfNode(masterCommons)
{
}

ProblemFullSetDownAlgorithm::~ProblemFullSetDownAlgorithm()
{
}

void ProblemFullSetDownAlgorithm::run()
{
  /// (Issam) I added this since we had a segFault when puting val(0)
  /// On deleted branching constr that were still in _inDualSol
  for (std::list<Problem *>::const_iterator probPtrIt =
          _masterCommons.problemList().begin();
          probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    (*probPtrIt)->resetSolution('d');

  for (std::list<Problem *>::const_iterator probPtrIt = _masterCommons.problemList().begin();
       probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    (*probPtrIt)->removeUnusedDynamicConstrsFromMemory();

  /// parameter true means that we also delete unused artificial variables
  /// (which are local and non related to stabilization)
  /// so we should delete dynamic variables after dynamic constraints
  _masterProbPtr->removeUnusedDynamicVarsFromMemory(true);

  if (printL(7))
    {
      std::cout << "_masterProbPtr->probVarSet().size(Unsuitable, 'd') = "
                << _masterProbPtr->probVarSet().size(Unsuitable, 'd') << std::endl;
      std::cout << "_masterProbPtr->probConstSet().size(Unsuitable, 'd') = "
                << _masterProbPtr->probConstrSet().size(Unsuitable, 'd') << std::endl;
    }

  Alg4ProblemSetDownOfNode::run();
}

ProblemSetupInfo * ProblemFullSetDownAlgorithm::recordProblemInfo(int globalTreatOrder)
{
  Problem * masterProbPtr = _masterCommons.problemList().front();
  ProblemSetupInfo * probInfoPtr = new ProblemSetupInfo(masterProbPtr->curNodePtr()->treatOrder());

  /// partial solution of the master
  for (VarPtr2DoubleMap::const_iterator mapIt = masterProbPtr->partialSolution().begin();
       mapIt != masterProbPtr->partialSolution().end(); ++mapIt)
    probInfoPtr->masterPartialSolutionInfo.push_back(VariableSolInfo(mapIt->first, mapIt->second));

  /// static variables of the master
  VarIndexManager & probVarPts = masterProbPtr->probVarSet();
  VarIndexManager::iterator varPtrIt;
  for (varPtrIt = probVarPts.begin(Active, 's'); varPtrIt != probVarPts.end(Active, 's'); ++varPtrIt)
    if (((*varPtrIt)->globalCurLb() != (*varPtrIt)->globalLb())
        || ((*varPtrIt)->globalCurUb() != (*varPtrIt)->globalUb())
        || ((*varPtrIt)->curCost() != (*varPtrIt)->costrhs()))
      probInfoPtr->modifiedStaticVarsInfo.push_back(new VariableInfo(*varPtrIt, Active));
  for (varPtrIt = probVarPts.begin(Inactive, 's'); varPtrIt != probVarPts.end(Inactive, 's'); ++varPtrIt)
    probInfoPtr->modifiedStaticVarsInfo.push_back(new VariableInfo(*varPtrIt, Inactive));

  for (varPtrIt = probVarPts.begin(Active, 'd'); varPtrIt != probVarPts.end(Active, 'd'); ++varPtrIt)
    if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
      {
        (*varPtrIt)->incrParticipation(3);
        if (printL(7))
            cout << "participation of " << hex << (long)(*varPtrIt) << dec << " incremeneted to : "
                 << (*varPtrIt)->participation() << endl;

        probInfoPtr->suitableMasterColumnsInfo.push_back(VariableSmallInfo(*varPtrIt, Active));
        if (printL(ProblemSetupPrintLevelPlus))
          std::cout << "Column " << (*varPtrIt)->name() << " is put to suitableMasterColumnsInfo as active "
                    << std::endl;
      }

  for (varPtrIt = probVarPts.begin(Inactive, 'd'); varPtrIt != probVarPts.end(Inactive, 'd'); ++varPtrIt)
    if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
      {
        (*varPtrIt)->incrParticipation(3);
        if(printL(7))
            cout << "participation of " << hex << (long)(*varPtrIt) << dec << " incremeneted to : "
                    << (*varPtrIt)->participation() << endl;

        probInfoPtr->suitableMasterColumnsInfo.push_back(VariableSmallInfo(*varPtrIt, Inactive));
        if (printL(ProblemSetupPrintLevelPlus))
          std::cout << "Column " << (*varPtrIt)->name() << " is put to suitableMasterColumnsInfo as inactive "
              << std::endl;
      }
  if (printL(1))
    std::cout << "Stored " << probVarPts.size(Active, 'd') << " active and "
              << probVarPts.size(Inactive, 'd') << " inactive" << std::endl;

  /// static constraints of the master
  ConstrIndexManager & probConstrPts = masterProbPtr->probConstrSet();
  ConstrIndexManager::iterator constrPtrIt;
  for (constrPtrIt = probConstrPts.begin(Active, 's');
       constrPtrIt != probConstrPts.end(Active, 's'); ++constrPtrIt)
    if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask)
        && (((*constrPtrIt)->curMinSlack() != (*constrPtrIt)->minSlack())
            || ((*constrPtrIt)->curMinSlack() != (*constrPtrIt)->maxSlack())))
      probInfoPtr->modifiedStaticConstrsInfo.push_back(ConstraintInfo(*constrPtrIt, Active));

  for (constrPtrIt = probConstrPts.begin(Inactive, 's');
       constrPtrIt != probConstrPts.end(Inactive, 's'); ++constrPtrIt)
    if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask))
      probInfoPtr->modifiedStaticConstrsInfo.push_back(ConstraintInfo(*constrPtrIt, Inactive));

  /// dynamic constraints of the master (cuts and branching constraints)
  for (constrPtrIt = probConstrPts.begin(Active, 'd'); constrPtrIt != probConstrPts.end(Active, 'd'); ++constrPtrIt)
    {
      if ((*constrPtrIt)->isTypeOf(VcId::InstMasterBranchingConstrMask))
        {
          (*constrPtrIt)->incrParticipation(1);
          if(printL(7))
            std::cout << "ProblemFullSetDownAlgorithm::recordProblemInfo() participation of brConstr "
                      << (*constrPtrIt)->name() <<  " at " << static_cast<void*>(*constrPtrIt)
                      << " was incremented to " << (*constrPtrIt)->participation() << std::endl;

          probInfoPtr->activeBranchingConstraintsInfo.push_back(ConstraintInfo(*constrPtrIt, Active));
          if (printL(ProblemSetupPrintLevelPlus))
            std::cout << "Master branching constraint " << (*constrPtrIt)->name()
                      << " is put to activeBranchingConstraintsInfo " << std::endl;
        }
      else if ((*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask))
        {
          (*constrPtrIt)->incrParticipation(1);
          if(printL(7))
            std::cout<<"ProblemFullSetDownAlgorithm::recordProblemInfo() participation of cut "
             << (*constrPtrIt)->name() <<  " at " << static_cast<void*>(*constrPtrIt)
             << " was incremented to " << (*constrPtrIt)->participation() << std::endl;

          probInfoPtr->suitableMasterCutsInfo.push_back(ConstraintInfo(*constrPtrIt, Active));
          if (printL(ProblemSetupPrintLevelPlus))
            std::cout << "Master cut " << (*constrPtrIt)->name() << " is put to suitableMasterCutsInfo " << std::endl;
        }
    }

  for (constrPtrIt = probConstrPts.begin(Inactive, 'd'); constrPtrIt != probConstrPts.end(Inactive, 'd'); ++constrPtrIt)
    probInfoPtr->suitableMasterCutsInfo.push_back(ConstraintInfo(*constrPtrIt, Inactive));

  for (std::vector<ColGenSpConf *>::const_iterator spConfIt = _masterCommons.colGenSubProbConfPts().begin();
       spConfIt != _masterCommons.colGenSubProbConfPts().end(); ++spConfIt)
    probInfoPtr->subProblemsInfo.push_back(SubProblemInfo(*spConfIt));

  for (std::list<Problem *>::const_iterator probPtrIt = ++_masterCommons.problemList().begin();
      probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    {
      /// variables of the subproblem
      VarIndexManager & probVarPts = (*probPtrIt)->probVarSet();
      for (varPtrIt = probVarPts.begin(Active, 's'); varPtrIt != probVarPts.end(Active, 's'); ++varPtrIt)
        {
          SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(*varPtrIt);

          if ((spVarPtr->globalCurLb() != spVarPtr->globalLb()) || (spVarPtr->globalCurUb() != spVarPtr->globalUb())
              || (spVarPtr->localCurLb() != spVarPtr->localLb()) || (spVarPtr->localCurUb() != spVarPtr->localUb())
              || (spVarPtr->curCost() != (*varPtrIt)->costrhs()))
            probInfoPtr->modifiedStaticVarsInfo.push_back(new SpVariableInfo(spVarPtr, Active));
        }
      for (varPtrIt = probVarPts.begin(Inactive, 's'); varPtrIt != probVarPts.end(Inactive, 's'); ++varPtrIt)
        {
          SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(*varPtrIt);
          probInfoPtr->modifiedStaticVarsInfo.push_back(new SpVariableInfo(spVarPtr, Inactive));
        }

      /// constraints of the subroblem
      ConstrIndexManager & probConstrPts = (*probPtrIt)->probConstrSet();
      for (constrPtrIt = probConstrPts.begin(Active, 's'); constrPtrIt != probConstrPts.end(Active, 's'); ++constrPtrIt)
        if (((*constrPtrIt)->curMinSlack() != (*constrPtrIt)->minSlack())
            || ((*constrPtrIt)->curMaxSlack() != (*constrPtrIt)->maxSlack())
            || ((*constrPtrIt)->curRhs() != (*constrPtrIt)->costrhs()))
          probInfoPtr->modifiedStaticConstrsInfo.push_back(ConstraintInfo(*constrPtrIt, Active));
      for (constrPtrIt = probConstrPts.begin(Inactive, 's');
           constrPtrIt != probConstrPts.end(Inactive, 's'); ++constrPtrIt)
        probInfoPtr->modifiedStaticConstrsInfo.push_back(ConstraintInfo(*constrPtrIt, Inactive));

      /// dynamic (branching) constraints of the subproblem
      for (constrPtrIt = probConstrPts.begin(Active, 'd'); constrPtrIt != probConstrPts.end(Active, 'd'); ++constrPtrIt)
        if ((*constrPtrIt)->isTypeOf(VcId::InstSubProbBranchingConstrMask))
          {
            (*constrPtrIt)->incrParticipation(1);
            if (printL(7))
              std::cout << "ProblemFullSetDownAlgorithm::recordProblemInfo() participation of brConstr "
                        << (*constrPtrIt)->name() <<  " at " << static_cast<void*>(*constrPtrIt)
                        << " was incremented to " << (*constrPtrIt)->participation() << std::endl;
            probInfoPtr->activeBranchingConstraintsInfo.push_back(ConstraintInfo(*constrPtrIt, Active));

            if (printL(ProblemSetupPrintLevelPlus))
              std::cout << "Subproblem branching constraint " << (*constrPtrIt)->name()
                        << " is put to activeBranchingConstraintsInfo " << std::endl;
          }
       probInfoPtr->solverOracleInfo.insert(std::make_pair(*probPtrIt, (*probPtrIt)->recordSolverOracleInfo()));
    }
  return probInfoPtr;
}

void Alg4ProblemSetupOfNode::resetPartialSolution(Problem * probPtr)
{
  if (_nodePtr->localFixedSolution() != NULL)
    {
      for (VarPtr2DoubleMap::const_iterator mapIt = _nodePtr->localFixedSolution()->solVarValMap().begin();
          mapIt != _nodePtr->localFixedSolution()->solVarValMap().end(); ++mapIt)
        probPtr->updatePartialSolution(mapIt->first, mapIt->second);
      _nodePtr->removeDebugSolutionFromThisNode(); /// we do not have yet a check to verify whether the debug
                                                   /// solution satisfies the local fixed solution
    }
}

Alg4ProblemSetupOfNode::Alg4ProblemSetupOfNode(MasterCommons4ProblemSetup & masterCommons):
  _nodePtr(NULL), _masterCommons(masterCommons),
  _masterProbPtr(_masterCommons.problemList().front())
{
}

Alg4ProblemSetupOfNode::~Alg4ProblemSetupOfNode()
{
}

bool Alg4ProblemSetupOfNode::run(Node * nodePtr)
{
  if (printL(ProblemSetupPrintLevel))
    std::cout << "ProblemSetupAlgorithm::run()" << std::endl;

  _nodePtr = nodePtr;
  for (std::list<Problem *>::const_iterator probPtrIt =
      _masterCommons.problemList().begin();
      probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    {
      bapcodInit().require((*probPtrIt)->curNodePtr() == NULL,
          "Error : current node of a problem should be null on setup");
      (*probPtrIt)->curNodePtr(nodePtr);
    }

  resetPartialSolution(_masterProbPtr);

  return false;
}

void Alg4ProblemSetupBase::deactivateVariable(Variable * varPtr,
    const VcStatus & status, bool removeFromForm)
{
  Problem * probPtr = varPtr->problemPtr();

  if (probPtr != NULL)
    {
      probPtr->probVarSet().insert(varPtr, status);
      if (removeFromForm)
        {
          varPtr->desactivate();
          _varsToRemoveFromForm.push_back(varPtr);
        }
      if (printL(ProblemSetupPrintLevel))
        {
          std::cout << "Variable " << varPtr->name()
              << " is deactivated and move to";
          if (status == Inactive)
            std::cout << " inactive var class" << std::endl;
          else
            std::cout << " unsuitable var class" << std::endl;
        }
    }
}

void Alg4ProblemSetupBase::activateVariable(Variable * varPtr)
{
    if (varPtr->activateVariable(false))
    {
        if (printL(ProblemSetupPrintLevel))
            std::cout << "Variable " << varPtr->name() << " is activated"
            << std::endl;
        _varsToAddToForm.push_back(varPtr);
    }
}

void Alg4ProblemSetupBase::deactivateConstraint(Constraint * constrPtr,
						const VcStatus & status,
						bool removeFromForm)
{
  Problem * probPtr = constrPtr->problemPtr();

  if (probPtr != NULL)
    {
      probPtr->probConstrSet().insert(constrPtr, status);
      if (removeFromForm)
        {
          constrPtr->desactivate();
          _constrsToRemoveFromForm.push_back(constrPtr);
          if(printL(7))
            cout<<"adding to _constrsToRemoveFromForm "<<constrPtr->name()
                <<" at "<<hex<<(long)constrPtr<<dec<<endl;
        }
      if (printL(ProblemSetupPrintLevel))
        {
          std::cout << "Constraint " << constrPtr->name()
              << " is deactivated and move to";
          if (status == Inactive)
            std::cout << " inactive constr class" << std::endl;
          else
            std::cout << " unsuitable constr class" << std::endl;
        }

      Variable * varPtr = constrPtr->posLocalArtVarPtr();

      if (varPtr != NULL)
        deactivateVariable(varPtr, status, removeFromForm);
      varPtr = constrPtr->negLocalArtVarPtr();

      if (varPtr != NULL)
        deactivateVariable(varPtr, status, removeFromForm);

      if (constrPtr->stabInfoPtr() != NULL)
        {
          varPtr = constrPtr->stabInfoPtr()->negInnerArtVarPtr();
          if (varPtr != NULL)
            deactivateVariable(varPtr, status, removeFromForm);
          varPtr = constrPtr->stabInfoPtr()->negOuterArtVarPtr();
          if (varPtr != NULL)
            deactivateVariable(varPtr, status, removeFromForm);
          varPtr = constrPtr->stabInfoPtr()->posInnerArtVarPtr();
          if (varPtr != NULL)
            deactivateVariable(varPtr, status, removeFromForm);
          varPtr = constrPtr->stabInfoPtr()->posOuterArtVarPtr();
          if (varPtr != NULL)
            deactivateVariable(varPtr, status, removeFromForm);
        }
    }
}

void Alg4ProblemSetupBase::activateConstraint(Constraint * constrPtr)
{
  if (constrPtr->activateConstraint(false))
    {
      if (printL(ProblemSetupPrintLevel))
        std::cout << "Constraint " << constrPtr->name() << " is activated" << std::endl;
      _constrsToAddToForm.push_back(constrPtr);

      Variable * varPtr = constrPtr->posLocalArtVarPtr();
      if (varPtr != NULL)
        _varsToAddToForm.push_back(varPtr);
      varPtr = constrPtr->negLocalArtVarPtr();
      if (varPtr != NULL)
        _varsToAddToForm.push_back(varPtr);

      if (constrPtr->stabInfoPtr() != NULL)
        {
          varPtr = constrPtr->stabInfoPtr()->negInnerArtVarPtr();
          if (varPtr != NULL)
              _varsToAddToForm.push_back(varPtr);
          varPtr = constrPtr->stabInfoPtr()->negOuterArtVarPtr();
          if (varPtr != NULL)
              _varsToAddToForm.push_back(varPtr);
          varPtr = constrPtr->stabInfoPtr()->posInnerArtVarPtr();
          if (varPtr != NULL)
              _varsToAddToForm.push_back(varPtr);
          varPtr = constrPtr->stabInfoPtr()->posOuterArtVarPtr();
          if (varPtr != NULL)
              _varsToAddToForm.push_back(varPtr);
        }
    }
}


void Alg4ProblemSetupBase::resetStaticVariables(Problem * probPtr,
    std::list<VariableInfo*>::const_iterator & varInfoIt)
{

  for (; (varInfoIt != _probSetupInfoPtr->modifiedStaticVarsInfo.end())
          && ((*varInfoIt)->varPtr->problemPtr() == probPtr); ++varInfoIt)
    {
      if ((*varInfoIt)->status == Active)
        {
          if ((*varInfoIt)->varPtr->vcIndexStatus() == Active)
            {
              if ((*varInfoIt)->needToChangeBounds())
                _varsToChangeBounds.push_back((*varInfoIt)->varPtr);
              if ((*varInfoIt)->cost != (*varInfoIt)->varPtr->curCost())
                _varsToChangeCost.push_back((*varInfoIt)->varPtr);
            }
          if ((*varInfoIt)->varPtr->vcIndexStatus() == Inactive)
            activateVariable((*varInfoIt)->varPtr);
          (*varInfoIt)->applyVarInfo();
        }
      if (((*varInfoIt)->status == Inactive)
          && ((*varInfoIt)->varPtr->vcIndexStatus() == Active))
        deactivateVariable((*varInfoIt)->varPtr, Inactive);
      (*varInfoIt)->varPtr->infoIsUpdated(true);
    }

  for (VarIndexManager::iterator varPtrIt = probPtr->probVarSet().begin(Active,'s');
       varPtrIt != probPtr->probVarSet().end(Active, 's'); ++varPtrIt)
    {
      if (!(*varPtrIt)->infoIsUpdated())
        {
          if ((*varPtrIt)->boundsAreDifferentFromDefaultInForm())
            {
              if (printL(ProblemSetupPrintLevel))
                std::cout << "Bounds of variable " << (*varPtrIt)->name() << " are reset to default" << std::endl;
              _varsToChangeBounds.push_back(*varPtrIt);
            }
          if ((*varPtrIt)->curCost() != (*varPtrIt)->costrhs())
            {
              if (printL(ProblemSetupPrintLevel))
                std::cout << "Cost of variable " << (*varPtrIt)->name() << " is reset to default " << std::endl;
              _varsToChangeCost.push_back(*varPtrIt);
            }
          (*varPtrIt)->resetBoundsAndCostToDefaults();
        }
      else
        (*varPtrIt)->infoIsUpdated(false);
    }

  for (VarIndexManager::iterator varPtrIt = probPtr->probVarSet().begin(
      Inactive, 's'); varPtrIt != probPtr->probVarSet().end(Inactive, 's');)
    {
      if (!(*varPtrIt)->infoIsUpdated())
        {
          Variable * varPtr = *varPtrIt;
          ++varPtrIt;
          varPtr->resetBoundsAndCostToDefaults();
          activateVariable(varPtr);
        }
      else
        {
          (*varPtrIt)->infoIsUpdated(false);
          ++varPtrIt;
        }
    }
}

void Alg4ProblemSetupBase::resetStaticConstraints(Problem * probPtr,
    std::list<ConstraintInfo>::const_iterator & constrInfoIt)
{

  for (; //Issam changed || to && in the following line
      (constrInfoIt != _probSetupInfoPtr->modifiedStaticConstrsInfo.end())
          && (constrInfoIt->constrPtr->problemPtr() == probPtr); ++constrInfoIt)
    {
      if (constrInfoIt->status == Active)
        {
          if (constrInfoIt->constrPtr->vcIndexStatus() == Active)
            {
              if (constrInfoIt->rhs != constrInfoIt->constrPtr->curRhs())
                _constrsToChangeRhs.push_back(constrInfoIt->constrPtr);
            }
          if (constrInfoIt->constrPtr->vcIndexStatus() == Inactive)
            activateConstraint(constrInfoIt->constrPtr);
          constrInfoIt->applyInfo();
        }

      if ((constrInfoIt->status == Inactive)
          && (constrInfoIt->constrPtr->vcIndexStatus() == Active))
        deactivateConstraint(constrInfoIt->constrPtr, Inactive);

      constrInfoIt->constrPtr->infoIsUpdated(true);
    }

  for (ConstrIndexManager::iterator constrPtrIt =
      probPtr->probConstrSet().begin(Active, 's');
      constrPtrIt != probPtr->probConstrSet().end(Active, 's'); ++constrPtrIt)
    if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask))
      {
        if (!(*constrPtrIt)->infoIsUpdated())
          {
            if ((*constrPtrIt)->curRhs() != (*constrPtrIt)->costrhs())
              {
                (*constrPtrIt)->curRhs((*constrPtrIt)->costrhs());
                _constrsToChangeRhs.push_back(*constrPtrIt);
              }
            (*constrPtrIt)->resetSlacksAndRhsToDefaults();
          }
        else
          (*constrPtrIt)->infoIsUpdated(false);
      }

  for (ConstrIndexManager::iterator constrPtrIt =
      probPtr->probConstrSet().begin(Inactive, 's');
      constrPtrIt != probPtr->probConstrSet().end(Inactive, 's');)
    if (!(*constrPtrIt)->isTypeOf(VcId::InstMastConvexityConstrMask))
      {
        if (!(*constrPtrIt)->infoIsUpdated())
          {
            Constraint * constrPtr = *constrPtrIt;
            ++constrPtrIt;
            constrPtr->resetSlacksAndRhsToDefaults();
            activateConstraint(constrPtr);
          }
        else
          {
            (*constrPtrIt)->infoIsUpdated(false);
            ++constrPtrIt;
          }
      }
    else
      ++constrPtrIt;
}

bool Alg4ProblemSetupBase::columnIsUnsuitable(MastColumn * colPtr)
{
  Solution * solPtr = colPtr->spSol();
  if (solPtr != NULL)
  {
    for (VarPtr2DoubleMap::const_iterator it = solPtr->solVarValMap().begin();
         it != solPtr->solVarValMap().end(); ++it)
    {
      SubProbVariable * spVarPtr = dynamic_cast<SubProbVariable *> (it->first);
      if ((it->second < spVarPtr->localCurLb()) || (it->second > spVarPtr->localCurUb()))
        return true;
    }
  }

  /// check now whether colPtr satisfies subproblem branching constraints
  /// note that subproblem it not setup at this point, so we use information
  /// from _newSpBranchingConstraints and _oldSpBranchingConstraints

  for (std::list<Constraint *>::iterator constrPtrIt =
      _newSpBranchingConstraints.begin();
      constrPtrIt != _newSpBranchingConstraints.end(); ++constrPtrIt)
    if (colPtr->cgSpConfPtr()->probPtr() == (*constrPtrIt)->problemPtr())
      {
        if ((*constrPtrIt)->violated(colPtr, 1.0))
          return true;
      }

  for (std::list<Constraint *>::iterator constrPtrIt =
      _oldSpBranchingConstraints.begin();
      constrPtrIt != _oldSpBranchingConstraints.end(); ++constrPtrIt)
    if (colPtr->cgSpConfPtr()->probPtr() == (*constrPtrIt)->problemPtr())
      {
        if ((*constrPtrIt)->violated(colPtr, 1.0))
          return true;
      }

  return false;
}

bool Alg4ProblemSetupBase::columnViolatesNewSpBranchingConstriants(MastColumn * colPtr)
{
  for (std::list<Constraint *>::iterator constrPtrIt =
      _newSpBranchingConstraints.begin();
      constrPtrIt != _newSpBranchingConstraints.end(); ++constrPtrIt)
    {
      if (((*constrPtrIt)->problemPtr() == colPtr->cgSpConfPtr()->probPtr())
          && ((*constrPtrIt)->violated(colPtr, 1.0)))
        {
          if (printL(ProblemSetupPrintLevelPlus))
            std::cout << "Column " << colPtr->name()
                      << " does not satisfy subproblem branching constraints and it is made unsuitable " << std::endl;
          return true;
        }
    }
  return false;
}

void Alg4ProblemSetupBase::resetMasterColumns()
{
    std::vector<Variable *> colsToActivate;
    colsToActivate.reserve(_probSetupInfoPtr->suitableMasterColumnsInfo.size());
  for (std::list<VariableSmallInfo>::const_iterator varInfoIt =
      _probSetupInfoPtr->suitableMasterColumnsInfo.begin();
      varInfoIt != _probSetupInfoPtr->suitableMasterColumnsInfo.end();
      ++varInfoIt)
    {
      if (columnViolatesNewSpBranchingConstriants(
          static_cast<MastColumn *>(varInfoIt->varPtr)))
        continue;

      if ((varInfoIt->status == Active) || _makeAllColumnsActive)
        {
          if (printL(7))
            std::cout << "accessing column at " << hex << (long)(varInfoIt->varPtr) << dec << std::endl;
          if (varInfoIt->varPtr->vcIndexStatus() == Active)
            {
              if (varInfoIt->cost != varInfoIt->varPtr->curCost())
                _varsToChangeCost.push_back(varInfoIt->varPtr);

              if (printL(ProblemSetupPrintLevelPlus))
                std::cout << "Column " << varInfoIt->varPtr->name() << " is already active " << std::endl;
            }
          if (varInfoIt->varPtr->vcIndexStatus() == Inactive || varInfoIt->varPtr->vcIndexStatus() == Unsuitable)
          {
              /// activation of columns should be done after deactivation of columns
              /// (otherwise there might be difference between the set of active columns
              ///  and the column pool in the prob.var.manager : activated column may not go to the pool
              ///  because a similar column is already there, and this similar column will be deactivated)
              colsToActivate.push_back(varInfoIt->varPtr);
          }
          varInfoIt->applyVarInfo();
        }
      else if (varInfoIt->status == Inactive)
        {
          if (varInfoIt->varPtr->vcIndexStatus() == Active)
            deactivateVariable(varInfoIt->varPtr, Inactive);
          else if (varInfoIt->varPtr->vcIndexStatus() == Unsuitable)
            deactivateVariable(varInfoIt->varPtr, Inactive, false);
        }
      varInfoIt->varPtr->infoIsUpdated(true);
    }

  /// For both active and inactive columns that are not in suitableMasterColumnsInfo,
  /// we treat them in the same way:
  /// - If they are new they can become unsuitable or inactive (unless param().UseColumnsPool() is false)
  /// - Otherwise, they become for sure unsuitable.

  VarIndexManager::iterator varPtrIt = _masterProbPtr->probVarSet().begin(Active, 'd');
  bool loopOnActiveList = true;

  while (true)
    {
      if (loopOnActiveList
          && (varPtrIt == _masterProbPtr->probVarSet().end(Active, 'd')))
        {
          varPtrIt = _masterProbPtr->probVarSet().begin(Inactive, 'd');
          loopOnActiveList = false;
        }
      if (!loopOnActiveList
          && (varPtrIt == _masterProbPtr->probVarSet().end(Inactive, 'd')))
        {
          break;
        }

      if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask)
          && !(*varPtrIt)->infoIsUpdated())
        {
          MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
          ++varPtrIt;

          /// column which is not in suitableMasterColumnsInfo can be suitable
          /// only if it is "new " (has been generated after creation of _probSetupInfoPtr)
          bool isUnsuitableColumn
            = ( ( (!bapcodInit().param().Search4NegRedCostColInInactivePool() || !bapcodInit().param().UseColumnsPool() )
                  && !_makeAllColumnsActive)
                || (colPtr->treatOrderId() <= _probSetupInfoPtr->treatOrderId)
                || columnIsUnsuitable(colPtr) );

          if (isUnsuitableColumn)
            deactivateVariable(colPtr, Unsuitable, loopOnActiveList);
          else
          {
            if (loopOnActiveList && !_makeAllColumnsActive)
              deactivateVariable(colPtr, Inactive); /// we will check this column for suitability second time
                                                    /// when we loop over inactive columns (to review it)
            if (!loopOnActiveList && _makeAllColumnsActive)
                activateVariable(colPtr);
          }
        }
      else
        {
          (*varPtrIt)->infoIsUpdated(false);
          varPtrIt++;
        }
    }

  for (auto varPtr : colsToActivate)
  {
      activateVariable(varPtr);
      varPtr->infoIsUpdated(false);
  }

  /// TO DO : review it, we will check suitability of new unsuitable columns two times
  if (bapcodInit().param().Search4NegRedCostColInInactivePool() || _makeAllColumnsActive)
    {
      for (varPtrIt = _masterProbPtr->probVarSet().begin(Unsuitable, 'd');
          varPtrIt != _masterProbPtr->probVarSet().end(Unsuitable, 'd'); )
        if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
          {
            MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
            ++varPtrIt;
            bool isSuitableColumn = (colPtr->treatOrderId() > _probSetupInfoPtr->treatOrderId)
                && !columnIsUnsuitable(colPtr);
            if (isSuitableColumn)
              {
                if (_makeAllColumnsActive)
                    activateVariable(colPtr);
                else
                    deactivateVariable(colPtr, Inactive, false);
              }
          }
        else
          varPtrIt++;
    }


}

void Alg4ProblemSetupBase::prepareBranchingConstraints()
{
  _newSpBranchingConstraints.clear();
  _oldSpBranchingConstraints.clear();
  _activeSpBranchingConstraintsExist = false;
  _nonLinearBranchingConstraints = false;
  /// first, we add each undefined constraints to its problem
  /// at the same time, we isolate component set branching constraints
  /// and subproblem branching constraints (we need them to filter active columns)
  std::vector<CompSetInstMastBranchConstr *> csBrConstrVector;
  for (std::list<Constraint*>::iterator brConstrPtrIt = _localBranchingConstraints.begin();
        brConstrPtrIt != _localBranchingConstraints.end(); ++brConstrPtrIt)
    {
      if ((*brConstrPtrIt)->vcIndexStatus() == Undefined)
        {
          if ((*brConstrPtrIt)->isTypeOf(VcId::MasterConstrMask))
            {
              _masterProbPtr->addConstraintToProblem(*brConstrPtrIt);
              InstMasterConstr * instMastConstrPtr = dynamic_cast<InstMasterConstr *>(*brConstrPtrIt);
              if (instMastConstrPtr != NULL)
                instMastConstrPtr->treatOrderId(_nodePtr->treatOrder());
            }
          else
            (*brConstrPtrIt)->probConfPtr()->probPtr()->addConstraintToProblem(*brConstrPtrIt);
        }

        if ((*brConstrPtrIt)->isTypeOf(VcId::InstSubProbBranchingConstrMask))
          _newSpBranchingConstraints.push_back(*brConstrPtrIt);

        if ((*brConstrPtrIt)->isTypeOf(VcId::Base4NonLinearConstraintMask))
          _nonLinearBranchingConstraints = true;

      if ((*brConstrPtrIt)->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
        csBrConstrVector.push_back(
            static_cast<CompSetInstMastBranchConstr *>(*brConstrPtrIt));
    }

  for (std::list<ConstraintInfo>::const_iterator constrInfoIt =
      _probSetupInfoPtr->activeBranchingConstraintsInfo.begin();
      constrInfoIt != _probSetupInfoPtr->activeBranchingConstraintsInfo.end();
      ++constrInfoIt)
    {
      /// we first verify whether constrInfoIt->constrPtr is in _localBranchingConstraints
      /// (possible in strong branching), and do not include it in this case
      bool alreadyPushed = false;
      for (std::list<Constraint*>::iterator brConstrPtrIt = _localBranchingConstraints.begin();
           brConstrPtrIt != _localBranchingConstraints.end(); ++brConstrPtrIt)
        if (constrInfoIt->constrPtr == *brConstrPtrIt)
          alreadyPushed = true;
      if (alreadyPushed)
        continue;

      if (constrInfoIt->constrPtr->isTypeOf(VcId::CompSetInstMastBranchConstrMask))
         csBrConstrVector.push_back(static_cast<CompSetInstMastBranchConstr *>(constrInfoIt->constrPtr));

      if (constrInfoIt->constrPtr->isTypeOf(VcId::InstSubProbBranchingConstrMask))
        _oldSpBranchingConstraints.push_back(constrInfoIt->constrPtr);

      if (constrInfoIt->constrPtr->isTypeOf(VcId::Base4NonLinearConstraintMask))
        _nonLinearBranchingConstraints = true;
    }

  for (std::list<Problem *>::const_iterator probPtrIt = ++_masterCommons.problemList().begin();
       probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
  {
      for (ConstrIndexManager::iterator constrPtrIt = (*probPtrIt)->probConstrSet().begin(Active, 'd');
           constrPtrIt != (*probPtrIt)->probConstrSet().end(Active, 'd'); )
      {
          if ((*constrPtrIt)->isTypeOf(VcId::BranchingConstrBaseTypeMask))
          {
              _activeSpBranchingConstraintsExist = true;
              break;
          }
      }
      if (_activeSpBranchingConstraintsExist)
          break;
  }


  /// we need to add constraints to csBrConstrList in the order such that
  /// constraints generated deeper in the search tree appear earlier in this list
  std::stable_sort(csBrConstrVector.begin(), csBrConstrVector.end(),
      CompSetInstMastBranchConstr::BCDepthWhenGeneratedComparator());
  std::list<CompSetInstMastBranchConstr *> csBrConstrList(
      csBrConstrVector.begin(), csBrConstrVector.end());

  if(!csBrConstrList.empty())
    {
      /// redundunt component set branching constraints will be put to the preprocessing list
      _problemInfeasible = _problemInfeasible || buildTreeOfColClasses(csBrConstrList);

      std::list<Problem*>::const_iterator probIt = _masterCommons.problemList().begin();
      Problem * probPtr = NULL;

      // We initialize the maps mapCompSetBrConstr... only for sub problems.
      probIt++;
      while(probIt != _masterCommons.problemList().end())
        {
          probPtr = *probIt;

          for(VarIndexManager::iterator varPtrIt = probPtr->probVarSet().begin(Active, 's');
                  varPtrIt != probPtr->probVarSet().end(Active, 's'); ++varPtrIt)
            {
              //The following if block is needed for component set branching (added Issam to be checked by Ruslan or FV)
              if((*varPtrIt)->isTypeOf(VcId::SubProbVariableMask))
                {
                  (*varPtrIt)->mapCompSetBrConstr2lowerBd().clear();
                  (*varPtrIt)->mapCompSetBrConstr2upperBd().clear();
                }
            }

          probIt++;
        }

    }

  for (std::vector<CompSetInstMastBranchConstr *>::const_iterator bcitPt =
          _nodePtr->treeOfColClasses().begin(); bcitPt != _nodePtr->treeOfColClasses().end();
          bcitPt++)
    {
      if (printL(3))
        std::cout << "MasterConf::setupBranchingConstraint():  non redundant csBrConstr include "
                  << (*bcitPt)->name() << std::endl;

      (*bcitPt)->recordInducedVarBounds();
    }
}

void Alg4ProblemSetupBase::resetBranchingConstraints(Problem * probPtr,
                                                     std::list<ConstraintInfo>::const_iterator & constrInfoIt)
{
  for (std::list<Constraint*>::iterator brConstrPtrIt = _localBranchingConstraints.begin();
       brConstrPtrIt != _localBranchingConstraints.end(); ++brConstrPtrIt)
    {
      if (((*brConstrPtrIt)->problemPtr() == probPtr) && !(*brConstrPtrIt)->inPreprocessedList())
        {
          if ((*brConstrPtrIt)->vcIndexStatus() == Active)
            {
              if ((*brConstrPtrIt)->curRhs() != (*brConstrPtrIt)->costrhs())
                {
                  (*brConstrPtrIt)->curRhs((*brConstrPtrIt)->costrhs());
                  _constrsToChangeRhs.push_back(*brConstrPtrIt);
                }
              if (printL(ProblemSetupPrintLevel))
                std::cout << "New branching constraint " << (*brConstrPtrIt)->name() << " is already active "
                          << std::endl;
            }
          else
            {
              (*brConstrPtrIt)->resetSlacksAndRhsToDefaults();
              if (printL(7))
                std::cout << "brConstr " << (*brConstrPtrIt)->name() << " " << std::hex << (long)(*brConstrPtrIt)
                          << std::dec <<  " is being activated " << std::endl;
              activateConstraint(*brConstrPtrIt);
            }
          (*brConstrPtrIt)->infoIsUpdated(true);
        }
    }

  for (; (constrInfoIt != _probSetupInfoPtr->activeBranchingConstraintsInfo.end())
         && (constrInfoIt->constrPtr->problemPtr() == probPtr); ++constrInfoIt)
    if (!constrInfoIt->constrPtr->inPreprocessedList())
      {
        if (constrInfoIt->constrPtr->vcIndexStatus() == Active)
          {
            if (constrInfoIt->rhs != constrInfoIt->constrPtr->curRhs())
              _constrsToChangeRhs.push_back(constrInfoIt->constrPtr);
            constrInfoIt->applyInfo();
            if (printL(ProblemSetupPrintLevel))
              std::cout << "Branching constraint " << constrInfoIt->constrPtr->name() << " is already active "
                        << std::endl;
          }
        else
          {
            constrInfoIt->applyInfo();
            activateConstraint(constrInfoIt->constrPtr);
          }
        constrInfoIt->constrPtr->infoIsUpdated(true);
      }

  for (ConstrIndexManager::iterator constrPtrIt = probPtr->probConstrSet().begin(Active, 'd');
       constrPtrIt != probPtr->probConstrSet().end(Active, 'd'); )
    {
      if ((*constrPtrIt)->isTypeOf(VcId::BranchingConstrBaseTypeMask) && !(*constrPtrIt)->infoIsUpdated())
        {
          Constraint * constrPtr = *constrPtrIt;
          ++constrPtrIt;
          deactivateConstraint(constrPtr, Unsuitable);
        }
      else
        {
          (*constrPtrIt)->infoIsUpdated(false);
          constrPtrIt++;
        }
    }

  probPtr->clearPreprocessingLists();
}

void Alg4ProblemSetupBase::resetConvexityConstraints()
{
  for (std::vector<ColGenSpConf *>::const_iterator spConfIt =
      _masterCommons.colGenSubProbConfPts().begin();
      spConfIt != _masterCommons.colGenSubProbConfPts().end(); ++spConfIt)
    {
      InstMastConvexityConstr * convConstrPtr =
          (*spConfIt)->lowerBoundMastConstrPtr();
      Double newRhs = convConstrPtr->newLocalRhs();
      if (newRhs != convConstrPtr->curRhs())
        {
          convConstrPtr->curRhs(newRhs);
          (*spConfIt)->lowerBoundPtr(new Double(newRhs));
          if (printL(ProblemSetupPrintLevel))
            std::cout << "Rhs of convexity constraint " << convConstrPtr->name() << " is set to " << newRhs
                      << std::endl;
          if (newRhs == InstMastConvexityConstr::lowerBoundWhenInactive)
            {
              if (convConstrPtr->vcIndexStatus() == Active)
                deactivateConstraint(convConstrPtr, Inactive);
            }
          else
            {
              if (convConstrPtr->vcIndexStatus() == Inactive)
                activateConstraint(convConstrPtr);
              else
                _constrsToChangeRhs.push_back(convConstrPtr);
            }
        }

      convConstrPtr = (*spConfIt)->upperBoundMastConstrPtr();
      newRhs = convConstrPtr->newLocalRhs();
      if (newRhs != convConstrPtr->curRhs())
        {
          convConstrPtr->curRhs(newRhs);
          (*spConfIt)->upperBoundPtr(new Double(newRhs));
          if (printL(ProblemSetupPrintLevel))
            std::cout << "Rhs of convexity constraint " << convConstrPtr->name() << " is set to " << newRhs
                      << std::endl;
          if (newRhs == InstMastConvexityConstr::upperBoundWhenInactive)
            {
              if (convConstrPtr->vcIndexStatus() == Active)
                deactivateConstraint(convConstrPtr, Inactive);
            }
          else
            {
              if (convConstrPtr->vcIndexStatus() == Inactive)
                activateConstraint(convConstrPtr);
              else
                _constrsToChangeRhs.push_back(convConstrPtr);
            }
        }
    }
}

void Alg4ProblemSetupBase::resetNonStabArtificialVariables()
{
  for (VarIndexManager::iterator varPtrIt = _masterProbPtr->probVarSet().begin(Active, 'a');
       varPtrIt != _masterProbPtr->probVarSet().end(Active, 'a'); ++varPtrIt)
    {
      if ((*varPtrIt)->isTypeOf(VcId::LocalArtificialVarMask))
        {
          LocalArtificialVar * artVarPtr = static_cast<LocalArtificialVar *>(*varPtrIt);
          LocalArtificialVar::LocalArtClassId locClassId = artVarPtr->localClassId();
          if ((locClassId != LocalArtificialVar::NegLocalId) && (locClassId != LocalArtificialVar::PosLocalId))
            continue;
        }
      (*varPtrIt)->resetCostFromDefaultCost();
      _varsToChangeCost.push_back(*varPtrIt);
    }
}

void Alg4ProblemSetupBase::resetMasterCuts()
{
  for (std::list<ConstraintInfo>::const_iterator constrInfoIt = _probSetupInfoPtr->suitableMasterCutsInfo.begin();
       constrInfoIt != _probSetupInfoPtr->suitableMasterCutsInfo.end(); ++constrInfoIt)
    {
      if (constrInfoIt->status == Active)
        {
          if (constrInfoIt->constrPtr->vcIndexStatus() == Active)
            {
              if (constrInfoIt->rhs != constrInfoIt->constrPtr->curRhs())
                _constrsToChangeRhs.push_back(constrInfoIt->constrPtr);
              constrInfoIt->applyInfo();
            }
          else
            {
              constrInfoIt->applyInfo();
              if (constrInfoIt->constrPtr->vcIndexStatus() == Undefined)
                //addConstraintToProblem(constrInfoIt->constrPtr, _masterProbPtr);
                _masterProbPtr->addConstraintToProblem(constrInfoIt->constrPtr);
              activateConstraint(constrInfoIt->constrPtr);
            }
        }
      if (constrInfoIt->status == Inactive)
        {
          if (constrInfoIt->constrPtr->vcIndexStatus() == Active)
            deactivateConstraint(constrInfoIt->constrPtr, Inactive);
          if (constrInfoIt->constrPtr->vcIndexStatus() == Unsuitable)
            _masterProbPtr->probConstrSet().insert(constrInfoIt->constrPtr, Inactive);
        }
      constrInfoIt->constrPtr->infoIsUpdated(true);
    }

  /// For both active and inactive cuts that are not in suitableMasterCutsInfo,
  /// we treat them in the same way:
  /// They become for sure unsuitable.

  ConstrIndexManager::iterator constrPtrIt = _masterProbPtr->probConstrSet().begin(
      Active, 'd');
  bool loopOnActiveList = true;

  while (true)
    {
      if (loopOnActiveList
          && constrPtrIt == _masterProbPtr->probConstrSet().end(Active, 'd'))
        {
          constrPtrIt = _masterProbPtr->probConstrSet().begin(Inactive, 'd');
          loopOnActiveList = false;
        }
      if (!loopOnActiveList
          && constrPtrIt == _masterProbPtr->probConstrSet().end(Inactive, 'd'))
        {
          break;
        }

      if (!(*constrPtrIt)->isTypeOf(VcId::BranchingConstrBaseTypeMask)
          && !(*constrPtrIt)->infoIsUpdated())
        {
          Constraint * constrPtr = *constrPtrIt;
          ++constrPtrIt;
          deactivateConstraint(constrPtr, Unsuitable, loopOnActiveList);
        }
      else
        {
          (*constrPtrIt)->infoIsUpdated(false);
          constrPtrIt++;
        }
    }

}

void Alg4ProblemSetupBase::printVarsList(std::ostream & os, const std::list<Variable *> & varsList,
                                         const std::string & stringToPrint)
{
  if (!varsList.empty())
    {
      std::list<Variable *>::const_iterator varPtrIt = varsList.begin();
      os << stringToPrint << " : " << (*varPtrIt)->name();
      for (++varPtrIt; varPtrIt != varsList.end(); ++varPtrIt)
        os << ", " << (*varPtrIt)->name();
      os << std::endl;
    }
}

void Alg4ProblemSetupBase::printConstrsList(std::ostream & os, const std::list<Constraint *> & constrsList,
                                            const std::string & stringToPrint)
{
  if (!constrsList.empty())
    {
      std::list<Constraint *>::const_iterator constrPtrIt = constrsList.begin();
      os << stringToPrint << " : " << (*constrPtrIt)->name();
      for (++constrPtrIt; constrPtrIt != constrsList.end(); ++constrPtrIt)
        os << ", " << (*constrPtrIt)->name();
      os << std::endl;
    }
}

void Alg4ProblemSetupBase::clearVarConstrLists()
{
    _varsToAddToForm.clear();
    _varsToRemoveFromForm.clear();
    _varsToChangeBounds.clear();
    _varsToChangeCost.clear();

    _constrsToAddToForm.clear();
    _constrsToRemoveFromForm.clear();
    _constrsToChangeRhs.clear();
}

void Alg4ProblemSetupBase::updateFormulation(Problem * probPtr)
{
  if (printL(ProblemSetupPrintLevel))
    {
      printVarsList(std::cout, _varsToAddToForm, "Vars to add to the form");
      printVarsList(std::cout, _varsToRemoveFromForm, "Vars to remove from the form");
      printVarsList(std::cout, _varsToChangeBounds, "Vars to change bound(s) in the form");
      printVarsList(std::cout, _varsToChangeCost, "Vars to change the cost in the form");
      printConstrsList(std::cout, _constrsToAddToForm, "Constrs to add to the form");
      printConstrsList(std::cout, _constrsToRemoveFromForm, "Constrs to remove from the form");
      printConstrsList(std::cout, _constrsToChangeRhs, "Constrs to change the rhs in the form");
    }

  probPtr->addVarsSimplyInForm(_varsToAddToForm);
  probPtr->delVarsSimplyInForm(_varsToRemoveFromForm);
  probPtr->updateBoundsInForm(_varsToChangeBounds);
  probPtr->updateObjCoeffsInForm(_varsToChangeCost);

  probPtr->addConstrsSimplyInForm(_constrsToAddToForm);
  probPtr->delConstrsSimplyInForm(_constrsToRemoveFromForm);
  probPtr->updateConstrRhsInForm(_constrsToChangeRhs);

  _varsToAddToForm.clear();
  _varsToRemoveFromForm.clear();
  _varsToChangeBounds.clear();
  _varsToChangeCost.clear();

  _constrsToAddToForm.clear();
  _constrsToRemoveFromForm.clear();
  _constrsToChangeRhs.clear();
}

void Alg4ProblemSetupBase::applySubproblemsInfo()
{
  for (std::list<SubProblemInfo>::const_iterator spInfoIt =
      _probSetupInfoPtr->subProblemsInfo.begin();
      spInfoIt != _probSetupInfoPtr->subProblemsInfo.end(); ++spInfoIt)
    {
      /// we update _locallyValidRhs of the convexity constraints
      /// then they can be changed by component set branching constraints reset
      spInfoIt->spConfPtr->lowerBoundMastConstrPtr()->newLocalRhs(spInfoIt->lb);
      spInfoIt->spConfPtr->upperBoundMastConstrPtr()->newLocalRhs(spInfoIt->ub);
    }
}

void Alg4ProblemSetupBase::fillLocalBranchingConstraints()
{
  for (std::list<BranchingConstrBaseType *>::const_iterator brConstrIt =
      _nodePtr->localNodeBrConstrList().begin();
      brConstrIt != _nodePtr->localNodeBrConstrList().end(); ++brConstrIt)
    if ((*brConstrIt)->isTypeOf(VcId::ConstraintMask))
      {
        Constraint * constrPtr = dynamic_cast<Constraint *>(*brConstrIt);
        if (constrPtr != NULL)
          _localBranchingConstraints.push_back(constrPtr);
      }
}

bool Alg4ProblemSetupBase::debugSolSatisfiesLocalBranchConstrs()
{
    if (_masterCommons.debugSolution() == NULL)
        return false;

    bool satisfiesAllConstraints = true;
    for (auto * brConstrPtr : _localBranchingConstraints)
    {
        Constraint * constrPtr = dynamic_cast<InstMasterBranchingConstr *>(brConstrPtr);

        if (constrPtr == NULL)
            continue;

        Double lhsValue = constrPtr->computeLhs(_masterCommons.debugSolution()->solVarValMap());
        if (((constrPtr->sense() == 'G') && (lhsValue + Double::precision < constrPtr->costrhs()))
            || ((constrPtr->sense() == 'L') && (lhsValue - Double::precision > constrPtr->costrhs())))
        {
            satisfiesAllConstraints = false;
        }
    }
    return satisfiesAllConstraints;
}

bool Alg4ProblemSetupBase::buildTreeOfColClasses(const std::list<CompSetInstMastBranchConstr*>& CSbrcList)
{
  bapcodInit().require(_nodePtr->treeOfColClasses().empty(),
                       "buildTreeOfColClasses(): treeOfColClasses should be empty");

  /**
   * We assume that each CompSetInstMastBranchConstr concerns a separate ColGenSpConf,
   * TODO: consider the case where they may concern several SP
   */

  if (CSbrcList.empty())
    return false;

  /// Split per ColGenSpConf
  std::map<ColGenSpConf *, std::list<CompSetInstMastBranchConstr *> > CSlistMap;

  for (std::list<CompSetInstMastBranchConstr *>::const_iterator curBcPt =
      CSbrcList.begin(); curBcPt != CSbrcList.end(); curBcPt++)
    {
      (*curBcPt)->reset();

      if (printL(5))
        std::cout << "CompSetInstMastBranchConstr::buildTreeOfColClasses(): candidat brConstr "
                  << (*curBcPt)->name() << " of spConf "
                  << ((*curBcPt)->ColGenSpConfPtr() != NULL ? (*curBcPt)->ColGenSpConfPtr()->name() : "undefined")
                  << std::endl;

      CSlistMap[(*curBcPt)->ColGenSpConfPtr()].push_back(*curBcPt);

    }

  bool masterInfeasible(false);

  for (std::map<ColGenSpConf *, std::list<CompSetInstMastBranchConstr *> >::const_iterator curBcPcPairPt =
      CSlistMap.begin(); curBcPcPairPt != CSlistMap.end(); curBcPcPairPt++)
    {
      ColGenSpConf * spcPtr = curBcPcPairPt->first;
      if (buildColGenSpConfTreeOfColClasses(spcPtr, curBcPcPairPt->second))
        {
          if (printL(5))
            std::cout << "ProblemBaseSetupAlgorithm::buildTreeOfColClasses(): ColGenSpConf infeasible "
                      << std::endl;

          masterInfeasible = true;
        }

      /// Record locTreeOfColClasses in ColGenSpConf
      if (spcPtr != NULL)
        {
          spcPtr->spConfHasClassInducedSpVarBounds(true); /// added by Ruslan
          ColClassesVector & cgSpTreeOfColClasses =
              _nodePtr->cgSpConfTreeOfColClassesMap()[spcPtr];
          _nodePtr->treeOfColClasses().insert(
              _nodePtr->treeOfColClasses().end(), cgSpTreeOfColClasses.begin(),
              cgSpTreeOfColClasses.end());

          Double locallyValidSpLowerBound(0);
          for (std::vector<CompSetInstMastBranchConstr *>::const_iterator bcPt =
              cgSpTreeOfColClasses.begin(); bcPt != cgSpTreeOfColClasses.end();
              bcPt++)
            {
              locallyValidSpLowerBound += (*bcPt)->marginLvalue();
            }
          spcPtr->lowerBoundMastConstrPtr()->defineLocalRhs(
              locallyValidSpLowerBound);
        }
    }

  return (masterInfeasible);
}

bool Alg4ProblemSetupBase::buildColGenSpConfTreeOfColClasses(ColGenSpConf * spcPtr,
                                                             const std::list<CompSetInstMastBranchConstr*>& CSbrcList)
{
  int printLevel = 5;
  ColClassesVector & treeOfColClasses = _nodePtr->cgSpConfTreeOfColClassesMap()[spcPtr];
  bapcodInit().require(treeOfColClasses.empty(),
                       "buildColGenSpConfTreeOfColClasses(): ColGenSpConfTreeOfColClasses should be empty");

  std::list<CompSetInstMastBranchConstr *> tempCSinTree;
  for (std::list<CompSetInstMastBranchConstr *>::const_iterator curBcPt = CSbrcList.begin();
       curBcPt != CSbrcList.end(); curBcPt++)
    {
      CompSetInstMastBranchConstr * removeBcPtr(NULL);

      for (std::list<CompSetInstMastBranchConstr *>::iterator existingBcPt = tempCSinTree.begin();
           existingBcPt != tempCSinTree.end(); existingBcPt++)
        {
          ComponentSequence::InclusionStatus res;
          res = compareCbS((*curBcPt)->compBoundSet(), (*existingBcPt)->compBoundSet());
          switch (res)
            {
          case ComponentSequence::different:
            {
              if (printL(printLevel))
                std::cout << " Column class " << (*curBcPt)->name() << " is independant of " << (*existingBcPt)->name()
                          << std::endl;

              /// Not comparable, continue;
              break;
            }
          case ComponentSequence::identical:
            {
              /// Keep only the tightest constraint
              if (printL(printLevel))
                std::cout << " Column class " << (*curBcPt)->name() << " is identical to " << (*existingBcPt)->name()
                          << std::endl;

              if ((*curBcPt)->sense() == 'G')
                {
                  if ((*curBcPt)->curRhs() > (*existingBcPt)->curRhs())
                    removeBcPtr = *existingBcPt;
                  else
                    /// That currently recorded must be removed
                    removeBcPtr = *curBcPt;
                }
              else
                {
                  if ((*curBcPt)->curRhs() < (*existingBcPt)->curRhs())
                    removeBcPtr = *existingBcPt;
                  else
                    /// That currently recorded must be removed
                    removeBcPtr = *curBcPt;
                }
              break;
            }
          case ComponentSequence::superclass:
            {
              /// (*curBcPt) is an ancestor of (*existingBcPt)
              if (printL(printLevel))
                std::cout << " Column class " << (*existingBcPt)->name() << " is a superclass of " << (*curBcPt)->name()
                          << std::endl;

              CompSetInstMastBranchConstr *unconst_curBcPt = *curBcPt;
              (*existingBcPt)->setOfPredCSconstrPtr().insert(unconst_curBcPt);
              break;
            }
          case ComponentSequence::subclass:
            {
              /// (*curBcPt) is child of (*existingBcPt)
              if (printL(printLevel))
                std::cout << " Column class " << (*curBcPt)->name() << " is a subclass of " << (*existingBcPt)->name()
                          << std::endl;

              (*curBcPt)->setOfPredCSconstrPtr().insert(*existingBcPt);
              break;
            }
            } /// Switch

        } /// List (*existingBcPt)

      if (printL(printLevel))
        std::cout << " add Column class to temp " << (*curBcPt)->name() << std::endl;

      tempCSinTree.push_back(*curBcPt);

      ///ReplaceColClass:
      if (removeBcPtr != NULL)
        {
          removeBcPtr->addToPreprocessedList();

          std::list<CompSetInstMastBranchConstr *>::iterator bcObsoletePt =
              tempCSinTree.end();

          if (printL(printLevel))
            std::cout << " remove Column class " << removeBcPtr->name() << std::endl;

          for (std::list<CompSetInstMastBranchConstr *>::iterator bcInPt =
              tempCSinTree.begin(); bcInPt != tempCSinTree.end(); bcInPt++)
            {
              if ((*bcInPt)->setOfPredCSconstrPtr().count(removeBcPtr))
                (*bcInPt)->setOfPredCSconstrPtr().erase(removeBcPtr);

              if ((*bcInPt) == removeBcPtr)
                bcObsoletePt = bcInPt;
            }
          bapcodInit().require(bcObsoletePt != tempCSinTree.end(),
              "buildColGenSpConfTreeOfColClasses should have found class to remove");

          if (printL(printLevel))
            std::cout << " remove Column class from temp " << (*bcObsoletePt)->name() << std::endl;

          tempCSinTree.erase(bcObsoletePt);
        }
    }

  /// Setup direct predecessor and compute marginLvalues
  int depth1(0);
  int depth2(0);
  for (std::list<CompSetInstMastBranchConstr *>::const_iterator bcitPt = tempCSinTree.begin();
       bcitPt != tempCSinTree.end(); bcitPt++)
    {
      if ((*bcitPt)->setOfPredCSconstrPtr().empty())
      /// predecessor is the root
        {
          continue;
        }

      depth1 = (*bcitPt)->setOfPredCSconstrPtr().size();
      for (CompSetInstBrConstrPtrSet::const_iterator bcPt =
          (*bcitPt)->setOfPredCSconstrPtr().begin();
          bcPt != (*bcitPt)->setOfPredCSconstrPtr().end(); bcPt++)
        {
          depth2 = (*bcPt)->setOfPredCSconstrPtr().size();
          if (printL(printLevel))
            std::cout << (*bcPt)->name() << " of depth " << depth2 << " is a predecessor of " << (*bcitPt)->name()
                      << " of depth " << depth1 << std::endl;

          if ((depth2 + 1) == depth1)
            {
              if (printL(printLevel))
                std::cout << (*bcPt)->name() << " is the direct predecessor of " << (*bcitPt)->name() << std::endl;

              (*bcitPt)->setDirPredCSconstrPtr(*bcPt);

              (*bcPt)->marginLvalue() -= (*bcitPt)->curRhs();
              break;
            }
        }
    }

  /**
   * Remove redundant classes, do not count global cardinality constraints
   * that are implemented directly as in the subproblem convexity constraint
   */
  bool thereAreRedundantClasses(false);
  for (std::list<CompSetInstMastBranchConstr *>::const_iterator bcitPt =
      tempCSinTree.begin(); bcitPt != tempCSinTree.end(); bcitPt++)
    {
      /// Redundant
      if ((!((*bcitPt)->marginLvalue()).positive(bapcodInit().param().BapCodReducedCostTolerance()))
          || ((*bcitPt)->compBoundSet().empty()))
        {
          thereAreRedundantClasses = true;
          for (std::list<CompSetInstMastBranchConstr *>::iterator bcPt =
              tempCSinTree.begin(); bcPt != tempCSinTree.end(); bcPt++)
            if ((*bcPt)->setOfPredCSconstrPtr().count(*bcitPt))
              (*bcPt)->setOfPredCSconstrPtr().erase(*bcitPt);

          (*bcitPt)->addToPreprocessedList();

          /// Memorise in cardinality constraints
          if ((*bcitPt)->compBoundSet().empty())
            {
              /// compBoundSet concern one SP only
              ColGenSpConf * spcPtr = (*bcitPt)->compBoundSet().cgSpConfPtr();
              if (spcPtr != NULL)
                {
                  if ((*bcitPt)->compBoundSet().activeSense() == 'G')
                    {
                      if (printL(printLevel))
                        std::cout << " Add Cardinality lower bound " << (*bcitPt)->name() << std::endl;

                      bapcodInit().require(spcPtr->lowerBoundMastConstrPtr() != NULL,
                                           "should inforce lower cardinality constraint");

                      if (spcPtr->upperBoundMastConstrPtr() != NULL)
                        {
                          /// Master problem infeasible
                          if ((*bcitPt)->curRhs()
                              /// we use newLocalRhs as curRhs() is out of date
                              > spcPtr->upperBoundMastConstrPtr()->newLocalRhs())
                            return (true);
                        }

                      if (printL(printLevel))
                        std::cout << " Set Cardinality lower bound in  " << spcPtr->lowerBoundMastConstrPtr()->name()
                                  << " to bound value = " << (*bcitPt)->costrhs() << std::endl;

                      spcPtr->lowerBoundMastConstrPtr()->defineLocalRhs(
                          (*bcitPt)->costrhs());
                    } /// Upper bound constraint
                  else
                    {
                      bapcodInit().require(
                          spcPtr->upperBoundMastConstrPtr() != NULL,
                          "should inforce upper cardinality constraint");

                      if (spcPtr->lowerBoundMastConstrPtr() != NULL)
                        {
                          /// Master problem infeasible
                          if ((*bcitPt)->curRhs()
                              /// we use newLocalRhs as curRhs() is out of date
                              < spcPtr->lowerBoundMastConstrPtr()->newLocalRhs())
                            return (true);
                        }

                      if (printL(printLevel))
                        std::cout << " Set Cardinality upper bound in  " << spcPtr->lowerBoundMastConstrPtr()->name()
                                  << std::endl;

                      spcPtr->upperBoundMastConstrPtr()->defineLocalRhs(
                          (*bcitPt)->costrhs());
                    }
                }
            }
        } /// Redundant
      else
        {
          if (printL(printLevel))
            std::cout << " add Column class to treeOfColClasses " << (*bcitPt)->name() << std::endl;

          treeOfColClasses.push_back(*bcitPt);
        }
    }

  /// Recompute direct predecessor if some class have been removed
  if (thereAreRedundantClasses)
    for (std::vector<CompSetInstMastBranchConstr *>::const_iterator bcitPt =
        treeOfColClasses.begin(); bcitPt != treeOfColClasses.end(); bcitPt++)
      {
        (*bcitPt)->setDirPredCSconstrPtr(NULL);
        if ((*bcitPt)->setOfPredCSconstrPtr().empty())
          continue;

        depth1 = (*bcitPt)->setOfPredCSconstrPtr().size();
        for (CompSetInstBrConstrPtrSet::const_iterator bcPt =
            (*bcitPt)->setOfPredCSconstrPtr().begin();
            bcPt != (*bcitPt)->setOfPredCSconstrPtr().end(); bcPt++)
          {
            depth2 = (*bcPt)->setOfPredCSconstrPtr().size();
            if (printL(printLevel))
              std::cout << (*bcPt)->name() << " of depth " << depth2 << " is a predecessor of " << (*bcitPt)->name()
                        << " of depth " << depth1 << std::endl;

            if ((depth2 + 1) == depth1)
              {
                if (printL(printLevel))
                  std::cout << (*bcPt)->name() << " is the direct predecessor of " << (*bcitPt)->name() << std::endl;

                (*bcitPt)->setDirPredCSconstrPtr(*bcPt);
                break;
              }
          }
      }

  return (false);
}

Alg4ProblemSetupBase::Alg4ProblemSetupBase(MasterCommons4ProblemSetup & masterCommons) :
  Alg4ProblemSetupOfNode(masterCommons), _probSetupInfoPtr(NULL), _problemInfeasible(false),
  _makeAllColumnsActive(false), _doSubproblemsSetup(true), _nonLinearBranchingConstraints(false),
  _newSpBranchingConstraints(), _oldSpBranchingConstraints(), _activeSpBranchingConstraintsExist(false),
  _varsToAddToForm(), _varsToRemoveFromForm(), _varsToChangeBounds(), _varsToChangeCost(), _constrsToAddToForm(),
  _constrsToRemoveFromForm(), _constrsToChangeRhs(), _localBranchingConstraints()
{
}

Alg4ProblemSetupBase::~Alg4ProblemSetupBase()
{
}

bool Alg4ProblemSetupBase::run(Node * nodePtr)
{
  _probSetupInfoPtr = nodePtr->probSetupInfoPtr();
  _problemInfeasible = Alg4ProblemSetupOfNode::run(nodePtr);
  return _problemInfeasible;
}

void Alg4ProblemSetupBase::setOptionDoSubproblemsSetup(const bool value)
{
  _doSubproblemsSetup = value;
}

void Alg4ProblemSetupBase::setOptionMakeAllColumnsActive(const bool value)
{
  _makeAllColumnsActive = value;
}

void Alg4ProblemSetupRootNode::resetConvexityConstraintsAtRoot()
{
  for (std::vector<ColGenSpConf *>::const_iterator spConfIt = _masterCommons.colGenSubProbConfPts().begin();
      spConfIt != _masterCommons.colGenSubProbConfPts().end(); ++spConfIt)
    {
      InstMastConvexityConstr * convConstrPtr = (*spConfIt)->lowerBoundMastConstrPtr();
      if ((convConstrPtr->curRhs() == InstMastConvexityConstr::lowerBoundWhenInactive)
           && (convConstrPtr->vcIndexStatus() == Active))
        deactivateConstraint(convConstrPtr, Inactive);
      convConstrPtr = (*spConfIt)->upperBoundMastConstrPtr();
      if ((convConstrPtr->curRhs() == InstMastConvexityConstr::upperBoundWhenInactive)
          && (convConstrPtr->vcIndexStatus() == Active))
        deactivateConstraint(convConstrPtr, Inactive);
    }
}

Alg4ProblemSetupRootNode::Alg4ProblemSetupRootNode(MasterCommons4ProblemSetup & masterCommons) :
    Alg4ProblemSetupBase(masterCommons)
{
}

Alg4ProblemSetupRootNode::~Alg4ProblemSetupRootNode()
{
}

bool Alg4ProblemSetupRootNode::run(Node * nodePtr)
{
  if (printL(ProblemSetupPrintLevel))
    std::cout << "ProblemRootSetupAlgorithm::run()" << std::endl;

  _problemInfeasible = Alg4ProblemSetupBase::run(nodePtr);

  resetConvexityConstraintsAtRoot();
  resetMasterColumns();
  resetNonStabArtificialVariables();

  updateFormulation(_masterProbPtr);

  _nodePtr = NULL;
  return _problemInfeasible;
}

Alg4ProblemSetupBranchingOnly::Alg4ProblemSetupBranchingOnly(MasterCommons4ProblemSetup & masterCommons) :
        Alg4ProblemSetupBase(masterCommons)
{
}

Alg4ProblemSetupBranchingOnly::~Alg4ProblemSetupBranchingOnly()
{
}

bool Alg4ProblemSetupBranchingOnly::run(Node * nodePtr)
{
    if (printL(ProblemSetupPrintLevel))
        std::cout << "ProblemBranchingOnlySetupAlgorithm::run()" << std::endl;

    _problemInfeasible = Alg4ProblemSetupBase::run(nodePtr);

    applySubproblemsInfo();
    fillLocalBranchingConstraints();
    prepareBranchingConstraints();
    if (nodePtr->debugSolutionAtThisNode() && !debugSolSatisfiesLocalBranchConstrs())
        nodePtr->removeDebugSolutionFromThisNode();

    std::list<ConstraintInfo>::const_iterator branchingConstrPtrIt =
            _probSetupInfoPtr->activeBranchingConstraintsInfo.begin();

    resetBranchingConstraints(_masterProbPtr, branchingConstrPtrIt);
    /// should be done after resetBranchingConstraints as component set branching constraints
    /// can modify convexity constraints
    resetConvexityConstraints();
    resetNonStabArtificialVariables();

    bool needToUpdateSpBranchConstrs = !_newSpBranchingConstraints.empty() || !_oldSpBranchingConstraints.empty()
                                       || _activeSpBranchingConstraintsExist;
    if (needToUpdateSpBranchConstrs)
        resetMasterColumns();

    updateFormulation(_masterProbPtr);

    bool considerSubproblems = needToUpdateSpBranchConstrs
                               || (_nonLinearBranchingConstraints && _doSubproblemsSetup);

    if (considerSubproblems)
    {
        for (std::list<Problem *>::const_iterator probPtrIt = ++_masterCommons.problemList().begin();
             probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
        {
            if (needToUpdateSpBranchConstrs)
            {
                resetBranchingConstraints(*probPtrIt, branchingConstrPtrIt);
                if (bapcodInit().param().colGenSubProbSolMode().status() != SolutionMethod::customSolver)
                    updateFormulation(*probPtrIt);
                else
                    clearVarConstrLists();
            }
            /// we need always to pass the state of the subproblem, otherwise taking into account branching
            /// constraints may not be correct. This may happen just after the first phase of strong branching,
            /// as full setup is done for the first node of the first phase of strong branching
            std::map<const Problem *, BcSolverOracleInfo *>::const_iterator mapIt
                    = _probSetupInfoPtr->solverOracleInfo.find(*probPtrIt);
            const BcSolverOracleInfo *solverOracleInfoPtr = NULL;
            if (mapIt != _probSetupInfoPtr->solverOracleInfo.end())
                solverOracleInfoPtr = mapIt->second;
            _problemInfeasible = _problemInfeasible || (*probPtrIt)->setupCustomizedSolver(solverOracleInfoPtr);
        }
    }

    _localBranchingConstraints.clear();
    _nodePtr = NULL;
    return _problemInfeasible;
}

Alg4ProblemSetupFull::Alg4ProblemSetupFull(MasterCommons4ProblemSetup & masterCommons) :
        Alg4ProblemSetupBase(masterCommons)
{
}

Alg4ProblemSetupFull::~Alg4ProblemSetupFull()
{
}

/// this only for the master for the moment
void Alg4ProblemSetupFull::resetPartialSolution(Problem * probPtr)
{
    probPtr->resetPartialSolution();

    Alg4ProblemSetupOfNode::resetPartialSolution(probPtr);

    for (std::list<VariableSolInfo>::const_iterator varInfoIt = _probSetupInfoPtr->masterPartialSolutionInfo.begin();
         varInfoIt != _probSetupInfoPtr->masterPartialSolutionInfo.end(); ++varInfoIt)
    {
        varInfoIt->applyVarInfo();
    }
}

void Alg4ProblemSetupFull::resetMaster(std::list<VariableInfo*>::const_iterator & varPtrIt,
                                       std::list<ConstraintInfo>::const_iterator & constrPtrIt,
                                       std::list<ConstraintInfo>::const_iterator & branchConstrPtrIt)
{
    resetBranchingConstraints(_masterProbPtr, branchConstrPtrIt);
    /// should be done after resetBranchingConstraints as component set branching constraints
    /// can modify convexity constraints
    resetConvexityConstraints();
    resetStaticVariables(_masterProbPtr, varPtrIt);
    resetStaticConstraints(_masterProbPtr, constrPtrIt);
    resetMasterCuts();
    resetNonStabArtificialVariables();

    updateFormulation(_masterProbPtr);
}

void Alg4ProblemSetupFull::resetSubproblem(Problem * probPtr,
                                           std::list<VariableInfo*>::const_iterator & varPtrIt,
                                           std::list<ConstraintInfo>::const_iterator & constrPtrIt,
                                           std::list<ConstraintInfo>::const_iterator & branchConstrPtrIt,
                                           const BcSolverOracleInfo * bcSolverOracleInfoPtr)
{
    resetBranchingConstraints(probPtr, branchConstrPtrIt);
    resetStaticVariables(probPtr, varPtrIt);
    resetStaticConstraints(probPtr, constrPtrIt);

    bool testOracles = bapcodInit().param().CheckOracleOptimality() || bapcodInit().param().CheckSpOracleFeasibility();
    if(testOracles || bapcodInit().param().colGenSubProbSolMode().status() != SolutionMethod::customSolver)
        updateFormulation(probPtr);
    else
        clearVarConstrLists();

    _problemInfeasible = _problemInfeasible || probPtr->setupCustomizedSolver(bcSolverOracleInfoPtr);
}

/// TO DO: for the moment, only master branching constraints are supported
bool Alg4ProblemSetupFull::run(Node * nodePtr)
{
    if (printL(ProblemSetupPrintLevel))
        std::cout << "ProblemFullSetupAlgorithm::run()" << std::endl;

    _problemInfeasible = Alg4ProblemSetupBase::run(nodePtr);

    std::list<VariableInfo*>::const_iterator staticVarPtrIt = _probSetupInfoPtr->modifiedStaticVarsInfo.begin();
    std::list<ConstraintInfo>::const_iterator staticConstrPtrIt = _probSetupInfoPtr->modifiedStaticConstrsInfo.begin();
    std::list<ConstraintInfo>::const_iterator branchingConstrPtrIt =
            _probSetupInfoPtr->activeBranchingConstraintsInfo.begin();

    applySubproblemsInfo();
    fillLocalBranchingConstraints();
    prepareBranchingConstraints();
    if (nodePtr->debugSolutionAtThisNode() && !debugSolSatisfiesLocalBranchConstrs())
        nodePtr->removeDebugSolutionFromThisNode();

    resetMaster(staticVarPtrIt,staticConstrPtrIt, branchingConstrPtrIt);

    for (std::list<Problem *>::const_iterator probPtrIt = ++_masterCommons.problemList().begin();
         probPtrIt != _masterCommons.problemList().end(); ++probPtrIt)
    {
        std::map<const Problem *, BcSolverOracleInfo *>::const_iterator mapIt
                = _probSetupInfoPtr->solverOracleInfo.find(*probPtrIt);
        const BcSolverOracleInfo * solverOracleInfoPtr = NULL;
        if (mapIt != _probSetupInfoPtr->solverOracleInfo.end())
            solverOracleInfoPtr = mapIt->second;
        resetSubproblem(*probPtrIt, staticVarPtrIt, staticConstrPtrIt,
                        branchingConstrPtrIt, solverOracleInfoPtr);
    }

    /// master columns should be reset after subproblems
    /// this implementation is not the best (updateFormulation is called two times for the master)
    /// one should de complete master reset here, but this is not compatible for now
    /// with the front position of the master in the problemList()
    /// TO DO: change the position of the master in the problemList() to the back one and update
    /// the code (ProblemSetup and Preprocessing)
    resetMasterColumns();
    updateFormulation(_masterProbPtr);

    _localBranchingConstraints.clear();
    _nodePtr->removeProbSetupObligationForFullSetup();
    _nodePtr = NULL;
    return _problemInfeasible;
}
