/**
 *
 * This file bcAlg4ProblemSetup.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPROBSETUPROOT_HPP_
#define BCPROBSETUPROOT_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcVarConstrC.hpp"
#include "bcProblemC.hpp"
#include "bcMasterConfC.hpp"
#include "bcSpVarConstrC.hpp"

using namespace VcIndexStatus;


struct VariableSmallInfo
{
public:
  Variable * const varPtr;
  VcStatus status;
  Double cost;

  VariableSmallInfo(Variable * varPtrV, const VcStatus statusV = Active) :
      varPtr(varPtrV), status(statusV), cost(varPtrV->curCost())
  {
  }

  virtual ~VariableSmallInfo()
  {
  }

  virtual void applyVarInfo() const
  {
    varPtr->resetCurCostByValue(cost);
  }
};

struct VariableInfo : public VariableSmallInfo
{
public:
  Double lb;
  Double ub;

  VariableInfo(Variable * varPtrV, const VcStatus statusV = Active) :
      VariableSmallInfo(varPtrV, statusV), lb(varPtrV->globalCurLb()),
          ub(varPtrV->globalCurUb())
  {
  }

  virtual void applyVarInfo() const
  {
    varPtr->resetCurCostByValue(cost);
    varPtr->globalCurLb(lb);
    varPtr->globalCurUb(ub);
  }

  virtual bool needToChangeBounds() const
  {
    return (varPtr->inCurForm() && ((lb != varPtr->globalCurLb()) || (ub != varPtr->globalCurUb())));
  }

};

struct SpVariableInfo : public VariableInfo
{
public:
  Double localLb;
  Double localUb;

  SpVariableInfo(SubProbVariable * varPtrV, const VcStatus statusV = Active) :
      VariableInfo(varPtrV, statusV)
  {
    SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(varPtrV);
    localLb = spVarPtr->localCurLb();
    localUb = spVarPtr->localCurUb();
  }

  virtual void applyVarInfo() const
  {
    SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(varPtr);
    spVarPtr->resetCurCostByValue(cost);
    spVarPtr->globalCurLb(lb);
    spVarPtr->globalCurUb(ub);
    spVarPtr->localCurLb(localLb);
    spVarPtr->localCurUb(localUb);
  }

  virtual bool needToChangeBounds() const
  {
    SubProbVariable * spVarPtr = static_cast<SubProbVariable *>(varPtr);
    return (spVarPtr->inCurForm() && ((localLb != spVarPtr->localCurLb()) || (localUb != spVarPtr->localCurUb())));
  }

};

struct ConstraintInfo
{
  Constraint * const constrPtr;
  VcStatus status;
  Double minSlack;
  Double maxSlack;
  Double rhs;

  ConstraintInfo(Constraint * constrPtrV, const VcStatus statusV = Active) :
      constrPtr(constrPtrV), status(statusV),
          minSlack(constrPtrV->curMinSlack()),
          maxSlack(constrPtrV->curMaxSlack()), rhs(constrPtrV->curRhs())
  {
  }

  void applyInfo() const
  {
    constrPtr->curMinSlack(minSlack);
    constrPtr->curMaxSlack(maxSlack);
    constrPtr->curRhs(rhs);
  }
};

struct SubProblemInfo
{
  ColGenSpConf * spConfPtr;
  Double lb;
  Double ub;

  SubProblemInfo(ColGenSpConf * spConfPtrV);
};

class BcSolverOracleInfo;

struct ProblemSetupInfo
{
  int treatOrderId;
  int numberOfNodes;
  bool fullSetupIsObligatory;

  std::list<VariableSmallInfo> suitableMasterColumnsInfo;
  std::list<ConstraintInfo> suitableMasterCutsInfo;
  std::list<ConstraintInfo> activeBranchingConstraintsInfo;
  std::list<SubProblemInfo> subProblemsInfo;
  std::list<VariableSolInfo> masterPartialSolutionInfo;


  /// - In these two lists we keep only static variables and constraints for
  ///   which at least one of the attributes in VariableInfo and ConstraintInfo is different from the default.
  ///   Default values are set by the user and can be changed by the preprocessing at the root
  /// - Unsuitable static variables or constraints are ignored: they are eliminated by the preprocessed at the root
  /// - We keep variables and constraints in the strict order: master -> subprob 1 -> subprob 2 -> ...

  std::list<VariableInfo*> modifiedStaticVarsInfo; // Pointers are needed for polymorphism. A cleaner solution
                                                   // would be : http://kremer.cpsc.ucalgary.ca/STL/1024x768/ref2.html
  std::list<ConstraintInfo> modifiedStaticConstrsInfo;

  std::map<const Problem *, BcSolverOracleInfo *> solverOracleInfo;

  ProblemSetupInfo(const int treatOrder);
  ProblemSetupInfo(const int treatOrder, const std::vector<Variable *> & initialSetOfActiveColumns,
                   const std::vector<Variable *> & initialSetOfInactiveColumns);
  void setFullSetupIsObligatory(const bool & value);

  ~ProblemSetupInfo();
};

class Alg4ProblemSetDownOfNode
{
protected:
  MasterCommons4ProblemSetup & _masterCommons;
  Problem * _masterProbPtr;
  BapcodInit & bapcodInit() const
  {
    return _masterProbPtr->bapcodInit();
  }

public:
  Alg4ProblemSetDownOfNode(MasterCommons4ProblemSetup & masterCommons);
  virtual ~Alg4ProblemSetDownOfNode();
  virtual void run();
  virtual ProblemSetupInfo * recordProblemInfo(int globalTreatOrder = -1);
};

class ProblemFullSetDownAlgorithm: public Alg4ProblemSetDownOfNode
{
public:
  ProblemFullSetDownAlgorithm(MasterCommons4ProblemSetup & masterCommons);
  virtual ~ProblemFullSetDownAlgorithm();
  void run();
  virtual ProblemSetupInfo * recordProblemInfo(int globalTreatOrder = -1);
};

class Alg4ProblemSetupOfNode
{
protected:
  Node * _nodePtr;
  MasterCommons4ProblemSetup & _masterCommons;
  Problem * _masterProbPtr;
  BapcodInit & bapcodInit() const
  {
    return _nodePtr->bapcodInit();
  }
  virtual void resetPartialSolution(Problem * probPtr);
public:
  Alg4ProblemSetupOfNode(MasterCommons4ProblemSetup & masterCommons);
  virtual ~Alg4ProblemSetupOfNode();
  virtual bool run(Node * nodePtr);
};

/// base class, does not have public methods
class Alg4ProblemSetupBase : public Alg4ProblemSetupOfNode
{
protected:
  const ProblemSetupInfo * _probSetupInfoPtr;
  bool _problemInfeasible;
  bool _makeAllColumnsActive;
  bool _doSubproblemsSetup;
  bool _nonLinearBranchingConstraints;
  std::list<Constraint *> _newSpBranchingConstraints;
  std::list<Constraint *> _oldSpBranchingConstraints;
  bool _activeSpBranchingConstraintsExist;

  std::list<Variable *> _varsToAddToForm;
  std::list<Variable *> _varsToRemoveFromForm;
  std::list<Variable *> _varsToChangeBounds;
  std::list<Variable *> _varsToChangeCost;
  std::list<Constraint *> _constrsToAddToForm;
  std::list<Constraint *> _constrsToRemoveFromForm;
  std::list<Constraint *> _constrsToChangeRhs;
  std::list<Constraint *> _localBranchingConstraints;

  void deactivateVariable(Variable * varPtr, const VcStatus & status,
      bool removeFromForm = true);
  void activateVariable(Variable * varPtr);
  void deactivateConstraint(Constraint * constrPtr, const VcStatus & status,
      bool removeFromForm = true);
  void activateConstraint(Constraint * constrPtr);

  void resetStaticVariables(Problem * probPtr,
      std::list<VariableInfo*>::const_iterator & varInfoIt);

  void resetStaticConstraints(Problem * probPtr,
      std::list<ConstraintInfo>::const_iterator & constrInfoIt);

  bool columnIsUnsuitable(MastColumn * colPtr);
  bool columnViolatesNewSpBranchingConstriants(MastColumn * colPtr);
  void resetMasterColumns();

  void prepareBranchingConstraints();  /// @todo must add  SpVar affected by CompSet Branching constraints in Problem

  void resetBranchingConstraints(Problem * probPtr,
      std::list<ConstraintInfo>::const_iterator & constrInfoIt);
  void resetConvexityConstraints();
  void resetNonStabArtificialVariables();

  void resetMasterCuts();

  void printVarsList(std::ostream & os, const std::list<Variable *> & varsList, const std::string & stringToPrint);
  void printConstrsList(std::ostream & os, const std::list<Constraint *> & varsList, const std::string & stringToPrint);
  void clearVarConstrLists();
  void updateFormulation(Problem * probPtr);

  void applySubproblemsInfo();

  void fillLocalBranchingConstraints();

  bool debugSolSatisfiesLocalBranchConstrs();

  bool buildTreeOfColClasses(const std::list<CompSetInstMastBranchConstr*>& CSbrcList);

  bool buildColGenSpConfTreeOfColClasses(ColGenSpConf * spcPtr,
                                         const std::list<CompSetInstMastBranchConstr*>& CSbrcList);


  Alg4ProblemSetupBase(MasterCommons4ProblemSetup & masterCommons);
  virtual ~Alg4ProblemSetupBase();
  virtual bool run(Node * nodePtr);
public:
  void setOptionMakeAllColumnsActive(const bool value);
  void setOptionDoSubproblemsSetup(const bool value);
};

class Alg4ProblemSetupRootNode : public Alg4ProblemSetupBase
{
  void resetConvexityConstraintsAtRoot();
public:
  Alg4ProblemSetupRootNode(MasterCommons4ProblemSetup & masterCommons);
  virtual ~Alg4ProblemSetupRootNode();
  virtual bool run(Node * nodePtr);
};

class Alg4ProblemSetupBranchingOnly : public Alg4ProblemSetupBase
{
public:
    Alg4ProblemSetupBranchingOnly(MasterCommons4ProblemSetup & masterCommons);
    virtual ~Alg4ProblemSetupBranchingOnly();
    virtual bool run(Node * nodePtr);
};

class Alg4ProblemSetupFull : public Alg4ProblemSetupBase
{
    virtual void resetPartialSolution(Problem * probPtr);
    void resetMaster(std::list<VariableInfo*>::const_iterator & varPtrIt,
                     std::list<ConstraintInfo>::const_iterator & constrPtrIt,
                     std::list<ConstraintInfo>::const_iterator & branchConstrPtrIt);
    void resetSubproblem(Problem * probPtr,
                         std::list<VariableInfo*>::const_iterator & varPtrIt,
                         std::list<ConstraintInfo>::const_iterator & constrPtrIt,
                         std::list<ConstraintInfo>::const_iterator & branchConstrPtrIt,
                         const BcSolverOracleInfo * bcSolverOracleInfoPtr);
public:
    Alg4ProblemSetupFull(MasterCommons4ProblemSetup & masterCommons);
    virtual ~Alg4ProblemSetupFull();
    virtual bool run(Node * nodePtr);
};

#endif
