/**
 *
 * This file bcAlg4PreprocessingOfNode.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPREPROCESSING_HPP__
#define BCPREPROCESSING_HPP__

#include "bcUsefulHeadFil.hpp"
#include "bcVarConstrC.hpp"
#include "bcPrintC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcStabilizationInfo.hpp"
#include "bcNodeC.hpp"

using namespace VcIndexStatus;


class Alg4PreprocessingOfNode
{
  ConstrPtrList _constrsListToPropagate;
  std::list<Constraint *> _mastConstrsToChangeRhs;
  std::map<ProbConfig*, std::list<Variable*> > _pcPtrToNonZeroVarPtrMap;
  std::set<ColGenSpConf *> _colGenSpPtsWithZeroUb;
  bool _needToPreprocessCompSetBrConstr;

protected:
  std::list <Problem *> const & _problemPts;
private:

  inline bool lessConstraint(Constraint * constrPtr)
  {
    return (!constrPtr->considerAsEqualityInPreprocessing() && (constrPtr->sense() == 'L'));
  }
  inline bool greaterConstraint(Constraint * constrPtr)
  {
    return (!constrPtr->considerAsEqualityInPreprocessing() && (constrPtr->sense() == 'G'));
  }
  inline bool lessOrEqualConstraint(Constraint * constrPtr)
  {
    return (constrPtr->considerAsEqualityInPreprocessing() || (constrPtr->sense() != 'G'));
  }
  inline bool greaterOrEqualConstraint(Constraint * constrPtr)
  {
    return (constrPtr->considerAsEqualityInPreprocessing() || (constrPtr->sense() != 'L'));
  }
  bool updateMaxSlack (VarConstr * varConstrPtr, Double const delta);
  bool updateMinSlack (VarConstr * varConstrPtr, Double const delta);
  bool updateLowerBound (Variable * varPtr,  Double newBound, Constraint * modifyingConsPtr = NULL, bool const isSpVar = false);
  bool updateUpperBound (Variable * varPtr,  Double newBound, Constraint * modifyingConsPtr = NULL, bool const isSpVar = false);
protected:
  bool updateGlobalLowerBound (SubProbVariable * varPtr, const Double & newBound, Constraint * modifyingConsPtr = NULL);
  bool updateGlobalUpperBound (SubProbVariable * varPtr, const Double & newBound, Constraint * modifyingConsPtr = NULL);
  bool updateLocalUpperBound (SubProbVariable * varPtr);
  bool updateLocalUpperBound (SubProbVariable * varPtr, Double newBound, Constraint * modifyingConsPtr = NULL);
  bool updateLocalLowerBound (SubProbVariable * varPtr);
  bool updateLocalLowerBound (SubProbVariable * varPtr, Double newBound, Constraint * modifyingConsPtr = NULL);
  bool computeCompSetBrConstrInducedGlobalBdOnSpVar();
private:
  bool propagate ();
  bool columnBecameUnsuitable_SolValNotWithinBounds(MastColumn * colPtr);
  bool columnBecameUnsuitable_InexistantNonZeroVar(MastColumn * colPtr);
  void applyPreprocessingListsInProbAndForm (bool const initialPreprocessing);
  void clearPreprocessingLists ();
  bool computeInitialConstrsSlacks();
  bool exitWhenInfeasible ();
  bool fixVariableValue(Variable * varPtr, const Double & value);
  void changeSubProblemBounds(ColGenSpConf * spConfPtr, const Double & value);
  bool propagateNonLinearMasterConstraints(const MastColumn * const colPtr, const Double & value);
  bool fixPartialSolution (Solution const * solPtr);
  bool preprocess (Solution const * solPtr, bool const initialPreprocessing);

  void deactivateLocalArtVarsOfConstr(Problem * probPtr, Constraint * constrPtr,
                                      const VcStatus & status, VarPtrList & varsToRemoveFromForm);
  void deactivateLocalArtVar(Problem * probPtr, LocalArtificialVar * artVarPtr,
                                                 const VcStatus & status, VarPtrList & varsToRemoveFromForm);

  bool initialUpdateOfSpVarBounds();

  const ControlParameters & param()
  {
    return _problemPts.front()->param();
  }

protected:
  bool preprocessRoot();
  bool preprocessConstraints (ConstrPtrList constrsToPropagate);
  bool preprocessAfterFixingPartialSolution(Solution const * solPtr);

public:
  Alg4PreprocessingOfNode (Alg4PreprocessingOfNode & that);
  Alg4PreprocessingOfNode (std::list <Problem *> const & problemPts);
  virtual ~Alg4PreprocessingOfNode(){}
  virtual bool run(Node* node) = 0;
} ;

class Algorithm4PreprocessingAtRoot : public Alg4PreprocessingOfNode
{
  public:
    Algorithm4PreprocessingAtRoot (std::list <Problem *> const & problemPts);
    virtual ~Algorithm4PreprocessingAtRoot(){}
    virtual bool run(Node* node){
      return preprocessRoot();
    }
};

class Algorithm4PreprocessingAtNodeOtherThanRoot : public Alg4PreprocessingOfNode
{
  public:
    Algorithm4PreprocessingAtNodeOtherThanRoot (std::list <Problem *> const & problemPts);
    virtual ~Algorithm4PreprocessingAtNodeOtherThanRoot(){}

    virtual bool run(Node* node){
      return preprocessConstraints(node->upCastedBranchingConstrList());
    }
};

class Algorithm4PreprocessingInDive : public Alg4PreprocessingOfNode
{
public:
  Algorithm4PreprocessingInDive(std::list<Problem *> const & problemPts) :
      Alg4PreprocessingOfNode(problemPts)
  {
  }
  virtual ~Algorithm4PreprocessingInDive(){}
  virtual bool run(Node * nodePtr)
  {
    return preprocessAfterFixingPartialSolution(nodePtr->localFixedSolution());
  }
};


#endif // BCPREPROCESSING_HPP__
