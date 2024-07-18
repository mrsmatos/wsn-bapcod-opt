/**
 *
 * This file bcAlg4Master.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCALGFORMASTER_HPP_
#define BCALGFORMASTER_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcProbConfigC.hpp"
#include "bcBoundLevC.hpp"


struct Alg4MasterSolAndBounds
{
  //Within a method in a subclass of Alg4Master , all _alg* attributes are considered to have
  //the meaning of the class that defines the called method.
  //As a consequence of that, in case of inheritance, an algorithm needs to copy its bounds before
  //calling the eval of a super class, update those copied bounds afterwards.

  Bound algIncStageLpPrimalBound; /// incumbent primal bound for the master LP at the current stage;
  Bound algIncStageLpDualBound; /// incumbent dual bound for the master LP at the current stage;

  Bound algIncLpPrimalBound; /// incumbent valid dual bound for the master LP
  Bound algIncLpDualBound;  /// incumbent valid dual bound for the master LP // _nodeFractionalDualBound;
  /// IMPORTANT : algIncLpPrimalBound IS FOR THE CURRENT PROBLEM WITHOUT CORE CUTS THAT ARE
  /// CHECKED ONLY IF AND INTEGER SOLUTION IS FOUND

  Bound algIncIpPrimalBound; /// incumbent valid dual bound for the master IP // stored also in Node
  Bound algIncIpDualBound; /// incumbent valid dual bound for the master IP //_curIpDualVal;

  //These values are tempory to the current Alg4Eval and might not be fully validated.
  //for example algIncIpPrimalBound might not have validated core cuts or benders cuts.

  //Issam: All algorithms should use these solutions but for now some of the algithms
  //still store the solutions in Problem.
  VarPtr2DoubleMap algIncLpPrimalSolMap;
  ConstrPtr2DoubleMap algIncLpDualSolMap;
  VarPtr2DoubleMap algIncIpPrimalSolMap;

  bool algIncIpPrimalBoundUpdated;

  Alg4MasterSolAndBounds(BcObjStatus::MinMaxIntFloat objStatus) :
    algIncStageLpPrimalBound(Bound::infPrimalBound(objStatus)),
    algIncStageLpDualBound(Bound::infDualBound(objStatus)),
    algIncLpPrimalBound(Bound::infPrimalBound(objStatus)),
    algIncLpDualBound(Bound::infDualBound(objStatus)),
    algIncIpPrimalBound(Bound::infPrimalBound(objStatus)),
    algIncIpDualBound(Bound::infDualBound(objStatus)),
    algIncIpPrimalBoundUpdated(false)
  {}

};

class Alg4Master
{
  friend class Node;

private:
  /// Issam : setters and getters of the elements in _solAndBnds are only definied for cur values.
  /// setting the other bounds should be done using setCurLpDualBoundAndUpdateDualBounds.
  /// This private encapsulation is crucial to keep the code tractable. If you feel like
  /// switching _solAndBnds to protected or public, or add public getters/setters,
  /// please think twice before doing it!
  Alg4MasterSolAndBounds _solAndBnds;

protected:
  /// temporary master value that can be a primal or a dual bound depending on the method //MARKFV
  /// Sometime we can have an estimation of a current bound. That s why we have the following
  /// 2 curBounds. But often only one of the current bounds is actually the mastCurObjVal.
  Bound _algCurLpPrimalBound;
  Bound _algCurLpDualBound;

  bool _solIsMasterLpFeasible;

  Problem * const _probPtr;
  Node * _currentNodePtr;

  VarPtr2DoubleMap _algCurLpPrimalSolMap;
  ConstrPtr2DoubleMap _algCurLpDualSolMap;

  int _maxLevelOfSubProbRestriction;

  void backupBoundsExceptIncIpPrimalBound(Alg4MasterSolAndBounds & backupBounds);

  void restoreBoundsExceptIncIpPrimalBound(const Alg4MasterSolAndBounds & backupBounds);

  Double totalObjVal() const
  {return _probPtr->objVal() + _probPtr->partialSolutionValue();}

  /// Issam : setters and getters of the elements in _solAndBnds are only definied for cur values.
  /// setting the other bounds should be done using setCurLpDualBoundAndUpdateDualBounds.
  /// This private encapsulation is crucial to keep the code tractable. If you feel like
  /// switching _solAndBnds to protected or public, or add public getters/setters,
  /// please think twice before doing it!
  const Alg4MasterSolAndBounds & algSolAndBnds() const
  {return _solAndBnds; }

  ///solution getters
  const VarPtr2DoubleMap & algCurLpPrimalSolMap() const
  { return _algCurLpPrimalSolMap;}

  const ConstrPtr2DoubleMap & algCurLpDualSolMap() const

  { return   _algCurLpDualSolMap;}

  const VarPtr2DoubleMap & algIncLpPrimalSolMap() const
  { return _solAndBnds.algIncLpPrimalSolMap;}

  const ConstrPtr2DoubleMap & algIncLpDualSolMap()const
  { return _solAndBnds.algIncLpDualSolMap;}

public:

  ///this only getter is public because the node needs to use it at the end of the algorithm.
  const VarPtr2DoubleMap & algIncIpPrimalSolMap()
  { return _solAndBnds.algIncIpPrimalSolMap;}

  bool algIncIpPrimalBoundUpdated() const
  {
    return _solAndBnds.algIncIpPrimalBoundUpdated;
  }

private:

  void clearAlgIncIpPrimalSolMap();

protected:

  /// bounds getters
  const Bound & algCurLpPrimalBound() const
  { return  _algCurLpPrimalBound; }

  const Bound & algCurLpDualBound() const
  { return  _algCurLpDualBound; }

  const Bound & algIncStageLpPrimalBound()const
  { return  _solAndBnds.algIncStageLpPrimalBound; }

  const Bound & algIncStageLpDualBound()  const
  { return  _solAndBnds.algIncStageLpDualBound; }

  const Bound & algIncLpPrimalBound()     const
  { return  _solAndBnds.algIncLpPrimalBound; }

  const Bound & algIncLpDualBound()       const
  { return  _solAndBnds.algIncLpDualBound; }

  const Bound & algIncIpPrimalBound()     const
  { return  _solAndBnds.algIncIpPrimalBound; }

  const Bound & algIncIpDualBound()       const
  { return  _solAndBnds.algIncIpDualBound; }

  void resetAlgIncLpPrimalBound(BcObjStatus::MinMaxIntFloat objStatus)
  {
    _solAndBnds.algIncLpPrimalBound = Bound::infPrimalBound(objStatus);
    _solAndBnds.algIncStageLpPrimalBound = Bound::infPrimalBound(objStatus);
  }

  void resetAlgIncStageLpDualBound()
  { _solAndBnds.algIncStageLpDualBound = _currentNodePtr->nodeIncLpDualBound(); }

  /// The two following update methods assumes that the curLpDualBounds/sols correspond to a
  /// feasible solution of the current problem.
  bool updateAlgPrimalLpBounds();
  bool updateAlgDualBounds();

  void updatePrimalIpSolAndBnds(const VarPtrSet & primalSol, const VarPtr2DoubleMap & partialSol);

  void rectifyIncumbentLpValue();

  inline const Double computeOptimGap() const
  {
    return computeOptimalityGap(_solAndBnds.algIncLpDualBound,_solAndBnds.algIncLpPrimalBound);
  }

  void markInfeasible();


  BapcodInit& bapcodInit() const
  {
    return _probPtr->bapcodInit();
  }

  ControlParameters & param() const
  {
    return bapcodInit().param();
  }

  ProgStatus& progStatus() const
  {
    return bapcodInit().progStatus();
  }

public:

  Alg4Master(Problem * const probPtr);

  virtual ~Alg4Master()
  {
  }

  virtual bool setupAlgo(Node * nodePtr);

  virtual void setDownAlgo();
};

#endif /* BCALGFORMASTER_HPP_ */
