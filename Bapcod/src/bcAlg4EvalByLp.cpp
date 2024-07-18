/**
 *
 * This file bcAlg4EvalByLp.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <vector>

#include "bcUsefulHeadFil.hpp"
#include "bcAlg4EvalByLp.hpp"
#include "bcStabilizationInfo.hpp"

using namespace std;

Alg4EvalByLp::Alg4EvalByLp(Problem * const probPtr, MasterCommons4EvalAlg & masterCommons) :
  Alg4EvalOfNode(probPtr, masterCommons), _updateIncDualBound(true), _keptStabInfoPtr(NULL),
  _keptLatestRedCostFixingGap(BapcodInfinity)
{                 
}


int Alg4EvalByLp::solveRestrictedMastLP()
{
  int formPrintLevel = 2;

  Time startbcTimeMastMPsol;
  int maxLevelOfRestriction(0);
  /// 'd' means do record the dual solution
  int solverReturnStatus(_probPtr->Problem::solveProb(maxLevelOfRestriction, 'd', printL(formPrintLevel)));
  
  _algCurLpPrimalBound = _algCurLpDualBound = Bound(totalObjVal(), _masterCommons.objStatus());
  updateAlgPrimalLpBounds();
  if (_updateIncDualBound)
      updateAlgDualBounds();
    
  bapcodInit().statistics().incrTimer("bcTimeMastMPsol", startbcTimeMastMPsol.getElapsedTime_dbl());
  if (printL(1))
    std::cout << " Restricted master LP is solved in " << startbcTimeMastMPsol.getElapsedTime_dbl()/100.0
              << " seconds" << std::endl;

  bapcodInit().statistics().incrRecord("bcAverageDualSolSize", int(_probPtr->inDualSol().size()));
  bapcodInit().statistics().incrCounter("bcCountMastSol");
  bapcodInit().statistics().setCounter("bcCountPrimalSolSize", int(_probPtr->inPrimalLpSol().size()));
  bapcodInit().statistics().incrRecord("bcAveragePrimalSolSize", int(_probPtr->inPrimalLpSol().size()));

  return (solverReturnStatus);
}

bool Alg4EvalByLp::setupAlgo(Node * nodePtr)
{
  if (Alg4EvalOfNode::setupAlgo(nodePtr))
    return true;

  ColGenEvalInfo * colGenEvalInfoPtr = dynamic_cast<ColGenEvalInfo *> (nodePtr->nodeEvalInfoPtr());

  bapcodInit().require(colGenEvalInfoPtr != NULL,
                       "BaPCod error: NodeEvalInfo for ColGenEval is not of type colGenEvalInfo.");

  /// we do not update basis if this node is treated immediately after its parent
  /// (so the basis in the LP solver remains unchanged)
  if ((_currentNodePtr->treatOrder() != colGenEvalInfoPtr->treatOrderId + 1)
      && (colGenEvalInfoPtr->masterLpBasisPtr != NULL))
    {
      /// we need to add to the basis the local branching constraints, as there were not in the problem
      /// when the basis was retrieved, TO DO : move basis reloading to the problem setup
      LpBasisRecord * basisPtrToReload = new LpBasisRecord(*(colGenEvalInfoPtr->masterLpBasisPtr));
      for (std::list<BranchingConstrBaseType *>::const_iterator
           constrIt = _currentNodePtr->localNodeBrConstrList().begin();
           constrIt != _currentNodePtr->localNodeBrConstrList().end(); ++constrIt)
        {
          if ((*constrIt)->isTypeOf(VcId::InstMasterBranchingConstrMask))
          {
            Constraint * constrPtr = dynamic_cast<Constraint *>(*constrIt);
            if (constrPtr != NULL)
              basisPtrToReload->_constrInBasis.push_back(
                  ConstrPtr_MpFormIndexStatus(constrPtr, MathProgSolverInterface::BasicVarConstr)
              );
          }
        }
      _probPtr->reloadMemorizedBasis(basisPtrToReload);
      delete basisPtrToReload;
    }

  if (colGenEvalInfoPtr->stabilizationInfoPtr != NULL)
    _keptStabInfoPtr = new StabilizationInfo(*(colGenEvalInfoPtr->stabilizationInfoPtr));
  _keptLatestRedCostFixingGap =  colGenEvalInfoPtr->latestReducedCostFixingGap;

  return false;
}

bool Alg4EvalByLp::eval()
{
  if (!progStatus().doRun())
    return false;

  if (solveRestrictedMastLP() <= 0)
    markInfeasible();
  else
    {
      _solIsMasterLpFeasible = true;

      /// for the moment, core cuts are not supported, TO DO : implement it
      /// if _updateIncDualBound = false, we are in the phase 1 of strong branching
      /// and checking solution feasibility can be too expensive
      if (_updateIncDualBound && checkIfCurSolIsInteger())
        updatePrimalIpSolAndBnds(_probPtr->inPrimalLpSol(), _probPtr->partialSolution());
    }
  
  if (printL(5))
    std::cout << "restrictedMasterLP :" << (_solIsMasterLpFeasible ? " " : " No ") << "LP solution found " << std::endl;

  return !_solIsMasterLpFeasible;
}

NodeEvalInfo * Alg4EvalByLp::recordNodeEvalInfo(int globalTreatOrder, NodeEvalInfo * nodeEvalInfoPtr)
{
  /// copy master solution, the copy will be used in branching
  _currentNodePtr->recordPrimalSol(_probPtr->inPrimalLpSol()); //Not needed for strong branching (Issam)

  LpBasisRecord * newLpBasisPtr = new LpBasisRecord(std::string("BasisN") + _currentNodePtr->ref());
    
  _probPtr->retrieveBasis(newLpBasisPtr);
    
  ColGenEvalInfo * colGenSolverInfoPtr = NULL;
  if (nodeEvalInfoPtr != NULL)
    {
      colGenSolverInfoPtr = dynamic_cast<ColGenEvalInfo *>(nodeEvalInfoPtr);

      bapcodInit().require(colGenSolverInfoPtr != NULL,
			   "BaPCod error: nodeEvalInfoPtr passed to LpEvalAlg::recordNodeEvalInfo is not of type ColGenEvalInfo");

      colGenSolverInfoPtr->masterLpBasisPtr = newLpBasisPtr;
      colGenSolverInfoPtr->stabilizationInfoPtr = _keptStabInfoPtr;
      colGenSolverInfoPtr->latestReducedCostFixingGap = _keptLatestRedCostFixingGap;
    }
  else
    colGenSolverInfoPtr = new ColGenEvalInfo(_keptStabInfoPtr , newLpBasisPtr, _keptLatestRedCostFixingGap);
  _keptStabInfoPtr = NULL;

  return Alg4EvalOfNode::recordNodeEvalInfo(globalTreatOrder, colGenSolverInfoPtr);
}

void Alg4EvalByLp::setDownAlgo()
{
  if (_keptStabInfoPtr != NULL)
    {
      delete _keptStabInfoPtr;
      _keptStabInfoPtr = NULL;
    }
  
  Alg4EvalOfNode::setDownAlgo();
}

void Alg4EvalByLp::setOptionUpdateIncDualBound(const bool & value)
{
  _updateIncDualBound = value;
}

