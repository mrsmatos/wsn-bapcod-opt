/**
 *
 * This file bcAlg4Master.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcAlg4Master.hpp"
#include "bcMathProgSolverInterfaceC.hpp"
#include "bcModelC.hpp"

using namespace std;

Alg4Master::Alg4Master(Problem * const probPtr) :
  _solAndBnds(probPtr->objStatus()),
  _algCurLpPrimalBound(Bound::infPrimalBound(probPtr->objStatus())),
  _algCurLpDualBound(Bound::infDualBound(probPtr->objStatus())),
  _solIsMasterLpFeasible(false), _probPtr(probPtr), _currentNodePtr(NULL), _algCurLpPrimalSolMap(),
  _algCurLpDualSolMap(), _maxLevelOfSubProbRestriction(0)
{
}

bool Alg4Master::setupAlgo(Node * nodePtr)
{
  _currentNodePtr = nodePtr;

  _solAndBnds.algIncLpDualBound = _currentNodePtr->nodeIncLpDualBound();
  _solAndBnds.algIncStageLpDualBound = _currentNodePtr->nodeIncLpDualBound();
  _solAndBnds.algIncIpDualBound = _currentNodePtr->nodeIncIpDualBound();

  _solAndBnds.algIncLpPrimalBound = Bound::infPrimalBound(_probPtr->objStatus());
  _solAndBnds.algIncStageLpPrimalBound = Bound::infPrimalBound(_probPtr->objStatus());
  _solAndBnds.algIncIpPrimalBound = _currentNodePtr->nodeIncIpPrimalBound();

  _solAndBnds.algIncIpPrimalBoundUpdated = false;
  
  if (printL(5))
    std::cout << "Alg4Master::setupAlgo() Solver setup of node " << nodePtr->ref() << " is done!" << std::endl;

  if (printL(6) && (_probPtr!=NULL))
    _probPtr->printForm();

  return false;
}

void Alg4Master::setDownAlgo()
{
  if (printL(5))
    std::cout << "evalAlg set down of node " << _currentNodePtr->ref() << " is done!" << std::endl;

  /// this is to "release" columns in this map (to decrease their participation)
  clearAlgIncIpPrimalSolMap();
  _currentNodePtr = NULL;
}

bool Alg4Master::updateAlgDualBounds()
{
  Bound roundedDualBd(_algCurLpDualBound);
  if (param().SafeDualBoundScaleFactor() > 0)
    roundedDualBd._val = ceil(roundedDualBd._val);
  else
    roundedDualBd.round();

  Bound fractionalDualBd(_algCurLpDualBound);

  bool status(false);

  if (_solAndBnds.algIncStageLpDualBound <= fractionalDualBd)
  {
    _solAndBnds.algIncStageLpDualBound = fractionalDualBd;

    status = true;

    if ( _maxLevelOfSubProbRestriction <= 0 ) // if we are in a relaxation/exact
    {

      if (_solAndBnds.algIncLpDualBound <= fractionalDualBd)
      {
        _solAndBnds.algIncLpDualBound = fractionalDualBd;
      }

      if (_solAndBnds.algIncIpDualBound <= roundedDualBd)
      {
        _solAndBnds.algIncIpDualBound = roundedDualBd;
      }
    }
  }

  return status;
}

bool Alg4Master::updateAlgPrimalLpBounds()
{
  bool status(false);

  if (_solAndBnds.algIncStageLpPrimalBound > _algCurLpPrimalBound)
    {
      _solAndBnds.algIncStageLpPrimalBound = _algCurLpPrimalBound;

      status = true;

      if ((_maxLevelOfSubProbRestriction >= 0) && (_solAndBnds.algIncLpPrimalBound > _algCurLpPrimalBound))
        {
          _solAndBnds.algIncLpPrimalBound = _algCurLpPrimalBound;
        }
    }

  return status;
}

void Alg4Master::clearAlgIncIpPrimalSolMap()
{
  for (VarPtr2DoubleMap::const_iterator mapIt = _solAndBnds.algIncIpPrimalSolMap.begin();
       mapIt != _solAndBnds.algIncIpPrimalSolMap.end(); ++mapIt)
    mapIt->first->decrParticipation(12);
  _solAndBnds.algIncIpPrimalSolMap.clear();
}

/// this function is added by Ruslan to store the best solution directly from the current solution of Problem
/// if the idea is to store all the solutions in Alg4Master, then they should be removed from Problem,
/// otherwise we would store the same thing two times, one time in _algCurLpPrimalSolMap, second time in Problem,
/// and we would lose time copying very often the solution from Problem to _algCurLpPrimalSolMap
void Alg4Master::updatePrimalIpSolAndBnds(const VarPtrSet & primalSol, const VarPtr2DoubleMap & partialSol)
{
  if (_algCurLpPrimalBound >= algIncIpPrimalBound())
    return;

  _solAndBnds.algIncIpPrimalBoundUpdated = true;
  _solAndBnds.algIncIpPrimalBound = _algCurLpPrimalBound;

  if (printL(1))
    std::cout << "New eval. alg. incumbent solution with value " << _solAndBnds.algIncIpPrimalBound << std::endl;

  clearAlgIncIpPrimalSolMap();
  for (VarPtrSet::const_iterator varPtrIt = primalSol.begin(); varPtrIt != primalSol.end(); ++varPtrIt)
    {
      _solAndBnds.algIncIpPrimalSolMap.insert(std::make_pair(*varPtrIt, (*varPtrIt)->val()));
      if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
        (*varPtrIt)->incrParticipation(18);
    }

  for (VarPtr2DoubleMap::const_iterator partSolMapIt = partialSol.begin();
       partSolMapIt != partialSol.end(); ++partSolMapIt)
    {
      VarPtr2DoubleMap::iterator mapIt = _solAndBnds.algIncIpPrimalSolMap.find(partSolMapIt->first);
      if (mapIt == _solAndBnds.algIncIpPrimalSolMap.end())
        {
          _solAndBnds.algIncIpPrimalSolMap.insert(std::make_pair(partSolMapIt->first, partSolMapIt->second));
          if (partSolMapIt->first->isTypeOf(VcId::MastColumnMask))
            partSolMapIt->first->incrParticipation(18);
        }
      else
        {
          mapIt->second += partSolMapIt->second;
        }
    }
}

/// added by Ruslan : dual bound can be slightly more than master solution value
/// because of floating point calculation, we make the dual bound equal to lp value in this case
void Alg4Master::rectifyIncumbentLpValue()
{
    if (printL(1))
        std::cout << "ColGen bounds rectification called, LP bounds = [" << std::setprecision(12)
                  << algIncLpDualBound() << "," << algIncLpPrimalBound()<< "], IP bounds = ["
                  << algIncIpDualBound() << "," << algIncIpPrimalBound() << "]" << std::setprecision(6)
                  << std::endl;

    double diff = _solAndBnds.algIncLpDualBound - _solAndBnds.algIncLpPrimalBound;
    if (diff > 0)
    {
        if (!gapSmallerThanTol(_solAndBnds.algIncLpPrimalBound, _solAndBnds.algIncLpDualBound, param()))
        {
            if (printL(0))
                std::cout << "BaPCod WARNING : dual bound " << std::setprecision(10) << _solAndBnds.algIncLpDualBound
                          << " is greater than the best LP value " << _solAndBnds.algIncLpPrimalBound
                          << ", this may be caused by tolerance issues" << std::setprecision(6) << std::endl;
            std::cerr << "BaPCod WARNING : dual bound is greater than the best master LP value" << std::endl;
        }

        _solAndBnds.algIncLpPrimalBound = _solAndBnds.algIncLpDualBound;

        /// thin block is commented, as now we increase the best LP value instead of decreasing the best bound
        /// this is done in order to avoid the error
        /// "primal solution is integer after node evaluation but the node duality gap is non-zero"

//      _solAndBnds.algIncLpDualBound = Bound::infDualBound(_probPtr->objStatus());
//      _solAndBnds.algIncIpDualBound = Bound::infDualBound(_probPtr->objStatus());
//      _solAndBnds.algIncStageLpDualBound = Bound::infDualBound(_probPtr->objStatus());
//
//      _algCurLpDualBound = _solAndBnds.algIncLpPrimalBound;
//
//      /// we set _maxLevelOfSubProbRestriction to 0 so that _solAndBnds.algIncLpDualBound
//      /// and _solAndBnds.algIncIpDualBound could be updated in updateAlgDualBounds
//      /// at this point, we are already out of the "stage decreasing loop"
//      _maxLevelOfSubProbRestriction = 0;
//
//      updateAlgDualBounds();
    }
}


void Alg4Master::markInfeasible()
{
  _solIsMasterLpFeasible = false;

  _solAndBnds.algIncLpPrimalBound = _solAndBnds.algIncLpDualBound = _solAndBnds.algIncIpDualBound
                                  = Bound::infPrimalBound(_probPtr->objStatus());
    
  if (printL(1))
    std::cout << "Alg4Master :: EARLY TERMINATION of evaluation algorithm : infeasibility is detected" << std::endl;
}

void Alg4Master::backupBoundsExceptIncIpPrimalBound(Alg4MasterSolAndBounds & backupBounds)
{
  backupBounds.algIncLpPrimalBound = algIncLpPrimalBound();
  backupBounds.algIncLpDualBound = algIncLpDualBound();
  backupBounds.algIncIpDualBound = algIncIpDualBound();
}
  
void Alg4Master::restoreBoundsExceptIncIpPrimalBound(const Alg4MasterSolAndBounds & backupBounds)
{
  _solAndBnds.algIncLpPrimalBound = backupBounds.algIncLpPrimalBound;
  _solAndBnds.algIncLpDualBound = backupBounds.algIncLpDualBound;
  _solAndBnds.algIncIpDualBound = backupBounds.algIncIpDualBound;
}

