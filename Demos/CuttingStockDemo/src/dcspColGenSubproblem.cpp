/**
 *
 * This file dcspColGenSubproblem.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dcspColGenSubproblem.hpp"
#include <vector>

using std::vector;

IntKnapsackSubProb::IntKnapsackSubProb(DcspData * dataPtr) :
    _dataPtr(dataPtr), _firstCallInitializationDone(false)
{
}

bool IntKnapsackSubProb::operator()(BcFormulation spPtr,
                                    int colGenPhase,
                                    double & objVal,
                                    double & dualBound,
                                    BcSolution & primalSol)
{
  int nbOfOrder(_dataPtr->orderPtrVector().size());

  bool doPrint(false);

  if (!_firstCallInitializationDone)
    {
      _firstCallInitializationDone = true;

      _flagOrderInUseV = std::vector<bool>(nbOfOrder);
      _profitV = std::vector<double>(nbOfOrder);
      _weightV = std::vector<double>(nbOfOrder);
      _boundV = std::vector<int>(nbOfOrder);
      _maxciV = std::vector<int>(nbOfOrder);

      for (int orderIndex = 0; orderIndex < nbOfOrder; ++orderIndex)
        {
          _weightV[orderIndex] = _dataPtr->orderPtrVector()[orderIndex]->width();
          _profitV[orderIndex] = 0;
          _boundV[orderIndex] = 0;
          _maxciV[orderIndex] = 0;
          _flagOrderInUseV[orderIndex] = false;
        }
    }
  else
    {
      for (int orderIndex = 0; orderIndex < nbOfOrder; ++orderIndex)
        {
          _profitV[orderIndex] = 0;
          _boundV[orderIndex] = 0;
          _maxciV[orderIndex] = 0;
          _flagOrderInUseV[orderIndex] = false;
        }
    }

  _vitems.clear();
  double capacity = _dataPtr->stockSheet().width();

  _minweight = capacity;

  Double capacityAlreadyUsed(0);
  Double partialCost(0);
  bool nonTrivialSol(false);

  BcVarArray X(spPtr, "X");

  for (int orderIndex = 0; orderIndex < nbOfOrder; ++orderIndex)
    {

      if (X[orderIndex].curUb() <= 0)
        {
          continue;
        }

      if (X[orderIndex].inCurForm())
        {
          _profitV[orderIndex] = -X[orderIndex].curCost();
          _boundV[orderIndex] = X[orderIndex].curUb();
          _maxciV[orderIndex] = 0;
          _flagOrderInUseV[orderIndex] = true;

          if (X[orderIndex].curLb() >= 1)
            {
              nonTrivialSol = true;
              X[orderIndex] = X[orderIndex].curLb();
              _boundV[orderIndex] -= X[orderIndex].curLb();
              primalSol += X[orderIndex];
              capacityAlreadyUsed += _weightV[orderIndex] * X[orderIndex].curLb();
              partialCost += X[orderIndex].curCost() * X[orderIndex].curLb();

            }

          if (doPrint)
            {
              std::cout << " orderIndex " << orderIndex << " profit "
                  << _profitV[orderIndex] << " weight " << _weightV[orderIndex]
                  << " bound " << _boundV[orderIndex] << std::endl;
            }
        }
    }

  int nbPosProd(0);
  for (int orderIndex = 0; orderIndex < nbOfOrder; ++orderIndex)
    {
      if (_flagOrderInUseV[orderIndex] == false)
        continue;
      if (_boundV[orderIndex] <= 0)
        continue;
      if (_profitV[orderIndex] <= 0)
        continue;

      Double ratio = _profitV[orderIndex] / _weightV[orderIndex];
      if (_weightV[orderIndex] < _minweight)
        _minweight = _weightV[orderIndex];

      int b(0), mult(1);
      while (b + mult < _boundV[orderIndex])
        {
          Double weight = _weightV[orderIndex] * mult;

          if (capacityAlreadyUsed + weight > capacity)
            break;
          _vitems.push_back(KnpItem(NULL, nbPosProd++, _profitV[orderIndex] * mult, weight,
                                    mult, orderIndex, ratio, 0));
          if (doPrint)
            {
              std::cout << " _vitems " << _vitems[nbPosProd - 1];
            }
          b += mult;
          mult *= 2;
        }
      Double maxweight = _weightV[orderIndex] * mult;
      _maxciV[orderIndex] = (b + mult);

      if (capacityAlreadyUsed + maxweight > capacity)
        continue;

      _vitems.push_back(KnpItem(NULL, nbPosProd++, _profitV[orderIndex] * mult, maxweight,
                                mult, orderIndex, ratio, 0));
      if (doPrint)
        {
          std::cout << " last _vitems " << _vitems[nbPosProd - 1];
        }
    }

  double incumbent(0);
  Double cap(capacity - capacityAlreadyUsed);
  int countUBeval = 0;
  int countBBnode = 0;

  if (nbPosProd > 0)
    {
      stable_sort(_vitems.begin(), _vitems.end());

      if (colGenPhase == 0)
        {
          if (doPrint)
            {
              std::cout << "pricing by Customized Solver..." << std::endl;
            }
          // solving root sp problem
          MCbinknap(nbPosProd, nbOfOrder, _profitV, _weightV, _boundV, _maxciV,
                    _vitems, _minweight, cap, countUBeval, countBBnode, incumbent,
                    true);
        }
      else if (colGenPhase == 1)
        {
          /// apply greedy heuristic
          if (doPrint)
            {
              std::cout << "pricing by Greedy Heuristic..." << std::endl;
              std::cout << "Phase = " << colGenPhase << std::endl;
            }

          GreedyKnap(nbPosProd, nbOfOrder, _profitV, _weightV, _boundV, _vitems,
              _minweight, cap, incumbent);
        }
      else
      {
          return false;
      }
    }

  BcVarArray Y(spPtr, "Y");
  Y[0] = 1;
  primalSol += Y[0];
  partialCost += Y[0].curCost();

  objVal = dualBound = partialCost - incumbent;

  incumbent = partialCost;
  for (int binaryItemIndex = 0; binaryItemIndex < nbPosProd; ++binaryItemIndex)
    {
      if (_vitems[binaryItemIndex].x > 0)
        {
          nonTrivialSol = true;
          capacityAlreadyUsed += _vitems[binaryItemIndex].w * _vitems[binaryItemIndex].x;
          incumbent -= _vitems[binaryItemIndex].p * _vitems[binaryItemIndex].x;
          int orderIndex = _vitems[binaryItemIndex].icl;

          X[orderIndex] = _vitems[binaryItemIndex].x * _vitems[binaryItemIndex].m;

          if (doPrint)
            {
              std::cout << " binaryItemIndex " << binaryItemIndex
                  << " orderIndex " << orderIndex << " profit "
                  << _profitV[orderIndex] << " weight " << _weightV[orderIndex]
                  << " bound " << _boundV[orderIndex]
                  << " _vitems[binaryItemIndex].x = "
                  << _vitems[binaryItemIndex].x
                  << " _vitems[binaryItemIndex].m = "
                  << _vitems[binaryItemIndex].m << "  X[orderIndex] = "
                  << X[orderIndex].curVal() << std::endl;
            }

          primalSol += X[orderIndex];
        }
    }

  if (capacityAlreadyUsed > _dataPtr->stockSheet().width())
    {
      std::cout << "ERROR  IntKnapsackSubProb::customizedSolver() : capacity consumption exceeds knapsack rhs"
                << std::endl;
      exit(1);
    }

  if ((Double) objVal != (Double) incumbent)
    {
      std::cout << "!!ERROR IntKnapsackSubProb::customizedSolver() : objVal != incumbent "
                << objVal << " " << incumbent << std::endl;
      exit(1);
    }

  return (true);

}

bool GreedyKnap(const int & n, const int & ni, const std::vector<double> & p,
                const std::vector<double> & w, std::vector<int> & ci,
                std::vector<KnpItem> & vit, const double & minw, const double & cap, double &inc)
{

  int residualCapacity = cap;

  /// Build a greedy solution according to ratio
  for (int orderIndex = 0; orderIndex < n; orderIndex++)
    {
      if (residualCapacity < minw)
        {
          break;
        }
      if (vit[orderIndex].w <= residualCapacity)
        {
          vit[orderIndex].x = 1;
          residualCapacity -= vit[orderIndex].w;
          inc += vit[orderIndex].p;
          if (vit[orderIndex].m > 1)
            {
              orderIndex += vit[orderIndex].m;
            }
        }
    }

  return (true);

}
