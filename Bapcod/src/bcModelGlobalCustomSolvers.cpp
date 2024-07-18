/**
 *
 * This file bcModelGlobalCustomSolvers.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcNetworkFlowC.hpp"
#include "bcModelGlobalCustomSolvers.hpp"
#include "bcProbConfigC.hpp"

#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <lemon/bellman_ford.h>


BcShortestPathOracleFunctor::BcShortestPathOracleFunctor(const BcNetwork & network, int mode) :
_netPtr((NetworkFlow *) network), _mode(mode)
{
}

#define IN_ // Indicates what is used for input
#define OUT_ // Indicates what is used for output

//details on the different parameters can be found on the parent class.  

bool BcShortestPathOracleFunctor::operator() (IN_ BcFormulation spPtr,
                                              IN_ int colGenPhase,
                                              OUT_ double & objVal,
                                              OUT_ double & dualBound,
                                              OUT_ BcSolution & primalSol)
{
  const lemon::ListDigraph & digraph = _netPtr->digraph();

  const lemon::ListDigraph::Node & s = _netPtr->sourceList().front();
  const lemon::ListDigraph::Node & t = _netPtr->sinkList().front();

  lemon::ListDigraph::ArcMap<double> reducedCosts(digraph);

  for (lemon::ListDigraph::ArcIt curArc(digraph); curArc != lemon::INVALID; ++curArc)
    {
      NetworkArc * netArcPtr = _netPtr->netArcPtr(curArc);
      reducedCosts[curArc] = 0;
      for (auto & pair : netArcPtr->varToCoeffMaps().front())
        reducedCosts[curArc] += pair.second * pair.first->curCost();
    }

  if (_mode != 2)
  {
    lemon::ListDigraph::Node curNode = t;
    lemon::ListDigraph::Node predNode;
    lemon::ListDigraph::Arc predArc;

    if (_mode == 0)
    {
      lemon::BellmanFord<lemon::ListDigraph, lemon::ListDigraph::ArcMap<double> >
              bf(digraph, reducedCosts);
      bf.run(s);

      do
      {
        predArc = bf.predArc(curNode);
        NetworkArc * netArcPtr = _netPtr->netArcPtr(predArc);
        for (auto & pair : netArcPtr->varToCoeffMaps().front())
            primalSol.updateVarVal(BcVar(pair.first),pair.second);
        curNode = digraph.source(predArc);
      }
      while (curNode != s);
      
      objVal = dualBound = bf.dist(t);
    }
    else //if _mode == 1
    {
      lemon::Dijkstra<lemon::ListDigraph, lemon::ListDigraph::ArcMap<double> >
              dijkstra(digraph, reducedCosts);
      dijkstra.run(s, t);

      do
      {
        predArc = dijkstra.predArc(curNode);
        NetworkArc * netArcPtr = _netPtr->netArcPtr(predArc);
        for (auto & pair : netArcPtr->varToCoeffMaps().front())
          primalSol.updateVarVal(BcVar(pair.first), pair.second);
        curNode = digraph.source(predArc);
      } while (curNode != s);
      
      objVal = dualBound = dijkstra.dist(t);

    }

  } else
  {
    std::cout << "BcShortestPathOracleFunctor::operator() ERR: DAG NOT SUPPORTED YET" << std::endl;
  }

  return true;
}