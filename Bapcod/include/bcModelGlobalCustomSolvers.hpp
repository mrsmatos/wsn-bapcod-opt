/**
 *
 * This file bcModelGlobalCustomSolvers.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELGLOBALCUSTOMSOLVER_H_
#define BCMODELGLOBALCUSTOMSOLVER_H_

#include "bcModelConstrC.hpp"
#include "bcModelNetworkFlow.hpp"
#include "bcModelFormulationC.hpp"
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include "bcTimeC.hpp"

class NetworkVertex;
class NetworkSet;
class ScalableResource;

class BcShortestPathOracleFunctor : public BcSolverOracleFunctor
{
  NetworkFlow* _netPtr;
  int _mode; /// 0 for Belman-Ford (Default), 1 for Dijkstra, 2 for DirectedAcyclicGraphs shortest path

public:

  BcShortestPathOracleFunctor(const BcNetwork & network, int mode = 0);

  #define IN_ // Indicates what is used for input
  #define OUT_ // Indicates what is used for output

  bool operator() (IN_ BcFormulation spPtr,
                   IN_ int colGenPhase,
                   OUT_ double & objVal,
                   OUT_ double & dualBound,
                   OUT_ BcSolution & primalSol) override;
};

#endif //BCMODELGLOBALCUSTOMSOLVER_H_

