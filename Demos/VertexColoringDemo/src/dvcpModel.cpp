/**
 *
 * This file dvcpModel.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dvcpModel.hpp"

void registerInitialSolution(BcModel & model, const VCPsolution & initialSolution)
{
  BcMaster master(model);
  BcColGenSpArray colGenSp(model);
  BcFormulation spForm = colGenSp[0];

  BcVarArray xVar(spForm, "X");
  BcVarArray yVar(spForm, "Y");

  BcSolution prevSol(NULL);
  for (VCPsolution::const_iterator vectIt = initialSolution.begin(); vectIt != initialSolution.end(); ++vectIt)
    {
      BcSolution newSol(spForm);
      for (std::set<int>::iterator setIt = vectIt->begin(); setIt != vectIt->end(); ++setIt)
        {
          xVar[*setIt] = 1;
          newSol += xVar[*setIt];
        }
      yVar[0] = 1;
      newSol += yVar[0];
      if (prevSol.solutionPtr() != NULL)
        newSol += prevSol;
      prevSol = newSol;
    }
  master.initializeWithSolution(prevSol);

  return;
}


void buildVCPModel(DvcpData & data, const VCPsolution & initSolution, BcModel & model)
{
    BcObjective objective(model);
    objective.setStatus(BcObjStatus::minInt);

    BcMaster master(model);

    BcColGenSpArray colGenSp(model);

    BcVarArray yVar(colGenSp[0], "Y");
    yVar.priorityForMasterBranching(2);
    yVar.priorityForSubproblemBranching(-1);
    yVar.localLb(1);
    yVar.type('I');
    objective += yVar(0);

    BcVarArray xVar(colGenSp[0], "X");
    xVar.defineIndexNames(MultiIndexNames('i'));
    xVar.type('B');
    xVar.priorityForMasterBranching(-1);
    xVar.priorityForSubproblemBranching(-1);
    xVar.priorityForRyanFosterBranching(1);

    BcConstrArray coverConstr(master, "COV");
    coverConstr.defineIndexNames(MultiIndexNames('i'));
    coverConstr.considerAsEqualityInPreprocessing(true);
    coverConstr >= 1;

    for (int vertexIndex = 0; vertexIndex < data.numV; ++vertexIndex)
    {
      coverConstr(vertexIndex) += xVar(vertexIndex);
    }

    BcConstrArray colorUsageConstr(colGenSp[0], "CLR");
    colorUsageConstr.defineIndexNames(MultiIndexNames('i'));
    colorUsageConstr <= 0;
    for (int vertexIndex = 0; vertexIndex < data.numV; ++vertexIndex)
    {
        colorUsageConstr(vertexIndex) += xVar[vertexIndex];
        colorUsageConstr[vertexIndex] -= yVar[0];
    }

    BcConstrArray colorConflictConstr(colGenSp[0], "EDG");
    colorConflictConstr.defineIndexNames(MultiIndexNames('i', 'j'));

    for (int edgeIndex = 0; edgeIndex < data.numE; ++edgeIndex)
    {
        int tail = data.T_e[edgeIndex];
        int head = data.H_e[edgeIndex];
        colorConflictConstr(tail, head) <= 1;
        colorConflictConstr[tail][head] += xVar[tail];
        colorConflictConstr[tail][head] += xVar[head];
    }

    if (initSolution.empty())
    {
        objective.setArtCostValue(data.numV);
    } else
    {
        objective.setArtCostValue(initSolution.size());
        registerInitialSolution(model, initSolution);
    }

}

