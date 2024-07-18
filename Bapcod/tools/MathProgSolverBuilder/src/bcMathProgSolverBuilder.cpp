/**
 *
 * This file bcMathProgSolverBuilder.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMathProgSolverBuilder.hpp"
#include "bcMathProgSolverFactory.hpp"
#include "bcUserControlParameters.hpp"
class BapcodInit;

#if _CPLEX_FOUND
#include "bcCplexSolverC.hpp"
#endif

#if _XPRESSMP_FOUND
#include "bcXpressMPSolverC.hpp"
#endif

#if _GLPK_FOUND
#include "bcGlpkSolverC.hpp"
#endif

#if _CLP_FOUND
#include "bcClpSolverC.hpp"
#endif

#if _GUROBI_FOUND
#include "bcGRBSolverC.hpp"
#endif

using namespace std;

enum SolverTypeEnum
{
  CPLEX = 1, GLPK = 2, KNAPSACK = 3, XPRESSMP = 4, XPRESS = 5, CLP = 6, GRB = 7,
};

/**
 * Build the map containing the type of the solver.
 */
static std::map<std::string, SolverTypeEnum> buildSolverTypeMap()
{
  std::map<std::string, SolverTypeEnum> map;
  map[CPLEX_SOLVER] = CPLEX;
  map[XPRESSMP_SOLVER] = XPRESSMP;
  map[GLPK_SOLVER] = GLPK;
  map[CLP_SOLVER] = CLP;
  map[GUROBI_SOLVER] = GRB;
  return map;
}
static std::map<std::string, SolverTypeEnum> solverTypeMap = buildSolverTypeMap();

MathProgSolverBuilder::MathProgSolverBuilder()
{
  //Nothing to do
}

MathProgSolverBuilder::~MathProgSolverBuilder()
{
  //Nothing to do
}

MathProgSolverInterface * MathProgSolverBuilder
                          ::buildLpMathProgSolverInterface(BapcodInit* bapcodInit, std::string solverName,
                                                           const int & ref, const std::string & name)
{
  int solverType = solverTypeMap[solverName];

  MathProgSolverInterface* retSolverInterface = nullptr;

  switch (solverType)
    {
#if _CPLEX_FOUND
  case CPLEX:
  retSolverInterface = MathProgSolverFactory<LpCplexInterface>::buildLpMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_CPLEX_FOUND
#if _XPRESSMP_FOUND
  case XPRESSMP:
  retSolverInterface = MathProgSolverFactory<LpXpressMPInterface>::buildLpMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_XPRESSMP_FOUND
#if _GLPK_FOUND
  case GLPK:
  retSolverInterface = MathProgSolverFactory<LpGlpkInterface>::buildLpMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif
#if _CLP_FOUND
  case CLP:
  retSolverInterface = MathProgSolverFactory<LpClpInterface>::buildLpMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_CPLEX_FOUND
#if _GUROBI_FOUND
  case GRB:
  retSolverInterface = MathProgSolverFactory<LpGRBInterface>::buildLpMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif
          /// should make this compatible with C++11
//  default:
//    retSolverInterface = 0;
//    string message("the solver ");
//    message.append(solverName).append(" is not found");
//    throw MathProgSolverBuilderException(message);
    }

    return retSolverInterface;
}

MathProgSolverInterface * MathProgSolverBuilder
                          ::buildMipMathProgSolverInterface(BapcodInit* bapcodInit, std::string solverName,
                                                            const int & ref, const std::string & name)
{
  int solverType = solverTypeMap[solverName];

  MathProgSolverInterface* retSolverInterface = nullptr;

  switch (solverType)
    {
#if _CPLEX_FOUND
  case CPLEX:
  retSolverInterface = MathProgSolverFactory<MipCplexInterface>::buildMipMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_CPLEX_FOUND
#if _XPRESSMP_FOUND
  case XPRESSMP:
  retSolverInterface = MathProgSolverFactory<MipXpressMPInterface>::buildMipMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_XPRESSMP_FOUND
#if _GLPK_FOUND
  case GLPK:
  retSolverInterface = MathProgSolverFactory<MipGlpkInterface>::buildMipMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif
#if _CLP_FOUND
  case CLP:
  retSolverInterface = MathProgSolverFactory<MipClpInterface>::buildMipMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif //_CPLEX_FOUND
#if _GUROBI_FOUND
  case GRB:
  retSolverInterface = MathProgSolverFactory<MipGRBInterface>::buildMipMathProgSolverInterface(bapcodInit, ref, name);
  break;
#endif
  /// should make this compatible with C++11
//  default:
//    retSolverInterface = 0;
//    string message("the solver ");
//    message.append(solverName).append(" is not found");
//    throw MathProgSolverBuilderException(message);
    }
  return retSolverInterface;
}
