/**
 *
 * This file bcMathProgSolverBuilder.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMATHPROGSOLVERBUILDER_H_
#define BCMATHPROGSOLVERBUILDER_H_

#include "bcMathProgSolverInterfaceC.hpp"
#include "bcMathProgSolverBuilderException.hpp"

/**
 * This class build a solver with the name.
 *
 ******************************************
 * DEVELOPPER
 * If a new solver is implemented, you need to change the three methods:
 * buildMathProgSolverInit, buildLpMathProgSolverInterface and
 * buildMipMathProgSolverInterface
 * You need to declare too a new std::string constant for the name of the
 * solver and a new enum associate with it.
 ******************************************
 */
class MathProgSolverBuilder
{
public:
  MathProgSolverBuilder();
  virtual ~MathProgSolverBuilder();

  /**
   * Build an init solver (CPLEX, XPressMP...) based on its name.
   *
   * \param the name of the solver. Use one of this fellow constant
   * CLPEX_SOLVER, GLPK_SOLVER, KNAPSACK_SOLVER, XPRESSMP_SOLVER,
   * XPRESS_SOLVER, COIN_SOLVER;
   * \return the solver.
   * \throw a MathProgSolverBuilderException if the solver did not implemented
   * or not found.
   */
//  virtual MathProgSolverInit* buildMathProgSolverInit(BapcodInit* bapcodInit, std::string solverName)
//      throw (const MathProgSolverBuilderException&);

  /**
   * Build an interface solver (LpCplexInterface...) based on its solver name,
   * the ref and name problem
   *
   * \param{solverName} the name of the solver. Use one of this fellow constant
   * CLPEX_SOLVER, GLPK_SOLVER, KNAPSACK_SOLVER, XPRESSMP_SOLVER,
   * XPRESS_SOLVER, COIN_SOLVER;
   * \param{ref} problem number
   * \param{name} problem name
   * \return the solver interface.
   * \throw a MathProgSolverBuilderException if the solver did not implemented
   * or not found.
   */
   virtual MathProgSolverInterface* buildLpMathProgSolverInterface(BapcodInit* bapcodInit,
      std::string solverName,
      const int & ref, const std::string & name);
          //throw (const MathProgSolverBuilderException&);

  /**
   * Build an interface solver (MipCplexInterface...) based on its solver name,
   * the ref and name problem
   *
   * \param{solverName} the name of the solver. Use one of this fellow constant
   * CLPEX_SOLVER, GLPK_SOLVER, KNAPSACK_SOLVER, XPRESSMP_SOLVER,
   * XPRESS_SOLVER, COIN_SOLVER;
   * \param{ref} problem number
   * \param{name} problem name
   * \return the solver interface.
   * \throw a MathProgSolverBuilderException if the solver did not implemented
   * or not found.
   */
  virtual MathProgSolverInterface* buildMipMathProgSolverInterface(BapcodInit* bapcodInit,
      std::string solverName,
      const int & ref, const std::string & name);
          //throw (const MathProgSolverBuilderException&);
};


#endif /* BCMATHPROGSOLVERBUILDER_H_ */
