/**
 *
 * This file bcMathProgSolverFactory.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMATHPROGSOLVERFACTORY_H_
#define BCMATHPROGSOLVERFACTORY_H_

#include "bcMathProgSolverInterfaceC.hpp"

/**
 * This class allow to create an instance of a MathProgSolverInit,
 * LpMathProgSolverInterface and a MipMathProgSolverInterface.
 *
 *****************************************************************
 * DEVELOPPER
 *
 * If a new solver is implemented, you need to implement the three methods of this class in ../src/bcMathProgSolverFactory.hpp
 *****************************************************************
 */
template<typename T>
class MathProgSolverFactory
{
public:
  static T* buildMathProgSolverInit(BapcodInit* bapcodInit)
  {
    return new T(bapcodInit);
  }

  static T* buildLpMathProgSolverInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name)
  {
    return new T(bapcodInit, ref, name);
  }

  static T* buildMipMathProgSolverInterface(BapcodInit* bapcodInit, const int & ref, const std::string & name)
  {
    return new T(bapcodInit, ref, name);
  }
};

#endif /* BCMATHPROGSOLVERFACTORY_H_ */
