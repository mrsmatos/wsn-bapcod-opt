/**
 *
 * This file bcMathProgSolverBuilderException.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMathProgSolverBuilderException.hpp"

MathProgSolverBuilderException::MathProgSolverBuilderException(
    std::string message, bool printingForced, std::ostream & os) :
    GlobalException(message, printingForced, os)
{
}

MathProgSolverBuilderException::~MathProgSolverBuilderException() throw ()
{
}

std::ostream & operator<<(std::ostream & os, const MathProgSolverBuilderException & exception)
{
  os << "MathProgSolverBuilderException: " << (GlobalException) exception;
  return os;
}

