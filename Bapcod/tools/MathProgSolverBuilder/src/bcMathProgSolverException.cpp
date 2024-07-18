/**
 *
 * This file bcMathProgSolverException.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcMathProgSolverException.hpp"

MathProgSolverException::MathProgSolverException(std::string message,
    bool printingForced, std::ostream & os) :
    GlobalException(message, printingForced, os)
{
}

MathProgSolverException::~MathProgSolverException() throw ()
{
  // Nothing to do
}

std::ostream & operator <<(std::ostream & os,
    const MathProgSolverException & exception)
{
  os << "MathProgSolverException: " << (GlobalException) exception;
  return os;
}

