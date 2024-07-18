/**
 *
 * This file bcMathProgSolverException.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMATHPROGSOLVEREXCEPTION_H_
#define BCMATHPROGSOLVEREXCEPTION_H_

#include "bcGlobalException.hpp"

/**
 * Exception throw when a solver is not found.
 */
class MathProgSolverException : public GlobalException
{
public:
  MathProgSolverException(std::string message="", bool printingForced = false, std::ostream & os = std::cerr);
  virtual ~MathProgSolverException() throw();
};

std::ostream& operator<<(std::ostream& os, const MathProgSolverException& exception);

#endif /* BCMATHPROGSOLVEREXCEPTION_H_ */
