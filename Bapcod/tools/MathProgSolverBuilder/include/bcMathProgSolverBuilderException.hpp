/**
 *
 * This file bcMathProgSolverBuilderException.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMATHPROGSOLVERBUILDEREXCEPTION_H_
#define BCMATHPROGSOLVERBUILDEREXCEPTION_H_

#include "bcGlobalException.hpp"

/**
 * Exception throw when a solver is not found.
 */
class MathProgSolverBuilderException : public GlobalException
{
public:

  MathProgSolverBuilderException(std::string message="", bool printingForced = false, std::ostream & os = std::cerr);

  virtual ~MathProgSolverBuilderException() throw ();
};

std::ostream& operator<<(std::ostream& os, const MathProgSolverBuilderException& exception);



#endif /* BCMATHPROGSOLVERBUILDEREXCEPTION_H_ */
