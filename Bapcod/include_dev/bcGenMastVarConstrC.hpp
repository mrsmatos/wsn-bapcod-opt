/**
 *
 * This file bcGenMastVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef GenMastVarConstrClasses_h
#define GenMastVarConstrClasses_h

#include "bcVcIdentifierC.hpp"


#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcGenVarConstrC.hpp"

/**
 * A generic Master Constraints implementing 
 * lower and upper bound on the number of 
 * subproblme solutions that can be selected
 */
class ConvexityGenConstr: public GenericConstr
{
 public:
  ConvexityGenConstr(Model * modelPtr,  MasterConf * masterPtr, const std::vector<ColGenSpConf *> & colGenSubProbConfPts);
  virtual ~ConvexityGenConstr();
  /* virtual const Double & genericCostRhs(const InstanciatedVarConstr * const ivarconstrPtr) const; */
  virtual bool genericCount(const InstanciatedConstr * const iconstrPtr, 
			    const InstanciatedVar * const ivarPtr) const;
  
  virtual const LpCoef genericCoef(const InstanciatedConstr * const iconstrPtr, 
				   const InstanciatedVar * const ivarPtr) const;
  
  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

#endif // GenVarConstrClasses_h


