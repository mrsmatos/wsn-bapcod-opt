/**
 *
 * This file bcOvfConfC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCOVFCONFC_H_
#define BCOVFCONFC_H_

#include "bcProbConfigC.hpp"

class OvfConf: public ProbConfig
{
public:
  /**
   *
   * @param problemClassPtr Simulate virtual constructor of problem
   */
  OvfConf(Model * modelPtr, Problem * problemClassPtr);

  virtual ~OvfConf();

  virtual Solution * solvePC();
  virtual Solution * getDissagregatedSolution(Solution * solPtr);
  virtual Constraint * castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately = false);
  virtual InstanciatedConstr * castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately = false);
  virtual Variable * castAndAddVariable(Variable * varPtr, const bool & insertImmediately = false);
  virtual InstanciatedVar * castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately = false);

  /**
   * Set initial variables and constraints
   */
  virtual void prepareProbConfig();

  virtual std::ostream & print(std::ostream& os = std::cout) const;
};


#endif /* BCOVFCONFC_H_ */
