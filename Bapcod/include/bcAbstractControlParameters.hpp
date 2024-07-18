/**
 *
 * This file bcAbstractControlParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCABSTRACTCONTROLPARAMETERS_HPP_
#define BCABSTRACTCONTROLPARAMETERS_HPP_

#include <iostream>

class ParameterManager;

class AbstractControlParameters
{
public:
  AbstractControlParameters()
  {
  }
  virtual ~AbstractControlParameters()
  {
  }

  /**
   * Add all parameters to parameterManager to be parsed.
   * @param parameterManager
   */
  virtual void addParameters(ParameterManager& parameterManager) = 0;
  virtual std::ostream& printParameters(std::ostream& os) const = 0;
  
  /**
   * This method is used for post treatment after the parameter are parsed.
   *
   */
  virtual void postTreatment() = 0;
};

#endif /* BCABSTRACTPARAMETERS_HPP_ */
