/**
 *
 * This file bcControlParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCCONTROLPARAMETERS_HPP_
#define BCCONTROLPARAMETERS_HPP_

#include "bcUserControlParameters.hpp"
#include "bcDevControlParameters.hpp"

class ControlParameters : public UserControlParameters,
    public DevControlParameters
{
public:
  ControlParameters();
  virtual ~ControlParameters();

  void addParameters(ParameterManager& parameterManager);
  std::ostream& printDevParameters(std::ostream& os) const;
  std::ostream& printUserParameters(std::ostream& os) const;
  std::ostream& printVRPSolverParameters(std::ostream& os) const;
  void postTreatment();
};

#endif /* BCCONTROLPARAMETERS_HPP_ */
