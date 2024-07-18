/**
 *
 * This file bcModelMasterC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELMASTERC_H_
#define BCMODELMASTERC_H_

#include "bcMultiIndexC.hpp"
#include "bcModelFormulationC.hpp"

class BcModel;
class BcMasterArray;
class MasterConf;
class Model;

class BcMaster : public BcFormulation
{
public:
    BcMaster(BcModel & model, const std::string & name = "master");
    BcMaster(BcFormulation & bcForm); /// for backward compatibility
};

/// this class is for backwards compatibility
class BcMasterArray: public BcFormulationArray
{
 public:
  BcMasterArray(BcModel & model, const std::string & name = "master");
  virtual BcFormulation & getElement(const MultiIndex & indexArray);
  virtual BcFormulation & createElement(const MultiIndex & indexArray);
};


#endif /* BCMODELMASTERC_H_ */
