/**
 *
 * This file bcModelColGenSpC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGCOLGENSPC_H_
#define BCMODELINGCOLGENSPC_H_

#include "bcModelConstrC.hpp"
#include "bcModelFormulationC.hpp"

class Model;

class BcColGenSpArray: public BcFormulationArray
{
  Double _fixedCost;
  Double _defaultUb;
  Double _defaultLb;
 public:
  BcColGenSpArray(BcModel & modelPointer, const std::string & name = "colGenSp");

  BcFormulation & getElement(const MultiIndex & indexArray) override;
  BcFormulation & createElement(const MultiIndex & indexArray) override;
  const BcColGenSpArray & operator<=(const double & defaultUb);
  const BcColGenSpArray & operator>=(const double & defaultLb);
  void setFixedCost(const Double & fixedCost);

  virtual ~BcColGenSpArray(){}

};



#endif /* BCMODELINGCOLGENSPC_H_ */
