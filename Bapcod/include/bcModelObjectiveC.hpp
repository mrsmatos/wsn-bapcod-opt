/**
 *
 * This file bcModelObjectiveC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGOBJECTIVEC_H_
#define BCMODELINGOBJECTIVEC_H_

#include <string>
#include <list>
#include "bcMultiIndexC.hpp"
#include "bcModelVarC.hpp"


class Model;
class BcModel;
class BcObjective;
class BcVar;
struct BcVarIndex;
struct BcVarCoef;
class ProbConfig;
class BcFormulation;

class BcObjective
{
 public:
  friend class BcVar;

 protected:
  Model * _modelPtr;

 public:
  BcObjective();
  BcObjective(BcModel & model);
  BcObjective(BcFormulation & bcForm);
  void setMinMaxStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense); /// for backwards compatibility
  void setStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense);
  void setArtCostValue(const double & cost);
  const BcObjective & addCumulatively(const BcRowExpression & exp);
  const BcObjective & add(const BcVar & var, const double & coef);
  const BcObjective & operator=(const BcVar & var);
  const BcObjective & operator+(const BcVar & var);
  const BcObjective & operator-(const BcVar & var);
  const BcObjective & operator+=(const BcVar & var);
  const BcObjective & operator-=(const BcVar & var);
  const BcObjective & operator=(const BcRowExpression & exp);
  const BcObjective & operator+(const BcRowExpression & exp);
  const BcObjective & operator+=(const BcRowExpression & exp);
  const BcObjective & operator-(const BcRowExpression & exp);
  const BcObjective & operator-=(const BcRowExpression & exp);
  const BcObjective & operator=(const BcVarCoef & coef);
  const BcObjective & operator+(const BcVarCoef & coef);
  const BcObjective & operator-(const BcVarCoef & coef);
  const BcObjective & operator+=(const BcVarCoef & coef);
  const BcObjective & operator-=(const BcVarCoef & coef);
  const BcObjective & operator=(BcVarIndex & modVarInd);
  const BcObjective & operator+(BcVarIndex & modVarInd);
  const BcObjective & operator-(BcVarIndex & modVarInd);
  const BcObjective & operator+=(BcVarIndex & modVarInd);
  const BcObjective & operator-=(BcVarIndex & modVarInd);
  /// sets an upper bound on optimal solution
  const BcObjective & operator<=(const double & ub);
  /// sets a lower bound on optimal solution
  const BcObjective & operator>=(const double & lb);
  /// sets an intial gess value estimating the order of magnitude of the optimal solution value,
  /// that is used for artificial variable setting
  const BcObjective & operator==(const double & guessVal);
};

/// this class is for backward compatibility
class BcObjectiveArray
{
 protected:
  Model * _modelPtr;

  BcObjective _curObjectivePtr;

 public:
  BcObjectiveArray(BcModel & model, const std::string & name = "objective");
  void setMinMaxStatus(const BcObjStatus::MinMaxIntFloat & newObjectiveSense);

  virtual BcObjective & createElement(const MultiIndex & indexArray);
  inline BcObjective & createElement(int firstIndex = 0, int secondIndex = -1, int thirdIndex = -1,
                                     int fourthIndex = -1, int fifthIndex = -1, int sixthIndex = -1,
                                     int seventhIndex = -1)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }
  inline BcObjective & operator()(int firstIndex = 0, int secondIndex = -1, int thirdIndex = -1, int fourthIndex = -1,
                                  int fifthIndex = -1, int sixthIndex = -1, int seventhIndex = -1)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }
};



#endif /* BCMODELINGOBJECTIVEC_H_ */
