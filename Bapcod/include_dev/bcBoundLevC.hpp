/**
 *
 * This file bcBoundLevC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BOUNDCLASSES_H
#define BOUNDCLASSES_H
#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcModelParameterC.hpp"

class ControlParameters;

struct BcObjStatus
{
    enum MinMaxIntFloat
    {
        minInt = 1,
        maxInt = -1,
        minFloat = 2, /// default
        maxFloat = -2
    };
};

class Bound : public Double
{
public:
  BcObjStatus::MinMaxIntFloat objStatus;
  static const Bound staticOne;
  static const Bound staticZero;
  static const double epsilonBound;
  Bound(BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(0), objStatus(objStatus_)
  {
  }

  virtual ~Bound()
  {
    return;
  }

  Bound(const Double & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  Bound(const long double & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  Bound(const double & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  Bound(const float & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  Bound(const int & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  Bound(const long & i, BcObjStatus::MinMaxIntFloat objStatus_) :
      Double(i), objStatus(objStatus_)
  {
  }

  const Bound round()
  {
      bool integerValuedBound = (objStatus == BcObjStatus::maxInt || objStatus == BcObjStatus::minInt);
      bool minimizationProblem = (objStatus == BcObjStatus::minInt || objStatus == BcObjStatus::minFloat);
      if (printL(2))
          std::cout << "Bound round() integerValuedBound = " << integerValuedBound
                    << " minimizationProblem  = " << minimizationProblem << std::endl;

      if (integerValuedBound)
      {
          if (minimizationProblem)
          {
              if (printL(5))
                  printf("Bound round(): val = %.10f, Dceil(_val) = %.10f\n",_val, Dceil(_val));
              _val = Dceil(_val);
          }
          else
          {
              _val = Dfloor(_val);
          }
      }
      return (*this);
  }

  static const Bound infDualBound(BcObjStatus::MinMaxIntFloat objStatus)
  {
      bool minimizationProblem = (objStatus == BcObjStatus::minInt || objStatus == BcObjStatus::minFloat);
      if (minimizationProblem)
          return Bound(-BapcodInfinity, objStatus);
      else
          return Bound(BapcodInfinity, objStatus);
  }

  static const Bound infPrimalBound(BcObjStatus::MinMaxIntFloat objStatus)
  {
      bool minimizationProblem = (objStatus == BcObjStatus::minInt || objStatus == BcObjStatus::minFloat);
      if (minimizationProblem)
          return Bound(BapcodInfinity, objStatus);
      else
          return Bound(-BapcodInfinity, objStatus);
  }

  bool negative(const double & epsilon = epsilonBound) const
  {
      bool minimizationProblem = (objStatus == BcObjStatus::minInt || objStatus == BcObjStatus::minFloat);
      if (minimizationProblem)
          return Double::negative(epsilon);
      else
          return Double::positive(epsilon);
  }

  bool positive(const double & epsilon = epsilonBound) const
  {
      bool minimizationProblem = (objStatus == BcObjStatus::minInt || objStatus == BcObjStatus::minFloat);
      if (minimizationProblem)
          return Double::positive(epsilon);
      else
          return Double::negative(epsilon);
  }
};

inline bool operator<(const Bound & a, const Bound & b)
{
  bool minimizationProblem = (a.objStatus == BcObjStatus::minInt || a.objStatus == BcObjStatus::minFloat);
  if (minimizationProblem)
    {
      return ((Double) a < (Double) b);
    }
  else
    {
      return ((Double) b < (Double) a);
    }
}

bool gapSmallerThanTol(const Bound & dualBound, const Bound & primalBound, const ControlParameters & param);

inline bool operator>(const Bound & a, const Bound & b)
{
  return b < a;
}

inline bool operator>=(const Bound & a, const Bound & b)
{
  return !(a < b);
}

inline bool operator<=(const Bound & a, const Bound & b)
{
  return !(b < a);
}

inline bool operator==(const Bound & a, const Bound & b)
{
  if (a < b)
    return false;
  if (b < a)
    return false;
  return true;
}

inline bool operator!=(const Bound & a, const Bound & b)
{
  return !(a == b);
}

inline bool operator!=(const Bound & x, const Double & y)
{
  return !(operator==(x, Bound(y, x.objStatus)));
}

inline bool operator!=(const Double & x, const Bound & y)
{
  return !(operator==(Bound(x, y.objStatus), y));
}

inline bool operator<(const Bound & x, const Double & y)
{
  return (operator<(x, Bound(y, x.objStatus)));
}

inline bool operator<(const Double & x, const Bound & y)
{
  return (operator<(Bound(x, y.objStatus), y));
}

inline bool operator>(const Bound & x, const Double & y)
{
  return (operator<(Bound(y, x.objStatus), x));
}

inline bool operator>(const Double & x, const Bound & y)
{
  return (operator<(y, Bound(x, y.objStatus)));
}

inline bool operator<=(const Bound & x, const Double & y)
{
  return !(operator<(Bound(y, x.objStatus), x));
}

inline bool operator<=(const Double & x, const Bound & y)
{
  return !(operator<(y, Bound(x, y.objStatus)));
}

inline bool operator>=(const Bound & x, const Double & y)
{
  return !(operator<(x, Bound(y, x.objStatus)));
}

inline bool operator>=(const Double & x, const Bound & y)
{
  return !(operator<(Bound(x, y.objStatus), y));
}

inline bool operator!=(const Bound & x, const double & y)
{
  return !(operator==(x, Bound(y, x.objStatus)));
}

inline bool operator!=(const double & x, const Bound & y)
{
  return !(operator==(Bound(x, y.objStatus), y));
}

inline bool operator<(const Bound & x, const double & y)
{
  return (operator<(x, Bound(y, x.objStatus)));
}

inline bool operator<(const double & x, const Bound & y)
{
  return (operator<(Bound(x, y.objStatus), y));
}

inline bool operator>(const Bound & x, const double & y)
{
  return (operator<(Bound(y, x.objStatus), x));
}

inline bool operator>(const double & x, const Bound& y)
{
  return (operator<(y, Bound(x, y.objStatus)));
}

inline bool operator<=(const Bound & x, const double & y)
{
  return !(operator<(Bound(y, x.objStatus), x));
}

inline bool operator<=(const double & x, const Bound & y)
{
  return !(operator<(y, Bound(x, y.objStatus)));
}

inline bool operator>=(const Bound & x, const double & y)
{
  return !(operator<(x, Bound(y, x.objStatus)));
}

inline bool operator>=(const double & x, const Bound & y)
{
  return !(operator<(Bound(x, y.objStatus), y));
}

const Double computeOptimalityGap(const Double & primalBd, const Double & dualBd);

#endif							 // BOUNDCLASSES_H
