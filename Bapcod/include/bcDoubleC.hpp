/**
 *
 * This file bcDoubleC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCDOUBLEC_H
#define BCDOUBLEC_H

#include "bcUsefulHeadFil.hpp"

inline int intMin(const int & v1, const int & v2)
{
  if (v1 < v2)
    return v1;

  return v2;
}

/// Min needs not be numerically precise
inline const double & Dmin(const double & a, const double & b)
{
  if (a < b)
    return a;

  return b;
}

/// Max needs not be numerically precise
inline const double & Dmax(const double & a, const double & b)
{
  if (b < a)
    return a;

  return b;
}

inline double Dfloor(const double & x);
inline int Ifloor(const double & x);
inline double Dceil(const double & x);
inline int Iceil(const double & x);
inline double Dfrac(const double & x);
inline double Dfrac(const double & x, const double & target);
inline bool Dnegative(const double & val, const double & epsilon);
inline bool Dpositive(const double & val, const double & epsilon);
inline double zero(const double & that, const double & epsilon);
inline double zero(const double & that, const double & epsilon);

class Double
{
public:
  static const Double staticZero;
  static const Double staticOne;
  static const double precision;
  static const double highprecision;
  double _val;
  const double & val() const
  {
    return _val;
  }
  double & val()
  {
    return _val;
  }

  /// Operator conversion:
  operator double() const
  {
    return _val;
  }

  operator int() const
  {
    return Ifloor(_val);
  }

  operator long() const
  {
    return Ifloor(_val);
  }

  Double() :
      _val(0)
  {
  }

  Double(const Double & x) :
      _val(x)
  {
  }
  Double(const double & x) :
      _val(x)
  {
  }

  bool negative(const double & epsilon = Double::precision) const
  {
    return (_val < -epsilon);
  }

  bool positive(const double & epsilon = Double::precision) const
  {
    return (_val > epsilon);
  }

  bool isZero(const double & epsilon = Double::precision) const
  {
    if (positive(epsilon))
      return false;
    return !negative(epsilon);
  }

  void Cfloor()
  {
    double xminus(_val - 1);
    _val = floor(_val + Double::highprecision * _val + Double::precision);
    if (_val < xminus)
      _val += 1;
    Czero();
  }

  void Cceil()
  {
    double xplus(_val + 1);
    _val = ceil(
        (double) (_val - Double::highprecision * _val - Double::precision));
    if (_val >= xplus)
      _val -= 1;
    Czero();
  }

  void Czero(const double & epsilon = Double::precision)
  {
    if (!positive(epsilon) && !negative(epsilon))
      _val = 0;
  }

  bool fractional(const double & epsilon = Double::precision) const
  {
    return (Dpositive(Dfrac(_val), epsilon));
  }

  const Double operator-() const
  {
    return -_val;
  }
};

inline double Fabs(const Double & x)
{
  return fabs(x._val);
}

inline bool Dnegative(const double & val, const double & epsilon =
    Double::precision)
{
  return (val < -epsilon);
}

inline bool Dpositive(const double & val, const double & epsilon =
    Double::precision)
{
  return (val > epsilon);
}

inline bool zeroTest(const double & that, const double & epsilon =
    Double::precision)
{
  if (Dpositive(that, epsilon))
    return false;

  return !Dnegative(that, epsilon);
}

inline double zero(const double & that, const double & epsilon = Double::precision)
{
  if (Dpositive(that, epsilon))
    return that;

  if (Dnegative(that, epsilon))
    return that;
  else
    return 0;
}

inline bool smallerThan(const Double & a, const Double & b,
                        const double epsilon = Double::precision,
                        const double relEpsilon = Double::highprecision)
{
  double tolerance((std::max)(fabs(a._val), fabs(b._val)) * relEpsilon + epsilon);

  return (a._val < (b._val - tolerance));
}

inline bool operator<(const Double & a, const Double & b)
{
  return smallerThan(a, b);
}

inline bool operator==(const Double & a, const Double & b)
{
  double tolerance((std::max)(fabs(a._val), fabs(b._val)) * Double::precision + Double::precision);
  return ((a._val >= (b._val - tolerance)) && (b._val >= (a._val - tolerance)));
}

inline std::ostream& operator<<(std::ostream& os, const Double & a)
{
  return os << a._val;
}

inline const std::istream& operator>>(const std::istream& is, Double& a)
{
  double x = 0.0;
  const_cast<std::istream *>(&is)->operator >>(x);
  a._val = x;

  return is;
}

inline double Dfloor(const double & x)
{
  double round(floor(x + Double::highprecision * x + Double::precision));
  double xminus(x - 1 + Double::highprecision * x + Double::precision);
  if (round < xminus)
    round += 1;

  return zero(round);
}

inline int Ifloor(const double & x)
{
  int round((int) floor(x + Double::highprecision * x + Double::precision));
  if (round < (x - 1 + Double::precision))
    round += 1;

  return round;
}

inline int Iceil(const double & x)
{
  int round((int) ceil(x - Double::highprecision * x - Double::precision));
  if (round >= x + 1 + Double::precision)
    round -= 1;

  return round;
}

inline double Dceil(const double & x)
{
  double round(ceil((double) (x - Double::highprecision * x - Double::precision)));
  double xplus(x + 1);
  if (round >= xplus)
    round -= 1;
  return zero(round);
}

inline double Lfrac(const double & x)
{
  return zero(x - Dfloor(x));
}

inline double Ufrac(const double & x)
{
  return zero(Dceil(x) - x);
}

/// Dfrac(x) returns the min  between (x - Dfloor(x)) and (Dceil(x) - x) 
/// It is maximum at value 0.5 for x = (Dfloor(x) + Dceil(x))/2. 
inline double Dfrac(const double & x)
{
  return Dmin(Lfrac(x), Ufrac(x));
}

/// Dfrac(x,target) returns the min between the normalized distance target for any target in [0,1]. 
/// It is maximum at value 0.5 for frac = (x - Dfloor(x)) = target = 1 - (Dceil(x) -x)
inline double Dfrac(const double & x, const double & target)
{
  const double & frac = Lfrac(x);
  if (frac > target)
    return 0.5 * (1 - ((frac - target) / (1 - target)));
  /// else (l <= target)
  return 0.5 * (1 - ((target - frac) / target));
}

/// For other comparators

inline bool operator!=(const Double & x, const Double & y)
{
  return !(operator==(x, y));
}

inline bool operator<(const Double & x, const double & y)
{
  return (operator<(x, (Double) y));
}
inline bool operator<(const double & x, const Double & y)
{
  return (operator<((Double) x, y));
}

inline bool operator==(const Double & x, const double & y)
{
  return (operator==(x, (Double) y));
}
inline bool operator==(const double & x, const Double & y)
{
  return (operator==((Double) x, y));
}

inline bool operator!=(const Double & x, const double & y)
{
  return !(operator==(x, (Double) y));
}

inline bool operator!=(const double & x, const Double & y)
{
  return !(operator==((Double) x, y));
}

inline bool operator>(const Double & x, const Double & y)
{
  return (operator<(y, x));
}

inline bool operator>(const Double & x, const double & y)
{
  return (operator<((Double) y, x));
}

inline bool operator>(const double & x, const Double& y)
{
  return (operator<(y, (Double) x));
}

inline bool operator<=(const Double & x, const Double & y)
{
  return !(operator<(y, x));
}
inline bool operator<=(const Double & x, const double & y)
{
  return !(operator<((Double) y, x));
}

inline bool operator<=(const double & x, const Double & y)
{
  return !(operator<(y, (Double) x));
}

inline bool operator>=(const Double & x, const Double & y)
{
  return !(operator<(x, y));
}

inline bool operator>=(const Double & x, const double & y)
{
  return !(operator<(x, (Double) y));
}

inline bool operator>=(const double & x, const Double & y)
{
  return !(operator<((Double) x, y));
}

/// Other operators

inline const Double operator+(const Double & x, const Double & y)
{
  return (x._val + y._val);
}

inline const Double operator+(const Double & x, const double & y)
{
  return (x._val + y);
}

inline const Double operator+(const double & x, const Double & y)
{
  return (x + y._val);
}

inline const Double operator-(const Double & x, const Double & y)
{
  return (x._val - y._val);
}

inline const Double operator-(const Double & x, const double & y)
{
  return (x._val - y);
}

inline const Double operator-(const double & x, const Double & y)
{
  return (x - y._val);
}

inline const Double operator*(const Double & x, const Double & y)
{
  return (x._val * y._val);
}

inline const Double operator*(const Double & x, const double & y)
{
  return (x._val * y);
}

inline const Double operator*(const double & x, const Double & y)
{
  return (x * y._val);
}

inline const Double operator/(const Double & x, const Double & y)
{
  return (x._val / y._val);
}

inline const Double operator/(const Double & x, const double & y)
{
  return (x._val / y);
}

inline const Double operator/(const double & x, const Double & y)
{
  return (x / y._val);
}

inline const Double operator+=(Double& x, const Double & y)
{
  return ((x._val) += y._val);
}

inline const Double operator+=(Double& x, const double & y)
{
  return ((x._val) += y);
}

inline const Double operator+=(double & x, const Double& y)
{
  return (x += y._val);
}

inline const Double operator-=(Double& x, const Double & y)
{
  return ((x._val) -= y._val);
}

inline const Double operator-=(Double& x, const double & y)
{
  return ((x._val) -= y);
}

inline const Double operator-=(double & x, const Double& y)
{
  return (x -= y._val);
}

inline const Double operator*=(Double& x, const Double & y)
{
  return ((x._val) *= y._val);
}

inline const Double operator*=(Double& x, const double & y)
{
  return ((x._val) *= y);
}

inline const Double operator*=(double & x, const Double& y)
{
  return (x *= y._val);
}

inline const Double operator/=(Double& x, const Double & y)
{
  return ((x._val) /= y._val);
}

inline const Double operator/=(Double& x, const double & y)
{
  return ((x._val) /= y);
}

inline const Double operator/=(double & x, const Double& y)
{
  return (x /= y._val);
}

/// For other comparators float

inline bool operator<(const Double & x, const float & y)
{
  return (operator<(x, (Double) y));
}

inline bool operator<(const float & x, const Double & y)
{
  return (operator<((Double) x, y));
}

inline bool operator==(const Double & x, const float & y)
{
  return (operator==(x, (Double) y));
}

inline bool operator==(const float & x, const Double & y)
{
  return (operator==((Double) x, y));
}

inline bool operator!=(const Double & x, const float & y)
{
  return !(operator==(x, (Double) y));
}

inline bool operator!=(const float & x, const Double & y)
{
  return !(operator==((Double) x, y));
}

inline bool operator>(const Double & x, const float & y)
{
  return (operator<((Double) y, x));
}

inline bool operator>(const float & x, const Double& y)
{
  return (operator<(y, (Double) x));
}

inline bool operator<=(const Double & x, const float & y)
{
  return !(operator<((Double) y, x));
}

inline bool operator<=(const float & x, const Double & y)
{
  return !(operator<(y, (Double) x));
}

inline bool operator>=(const Double & x, const float & y)
{
  return !(operator<(x, (Double) y));
}

inline bool operator>=(const float & x, const Double & y)
{
  return !(operator<((Double) x, y));
}

/// Other operators float

inline const Double operator+(const Double & x, const float & y)
{
  return (x._val + y);
}

inline const Double operator+(const float & x, const Double & y)
{
  return (x + y._val);
}

inline const Double operator-(const Double & x, const float & y)
{
  return (x._val - y);
}

inline const Double operator-(const float & x, const Double & y)
{
  return (x - y._val);
}

inline const Double operator*(const Double & x, const float & y)
{
  return (x._val * y);
}

inline const Double operator*(const float & x, const Double & y)
{
  return (x * y._val);
}

inline const Double operator/(const Double & x, const float & y)
{
  return (x._val / y);
}

inline const Double operator/(const float & x, const Double & y)
{
  return (x / y._val);
}

inline const Double operator+=(Double & x, const float & y)
{
  return ((x._val) += y);
}

/**
 * Force the conversion of y._val into float,
 * with loss of precision,
 * then before cast it again into a Double.
 */
inline const Double operator+=(float & x, const Double & y)
{
    double xx = x;
  return (xx += y._val);
}

inline const Double operator-=(Double & x, const float & y)
{
  return ((x._val) -= y);
}

/**
 * Force the conversion of y._val into float,
 * with loss of precision,
 * then before cast it again into a Double.
 */
inline const Double operator-=(float & x, const Double & y)
{
    double xx = x;
  return (xx -= y._val);
}

inline const Double operator*=(Double & x, const float & y)
{
  return ((x._val) *= y);
}

/**
 * Force the conversion of y._val into float,
 * with loss of precision,
 * then before cast it again into a Double.
 */
inline const Double operator*=(float & x, const Double & y)
{
  double xx = x;
  return (xx *= y._val);
}

inline const Double operator/=(Double & x, const float & y)
{
  return ((x._val) /= y);
}

/**
 * Force the conversion of y._val into float,
 * with loss of precision,
 * then before cast it again into a Double.
 */
inline const Double operator/=(float & x, const Double & y)
{
  double xx = x;
  return (xx /= y._val);
}

/// For other comparators int

inline bool operator<(const Double & x, const int & y)
{
  return (operator<(x, (Double) y));
}

inline bool operator<(const int & x, const Double & y)
{
  return (operator<((Double) x, y));
}

inline bool operator==(const Double & x, const int & y)
{
  return (operator==(x, (Double) y));
}

inline bool operator==(const int & x, const Double & y)
{
  return (operator==((Double) x, y));
}

inline bool operator!=(const Double & x, const int & y)
{
  return !(operator==(x, (Double) y));
}

inline bool operator!=(const int & x, const Double & y)
{
  return !(operator==((Double) x, y));
}

inline bool operator>(const Double & x, const int & y)
{
  return (operator<((Double) y, x));
}

inline bool operator>(const int & x, const Double& y)
{
  return (operator<(y, (Double) x));
}

inline bool operator<=(const Double & x, const int & y)
{
  return !(operator<((Double) y, x));
}

inline bool operator<=(const int & x, const Double & y)
{
  return !(operator<(y, (Double) x));
}

inline bool operator>=(const Double & x, const int & y)
{
  return !(operator<(x, (Double) y));
}

inline bool operator>=(const int & x, const Double & y)
{
  return !(operator<((Double) x, y));
}

/// Other operators int

inline const Double operator+(const Double & x, const int & y)
{
  return (x._val + y);
}

inline const Double operator+(const int & x, const Double & y)
{
  return (x + y._val);
}

inline const Double operator-(const Double & x, const int & y)
{
  return (x._val - y);
}

inline const Double operator-(const int & x, const Double & y)
{
  return (x - y._val);
}

inline const Double operator*(const Double & x, const int & y)
{
  return (x._val * y);
}

inline const Double operator*(const int & x, const Double & y)
{
  return (x * y._val);
}

inline const Double operator/(const Double & x, const int & y)
{
  return (x._val / y);
}

inline const Double operator/(const int & x, const Double & y)
{
  return (x / y._val);
}

inline const Double operator+=(Double & x, const int & y)
{
  return ((x._val) += y);
}

inline int operator+=(int & x, const Double & y)
{
  return (x += (int) y);
}

inline const Double operator-=(Double & x, const int & y)
{
  return ((x._val) -= y);
}

inline int operator-=(int & x, const Double & y)
{
  return (x -= (int) y);
}

inline const Double operator*=(Double & x, const int & y)
{
  return ((x._val) *= y);
}

inline int operator*=(int & x, const Double & y)
{
  return (x = (int) (x * y));
}

inline const Double operator/=(Double & x, const int & y)
{
  return ((x._val) /= y);
}

inline int operator/=(int & x, const Double & y)
{
  return (x = (int) (x / y));
}

/// For other comparators long

inline bool operator<(const Double & x, const long & y)
{
  return (operator<(x, (Double) y));
}

inline bool operator<(const long & x, const Double & y)
{
  return (operator<((Double) x, y));
}

inline bool operator==(const Double & x, const long & y)
{
  return (operator==(x, (Double) y));
}

inline bool operator==(const long & x, const Double & y)
{
  return (operator==((Double) x, y));
}

inline bool operator!=(const Double & x, const long & y)
{
  return !(operator==(x, (Double) y));
}

inline bool operator!=(const long & x, const Double & y)
{
  return !(operator==((Double) x, y));
}

inline bool operator>(const Double & x, const long & y)
{
  return (operator<((Double) y, x));
}

inline bool operator>(const long & x, const Double& y)
{
  return (operator<(y, (Double) x));
}

inline bool operator<=(const Double & x, const long & y)
{
  return !(operator<((Double) y, x));
}

inline bool operator<=(const long & x, const Double & y)
{
  return !(operator<(y, (Double) x));
}

inline bool operator>=(const Double & x, const long & y)
{
  return !(operator<(x, (Double) y));
}

inline bool operator>=(const long & x, const Double & y)
{
  return !(operator<((Double) x, y));
}

inline const Double operator+(const Double & x, const long & y)
{
  return (x._val + y);
}

inline const Double operator+(const long & x, const Double & y)
{
  return (x + y._val);
}

inline const Double operator-(const Double & x, const long & y)
{
  return (x._val - y);
}

inline const Double operator-(const long & x, const Double & y)
{
  return (x - y._val);
}

inline const Double operator*(const Double & x, const long & y)
{
  return (x._val * y);
}

inline const Double operator*(const long & x, const Double & y)
{
  return (x * y._val);
}

inline const Double operator/(const Double & x, const long & y)
{
  return (x._val / y);
}

inline const Double operator/(const long & x, const Double & y)
{
  return (x / y._val);
}

inline const Double operator+=(Double & x, const long & y)
{
  return ((x._val) += y);
}

inline long operator+=(long & x, const Double & y)
{
  return (x += (long) y);
}

inline const Double operator-=(Double & x, const long & y)
{
  return ((x._val) -= y);
}

inline long operator-=(long & x, const Double & y)
{
  return (x -= (long) y);
}

inline const Double operator*=(Double & x, const long & y)
{
  return ((x._val) *= y);
}

inline long operator*=(long & x, const Double & y)
{
  return (x = (long) (x * y));
}

inline const Double operator/=(Double & x, const long & y)
{
  return ((x._val) /= y);
}

inline long operator/=(long & x, const Double & y)
{
  return (x = (long) (x / y));
}

struct LpCoef
{
  bool first; // takes value true if second is non-zero
  Double second;
public:
  static const LpCoef ZeroCoef;
  static const LpCoef UnitCoef;
  static const LpCoef MinusOneCoef;
  LpCoef(const bool & status, const Double & x) : first(status), second(x)  {}
  LpCoef(const LpCoef & a) : first(a.first), second(a.second) {}
  LpCoef(const Double & x) : first(false), second(0.0)
  {
      if (!zeroTest(x))
      {
          first = true;
          second = x;
      }
  }
};

inline const LpCoef operator+=(LpCoef & x, const LpCoef & y)
{
  if (!x.first)
    {
      x.second = y.second;
      return x;
    }
  if (!y.first)
   {
      return x;
    }
  x.second += y.second;
  return x;
}

inline std::ostream& operator<<(std::ostream& os, const LpCoef & a)
{
  return os << " LpCoef < " << a.first << ", " << " > " << a.second ;
}


#endif // DOUBLEC_H
