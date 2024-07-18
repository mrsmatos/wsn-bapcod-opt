/**
 *
 * This file bcTimeC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef TIMECLASS_H
#define TIMECLASS_H


#include "bcUsefulHeadFil.hpp"
#include <time.h>
#include <ctime>
#include <string>
#include <sstream>

#ifndef _MSC_VER
/// Only available on Unix
#include <sys/times.h>
#include <sys/time.h>
#else
#include <windows.h> // by Artur
#endif //_MSC_VER

#include <string.h>

const long TICKFREQU = 1e2L;

#include <boost/timer/timer.hpp>

class Time
{
  boost::timer::cpu_timer timer; // this is both wallclock time AND cpu time

public:
  Time() {};

  /// Returns the ellapsed time since start in ticks (more specifically in centiseconds)
  /// Wall time or CPU time depending on the parameter.
  long getElapsedTime() const;
  
  double getElapsedTime_dbl() const;

  bool exceeds(const int & day, const int & month, const int & year);
};

inline std::ostream & printTime(long ticks, std::ostream & os)
{
  long t = ticks;
  long sec = t / TICKFREQU;
  t %= TICKFREQU;
  long min = sec / 60;
  sec %= 60;
  long hour = min / 60;
  min %= 60;
  return os << "TIME = " << hour << "h" << min << "m" << sec << "s" << t << "t = " << ticks << std::endl;
}

#endif