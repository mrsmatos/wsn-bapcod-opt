/**
 *
 * This file bcTimeC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcTimeC.hpp"
#include <stdlib.h>

#include "bcPrintC.hpp"

long Time::getElapsedTime() const
{
  if (printL(5))
    std::cout << "Time::getElapsedTime(): " << timer.format();

  if (BapcodClockType == 0)
    return timer.elapsed().wall / (1e9L / TICKFREQU);
  else
    return (timer.elapsed().user + timer.elapsed().system) / (1e9L / TICKFREQU) ;
}

double Time::getElapsedTime_dbl() const
{
  if (printL(5))
    std::cout << "Time::getElapsedTime(): " << timer.format();

  if (BapcodClockType == 0)
    return timer.elapsed().wall / (1.0 * (1e9 / TICKFREQU));
  else
    return (timer.elapsed().user + timer.elapsed().system) / (1.0 * (1e9 / TICKFREQU)) ;
}

bool Time::exceeds(const int & day, const int & month, const int & year)
{
  std::time_t t = std::time(NULL);

  std::tm * threshold = localtime(&t);

  /// Seconds [0-60] (1 leap second)
  threshold->tm_sec = 0;

  /// Minutes [0-59]
  threshold->tm_min = 0;

  /// Hours [0-23]
  threshold->tm_hour = 0;

  /// Day [1-31]
  threshold->tm_mday = day;

  /// Month [0-11]
  threshold->tm_mon = month - 1;

  /// Year - 1900
  threshold->tm_year= year - 1900;

  return (t > mktime (threshold));
}

