/**
 *
 * This file bcErrorC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef ERRORC_H
#define ERRORC_H
#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcTimeC.hpp"

class Time;
class ProgStatistics;
class Time;
class ProgStatus
{
 public:
  enum Type 
  {
    quit = 1,
    run = 0,
    terminate = 2
  };
  
 protected:
  Type _status;
  std::string _msg;
  Time* _startTime;
  ProgStatistics* _statistics;
  long _timeLimit;

 public:
  ProgStatus(Time* startTime, ProgStatistics* statistics);
  const Type &  stat();
  const std::string & msg() const;
  bool doRun();
  Time * startTime() const {return  _startTime;}
  void startTime(Time * stPtr)  {_startTime = stPtr;}
  void timeLimit(long tl) {_timeLimit = tl;}
  const Type &  setStat(const Type & stat);
  void  resetStat();
  void pushMsg(const std::string & message);
  std::ostream & print(std::ostream & os = std::cout) const;
};

inline std::ostream & operator<<(std::ostream& os, const ProgStatus & a) 
{
  return a.print(os);
}

inline void assure(std::ifstream& in, const char* filename = "", std::ostream& os = std::cout)
{
  if(!in)
    {
      os << "Could not open file " << filename << std::endl;
      exit(1);
    }
}

inline void assure(std::ofstream& in, const char* filename = "", std::ostream& os = std::cout)
{
  if(!in)
    {
      os << "Could not open file " << filename << std::endl;
      exit(1);
    }
}

#endif // ERRORC_H
