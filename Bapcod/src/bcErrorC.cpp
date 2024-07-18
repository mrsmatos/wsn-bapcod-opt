/**
 *
 * This file bcErrorC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcDoubleC.hpp"
#include "bcErrorC.hpp"

#include "bcPrintC.hpp"
using namespace std;

ProgStatus::ProgStatus(Time* startTime, ProgStatistics* statistics):
    _status(run), _startTime(startTime), _statistics(statistics), _timeLimit(2147483645)
    {}

const ProgStatus::Type &  ProgStatus::setStat(const ProgStatus::Type & stat)
{
  if (printL(1))
    std::cout << " ProgStatus::setStat " <<  stat << std::endl;

  if (stat == quit)
    {
      if (printL(1)) 
	    _statistics->print();

      if (printL(1))
        print();

#ifdef ADLIB
      return(stat);
#else
      std::cerr << " Program exited prematurely " << std::endl;
      exit(1);
#endif
    }
  
  /// Should have  (_status != pStatus::exit) at this stage
  if (_status == run) 
    return(_status = stat);
  else  
    return(_status);
}

void  ProgStatus::resetStat()
{
  _status = run;
  _msg = std::string();

  return;
}
const ProgStatus::Type &  ProgStatus::stat()
{
  return(_status);
}

const std::string & ProgStatus::msg() const
{
  return(_msg);
}

bool ProgStatus::doRun()
{
  if (_status != run) return false;

  long elapsedTime = _startTime->getElapsedTime_dbl();

  if (printL(1))
    std::cout << "ProgStatus:: elapsedTime =  " << elapsedTime << " <? time limit = " << _timeLimit << std::endl;

  if (elapsedTime > _timeLimit)
    {
      setStat(ProgStatus::terminate);
      if (printL(-1))
        std::cout << "SEARCH IS INTERRUPTED as the time limit is reached. "  << std::endl;
      return false;
    }

  return true;
}

void ProgStatus::pushMsg(const std::string & message)
{
  std::ostringstream oss;
  oss << message << std::endl;
  _msg.append(oss.str());

  return;
}


std::ostream& ProgStatus::print(std::ostream & os) const
{
  os << "Program Status = " << _status << std::endl;
  if (_status == run) 
    os << "Program exited normaly. Messages = "  << _msg << std::endl;
  else 
    os << "Program exited because of " << _msg << std::endl;
  
  return(os);
}

