/**
 *
 * This file bcStatisticC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcPrintC.hpp"
#include "bcTimeC.hpp"

ProgStatistics::Record::Record(const std::string & key, 
			                   const long  & counter,
			                   const double  & time,
			                   const Double & value)
  : _key(key), 
    _counter(counter), 
    _time(time), 
    _value(value)
{
  return;
}

ProgStatistics::Record::Record(const Record & record)
  : _key(record._key),
    _counter(record._counter),
    _time(record._time),
    _value(record._value)
{
  return;
}

void ProgStatistics::Record::operator+=(const Record & that)
{
  if (that._counter != ProgStatisticsUndefinedStatus)
    {
      if (_counter  == ProgStatisticsUndefinedStatus) 
	    _counter = that._counter;
      else 
	    _counter += that._counter;
    }
  
  if (that._time != ProgStatisticsUndefinedStatus)
    {
      if (_time  == ProgStatisticsUndefinedStatus) 
	    _time = that._time;
      else 
	    _time += that._time;
    }
  
  if (that._value != ProgStatisticsUndefinedStatus)
    {
      if (_value  == ProgStatisticsUndefinedStatus) 
	    _value = that._value;
      else 
	    _value += that._value;
    }
  
  return;
}

void ProgStatistics::Record::operator/=(const int & den)
{
  if (_counter != ProgStatisticsUndefinedStatus)
    _counter = (long) Dfloor((double) _counter / (double) den) ;
  
  if (_time != ProgStatisticsUndefinedStatus) 
    _time  = (long) Dfloor((double) _time / (double) den) ;
  
  if (_value != ProgStatisticsUndefinedStatus) 
    _value /= (double) den ;
  
  return;
}
std::ostream& ProgStatistics::Record::print(std::ostream& os) const
{
  if (_counter != ProgStatisticsUndefinedStatus)
    os << _counter << " " ;
  
  if (_time != ProgStatisticsUndefinedStatus)
    {
      os << _time << " & "; /// added by Ruslan
      long t = _time;
      long sec = t / TICKFREQU;
      t %= TICKFREQU;
      long min = sec / 60;
      sec %= 60;
      long hour = min / 60;
      min %= 60;
      os <<  hour << "h" << min << "m" << sec << "s" << t << "t ";
    }
  
  if (_value != ProgStatisticsUndefinedStatus) 
    os << std::setprecision(12) << _value << std::setprecision(6) << " " ;
  
  return(os);
}
std::ostream& ProgStatistics::Record::plainPrint(std::ostream& os) const
{
  if (_counter != ProgStatisticsUndefinedStatus)
    os << _counter  ;
  
  if (_time != ProgStatisticsUndefinedStatus) 
    os <<  _time ;
  
  if (_value != ProgStatisticsUndefinedStatus) 
    os << std::setprecision(12) << _value << std::setprecision(6)  ;
  
  return(os);
}

std::ostream& ProgStatistics::Record::printAverage(std::ostream& os) const
{
  if (_counter != ProgStatisticsUndefinedStatus && _value != ProgStatisticsUndefinedStatus)
    os <<  (_value/_counter);
  
  return(os);
}

ProgStatistics::ProgStatistics(const std::vector< std::string> & selectedKeys): 
  _selectedKeys(selectedKeys)
{
  return;
}

ProgStatistics::ProgStatistics(const std::vector<ProgStatistics> & prStatVector)
{
  if (prStatVector.empty())
    return;
  
  _selectedKeys = prStatVector[0].selectedKeys();
  _prName = "average";
  
  for (std::vector< std::string>::const_iterator it =  _selectedKeys.begin(); it != _selectedKeys.end(); ++it)
    {
      Record aggregateRec(*it);
      int nbEntries(0);
      for (std::vector<ProgStatistics>::const_iterator psIt =  prStatVector.begin(); psIt != prStatVector.end(); psIt++)
          if (psIt->_recordMultiSet.count(*it))
          {
              aggregateRec += *(psIt->_recordMultiSet.lower_bound(*it));
              nbEntries++;
          }
      aggregateRec /= nbEntries;
      _recordMultiSet.insert(aggregateRec);
    }
  {
    std::string key("bcTimeMain");
    Record aggregateRec(key);
    int nbEntries(0);
    for (std::vector<ProgStatistics>::const_iterator psIt =  prStatVector.begin(); psIt != prStatVector.end(); ++psIt)
      if (psIt->_recordMultiSet.count(key))
      {
          aggregateRec += *(psIt->_recordMultiSet.lower_bound(key));
          nbEntries++;
      }
    aggregateRec /= nbEntries;
    _recordMultiSet.insert(aggregateRec);
  }
  
  return;
}

ProgStatistics::ProgStatistics()
{
}

ProgStatistics::~ProgStatistics()
{
}

std::ostream & ProgStatistics::printStat(std::ostream& os)
{
    for (std::vector< std::string>::const_iterator it =  _selectedKeys.begin(); it != _selectedKeys.end(); ++it)
    {
        if (_recordMultiSet.count(*it))
        {
            std::multiset<Record, smallerKey >::const_iterator msit = _recordMultiSet.lower_bound(*it);
            os << msit->_key;
            os << "=";
            msit->plainPrint(os);
            os << std::endl;
        }
        else
        {
            os << *it << "=0" << std::endl;
        }
    }
    std::multiset<Record, smallerKey >::const_iterator msit = _recordMultiSet.lower_bound(std::string("bcTimeMain"));
    os << msit->_key;
    os << "=";
    msit->plainPrint(os);
    os << std::endl;

    return(os);
}

std::ostream& ProgStatistics::print(std::ostream& os)
{
  os << _prName <<std::endl;
  for (std::multiset<Record, smallerKey >::const_iterator it =  _recordMultiSet.begin();
       it != _recordMultiSet.end(); ++it)
    {
      os << it->_key;
      os << " ";
      it->print(os);
      os << std::endl;
    }
  
  return(os);
}

std::ostream& ProgStatistics::selectPrint(std::ostream& os, const bool printEndLine)
{
  os << "BaPStat: "; // added a string for easy "grep"
  os << _prName;
  for (std::vector< std::string>::const_iterator it =  _selectedKeys.begin(); it != _selectedKeys.end(); ++it)
    {
      os <<  " & ";
      if (_recordMultiSet.count(*it)) 
	    _recordMultiSet.lower_bound(*it)->plainPrint(os);
      else  
	    os <<  0;
    }
    os <<  " & ";
    _recordMultiSet.lower_bound( std::string("bcTimeMain"))->print(os);
  
    if (printEndLine)
    {
      os <<  " \\\\ ";
      os << std::endl;
    }
  
  return(os);
}

std::ostream& ProgStatistics::titlePrint(std::ostream& os)
{
  os << std::endl;
  
  os << "prName";
  for (std::vector< std::string>::const_iterator it =  _selectedKeys.begin(); it != _selectedKeys.end(); ++it)
    {
      os  << " & " <<  *it ;
    }
  
  os <<  " & ";
  os << "bcTimeMain";
  os <<  " \\\\ " << std::endl;

  return os;
}

void ProgStatistics::replaceRecord(const Record  & oldRec, const Record  & newRec)
{
  _recordMultiSet.erase(oldRec);
  _recordMultiSet.insert(newRec);
  
  return;
}

void ProgStatistics::record(const std::string & key, const long & counter, const double & time, const Double & value)
{
  _recordMultiSet.insert(Record(key, counter, time, value));
  
  return;
}

void ProgStatistics::incrRecord(const std::string  & key, const Double & val)
{
  Record rec(key, 1, ProgStatisticsUndefinedStatus, val);
  if (_recordMultiSet.count( rec))
    {
      Record newRec(*_recordMultiSet.lower_bound(rec));
      newRec += rec;
      replaceRecord(*_recordMultiSet.lower_bound(rec), newRec);
    }
  else 
    _recordMultiSet.insert(rec);
  
  return;
}

void ProgStatistics::recCounter(const std::string  & key, const long  & counter)
{
  _recordMultiSet.insert(Record(key, counter));
  return;
}

void ProgStatistics::recTime(const std::string  & key, const double  & time)
{
  _recordMultiSet.insert(Record(key, ProgStatisticsUndefinedStatus, time, ProgStatisticsUndefinedStatus));
  return;
}

void ProgStatistics::recValue(const std::string  & key, const Double & value)
{
  _recordMultiSet.insert(Record(key,ProgStatisticsUndefinedStatus, ProgStatisticsUndefinedStatus, value));

  return;
}

void ProgStatistics::setCounter(const std::string  & key, const long  & counter)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    {
      std::multiset<Record, smallerKey >::iterator it = _recordMultiSet.lower_bound(rkey);
      replaceRecord(*it, Record(key,counter,it->_time,it->_value));
    }
  else 
    recCounter(key, counter);
  
  return;
}

void ProgStatistics::incrValue(const std::string  & key,const Double & value)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    {
      std::multiset<Record, smallerKey >::iterator it = _recordMultiSet.lower_bound(rkey);
      replaceRecord(*it, Record(key,it->_counter,it->_time,value));
    }
  else 
    recValue(key, value);
  
  return;
}

void ProgStatistics::incrCounter(const std::string  & key, const long  & incr)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    {
      std::multiset<Record, smallerKey >::iterator it = _recordMultiSet.lower_bound(rkey);
      replaceRecord(*it, Record(key,it->_counter + incr,it->_time,it->_value));
    }
  else 
    recCounter(key, incr);
  
  return;
}
void ProgStatistics::incrTimer(const std::string  & key, const double  & timer)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    {
      std::multiset<Record, smallerKey >::iterator it = _recordMultiSet.lower_bound(rkey);
      replaceRecord(*it, Record(key,it->_counter,it->_time + timer,it->_value));
    }
  else 
    recTime(key, timer);
 
  return;
}

long ProgStatistics::getCounter(const std::string  & key)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey)) 
    return( _recordMultiSet.lower_bound(rkey)->_counter);
  else 
    return(0);
}

double ProgStatistics::getValue(const std::string & key)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    return( _recordMultiSet.lower_bound(rkey)->_value);
  else
    return(0);
}

double ProgStatistics::getTime(const std::string & key)
{
  Record rkey(key);
  if (_recordMultiSet.count(rkey))
    return( _recordMultiSet.lower_bound(rkey)->_time);
  else
    return(0);
}

bool ProgStatistics::Record::operator<(const ProgStatistics::Record & that) const
{
  if (_key < that._key)
    return(true);
  
  if (_key > that._key) 
    return(false);

  if (_counter < that._counter) 
    return(true);

  if (_counter > that._counter) 
    return(false);

  if (_time < that._time) 
    return(true);
  
  if (_time > that._time) 
    return(false);

  return((_value < that._value));
}
