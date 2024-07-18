/**
 *
 * This file bcPrintC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef PRINTCLASS_H
#define PRINTCLASS_H

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcBapcodUtilities.hpp"

class PrintLevel
{
  static int printLevel;

  PrintLevel() = delete;
public:
  static void setPrintLevel(const int &i) {printLevel = i;}
  friend bool printL(const int& ifPlevel);
};

inline bool printL(const int& ifPlevel)
{
  return (PrintLevel::printLevel >= ifPlevel);
}

class ProgStatistics
{
protected:
  std::string _prName;
  struct Record
  {
    std::string _key;
    long _counter;
    double _time;
    Double _value;
    Record(const std::string & key = "NULL", const long & counter = ProgStatisticsUndefinedStatus,
           const double & time = ProgStatisticsUndefinedStatus,
           const Double & value = ProgStatisticsUndefinedStatus);

    Record(const Record& record);

    bool operator<(const Record & that) const;
    void operator+=(const Record & that);
    void operator/=(const int & den);
    std::ostream & print(std::ostream& os = std::cout) const;
    std::ostream& plainPrint(std::ostream& os = std::cout) const;
    std::ostream& printAverage(std::ostream& os) const;
  };

  struct smallerKey
  {
    bool operator()(const Record & a, const Record & b) const
    {
      return (a._key < b._key);
    }
  };

protected:
  Record _emptyRec;
  std::multiset<Record, smallerKey> _recordMultiSet;
  std::vector<std::string> _selectedKeys;
  void replaceRecord(const Record & oldRec, const Record & newRec);
public:
  ProgStatistics(const std::vector<std::string> & selectedKeys);
  ProgStatistics(const std::vector<ProgStatistics> & prStatVector);
  ProgStatistics(); /// compute the average of whatever is in programStatisticsVector

  virtual ~ProgStatistics();

  void recName(const std::string & name)
  {
    _prName = name;
  }

  void record(const std::string & key, const long & counter, const double & time,
              const Double & value);

  void incrRecord(const std::string & key, const Double & val);
  const Record & getRecord(const std::string & key)
  {
    if (_recordMultiSet.count(key))
      return *_recordMultiSet.lower_bound(key);
    else
      return _emptyRec;
  }

  void recCounter(const std::string & key, const long & counter);
  void recTime(const std::string & key, const double & time);
  void recValue(const std::string & key, const Double & value);
  void setCounter(const std::string & key, const long & counter);
  void incrCounter(const std::string & key, const long & incr = 1);
  void incrTimer(const std::string & key, const double & timer);
  void incrValue(const std::string & key, const Double & value);
  void incrSummedValue(const std::string & key, const Double & value);
  const std::vector<std::string> & selectedKeys() const
  {
    return _selectedKeys;
  }

  long getCounter(const std::string & key);
  double getValue(const std::string & key);
  double getTime(const std::string & key);
  std::ostream & printStat(std::ostream& os = std::cout);
  std::ostream & print(std::ostream& os = std::cout);

  std::ostream& selectPrint(std::ostream& os = std::cout, const bool printEndLine = true);
  std::ostream& titlePrint(std::ostream& os = std::cout);
};

/**
 * operator<<: Appends any item (preceeded by a white space)
 * To a std::string using overloaded << if (printLevel >= 5)
 */
template<typename T>
std::string& operator<<(std::string & s, const T & item);

template<typename T>
std::string& operator<<(std::string & s, const T & item)
{
  if (printL(6))
    return s;
  std::ostringstream oss;
  oss << " " << item;
  s.append(oss.str());

  /// We took over the control of memeory by using str()
  return s;
}

/**
 * Prints the contents of any sequence for which print() is defined if (printLevel >= ifPlevel)
 */
template<typename Container>
std::ostream& printSeq(const Container & c, const char * nm = "", const int& ifPlevel = 5,
                       std::ostream& os = std::cout);

template<typename Container>
std::ostream& printSeq(const Container & c, const char * nm, const int& ifPlevel, std::ostream& os)
{
  if (!printL(ifPlevel))
    return os;

  /// Only if you provide a std::string
  if (*nm != '\0')
    /// Is this printed
    os << nm << std::endl;

  if (c.empty())
  {
    os << "Container is empty" << std::endl;
    return os;
  }

  os << "Container size is " << c.size() << std::endl;

  for (typename Container::const_iterator it = c.begin(); it != c.end(); ++it)
    (*it).print(os);

  return os;
}

#endif // PRINTCLASS_H
