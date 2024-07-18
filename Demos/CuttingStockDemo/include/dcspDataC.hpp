/**
 *
 * This file dcspDataC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DcspDataClasses_h
#define DcspDataClasses_h

#include "bcDoubleC.hpp"
#include <vector>
#include <iostream>
#include "dcspParameters.hpp"

class Order
{
private:
  int _orderIndex;
  double _demand;
  double _width;
 public:
  Order(const int & orderIndex, 
	const double & demand,
	const double & width);
  virtual ~Order();

  const double & demand() const;
  const double & width() const;
  std::ostream& print(std::ostream& os = std::cout) const;
};

inline bool operator<(const Order &oa, const Order &ob)
{
  // higher priority to larger width 
  if (oa.width() > ob.width())
    {
      return true;
    }
  return false;
}

inline std::ostream& operator<<(std::ostream& os, const Order & m)
{
  return m.print(os);
}

class StockSheet
{
  int _stockSheetIndex;
  double _width;

public:
  StockSheet(const int & stockSheetIndex, 
	     const double & width);
  virtual ~StockSheet();

  const double & width() const
  {
    return _width;
  }
 
  std::ostream& print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const StockSheet & m)
{
  return m.print(os);
}

class DcspData
{
private:
  std::string _modelName;
  std::vector< Order * > _orderPtrVector;
  StockSheet * _stockSheetPtr; /// for use in case of a single stock sheet type
  std::vector< StockSheet * > _stockSheetPtrVector; /// for use in case of multiple stock sheet type
  ApplicationSpecificParam* _dcspParameterPtr;
public:
  DcspData(ApplicationSpecificParam* dcspParameterPtr);
  virtual ~DcspData();

  void readData(const std::string & inputFileName);
  void generateAtRandom(const int & genSeed);
  const std::vector<Order*>& orderPtrVector() const;
  const StockSheet& stockSheet() const;
  const ApplicationSpecificParam& dcspParameter() const;
  const std::string& modelName() const;
  void modelName(const std::string& modelName);

  std::ostream& print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const DcspData & data)
{
  return data.print(os);
}

#endif // DcspDataClasses_h
