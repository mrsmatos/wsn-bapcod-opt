/**
 *
 * This file dcspDataC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dcspDataC.hpp"
#include "bcRandGen.hpp"
#include "bcBapcodUtilities.hpp"

using namespace std;

/**
 * Exit the application if we can't read the data.
 * @param value wich value we can't read.
 */
inline void readingError(string value)
{
  cerr << "\nerror: could not read " << value << endl;
  exit(EXIT_FAILURE);
}

Order::Order(const int & orderIndex, 
	         const double & demand,
	         const double & width) :
  _orderIndex(orderIndex), _demand(demand), _width(width)
{
}

Order::~Order()
{
}

const double & Order::demand() const
{
  return _demand;
}

const double & Order::width() const
{
  return _width;
}

ostream& Order::print(ostream& os) const
{
  os << "Order: " << endl;
  os << "   orderIndex = " << _orderIndex << endl;
  os << "   demand = " << _demand << endl;
  os << "   width = " << _width << endl;
  return os;
}

StockSheet::StockSheet(const int & stockSheetIndex, const double & width) :
  _width(width)
{
}

StockSheet::~StockSheet()
{
}

ostream& StockSheet::print(ostream& os) const
{
  os << "StockSheet: " << endl;
  os << "   stockSheetIndex = " << _stockSheetIndex << endl;
  os << "   width = " << _width << endl;
  return os;
}


DcspData::DcspData(ApplicationSpecificParam* dcspParameterPtr): _dcspParameterPtr(dcspParameterPtr)
{
}

DcspData::~DcspData()
{
}

void DcspData::readData(const string& inputFileName)
{
  ifstream is(inputFileName.c_str(), ios::in);
  if (!is)
    {
      cout << "cannot find input file" << endl;
      exit(EXIT_FAILURE);
    }
  if (is.eof())
    {
      cout << "empty input file" << endl;
      exit(EXIT_FAILURE);
    }

  //Read of the model name.
  if (!(is >> _modelName))
    {
      readingError("modelName");
    }

  //Read of the stock sheet width.
  int stockSheetWidth = 0;
  if (!(is >> stockSheetWidth))
    {
      readingError("stockSheetWidth");
    }
  _stockSheetPtr = new StockSheet(0, stockSheetWidth);

  //Read of the stock sheet width.
  int nbOrders = 0;
  if (!(is >> nbOrders))
    {
      readingError("nbOrders");
    }

  //Read (till the end of file) of the width roll and number of demand for this roll.
  int orderWidth = 0;
  int orderDemand = 0;
  for (int currentOrder = 0; currentOrder < nbOrders; ++currentOrder)
    {
      if (!(is >> orderWidth))
        {
          readingError("orderWidth");
        }
      if (!(is >> orderDemand))
        {
          readingError("orderDemand");
        }
      _orderPtrVector.push_back(
          new Order(currentOrder, orderDemand, orderWidth));
    }
}

void DcspData::generateAtRandom(const int& genSeed)
{
  bool generateIntData = true;

  int nbOrders = _dcspParameterPtr->nbOrders();
  int maxDem = _dcspParameterPtr->maxDemand();
  double stockSheetWidth(10000);
  double maxOrderWidth(stockSheetWidth / 4.0);
  double minOrderWidth(stockSheetWidth / 20.0);

  int i;
  RandGen orderWidthGen((int) (maxOrderWidth - minOrderWidth), 0,
      (1 + genSeed) * 17 + 51);
  RandGen demGen(maxDem, 1, (1 + genSeed) * 17 + 51);
  double width;
  double demand;
  string outputFileName("randomDCSPdata");
  outputFileName = outputFileName + genSeed;
  ofstream os(outputFileName.c_str(), ios::out);

  _modelName = "RandInst";
  _modelName = _modelName + genSeed;
  os << _modelName << endl;
  os << stockSheetWidth << endl;
  os << nbOrders << endl;

  double totalOrderWidth(0);
  for (int orderIndex = 0; orderIndex < nbOrders; ++orderIndex)
    {
      width = minOrderWidth + orderWidthGen();
      if (generateIntData)
        width = Dfloor(width);
      demand = demGen(); // in [1,maxdem]; //
      _orderPtrVector.push_back(new Order(orderIndex, demand, width));
      os << width << " " << demand << endl;
      totalOrderWidth += width * demand;
    }
  _stockSheetPtr = new StockSheet(0, stockSheetWidth);
}

const vector<Order*>& DcspData::orderPtrVector() const
{
  return _orderPtrVector;
}

const StockSheet& DcspData::stockSheet() const
{
  return *_stockSheetPtr;
}

const ApplicationSpecificParam& DcspData::dcspParameter() const
{
  return *_dcspParameterPtr;
}

const std::string& DcspData::modelName() const
{
  return _modelName;
}

void DcspData::modelName(const std::string& modelName)
{
  _modelName = modelName;
}

ostream& DcspData::print(ostream& os) const
{
  os << stockSheet() << endl;

  os << "Nb orders: " << _orderPtrVector.size() << endl;

  for (vector<Order*>::const_iterator it=_orderPtrVector.begin(); it != _orderPtrVector.end(); ++it)
    {
      os << *it << endl;
    }

  return os;
}

