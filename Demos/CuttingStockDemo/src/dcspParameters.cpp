/**
 *
 * This file dcspParameters.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dcspParameters.hpp"

ApplicationSpecificParam::ApplicationSpecificParam(std::string parameterFileName, int argc, char* argv[], std::ostream& os) :
        ParameterParser(parameterFileName),
        _readProb(std::string("readProb"), false,
            "Reading problems from data or generate them."),
        _nbInstances(std::string("nbInstances"), 10,
            "Nb of Instance  on which we establish statistics."),
        _firstInstance(std::string("firstInstance"), 0,
            "Nb of the Instance  on which we start statistics."),
        _nbOrders(std::string("nbOrders"), 10,
            "Nb orders during random generation problem."),
        _maxDemand(std::string("maxDemand"), 10, "Max demand"),
        _inputFile(std::string("inputFile,i"), "inputList", "inputFile data")
{
  addApplicationParameter(_readProb);
  addApplicationParameter(_nbInstances);
  addApplicationParameter(_firstInstance);
  addApplicationParameter(_nbOrders);
  addApplicationParameter(_maxDemand);
  addApplicationParameter(_inputFile);

  parse(argc, argv);

  print(os);
}

std::ostream& ApplicationSpecificParam::print(std::ostream& os) const
{
  os << "Cutting Stock Parameters  " << std::endl;
  os <<  _readProb << std::endl;
  os << _nbInstances << std::endl;
  os << _firstInstance << std::endl;
  os << _nbOrders << std::endl;
  os << _maxDemand << std::endl;
  os << _inputFile << std::endl;
  return (os);
}

bool ApplicationSpecificParam::readProb() const
{
  return _readProb();
}

int ApplicationSpecificParam::nbInstances() const
{
  return _nbInstances();
}

int ApplicationSpecificParam::firstInstance() const
{
  return _firstInstance();
}

void ApplicationSpecificParam::firstInstance(int first)
{
  return _firstInstance(first);
}

int ApplicationSpecificParam::nbOrders() const
{
  return _nbOrders();
}

int ApplicationSpecificParam::maxDemand() const
{
  return _maxDemand();
}

const std::string& ApplicationSpecificParam::inputFile() const
{
  return _inputFile();
}

