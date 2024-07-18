/**
 *
 * This file dcspParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef PROBPARAMETERS_H
#define PROBPARAMETERS_H

#include "bcModelParameterC.hpp"
#include "bcParameterC.hpp"
#include "bcParameterParserC.hpp"
#include <iostream>
#include <string>

/* PROBLEM SPECIFIC PARAMETERS */
class ApplicationSpecificParam : public ParameterParser
{
private:
  ApplicationParameter<bool> _readProb;
  ApplicationParameter<int> _nbInstances;
  ApplicationParameter<int> _firstInstance;
  ApplicationParameter<int> _nbOrders;
  ApplicationParameter<int> _maxDemand;
  ApplicationParameter<std::string> _inputFile;
 
public:
  ApplicationSpecificParam(std::string parameterFileName, int argc, char* argv[], std::ostream& os = std::cout);

  virtual ~ApplicationSpecificParam() {}
  std::ostream& print(std::ostream& os = std::cout) const;

  bool readProb() const;
  int nbInstances() const;
  int firstInstance() const;
  void firstInstance(int firstInstance);
  int nbOrders() const;
  int maxDemand() const;
  const std::string& inputFile() const;
};

inline std::ostream& operator<<(std::ostream& os, const ApplicationSpecificParam & that)
{return that.print(os);}


#endif // PROBPARAMETERS_H
