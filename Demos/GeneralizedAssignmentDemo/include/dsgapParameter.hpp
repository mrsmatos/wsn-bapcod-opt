/**
 *
 * This file dsgapParameter.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DSGAPPARAMETER_HPP_
#define DSGAPPARAMETER_HPP_
#include "bcParameterParserC.hpp"

/**
 * This class permits to read and print some parameters for the application.
 * It uses the class ParameterParser to help its implementation.
 */
class ApplicationSpecificParam : public ParameterParser
{
private:
  ApplicationParameter<int> _testSimpleLp;

  ///Do we read from files or generate the problems.
  ApplicationParameter<bool> _readProb;

  ///If readProb is true, name of the input file.
  ApplicationParameter<std::string> _inputFile;

  ///If readProb is false, number of problems needed to be generate.
  ApplicationParameter<int> _nbProblems;

  ///If readProb is false, defines the first random instance that is being solved
  ApplicationParameter<int> _firstInstance;

  ///If readProb is false, number type of jobs.
  ApplicationParameter<int> _nbJobs;

  ///TODO: Comment this parameter.
  ApplicationParameter<double> _tightness;

  ///If readProb is false, max capacity per machine.
  ApplicationParameter<double> _machineMaxCapacity;

  ApplicationParameter<double> _maxJobCapUsage;
  ApplicationParameter<double> _minJobCapUsage;

  ApplicationParameter<bool> _justSolveMIP;

public:
  ApplicationSpecificParam(std::string parameterFileName);
  virtual ~ApplicationSpecificParam();

  /**
   * Print in the given \p os the parameter value.
   * @param os
   * @return output stream
   */
  std::ostream& print(std::ostream& os = std::cout) const;

  int testSimpleLp() const;
  const std::string& inputFile() const;
  double machineMaxCapacity() const;
  double maxJobCapUsage() const;
  double minJobCapUsage() const;
  int nbJobs() const;
  int nbProblems() const;
  int firstInstance() const;
  bool readProb() const;
  double tightness() const;
  bool justSolveMIP() const;
};

inline std::ostream& operator<<(std::ostream& os, const ApplicationSpecificParam & that)
{
  return that.print(os);
}

#endif /* DSGAPPARAMETER_HPP_ */
