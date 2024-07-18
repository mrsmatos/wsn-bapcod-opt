/**
 *
 * This file bcBapcodInit.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCBAPCODINIT_H_
#define BCBAPCODINIT_H_

#include <iostream>
#include <vector>
#include "bcDoubleC.hpp"
#include "bcErrorC.hpp"
#include "bcPrintC.hpp"
#include "bcControlParameters.hpp"

class Time;
class Model;

class BapcodInit
{
private:
  Time * _startTime;
  Model * _modelPtr;

  ProgStatistics _statistics;
  ProgStatus _progStatus;
  std::vector<ProgStatistics> _programStatisticsVector;

  int _testLevel;

  ControlParameters _controlParam;
  bool _createStatsAndBapTreeFiles; /// allows to disable the generation of the statistics and baptree files

  BapcodInit(const BapcodInit& bapcodInit) = delete;

public:

   BapcodInit();

   /**
   * Construct the object and read BapcodParameters. Reads bapcod parameters and application parameters
   * @param argc
   * @param argv
   * @param filename name of the config file. Not used for now.
   * @param printParameters if false, the parameters are not printed.
   * @param createFiles if false, the statistics and baptree files are not generated.
   * @param os output stream.
   */
  BapcodInit(int argc, char *argv[], std::string filename = "", bool toPrintParameters = false,
             bool printVrpSolverHeader = false, bool printBapcodHeader = true);

  virtual ~BapcodInit();

  /**
   * Reset all variables;
   */
  void reset(bool analyseParameters = true);

  void outputBaPCodStat(const std::string & modelName, std::ostream & os = std::cout);

  void printParameters(bool VRPSolverParametersOnly = false, bool allParameters = false,
                       std::ostream & os = std::cout) const;

  inline const ControlParameters & param() const
  {
    return _controlParam;
  }

  inline ControlParameters & param()
  {
    return _controlParam;
  }

  Model * modelPtr() const;
  void modelPtr(Model * modelPtr);

  const Time & startTime() const;
  const ProgStatus& progStatus() const;
  ProgStatus & progStatus();

  const ProgStatistics & statistics() const;
  ProgStatistics & statistics();

  void recordProgStatisticsInstance(const ProgStatistics & instance);

  const int & testLevel() const;

  /**
   * Expect requirement to be true or else generate error message and stop
   */
  template<typename T>
  T require(T requirement, const char *const msg = "Requirement is not statisfied",
            const ProgStatus::Type & terminationStatus = ProgStatus::quit,
            const int& ifTlevel = 1, std::ostream & os = std::cerr)
  {
    if ((testLevel() >= ifTlevel) && (requirement == 0))
      {
        os << msg << " error code = " << requirement << std::endl;

        _progStatus.setStat(terminationStatus);
        _progStatus.pushMsg(msg);
      }

    return requirement;
  }

  /**
   * Expect error to be false
   * or else generate error message and stop
   */
  template<typename T>
  T check(T error, const char * const msg = "Check has failed",
          const ProgStatus::Type & terminationStatus = ProgStatus::quit,
          const int& ifTlevel = 1, std::ostream & os = std::cerr)
  {

    if ((testLevel() >= ifTlevel) && (error != 0))
      {
        os << msg << " error code = " << error << std::endl;

        _progStatus.setStat(terminationStatus);
        _progStatus.pushMsg(msg);
      }
    return error;
  }

private:

  /**
   * Read the parameter.
   * @param argc
   * @param argv
   * @param filename name of the config file. Not used for now.
   * @param printParameters if false, the parameters are not printed.
   * @param os output stream.
   */
  void initializeParameters(int argc, char *argv[], std::string filename);
  void initializeParameters();

}; //BapcodInit


#endif /* BCBAPCODINIT_H_ */
