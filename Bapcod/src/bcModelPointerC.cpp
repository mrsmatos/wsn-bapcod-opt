/**
 *
 * This file bcModelPointerC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcBoundLevC.hpp"
#include "bcModelC.hpp"
#include "bcModelingLanguageC.hpp"
#include "bcPrintC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcBapcodInit.hpp"
#include "bcGlobalException.hpp"
#include "bcOutputRunSummary.hpp"

/**
 * Generic code
 */
using namespace std;

BcInitialisation::BcInitialisation():
   _bcpInitptr(new BapcodInit())
{
}

BcInitialisation::BcInitialisation(std::string filename, bool toPrintParameters, bool printVrpSolverHeader,
                                   bool printBapcodHeader) :
    BcInitialisation(0, {}, filename, toPrintParameters, printVrpSolverHeader, printBapcodHeader)
{
}


BcInitialisation::BcInitialisation(int argc, char *argv[], std::string filename, bool printParameters,
                                   bool printVrpSolverHeader, bool printBapcodHeader):
  _bcpInitptr(new BapcodInit(argc, argv, filename, printParameters, printVrpSolverHeader, printBapcodHeader))
{
  _bcpInitptr->reset(false);
}

BcInitialisation::~BcInitialisation()
{
  delete _bcpInitptr;
}

void BcInitialisation::bcReset()
{
  _bcpInitptr->reset(); 
}

ControlParameters & BcInitialisation::param()
{
  return _bcpInitptr->param();
}

void BcInitialisation::outputBaPCodStatistics(const std::string & instanceName, std::ostream & os)
{
    _bcpInitptr->outputBaPCodStat(instanceName, os);
}

void BcInitialisation::outputBaPCodStatistics(const std::string & instanceName, std::ostream & os, std::ostream & fs,
                                              std::string comFileName)
{
  _bcpInitptr->outputBaPCodStat(instanceName, os);
  _bcpInitptr->statistics().selectPrint(fs);
  if (printL(0))
    std::cout << "comFileName = " << comFileName << std::endl;
  if (comFileName != "")
    {
      std::ofstream cs(comFileName.c_str(), std::ios::out);
      _bcpInitptr->statistics().printStat(cs);
      cs << "statLine=\"";
      _bcpInitptr->statistics().selectPrint(cs, false);
      cs << "\"" << std::endl;
      cs.close();
    }
    _bcpInitptr->reset(); /// should we reset?
}

void BcInitialisation::addStatisticsValue(const std::string & statName, double value)
{
  _bcpInitptr->statistics().recValue(statName, value);
}

void BcInitialisation::incrStatisticsCounter(const std::string & statName, int value)
{
  _bcpInitptr->statistics().incrCounter(statName, value);
}

double BcInitialisation::getStatisticValue(const std::string & statName)
{
  return(_bcpInitptr->statistics().getValue(statName));
}

long BcInitialisation::getStatisticCounter(const std::string & statName)
{
  return(_bcpInitptr->statistics().getCounter(statName));
}

double BcInitialisation::getStatisticTime(const std::string & statName)
{
  return(_bcpInitptr->statistics().getTime(statName));
}

void BcInitialisation::outputResultsForTests()
{
    OutputRunSummary summary(_bcpInitptr, "");
}

BcModel::BcModel(const BcInitialisation & bapcodInitPtr,
		         const std::string & modelName,
		         const BcObjStatus::MinMaxIntFloat & objectiveSense) :
  _modelPtr(new Model(bapcodInitPtr._bcpInitptr, modelName, objectiveSense)), _deleteInDestructor(true)
{
  if (printL(5)) 
    std::cout << "NEW ModelPtr(" << modelName << ") " << std::endl;
}

BcModel::BcModel(Model * modelPtr) :
  _modelPtr(modelPtr), _deleteInDestructor(false)
{
}

std::ostream& operator<<(std::ostream& os, const BcModel & that)
{
  if (that._modelPtr == NULL)
    return os << "operator<<BcModel: undefined model" << std::endl;

  if (printL(5))
    std::cout << "operator<<(BcModel) " << std::endl;

  that._modelPtr->setup();

  return that._modelPtr->print(os);
}

BcModel::~BcModel()
{
  if (_deleteInDestructor && (_modelPtr != NULL))
    {
      delete _modelPtr;
      _modelPtr = NULL;
    }
}

BcModel::operator Model *() const
{
  return _modelPtr;
}

BcFormulation BcModel::master()
{
    if (_modelPtr == NULL)
    {
        throw GlobalException("ModelPtr::master() undefined pointer", true);
    }
    return BcFormulation(_modelPtr->masterConfPtr());
}

void BcModel::attach(BcSolutionFoundCallback * solutionFoundCallbackPtr)
{
  if(_modelPtr == NULL)
    {
      throw GlobalException("ModelPtr::objective() undefined pointer", true);
    }
  _modelPtr->setSolutionFoundCallback(solutionFoundCallbackPtr);
}

BcSolution BcModel::solve() 
{
  return BcSolution(_modelPtr->solve());
}

BcSolution BcModel::enumerateAllColumns(int & nbEnumColumns)
{
  return BcSolution(_modelPtr->enumerateAllColumns(nbEnumColumns));
}

/// returns triples (reduced cost, cost, solution pointer)
void BcModel::getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions)
{
    _modelPtr->getEnumeratedSolutions(enumSolutions);
}

const std::string& BcInitialisation::instanceFile() const
{
  return _bcpInitptr->param().instance_file();
}

const std::string& BcInitialisation::configFile() const
{
  return _bcpInitptr->param().configFileName();
}

const std::string& BcInitialisation::statFile() const
{
  return _bcpInitptr->param().statistics_file();
}
