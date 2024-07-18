/**
 *
 * This file bcBapcodInit.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcBapcodInit.hpp"
#include "bcFormC.hpp"
#include "bcVarConstrC.hpp"
#include "bcKnapsackSolverC.hpp"
#include "bcMathProgSolverBuilder.hpp"
#include "bcModelC.hpp"
#include "bcParameterManager.hpp"
#include "bcOutputRunSummary.hpp"
#include "bcAlg4GenChildrenInBranching.hpp"
#include "bcPrintC.hpp"
#include "bcBapcodVersion.hpp"

#if BCP_RCSP_IS_FOUND
#include "rcsp_version.hpp"
#endif

const Double Double::staticZero = 0.0;
const Double Double::staticOne = 1.0;
const double Double::precision = 1e-06;
const double Double::highprecision = 1e-10;

const double Bound::epsilonBound = 0.000001;

const LpCoef LpCoef::ZeroCoef = LpCoef(false,0.0);
const LpCoef LpCoef::UnitCoef = LpCoef(true,1.0);
const LpCoef LpCoef::MinusOneCoef = LpCoef(true,-1.0);

const double CandidateBranchGroup::treeDepthScoreNormConstant = 1;

/// this variable seems to be the only one global variable
/// it should not impact the thread-safety of BaPCod, as long as DEFAULTPRINTLEVEL parameter stays the same
int PrintLevel::printLevel = 0;

BapcodInit::BapcodInit() :
    _startTime(NULL), _modelPtr(NULL), _statistics(), _progStatus(NULL, &_statistics), _testLevel(1),
    _controlParam(), _createStatsAndBapTreeFiles(true)
{
    initializeParameters();
}

BapcodInit::BapcodInit(int argc, char *argv[], std::string filename, bool toPrintParameters, bool printVrpSolverHeader,
                       bool printBapcodHeader):
    _startTime(NULL), _modelPtr(NULL), _statistics(), _progStatus(NULL, &_statistics), _testLevel(1),
    _controlParam()
{
    if (filename == "NOT_SPECIFIED")
        initializeParameters();
    else
        initializeParameters(argc, argv, filename);

    if (printL(-1))
    {
#if BCP_RCSP_IS_FOUND
        if (printVrpSolverHeader)
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      << std::endl
                      << "VRPSolver v" << BCP_RCSP_VERSION << ", " << BCP_RCSP_VERSION_DATE
                      << ", © Inria Bordeaux, France, vrpsolver.math.u-bordeaux.fr" << std::endl
                      << "      Corresponds to the solver by Pessoa, Sadykov, Uchoa and Vanderbeck (2020)" << std::endl
                      << "                 Paper: dx.doi.org/10.1007/s10107-020-01523-z" << std::endl;
#endif
        if (printBapcodHeader)
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      << std::endl
                      << "   BaPCod v" << BAPCOD_VERSION << ", " << BAPCOD_VERSION_DATE
                      << ", © Inria Bordeaux, France, bapcod.math.u-bordeaux.fr " << std::endl
                      << "           THIS CODE IS PROVIDED AS IS, USE IT AT YOUR OWN RISK" << std::endl
                      << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                      << std::endl;
    }

    if (printL(1))
    {
        printParameters(false, true, std::cout);
    }
    else if (printL(0))
    {
        if (toPrintParameters)
        {
            printParameters(false, false, std::cout);
            if (_createStatsAndBapTreeFiles)
                _statistics.titlePrint();
        }
    }
    else if (printL(-1) && printVrpSolverHeader)
    {
        printParameters(true, false, std::cout);
    }

    if (param().statistics_file() != "")
    {
        std::ofstream ofs(param().statistics_file().c_str(), std::ios::out);
        _statistics = ProgStatistics(param().outputKeyVector());
        _statistics.titlePrint(ofs);
        ofs.close();
    }

    _startTime = new Time();
    _progStatus.startTime(_startTime);
    _progStatus.timeLimit(_controlParam.GlobalTimeLimitInTick());
}

BapcodInit::~BapcodInit()
{
	if (_programStatisticsVector.size() > 1)  // if added by Boris to avoid a bug occuring when
											  // no output of stats is done before the
											  // destructor is called
	{
		ProgStatistics average(_programStatisticsVector);
        if (printL(0))
		  std::cout << " average print" << std::endl;
		average.titlePrint();
		average.print();
		average.selectPrint();
        if (param().statistics_file() != "")
        {
            std::ofstream ofs(param().statistics_file().c_str(), std::ios::app);
            average.selectPrint(ofs);
            ofs.close();
        }
	}

  if (printL(1))
    progStatus().print(std::cout);

  delete _startTime;
}

/// initialize parameters with default values
void BapcodInit::initializeParameters()
{
    ParameterManager paramManager;
    _controlParam.addParameters(paramManager);
    _startTime = new Time();
    _progStatus.startTime(_startTime);
}

void BapcodInit::initializeParameters(int argc, char *argv[], std::string filename)
{
    //Initialize Bapcod parameters and application parameters.
    ParameterManager paramManager;

    //The method here should be normally named "addParametersTo"
    _controlParam.addParameters(paramManager);

    //The following line was added by Issam as filename was not used by
    //This method. There is a problem in the design of the paramManager
    //Fixing it will require a lot of modifications so i am just using
    //a workaround here and in the method ParameterManager::parse
    _controlParam.configFileName(filename);

    paramManager.parse(argc, argv, _controlParam.configFileName());

    _controlParam.postTreatment();

    PrintLevel::setPrintLevel(param().DEFAULTPRINTLEVEL());

    _testLevel = param().DEFAULTTESTLEVEL();
}

void BapcodInit::reset(bool analyseParameters)
{
    delete _startTime;
    _startTime = new Time();
    _progStatus.startTime(_startTime);

    if (analyseParameters)
    {
        _controlParam.postTreatment();
        PrintLevel::setPrintLevel(param().DEFAULTPRINTLEVEL());
    }
    _progStatus.timeLimit(param().GlobalTimeLimitInTick());
    _progStatus.resetStat();

    _statistics = ProgStatistics(param().outputKeyVector);
}

void BapcodInit::outputBaPCodStat(const std::string & modelName, std::ostream & os)
{
  _statistics.recName(modelName);
  _statistics.print(os);
  _statistics.selectPrint(os);
    if (param().statistics_file() != "")
    {
        std::ofstream ofs(param().statistics_file().c_str(), std::ios::app);
        _statistics.selectPrint(ofs);
        ofs.close();
    }

  OutputRunSummary pSol(this, modelName);
  recordProgStatisticsInstance(_statistics);

  printTime(_startTime->getElapsedTime(), os);
  if (printL(1))
    progStatus().print(os);
}

const Time & BapcodInit::startTime() const
{
  return *_startTime;
}

const ProgStatus& BapcodInit::progStatus() const
{
  return _progStatus;
}

ProgStatus& BapcodInit::progStatus()
{
  return _progStatus;
}

void BapcodInit::printParameters(bool VRPSolverParametersOnly, bool allParameters, std::ostream & os) const
{
  if (VRPSolverParametersOnly)
  {
      param().printVRPSolverParameters(os);
      return;
  }

  param().printUserParameters(os);
  if (allParameters)
      param().printDevParameters(os);
}

Model * BapcodInit::modelPtr() const
{
  return _modelPtr;
}

void BapcodInit::modelPtr(Model * modelPtr)
{
  _modelPtr = modelPtr;
}

const ProgStatistics & BapcodInit::statistics() const
{
  return _statistics;
}

ProgStatistics& BapcodInit::statistics()
{
  return _statistics;
}

void BapcodInit::recordProgStatisticsInstance(const ProgStatistics& instance)
{
  _programStatisticsVector.push_back(instance);
}

const int& BapcodInit::testLevel() const
{
  return _testLevel;
}


