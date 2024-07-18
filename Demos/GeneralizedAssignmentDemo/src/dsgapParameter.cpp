/**
 *
 * This file dsgapParameter.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dsgapParameter.hpp"

ApplicationSpecificParam::ApplicationSpecificParam(std::string parameterFileName) :
    ParameterParser(parameterFileName),
        _testSimpleLp(std::string("testSimpleLp"), 0, "Test a small Lp to test solveLpFunctionality"),
        _readProb(std::string("readProb"), false, "Reading problems from data or generate them."),
        _inputFile(std::string("inputFile,i"), "data/inputList", "inputFile data"),
        _nbProblems(std::string("nbProblems"), 2, "Nb Generate Problems."),
        _firstInstance(std::string("firstInstance"), 0, "First Problem to be solved."),
        _nbJobs(std::string("nbJobs"), 20, "Nb Jobs per problem."),
        _tightness("tightness", 0.95, "TODO"),
        _machineMaxCapacity("machineMaxCapacity", 10000, "Max load capacity per machine"),
        _maxJobCapUsage("maxJobCapUsage", 2500, "TODO"),
        _minJobCapUsage("minJobCapUsage", 500, "TODO"),
        _justSolveMIP("justSolveMIP", false, "Undocumented")
{
  addApplicationParameter(_testSimpleLp);
  addApplicationParameter(_readProb); //add _readProb to the parser
  addApplicationParameter(_inputFile);
  addApplicationParameter(_nbProblems);
  addApplicationParameter(_firstInstance);
  addApplicationParameter(_nbJobs);
  addApplicationParameter(_tightness);
  addApplicationParameter(_machineMaxCapacity);
  addApplicationParameter(_maxJobCapUsage);
  addApplicationParameter(_minJobCapUsage);
  addApplicationParameter(_justSolveMIP);
}

ApplicationSpecificParam::~ApplicationSpecificParam()
{
}

std::ostream& ApplicationSpecificParam::print(std::ostream& os) const
{
  os << std::endl;
  os << "##### GeneralizedAssignment Parameters #####" << std::endl;
  os << _readProb << std::endl;
  os << _inputFile << std::endl;
  os << _nbProblems << std::endl;
  os << _firstInstance << std::endl;
  os << _nbJobs << std::endl;
  os << _machineMaxCapacity << std::endl;
  os << _minJobCapUsage << std::endl;
  os << _maxJobCapUsage << std::endl;
  os << _tightness << std::endl;
  os << _justSolveMIP << std::endl;
  os << "#############################################" << std::endl << std::endl;

  return os;
}

const std::string& ApplicationSpecificParam::inputFile() const
{
  return _inputFile();
}

double ApplicationSpecificParam::machineMaxCapacity() const
{
  return _machineMaxCapacity();
}

double ApplicationSpecificParam::maxJobCapUsage() const
{
  return _maxJobCapUsage();
}

double ApplicationSpecificParam::minJobCapUsage() const
{
  return _minJobCapUsage();
}

int ApplicationSpecificParam::nbJobs() const
{
  return _nbJobs();
}

int ApplicationSpecificParam::nbProblems() const
{
  return _nbProblems();
}

int ApplicationSpecificParam::firstInstance() const
{
  return _firstInstance();
}

bool ApplicationSpecificParam::readProb() const
{
  return _readProb();
}

double ApplicationSpecificParam::tightness() const
{
  return _tightness();
}

bool ApplicationSpecificParam::justSolveMIP() const
{
    return _justSolveMIP();
}

int ApplicationSpecificParam::testSimpleLp() const
{
    return _testSimpleLp();
}
