/**
 *
 * This file dsgapDataC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dsgapDataC.hpp"
#include "bcBapcodUtilities.hpp"

using namespace std;

Machine::Machine(const int & idNb,
		 const Double & capacity):
  _idNb(idNb), _capacity(capacity)
{
}

std::ostream& Machine::print(std::ostream& os) const
{
  os << "Machine: "  << std::endl;
  os << "   idNb = " << _idNb  << std::endl;
  os << "   capacity = " << _capacity  << std::endl;
  return (os);
}


Job::Job(const int & idNb,
	 const Double & nominalCapacityUsage):
  _idNb(idNb), _nominalCapacityUsage(nominalCapacityUsage), _maxCostOnMachine(0)
{
}

Job::Job(const int & idNb,
	 const int & nbMachines,
	 const std::vector<int> & resConsumption, 
	 const std::vector<Double> & cost):
  _idNb(idNb)
{
  _maxCostOnMachine = 0;
  for (int i = 0; i < nbMachines; i++)
    {
      _capUsageOnMachine.push_back(resConsumption[i]);
      _costOnMachine.push_back(cost[i]);
      if (cost[i] > _maxCostOnMachine) 
	_maxCostOnMachine = cost[i];
    }
}

Job::~Job()
{}

void Job::randomWeightGeneration(const int & nbMachines)
{
  RandGen jobCapUsageOrCostGen( (int)(_nominalCapacityUsage * 1.2), 
			  (int)(_nominalCapacityUsage * 0.8), 
			  nbMachines * 17 + 51);

  Double cost(0);
  for(int k = 0; k <  nbMachines ; ++k)
    {
      _capUsageOnMachine.push_back(jobCapUsageOrCostGen());
      cost = jobCapUsageOrCostGen();
      _costOnMachine.push_back(cost);
      if (cost > _maxCostOnMachine) _maxCostOnMachine = cost;
    }
  return;
}

std::ostream& Job::print(std::ostream& os) const
{ 
  os << "Job: "  << std::endl;
  os << "   idNb = " << _idNb  << std::endl;
  os << "   nominalCapacityUsage = " << _nominalCapacityUsage << std::endl;
  os << "   maxCostOnMachine = " <<  _maxCostOnMachine << std::endl;
  int m = 0;
  for (std::vector<Double>::const_iterator itc = _costOnMachine.begin();
       itc != _costOnMachine.end(); ++m, ++itc)
           os << "   _costOnMachine[" << m << "] =  " << *itc << std::endl;
  m = 0;
  for (std::vector<Double>::const_iterator itc = _capUsageOnMachine.begin();
       itc != _capUsageOnMachine.end(); ++m, ++itc)
           os << "   _capUsageOnMachine[" << m << "] =  " << *itc << std::endl;
  return (os);
}

GapData::GapData(ApplicationSpecificParam* paramPtr) : _paramPtr(paramPtr)
{
}

GapData::~GapData()
{
  for(std::vector<Machine *>::iterator it = _machinePts.begin () ; it != _machinePts.end () ; it++) 
    {
      delete *it;      
    }
  
  for(std::vector<Job *>::iterator it = _jobPts.begin () ; it != _jobPts.end () ; it++) 
    {
      delete *it;      
    }
}

void GapData::readData(const  std::string & inputFileName)
{
  // read input data from a file whose name is inputFileName

  std::ifstream inps(inputFileName.c_str(), std::ios::in);
  if (!inps)
    {
      std::cout << "cannot find input file" << std::endl;
      exit(EXIT_FAILURE);
    }
  if (inps.eof())
    {
      std::cout << "empty input file" << std::endl;
      exit(EXIT_FAILURE);
    }

//  check(!inps,"cannot find input file"); /* make sure that the input file exists */
//  check(inps.eof(),"empty input file"); /* make sure that the input file is not empty */

  int machineIndex, jobIndex;
  int nbMachines, nbJobs;

  if (!(inps >> nbMachines))
    printf("\nerror: could not read nbMachines\n");

  if (!(inps >> nbJobs))
    printf("\nerror: could not read nbJobs\n");

  std::vector< std::vector<int> > resConsumption(nbJobs, std::vector<int>(nbMachines)); 
  std::vector< std::vector<Double> > cost(nbJobs, std::vector<Double>(nbMachines)); 

  for (machineIndex = 0; machineIndex < nbMachines; machineIndex++)
    {
      for (jobIndex = 0; jobIndex < nbJobs; jobIndex++)
	{
	  if (!(inps >> cost[jobIndex][machineIndex]))
	    printf("\nerror: could not read the cost for job %d and machine %d\n",jobIndex+1,machineIndex+1);
	}
    }
  for (machineIndex = 0; machineIndex < nbMachines; machineIndex++)
    {
      for (jobIndex = 0; jobIndex < nbJobs; jobIndex++)
	{
	  if (!(inps >> resConsumption[jobIndex][machineIndex]))
	    printf("\nerror: could not read the resource consumption for job %d and machine %d\n",jobIndex+1,machineIndex+1);
	}
    }
  int cap, maxCap = 0;
  for (machineIndex = 0; machineIndex < nbMachines; machineIndex++)
    {
      if (!(inps >> cap))
	printf("\nerror: could not read the capacity for machine %d\n",machineIndex+1);
      _machinePts.push_back( new Machine(machineIndex, cap) );
      if (cap > maxCap)
	maxCap = cap;
    }

  for (jobIndex = 0; jobIndex < nbJobs; jobIndex++)
    {
      _jobPts.push_back( new Job(jobIndex, nbMachines, resConsumption[jobIndex], cost[jobIndex]) );
    }
}

std::ostream&  GapData::print(std::ostream& os) const
{
  os << "GapData" << std::endl;
  os << "   modelName: "  << _modelName << std::endl;
  os << "   vector<Machine>"   << std::endl;
  for (int machineIndex = 0; machineIndex < machinePts().size(); ++machineIndex)
    { os << *machinePts()[machineIndex] << std::endl;}
  os << "   vector<Job>"   << std::endl;
  for (int jobIndex = 0; jobIndex < jobPts().size(); ++jobIndex)
    { os << *jobPts()[jobIndex] << std::endl;}
  return os;
}

void GapData::generateAtRandom(const int & genSeed)
{
  // int nbJobs(applicationSpecificParamPtr()->nbJobs());
  // double tightness(applicationSpecificParamPtr()->tightnessFactor())
  int nbJobs;
  double tightness;
  double machineMaxCapacity;
  double maxJobCapUsage;
  double minJobCapUsage;

  if (_paramPtr != NULL)
    {
      nbJobs = _paramPtr->nbJobs();
      tightness = _paramPtr->tightness();
      machineMaxCapacity = _paramPtr->machineMaxCapacity();
      maxJobCapUsage = _paramPtr->maxJobCapUsage();
      minJobCapUsage = _paramPtr->minJobCapUsage();
    }
  else
    {
      nbJobs = 20;
      tightness = 0.95;
      machineMaxCapacity = 10000;
      maxJobCapUsage = 2500;
      minJobCapUsage = 500;
    }

  RandGen jobCapUsageGen((int)(maxJobCapUsage-minJobCapUsage),0, (1 + genSeed) * 17 + 51);

  string outputFileName("randomDSGAPdata");
  outputFileName = outputFileName + genSeed;
  ofstream os(outputFileName.c_str(),std::ios::out);

  _modelName = "RandInst" ;
  _modelName = _modelName + genSeed;
  double totalJobCapUsage(0);
  double capUsage(0);
  for(int jobIndex = 0; jobIndex < nbJobs ; ++jobIndex)
    {
      std::string jobName("JOB");
      jobName = jobName + jobIndex;
      capUsage =  minJobCapUsage + jobCapUsageGen();
      _jobPts.push_back(new Job(jobIndex, capUsage));
      std::cout << "Job " <<  jobIndex << " capUsage = " << capUsage << std::endl;
      totalJobCapUsage += capUsage;
    }


  //int nbMachines = 1 + (int) ((tightness * totalJobCapUsage) / machineMaxCapacity);

  int nbMachines = 0;
  double totalMachineCapacity = 0;

  int machineIndex = 0;
  do 
    {
      ++nbMachines;
      _machinePts.push_back( new  Machine(machineIndex,  machineMaxCapacity) ); //
      totalMachineCapacity += machineMaxCapacity;
      std::cout << "Machine " <<  machineIndex << " machineMaxCapacity = " << machineMaxCapacity << " totalMachineCapacity = " << totalMachineCapacity  << " totalJobCapUsage = " << totalJobCapUsage << std::endl;
      machineMaxCapacity *= tightness;
      machineMaxCapacity = Dceil(machineMaxCapacity);
      ++machineIndex;
    } while (totalJobCapUsage > tightness * totalMachineCapacity);

  for(vector<Job *>::iterator ptr = _jobPts.begin();
      ptr != _jobPts.end(); ptr++)
    (*ptr)->randomWeightGeneration(nbMachines);

  return;
}

const std::vector<Job*>& GapData::jobPts() const
{
  return _jobPts;
}

const std::vector<Machine*>& GapData::machinePts() const
{
  return _machinePts;
}

const std::string& GapData::modelName() const
{
  return _modelName;
}

ApplicationSpecificParam* GapData::paramPtr() const
{
  return _paramPtr;
}
