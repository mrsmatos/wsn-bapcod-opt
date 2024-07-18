/**
 *
 * This file dsgapDataC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DgapDataClasses_h
#define DgapDataClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "dsgapParameter.hpp"

class Machine
{
  int _idNb;
  Double _capacity;

  //Machine(){} // should not be used
public:
  Machine(const int & idNb,
	  const Double & capacity);

  //  const int & idNb() const {return _idNb;}
  const Double & capacity()  const { return _capacity;}
  std::ostream& print(std::ostream& os = std::cout) const;
};


inline std::ostream& operator<<(std::ostream& os, const Machine & m)
{return m.print(os);}


class Job
{
  int _idNb;
  Double _nominalCapacityUsage;
  Double _maxCostOnMachine;
  std::vector<Double> _costOnMachine; 
  std::vector<Double> _capUsageOnMachine; 
  //Job(){} // should not be used
public:
  Job(const int & idNb,
      const Double & nominalCapacityUsage);
  Job(const int & idNb,
      const int & nbMachines,
      const std::vector<int> & resConsumption, 
      const std::vector<Double> & cost);
  ~Job();
  const Double & nominalCapUsage() const {return _nominalCapacityUsage;}
  const Double & maxCostOnMachine() const {return _maxCostOnMachine;}
  const Double & costOnMachine(const int & machineRef) const {return _costOnMachine[machineRef];}
  const Double & capUsageOnMachine(const int & machineRef) const {return _capUsageOnMachine[machineRef];}
  void randomWeightGeneration(const int & nbMachines);
  std::ostream& print(std::ostream& os = std::cout) const;
};



inline std::ostream& operator<<(std::ostream& os, const Job & m)
{return m.print(os);}


class GapData
{
  std::string _modelName;
  std::vector<Machine *> _machinePts; 
  std::vector<Job *> _jobPts; 
  ApplicationSpecificParam* _paramPtr;

public:

  GapData(ApplicationSpecificParam* paramPtr = NULL);
  virtual ~GapData();

  void readData(const std::string & inputFileName);
  void generateAtRandom(const int & genSeed);
//  std::ostream& print(std::ostream& os = std::cout);
  const std::vector<Job*>& jobPts() const;
  const std::vector<Machine*>& machinePts() const;
  const std::string& modelName() const;
  ApplicationSpecificParam* paramPtr() const;
  std::ostream& print(std::ostream& os = std::cout) const;
};


inline std::ostream& operator<<(std::ostream& os, const GapData & dStruct)
{return dStruct.print(os);}



struct GreedyCell
{
  int jobRef;
  int machineRef;
  Double ratio;
  GreedyCell(const int & j = -1, 
	     const int & m = -1, 
	     const Double & costOnMachine = 1, 
	     const Double & capUsageOnMachine = 1, 
	     const Double machineCap = 1):
    jobRef(j), 
    machineRef(m), 
    ratio(costOnMachine * capUsageOnMachine / machineCap)
  {}
  //bool operator<(const GreedyCell & that) {return (ratio < that.ratio);}
};

struct GreedyCellOrdering
{
  bool operator()(const GreedyCell & a, const GreedyCell & b) 
  {return (a.ratio < b.ratio);}
};


#endif // DgapDataClasses_h
