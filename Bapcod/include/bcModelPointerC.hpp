/**
 *
 * This file bcModelPointerC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef ModelPointerC_h
#define ModelPointerC_h

#include  "bcModelSolutionC.hpp"
#include  "bcModelObjectiveC.hpp"

#include "bcControlParameters.hpp"

class Solution;
class BcModel;
class BcMasterArray;
class BcColGenSpArray;
class BcObjectiveArray;
class BcObjective;
class BapcodInit;


class BcInitialisation // BcInitialisation
{
  friend class BcModel;

  BapcodInit * _bcpInitptr;
 public:
  BcInitialisation();
  explicit BcInitialisation(std::string filename, bool toPrintParameters = false, bool printVrpSolverHeader = false,
                            bool printBapcodHeader = true);

  BcInitialisation(int argc, char *argv[], std::string filename = "config/bcParameters.cfg",
				   bool printParameters = false, bool printVrpSolverHeader = false, bool printBapcodHeader = true);
  ~BcInitialisation();
  void bcReset();
  ControlParameters & param();
  void outputBaPCodStatistics(const std::string & instanceName, std::ostream & os = std::cout);
  void outputBaPCodStatistics(const std::string & instanceName, std::ostream & os, std::ostream & fs,
							  std::string comFileName = "");
  void addStatisticsValue(const std::string & statName, double value);
  void incrStatisticsCounter(const std::string & statName, int value = 1);
  double getStatisticValue(const std::string & statName);
  long getStatisticCounter(const std::string & statName);
  double getStatisticTime(const std::string & statName);

  void outputResultsForTests();

    /**
   * @return the input filename (this file contains some input data)
   */
  const std::string& instanceFile() const;
  const std::string& configFile() const; // added by Ruslan
  const std::string& statFile() const; // added by Ruslan
};

class BcSolutionFoundCallback
{
protected:
	bool _disagregateSolutionBeforeProcessing;	// added by Boris: if false, the non disaggregated solution is passed to the callback
public :
  BcSolutionFoundCallback() : _disagregateSolutionBeforeProcessing(true){}
  virtual ~BcSolutionFoundCallback(){}
  virtual bool operator()(BcSolution newSolution) const //returns false if solution is not feasible
  {
    return true;
  }
  bool disagregateSolutionBeforeProcessing() {return _disagregateSolutionBeforeProcessing;}
  void disagregateSolutionBeforeProcessing(bool b) {_disagregateSolutionBeforeProcessing=b;}
};

class BcModel  // BcModel
{
  friend class BcObjectiveArray;
  friend class BcMasterArray;
  friend class BcColGenSpArray;
  friend class BcVarArray;
  friend class BcConstrArray;
  friend class Model;
  Model * _modelPtr;
  bool _deleteInDestructor;

  inline void clearModel() { _modelPtr = NULL; }
public:
  BcModel(const BcInitialisation &  bapcodInitPtr,
	      const std::string & modelName = "Model",
	      const BcObjStatus::MinMaxIntFloat & objectiveSense = BcObjStatus::minInt);
  BcModel(Model * modelPtr);
  ~BcModel();

  BcFormulation master();
  void attach(BcSolutionFoundCallback * solutionFoundCallbackPtr);

  BcSolution solve();
  BcSolution enumerateAllColumns(int & nbEnumColumns);
  operator Model *() const;

  friend std::ostream& operator<<(std::ostream& os, const BcModel & that);

  /// returns triples (reduced cost, cost, solution pointer)
  void getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions);

};



#endif // ModelPointerC_h
