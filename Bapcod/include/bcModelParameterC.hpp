/**
 *
 * This file bcModelParameterC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCAPPLICATIONPARAMETERSC_H_
#define BCAPPLICATIONPARAMETERSC_H_

#include "bcParameterC.hpp"
#include "bcDoubleC.hpp"

class VarConstr;
class Variable;
struct SolutionVarInfo;

/*** Global constants  ***/

const int ProgStatisticsUndefinedStatus = -1;
const double BapcodInfinity = 1e12;
const double BapcodClockType = 0; /// set 0 for wall time, 1 for user + system time


/**
 * This class describes an application parameter inside a
 * configuration file read by the class ProgramParameter.
 */
template<typename T, typename U=T>
class ApplicationParameter : public Parameter<T,U>
{
private:
  /**
   * The name the parameter to be found inside the configuration file
   */
  std::string _name;

  /**
   * The variable to stock the value of the parameter
   */
  T _value;

  /**
   * The defaultValue of the parameter.
   */
  T _defaultValue;

  /**
   * The description of the parameter.
   */
  std::string _description;

public:
  /**
   * Add a parameter to be read inside the file parameterFileName.
   *
   * @param name: parameter name.
   * @param defaultValue: defaultValue of the parameter.
   * @param description: description of the parameter.
   * @param configOption: is it a config Option?
   * @param hideOption: does this parameter is hidden when the application is launch with help command?
   */
  ApplicationParameter(std::string name, const T& defaultValue,
      std::string description = "Undocumented",
      bool configOption = true,
      bool hideOption = false) :
        _name(name), _value(defaultValue),
        _defaultValue(defaultValue),
        _description(description)
  {
  }

  /**
   * Add a parameter without default value.
   * Use this constructor if the parameters is a multi-token
   * (can have multiple values).
   *
   * @param name: parameter's name.
   * @param description: parameter's description.
   */
  ApplicationParameter(const std::string & name,
      const std::string & description) :
        _name(name),
        _description(description)
  {
  }

  virtual ~ApplicationParameter() {}

  inline const std::string& name() const {
    return _name;
  }

  inline const T& value() const {
    return _value;
  }

  inline T& value() {
    return _value;
  }

  inline T* valuePtr() {
    return &_value;
  }

  inline const T& defaultValue() const {
    return _defaultValue;
  }

  inline const std::string& description() const {
    return _description;
  }

  inline const T& operator()() const {
    return value();
  }

  inline T& operator()() {
    return value();
  }

  inline void operator()(const T& val) {
    _value = val;
  }

  inline void operator=(const T& val)
  {
    operator()(val);
  }

  inline operator T() {
    return value();
  }
};

template<typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const ApplicationParameter<T,U>& option)
{
  os << option.name() << " = " << option.value();
  return os;
}

template<>
inline std::ostream& operator<<(std::ostream& os, const ApplicationParameter<bool>& option)
{
  os << option.name() << " = " << std::boolalpha << option.value();
  return os;
}

/**
 * This class describes some other kind of parameters (Multitoken or some
 * specific classes that already inherits from Parameter).
 */
template<typename T, typename D>
class ApplicationAdvancedParameter
{
private:
  /**
   * The name the parameter to be found inside the configuration file
   */
  std::string _name;

  /**
   * The variable to stock the value of the parameter
   */
  T _value;

  /**
   * The defaultValue of the parameter.
   */
  D _defaultValue;

  /**
   * The description of the parameter.
   */
  std::string _description;

public:
  /**
   * Add a parameter to be read inside the file parameterFileName.
   *
   * @param name: parameter name.
   * @param defaultValue: defaultValue of the parameter.
   * @param description: description of the parameter.
   * @param configOption: is it a config Option?
   * @param hideOption: does this parameter is hidden when the application is launch with help command?
   */
  ApplicationAdvancedParameter(std::string name, const D& defaultValue,
      std::string description,
      bool configOption = true,
      bool hideOption = false) :
        _name(name), _value(defaultValue),
        _defaultValue(defaultValue),
        _description(description)
  {
  }

  /**
   * Add a parameter without default value.
   * Use this constructor if the parameters is a multi-token
   * (can have multiple values).
   *
   * @param name: parameter's name.
   * @param description: parameter's description.
   */
  ApplicationAdvancedParameter(const std::string & name,
      const std::string & description) :
        _name(name),
        _description(description)
  {
  }

  virtual ~ApplicationAdvancedParameter() {}

  inline const std::string& name() const {
    return _name;
  }

  inline const T& value() const {
    return _value;
  }

  inline T& value() {
    return _value;
  }

  inline T* valuePtr() {
    return &_value;
  }

  inline const D& defaultValue() const {
    return _defaultValue;
  }

  inline const std::string& description() const {
    return _description;
  }

  inline const T& operator()() const {
    return value();
  }

  inline T& operator()() {
    return value();
  }

  inline void operator()(const T& val) {
    _value = val;
  }

  inline void operator=(const T& val)
  {
    operator()(val);
  }

  inline operator T() {
    return value();
  }

};

template<typename T, typename D>
inline std::ostream& operator<<(std::ostream& os,
    const ApplicationAdvancedParameter<T, D>& option)
{
  os << option.name() << " = " << option.value();
  return os;
}

template<bool, typename D>
inline std::ostream& operator<<(std::ostream& os,
    const ApplicationAdvancedParameter<bool, D>& option)
{
  os << option.name() << " = " << std::boolalpha << option.value();
  return os;
}

class PricingSolverCutsMessage
{
public:
  enum PricingSolverStatusMenum
  {
    noMessage = 0,
    stopCutGeneration = 1,
    doCutsRollback = 2,
    interruptSolution = 3
  };
};

class SolutionMethod : public Parameter<SolutionMethod>
{
public:
  enum SMenum
  {
    undefined = -1,
    none = 0,
    lpSolver = 1,
    mipSolver = 2,
    customSolver = 3,
    custom2mipSolver = 4,
  }; // If a new method is added, the check in set() method below should be corrected.
private:
  SMenum _status;
public:
  SolutionMethod(const int & stat = 0);
  virtual ~SolutionMethod(){}
  inline bool set(int status)
  {
    /**
     * Don't forget to this bounds (none=-1, kelleyInit2Bundle=9),
     * if a new method is added in SMenum
     */
    if (status >= -1 && status <= 10)
      {
        _status = static_cast<SMenum>(status);
        return true;
      }
    else
      {
        _status = undefined;
        return false;
      }
  }
  const SMenum & status() const;
  int getStatusAsInteger() const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream & operator<<(std::ostream& os, const SolutionMethod & that)
{
  return that.print(os);
}

class ColGenProximalStabilizationMode : public Parameter<ColGenProximalStabilizationMode>
{
public:
  enum ColGenProximalStabilizationModeMenum
  {
    undefined = -1,
    curvatureMode = 0,
    explicitMode = 1,
    multiPointMode = 2
   };
private:
  ColGenProximalStabilizationModeMenum _status;
public:
  ColGenProximalStabilizationMode(const int & stat = 0);
  virtual ~ColGenProximalStabilizationMode(){}
  const ColGenProximalStabilizationModeMenum & status() const;
  int getStatusAsInteger() const;
  virtual std::ostream&
  print(std::ostream& os = std::cout) const;
  bool set(int);
};

inline std::ostream& operator<<(std::ostream& os, const ColGenProximalStabilizationMode & that)
{
  return that.print(os);
}

class StabilizationFunctionType : public Parameter<StabilizationFunctionType>
{
public:
  enum StabilizationFunctionTypeMenum
  {
    undefined = -1,
    none = 0,
    boxStep = 1,
    threePiece = 2,
    fivePiece = 3,
    bundle = 4
  };
private:
  StabilizationFunctionTypeMenum _status;
public:
  StabilizationFunctionType(const int & stat = 0);
  virtual ~StabilizationFunctionType(){}
  const StabilizationFunctionTypeMenum & status() const;
  int getStatusAsInteger() const;
  virtual std::ostream&
  print(std::ostream& os = std::cout) const;
  bool set(int);
};

inline std::ostream& operator<<(std::ostream& os, const StabilizationFunctionType & that)
{
  return that.print(os);
}

class SmoothingMode : public Parameter<SmoothingMode>
{
public:
  enum SmoothingMenum
  {
    undefined = -1, none = 0, Wentges= 1, Neame = 2, ActiveColDs = 3
  };
private:
  SmoothingMenum _status;
public:
  SmoothingMode(const int & stat = 0);
  virtual ~SmoothingMode() {}
  const SmoothingMenum & status() const;
  int getStatusAsInteger() const;
  virtual std::ostream& print(std::ostream& os = std::cout) const;
  bool set(int);
};

inline std::ostream& operator<<(std::ostream& os, const SmoothingMode & that)
{
  return that.print(os);
}

class MasterInitMode : public Parameter<MasterInitMode>
{
public:
  enum MIMenum
  {
    defaultInit = -1,
    noArtCol = 0,
    globalArtCol = 1,
    subProbArtCol = 2,
    localArtCol = 3,
    incSolCol = 4,
    incSolColAndGac = 5,
    incSolColAndLac = 6,
    localAndGlobAc = 7
  };
private:
  MIMenum _status;
public:
  MasterInitMode(const int & stat = -1);
  const MIMenum & status() const
  {
    return _status;
  }
  bool set(const int & stat);
  virtual ~MasterInitMode(){}
  virtual std::ostream& print(std::ostream& os = std::cout) const;
};

inline std::ostream&
operator<<(std::ostream& os, const MasterInitMode & that)
{
  return that.print(os);
}

class SolutionStatus : public MultitokenParameter<int, std::set<int>, SolutionStatus>
{
public:
  enum
  {
    Undefined = -1,
    Optimum = 0,
    Infeasible = 1,
    Unbounded = 2,
    UnSolved = 3,
    PrimalFeasSolFound = 4,
    DualFeasSolFound = 5,
    OptimumUnscalInfeas = 6
  };
public:
  SolutionStatus();
  SolutionStatus(const int & stat);
  SolutionStatus(const int & stat1, const int & stat2);
  SolutionStatus(const SolutionStatus & stat);
  SolutionStatus(const std::vector<int> & statVect);
  SolutionStatus(const int * first, const int * last);
  virtual ~SolutionStatus(){}
  virtual std::pair<std::set<int>::iterator, bool> insert(const int & stat);
  template<typename InputIter>
  inline void insert(InputIter first, InputIter last)
  {
    int stat;
    while (first != last)
      {
        stat = *first++;
        std::set<int>::insert(stat);
      }
  }

  virtual SolutionStatus & operator=(const SolutionStatus & right);
  virtual SolutionStatus & operator=(const int & stat);
  virtual bool intersects(const SolutionStatus & right) const;
  virtual std::string stat2string(const int & stat) const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;
private:
  void validate_one(std::string token);
};

class SelectionStrategy : public Parameter<SelectionStrategy>
{
public:
  enum PriorityEnum
  {
    Undefined = -1, 
    NotConsideredForSelection = 0, // Both Branching and Cut Separation Priority
    FirstFound = 1, // Both Branching and Cut Separation Priority
    HighestPriority = 2, // Both Branching and Cut Separation Priority
    MostFractional = 3, // Branching Priority
    LeastFractional = 4, // Branching Priority
    FracWeightedPriority = 5, // Branching Priority is weighted by the fractional part
    Closest2RoundUp = 6, // Branching Priority
    Closest2RoundDown = 7, // Branching Priority
    GuidedSearch = 8, //
    LeastCost = 9, //
    LeastReducedCost = 10, //
    LeastGreedyCost = 11, //
    LeastSteepestEdgeCost = 12, //
    LeastPseudoCost = 13, //
    LeastInfeasibility = 15, //
    MostViolated = 16, //
    NotConsideredForIntegralityCheck = 17
  };
private:
  PriorityEnum _selectedRule;
  Double _incumbentCriteria;
public:
  SelectionStrategy(const int & stat = -1);
  SelectionStrategy(const PriorityEnum & stat);
  virtual ~SelectionStrategy() {}
  const PriorityEnum & selectedRule() const
  {
    return _selectedRule;
  }
  bool set(const int & stat);
  void initializeIncumbent();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  const int computeCriteriaAndValToCompareToIncumbent(const SolutionVarInfo * varInfo, Double & challengerRoundedValue,
                                                      const int & printLevel = 5);
  /// returns 1 if challenger is better than incumbent, -1 if if challenger is worse than incumbent, and 0 otherwise. */
  const int computeCriteriaToCompareToIncumbent(const SolutionVarInfo * varInfo, const Double & RoundedValue,
                                                const int & printLevel = 5);
  /// returns 1 if challenger is better than incumbent, -1 if if challenger is worse than incumbent, and 0 otherwise. */
};

/// replacing a typedef
class SelectionStrategyVector : public std::vector<SelectionStrategy>
{
};

class MultitokenSelectionStrategyVector : public MultitokenParameter<SelectionStrategy, SelectionStrategyVector,
                                                                     MultitokenSelectionStrategyVector>
{
public:
  void validate_one(std::string token)
  {
    std::istringstream iss(token);
    int element = 0;
    iss >> element;
    SelectionStrategy s(element);
    push_back(s);
  }
};

inline std::ostream & operator<<(std::ostream& os, const SelectionStrategy & that)
{
  return that.print(os);
}

inline std::ostream & operator<<(std::ostream& os, const SolutionStatus & that)
{
  return that.print(os);
}


#endif /* BCAPPLICATIONPARAMETERSC_H_ */

