/**
 *
 * This file bcModelFormulationC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELformulation_H_
#define BCMODELformulation_H_

#include "bcUsefulHeadFil.hpp"

#include "bcModelPointerC.hpp"
#include "bcMultiIndexC.hpp"
#include "bcModelSolutionC.hpp"

class ProbConfig;
class BcSolverOracleFunctor;
class BcFracSolBasedHeuristicFunctor;
class BcMasterHeuristicFunctor;
class BcDivingFixingFunctor;
class BcEnumSolBasedHeuristicFunctor;
class Model;
class BcFormIndex;
class BcSolution;
class ProbConfig;
class BcRyanAndFosterBranchConstr;
class BcSoftConflictsCut;
class BcCustomNonLinearCut;
class BcExtendedArcCut;
class BcNetwork;
class BcArcInfo;

#ifdef BCP_RCSP_IS_FOUND
namespace bcp_rcsp
{
    struct GraphData;
}
#endif

class BcFormulation
{
  friend class BcColGenSpArray;
  friend class BcMasterArray;
  friend class BcFormulationArray;
  friend class BcVarArray;
  friend class BcConstrArray;

protected:
  ProbConfig * _probConfPtr;
public:
  BcFormulation(ProbConfig * pConfPointer = NULL);
  BcFormulation(const BcFormulation & that);
  BcFormulation(BcModel & model, const std::string & name = "directForm");
  virtual ~BcFormulation();

  BcSolution solve(bool printForm = false, bool printOutput = false);
  SolutionStatus getStatus();
  void setFixedCost(const double & value);

  ProbConfig * probConfPtr() const;

  bool isDefined() const;
  bool isMaster() const;
  bool isColGenSp() const;
  BcFormulation master() const;
  const std::list< BcFormulation > & colGenSubProblemList() const;
  void initializeWithColumns(BcSolution & sol, bool activeColumns = true);
  void initializeWithSolution(BcSolution & sol);
  const BcFormulation & operator+=(BcSolution & sol); /// for backwards compatibility

  void resetObjective(); /// set objective to zero, useful before modifying objective between two calls of solve()
  void resetConstrRHS(); /// set objective to zero, useful before modifying RHS of constraints between two calls of solve()
  void update(); /// apply to the formulation the latest modifications (variable bounds, objective coefficients,
                 /// and constraint right-hand sides) that have been done since previous call to update() or solve()

  /// record and append a solution
  virtual void priorityLevel(const double & prLevel);
  virtual const MultiIndex & id() const;
  virtual const std::string & name() const;
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  const BcFormulation & attach(BcSolverOracleFunctor * oraclePtr);
  const BcFormulation & attach(BcFracSolBasedHeuristicFunctor * functorPtr);
  const BcFormulation & attach(BcMasterHeuristicFunctor * functorPtr);
  const BcFormulation & attach(BcDivingFixingFunctor * functorPtr);
  const BcFormulation & attach(BcEnumSolBasedHeuristicFunctor * functorPtr);
  virtual const BcFormulation & operator<=(const double & ubOnNbOfUsedSol);
  virtual const BcFormulation & operator>=(const double & lbOnNbOfUsedSol);
  virtual const BcFormulation & operator==(const double & nbOfUsedSol);
  friend std::ostream& operator<<(std::ostream& os, const BcFormulation & that);

  double zeroReducedCostThreshold() const;
  bool rollbackPointSavedStatus() const;
  bool enumeratedStatus() const;

  const BcNetwork network();

#ifdef BCP_RCSP_IS_FOUND
  const bcp_rcsp::GraphData * RCSPGraph();
#endif

  void getCustomNonLinearActiveCutsList(std::list<BcCustomNonLinearCut> & cutsList) const;
  void getRyanAndFosterBranchingConstrsList(std::list<BcRyanAndFosterBranchConstr> & ryanAndFosterBranchConstrList) const;
  void getActiveSoftConflictCutsList(std::list<BcSoftConflictsCut> & softConflictsCutsList) const;
  void getColumnsInPrimalLpSolution(std::vector<std::pair<double, BcSolution>> & columnsInSol) const;

  bool currentNodeIsRoot();
  bool debugSolutionIsValidAtThisNode();
};


class BcFormulationArray
{
 protected:
  Model * _modelPtr;
  std::string _name;
  BcFormulation _curForm; /// cache for operator()

 public:
  BcFormulationArray(BcModel & modelPointer, std::string  name = "directForm");
  virtual BcFormIndex operator[](const int & index);
  virtual BcFormulation & getElement(const MultiIndex & indexArray);
  virtual BcFormulation & createElement(const MultiIndex & indexArray);
  inline BcFormulation & createElement(int firstIndex = 0, int secondIndex = -1, int thirdIndex = -1,
                                       int fourthIndex = -1, int fifthIndex = -1, int sixthIndex = -1,
                                       int seventhIndex = -1)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }

    inline BcFormulation & operator()()
    {
        return createElement(MultiIndex(0));
    }

    inline BcFormulation & operator()(int firstIndex)
  {
    return createElement(MultiIndex(firstIndex));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex, int thirdIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex, int thirdIndex, int fourthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex ));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex, int thirdIndex, int fourthIndex, int fifthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
				    fifthIndex));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex, int thirdIndex, int fourthIndex,
                                    int fifthIndex, int sixthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex));
  }

  inline BcFormulation & operator()(int firstIndex, int secondIndex ,
                                    int thirdIndex , int fourthIndex ,
                                    int fifthIndex , int sixthIndex ,
                                    int seventhIndex )
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }
  inline BcFormulation & operator()(const MultiIndex & indexArray)
  {
    return createElement(indexArray);
  }

  virtual ~BcFormulationArray();
};


class BcFormIndex: public std::pair < BcFormulationArray * , MultiIndex >
{
 public:
  BcFormIndex(BcFormulationArray * modFormPtr, const MultiIndex & indexArray);
  BcFormIndex(BcFormulationArray * modFormPtr);
  virtual ~BcFormIndex();
  BcFormIndex & operator[](const int & index);
  virtual const BcFormulation & attach(BcFracSolBasedHeuristicFunctor * functorPtr);
  virtual const BcFormulation & attach(BcMasterHeuristicFunctor * functorPtr);
  virtual const BcFormulation & attach(BcDivingFixingFunctor * functorPtr);
  virtual const BcFormulation & attach(BcSolverOracleFunctor * oraclePtr);
  virtual const BcFormulation & operator<<(BcMasterHeuristicFunctor * functorPtr);
  virtual const BcFormulation & operator<<(BcFracSolBasedHeuristicFunctor * functorPtr);
  virtual const BcFormulation & operator<<(BcDivingFixingFunctor * functorPtr);
  virtual const BcFormulation & operator<<(BcSolverOracleFunctor * oraclePtr);
  virtual const BcFormulation & operator<=(const double & ubOnNbOfUsedSol);
  virtual const BcFormulation & operator>=(const double & lbOnNbOfUsedSol);
  virtual const BcFormulation & operator==(const double & nbOfUsedSol);
  virtual operator BcFormulation ();
};

#define IN_ // Indicates what is used for input
#define OUT_ // Indicates what is used for output

class BcMasterHeuristicFunctor
{
  public:
   BcMasterHeuristicFunctor(){}
   virtual ~BcMasterHeuristicFunctor(){}

  virtual bool operator() (IN_ BcFormulation spPtr,
                           IN_ std::vector<std::pair<BcSolution, double> > & colsInFixedSolution,
                           OUT_ BcSolution & primalSol);
};

struct BcFracSolBasedHeuristicFunctorInput
{
  BcFormulation spPtr;
  std::vector<std::pair<BcVar, double> > pureMastVarsInFixedSolution;
  std::vector<std::pair<BcSolution, double> >    colsInFixedSolution;

  std::vector<std::pair<BcVar, double> > pureMastVarsInMasterSolution;
  std::vector<std::pair<BcSolution, double> >    colsInMasterSolution;
};

class BcFracSolBasedHeuristicFunctor
{
  public:
   BcFracSolBasedHeuristicFunctor(){}
   virtual ~BcFracSolBasedHeuristicFunctor(){}

  virtual bool operator() (IN_ const BcFracSolBasedHeuristicFunctorInput & in,
                           OUT_ BcSolution & primalSol);
};


struct BcDivingFixingFunctorInput
{
  BcFormulation spPtr;
  std::vector<std::pair<BcVar, double> > pureMastVarsInFixedSolution;
  std::vector<std::pair<BcSolution, double> >    colsInFixedSolution;

  std::vector<std::pair<BcVar, double> > pureMastVarsInMasterSolution; // Only those not in tabou list
  std::vector<std::pair<BcSolution, double> >    colsInMasterSolution; // Only those not in tabou list

  std::vector<std::pair<BcVar, double> > pureMastVarsInMasterSolAndTabouList;
  std::vector<std::pair<BcSolution, double> >    colsInMasterSolAndTabouList;
};

struct BcDivingFixingFunctorOutput
{
  BcFormulation spPtr;
  std::vector<std::pair<BcVar, int> > pureMastVarsToFix;
  std::vector<std::pair<int, int> >           colsToFix;
};

class BcDivingFixingFunctor
{
  public:

  BcDivingFixingFunctor(){}
  virtual ~BcDivingFixingFunctor(){}

  void colCutGenTerminated(IN_ BcFormulation spPtr,
                           IN_ std::vector<std::pair<BcSolution, double> > & colsInFixedSolution,
                           IN_ std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                           OUT_ BcSolution & primalSol);

  virtual void operator() (const BcDivingFixingFunctorInput & in, BcDivingFixingFunctorOutput & out);
};

class BcEnumSolBasedHeuristicFunctor
{
public:

    BcEnumSolBasedHeuristicFunctor(){}
    virtual ~BcEnumSolBasedHeuristicFunctor(){}

    virtual bool operator() (IN_ BcFormulation spPtr,
                             IN_ double incPrimalIpBound,
                             const std::vector<BcSolution> & enumSolution,
                             OUT_ std::vector<std::pair<int, double> > & solution,
                             OUT_ BcSolution & addSolution);
};

class BcSolverOracleInfo
{
  public:
  BcSolverOracleInfo(){}
  virtual ~BcSolverOracleInfo(){}
};

class BcSolverOracleFunctor // BcSolverOracleFunctor
{
 public:
  BcSolverOracleFunctor(){}
  virtual ~BcSolverOracleFunctor(){}

  /// should return false if solver cannot be prepared (BapCod quites with an error status in this case)
  virtual bool prepareSolver(){ return true;}

  /// return true if infeasible and false if feasible
  virtual bool setupNode(BcFormulation spPtr, const BcSolverOracleInfo * infoPtr);

  virtual BcSolverOracleInfo * recordSolverOracleInfo(const BcFormulation spPtr);

  /// returns -1 if not enumerated, otherwise returns the number of enumerated solutions
  virtual void reducedCostFixingAndEnumeration(IN_ BcFormulation spPtr,
                                               IN_ const int & enumerationMode,
                                               IN_ const double & threshold);
  /// if maxNumberOfSolutions is positive, then at most this number of enum. solutions with smallest reduced cost
  /// will be returned in enumeratedSol, reduced costs of them will be returned in reducedCosts
  /// if maxNumberOfSolutions is not positive, all enum. solutions will be returned without calculating reduced costs
  virtual bool getEnumeratedStatus() const;
  virtual void getEnumeratedSolutions(IN_ BcFormulation spPtr,
                                      IN_ const int & maxNumberOfSolutions,
                                      OUT_ BcSolution & enumeratedSol,
                                      OUT_ std::vector<double> & reducedCosts);
  virtual void checkEnumeratedSolutions(BcFormulation spPtr, const std::vector<Solution *> & solPts,
                                        std::vector<bool> & solIsEnumerated);
  virtual bool isProperSolution(IN_ BcFormulation spPtr,
                                IN_ BcSolution & solution);
  virtual bool solSatisfiesCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                               IN_ const BcSolution & solution);
  virtual void columnGenerationTerminated(IN_ BcFormulation spPtr, bool afterRedCostFixing, int nodeOrder,
                                          int nodeDepth, int cutSepRound, double dualBound, double elapsedTime,
                                          bool masterConverged);
  virtual bool improveCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                          IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                          IN_ const bool & masterConverged);
  virtual bool drawPrimalSolutionToDotFile(IN_ BcFormulation spPtr,
                                           IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                           IN_ const std::string & filename);
  virtual bool lightenCurrentSpRelaxation(IN_ BcFormulation spPtr, const int & masterConverged, const int & callMode);

  virtual bool operator() (IN_ BcFormulation spPtr,
                           IN_ int colGenPhase,
                           OUT_ double & objVal,
                           OUT_ double & dualBound,
                           OUT_ BcSolution & primalSol);

  /// left for backward compatibility
    virtual bool operator() (IN_ BcFormulation spPtr,
                             OUT_ double & objVal,
                             IN_ OUT_ double & primalBound,
                             IN_ OUT_ double & dualBound,
                             OUT_ BcSolution & primalSol,
                             OUT_ BcDualSolution & dualSol,
                             IN_ OUT_ int & phaseOfStageApproach);

   virtual int getMessageIdToCutGeneration() const;
   virtual int getNumberOfEnumeratedSolutions() const;
   virtual bool getDebugSolution(BcFormulation spPtr, BcSolution & primalSol);
   virtual bool setDebugSolution(const std::vector<std::vector<int> > & ids, bool vertexBased);
};


inline const BcFormulation & BcFormIndex::attach(BcFracSolBasedHeuristicFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::attach(BcMasterHeuristicFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::attach(BcDivingFixingFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::attach(BcSolverOracleFunctor * oraclePtr)
{
  return (first->operator()(second)).attach(oraclePtr);
}

inline const BcFormulation & BcFormIndex::operator<<(BcFracSolBasedHeuristicFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::operator<<(BcMasterHeuristicFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::operator<<(BcDivingFixingFunctor * functorPtr)
{
  return (first->operator()(second)).attach(functorPtr);
}

inline const BcFormulation & BcFormIndex::operator<<(BcSolverOracleFunctor * oraclePtr)
{
  return (first->operator()(second)).attach(oraclePtr);
}

inline const BcFormulation & BcFormIndex::operator<=(const double & ubOnNbOfUsedSol)
{
  return (first->operator()(second)).operator<=(ubOnNbOfUsedSol);
}

inline const BcFormulation & BcFormIndex::operator>=(const double & lbOnNbOfUsedSol)
{
  return (first->operator()(second)).operator>=(lbOnNbOfUsedSol);
}

inline const BcFormulation & BcFormIndex::operator==(const double & nbOfUsedSol)
{
  return (first->operator()(second)).operator==(nbOfUsedSol);
}

inline BcFormIndex::operator BcFormulation()
{
  return first->operator()(second);
}



#endif /* BCMODELformulation_H_ */
