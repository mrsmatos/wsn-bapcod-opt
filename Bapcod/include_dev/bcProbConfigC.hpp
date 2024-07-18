/**
 *
 * This file bcProbConfigC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef ProbConfigClasses_h
#define ProbConfigClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcMastVarConstrC.hpp"
#include "bcMastColumnC.hpp"
#include "bcNodeC.hpp"
#include "bcProblemC.hpp"
#include "bcSolutionC.hpp"
#include "bcGenBranchingConstrC.hpp"
#include <unordered_map>

#if(BOOST_VERSION >= 104000)
#include <boost/unordered_set.hpp>
#else
#include <ext/hash_set>
#endif //BOOST_LIB_VERSION
class ColGenSpConf;
class OvfConf;
class GenericVarConstr;
class GenericVar;
class Model;
class GenericConstr;
class CompSetInstMastBranchConstr;
class NetworkFlow;
class BcFormulation;

#ifdef BCP_RCSP_IS_FOUND
namespace bcp_rcsp
{
    struct GraphData;
}
#endif /* BCP_RCSP_IS_FOUND */

/**
   @brief ProbConfig is a object that represents Mathematical Program along Side Solver Method Configuration
   @details It is the common base to special configuration such as Benders or Dantzig-Wolfe Master, Pricing or Sepration suproblems
*/
class ProbConfig
{
 public:
  enum ProbConfType {undefined = -1, ovf = 0, master = 1, colGenSp = 2, sideProb = 4};
 
 protected:
  ProbConfType _configType;
  Model * _modelPtr;
  std::string _genericName;
  int _ref;
  int _pcSolCount;
  int _pcNodeCount;
  IndexCell _id; 

  Bound _primalIncBound; /// the value of the best primal solution encountered so far through the whole BaB tree.
  Bound _dualIncBound; /// the value of the best dual bounded computed so far by collecting subtreebounds from the whole BaB tree.
  Double _cutOffValue; /// if one finds a dual bound above this cutoff value, the solver stops.

  /// Only defined for master: subProblems belong to MasterConf
  /// Original formulation belongs to master
  OvfConf * _ovfConfPtr;
  /// vector of col gen sp
  std::vector<ColGenSpConf *> _colGenSubProbConfPts;
  std::list<Problem *> _problemList; /// added by Ruslan: list of problems for the master and the subproblems

  NetworkFlow * _networkFlowPtr; /// added by Ruslan (needed for communication between the RCSP pricing solver
                                 ///                  and the rank-1 cuts separation)

#ifdef BCP_RCSP_IS_FOUND
  bcp_rcsp::GraphData * _rcspGraphPtr;
#endif

  std::unordered_map<InstanciatedVar *, int> _instVarToIdMap;
  std::vector<InstanciatedVar *> _instVarPts;
  std::vector<double> _instVarRedCosts;
  std::unordered_map<int, InstanciatedVar *> _resIdToVarIdMap;

  GenericVar * _defaultGenericVar;
  GenericConstr * _defaultGenericConstr;
  MultiIndex _multiIndex;
  MultiIndexNames _multiIndexNames;
  std::map< std::string , GenericVar * > _name2GenericVarPtrMap;
  std::map< std::string , GenericConstr * > _name2GenericConstrPtrMap;
  std::map< std::string , GenericCutConstr * > _name2GenericCutConstrPtrMap;
  std::map< std::string , GenericBranchingConstr * > _name2GenericBranchingConstrPtrMap;

  /**
   * To hold generated InstanciatedVar until calling prepareConfig()
   */
  std::list< InstanciatedVar *> _iVarPts;
  
  /**
   * To hold generated InstanciatedConstr  until calling prepareConfig()
   */
  std::list< InstanciatedConstr *> _iConstrPts;

  /// Variables do not belong to ProbConfig
  std::list< Variable *> _pcVarPtrList;

  /// Constraints do not belong to ProbConfig
  std::list< Constraint *> _pcConstrPtrList;
  
  /// Constraints belong to ProbConfig
  std::list< Constraint *> _pcDelayedConstrPtrList;

  std::set<GenericBranchingConstr *, DynamicGenConstrSort> _candidateBranchingGenericConstr;
  std::set<GenericCutConstr *, DynamicGenConstrSort> _candidateCutGenericConstr;
  
  /// To hold CS branching constraints
  ColClassesVector _treeOfColClasses;

  Bound _dualBoundContrib;
  Problem * _probPtr;
  bool _curFormIsInfeasible;
  Solution * _primalSolPtr;
  bool _isPrepared;

  std::list< BcFormulation > _colGenSubProbFormList;

 public:

  int ref() const {return _probPtr->ref();}
  const std::string name() const {return _probPtr->name();}
  int pcSolCount() {return _pcSolCount;}
  int pcNodeCount() {return _pcNodeCount;}
  void increasePCSolCount();
  void increasePCNodeCount();
  virtual const bool & isPrepared() const {return _isPrepared;}
  
  virtual void addVariablesToForm(){}

  ProbConfig(ProbConfType configType,
	         Model * modelPtr,
	         std::string genericName,
	         const IndexCell& id,
	         const Bound & primalIncBound,
	         const Bound & dualIncBound,
	         Problem * problemPtr);
  virtual ~ProbConfig();
  const ProbConfType & configType() const {return _configType;}
  Solution * primalSolutionPtr() const {return _primalSolPtr;}
  virtual MasterConf * mastConfPtr() const {return NULL;}
  const std::string & genericName() const {return _genericName;}

  Model * modelPtr() const {return  _modelPtr;}
  void ovfConfPtr(OvfConf * OVFptr){_ovfConfPtr = OVFptr;}
  OvfConf * ovfConfPtr(){return _ovfConfPtr;}

  GenericVar * defaultGenericVarPtr() const;
  GenericConstr * defaultGenericConstrPtr() const;
  virtual void insertPureVar(InstanciatedVar * ivarPtr){}
  virtual void insertPureConstr(InstanciatedConstr * iconstrPtr){}

  virtual GenericVar * getGenericVar(const std::string & name) const;
  virtual GenericConstr * getGenericConstr(const std::string & name) const;
  virtual GenericCutConstr * getGenericCutConstr(const std::string & name) const;
  virtual GenericBranchingConstr * getGenericBranchingConstr(const std::string & name) const;
  virtual Solution * getDissagregatedSolution(Solution * curSolPtr); /// project solution in original space if need be
  virtual Solution * getAggregatedSolution(Solution * curSolPtr);
  Solution * getSolution(const VarPtrSet & curSol) ;
  Solution * getSolution(const MasterVarSolution & varPtrList);
  const Double & cutOffValue() const {return _cutOffValue;} /// if one finds a dual bound above this custoff value, the solver stops.
  virtual void resetCutOffValue(const Bound & incumbentVal); /// if one finds a dual bound above this custoff value, the solver stops.
  
  virtual void insertGenericVar(GenericVar * gvPtr);
  virtual void insertGenericConstr(GenericConstr * gcPtr);
  virtual void insertGenericCutConstr(GenericCutConstr * gcPtr);
  virtual void insertGenericBranchingConstr(GenericBranchingConstr * gbcPtr);

  virtual InstanciatedConstr * checkConstraint4Insertion(Constraint * constrPtr, const int & insertionLevel = 2)
  {return NULL;}
  virtual InstanciatedConstr * checkConstraint4Insertion(InstanciatedConstr * iconstrPtr,const int & insertionLevel = 2)
  {return NULL;}

  virtual void createDefaultGenericVarConstr();
  virtual const Bound & primalIncBound() const {return _primalIncBound;}
  const Bound & dualIncBound() const {return _dualIncBound;}
  virtual bool updateDualIncBound(const Bound & newDualBound);
  virtual bool updatePrimalIncBound(const Bound & newIncVal);
  virtual void resetDualIncBound();
  virtual void resetPrimalIncBound();
  virtual const Bound & dualBoundContrib() const;

  virtual const bool enumeratedStatus() const;

  virtual void insertInstVar(InstanciatedVar * iVarPtr);
  virtual void insertInstConstr(InstanciatedConstr * iConstrPtr);
  virtual Constraint * castAndAddConstraint(Constraint * constrPtr, const bool & insertImmediately = false);
  virtual InstanciatedConstr * castAndAddConstraint(InstanciatedConstr * iconstrPtr, const bool & insertImmediately = false);
  virtual Variable * castAndAddVariable(Variable * varPtr, const bool & insertImmediately = false);
  virtual InstanciatedVar * castAndAddVariable(InstanciatedVar * ivarPtr, const bool & insertImmediately = false);
  virtual Problem * probPtr() const;
  virtual Solution * solvePC(bool showOutput); /// solves the Problem configuration according to the method selected in _solMode
  virtual void recordCurrentPrimalSol();

  const MultiIndex & multiIndex() const {return _multiIndex;}

  /**
   * Return true if preprocessing has detected that the current formulation admits not feasible solution, false otherwise
   */
  virtual const bool & curFormIsInfeasible() const 
  {
    return _curFormIsInfeasible;
  }

  virtual Solution * extractCurrentSol();

  /**
   * Only active for Master
   */
  virtual void updatePrimalIncSolution(Solution * solPtr){}
  virtual void recordInitialActiveSetOfColumns(Solution * activeSolPtr){}
  virtual void recordInitialInactiveSetOfColumns(Solution * activeSolPtr){}

  /**
   * Set initial variables and constraints
   */
  virtual void prepareProbConfig();

  /// Variables belong to ProbConfig
  virtual const std::list< Variable *> & pcVarPtrList() const 
  {
    return _pcVarPtrList; 
  }
  
  //added by Issam.
  virtual const std::list< InstanciatedVar *> & iVarPtrList() const 
  {
    return _iVarPts; 
  }
  
  //added by Issam.
  virtual const std::list< InstanciatedConstr *> & iConstrPtrList() const 
  {
    return _iConstrPts;
  }
  

  /// Constraints belong to ProbConfig
  virtual const std::list< Constraint *> & pcConstrPtrList() const 
  {
    return _pcConstrPtrList;
  }
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual std::ostream & printForm(std::ostream& os = std::cout) const;
  virtual void nicePrintAllConstraints(std::ostream& os = std::cout) const;

  virtual std::set<GenericCutConstr *, DynamicGenConstrSort> & candidateCutGenericConstr()
    {
      return _candidateCutGenericConstr;
    }

  virtual std::set<GenericBranchingConstr * , DynamicGenConstrSort> & candidateBranchingGenericConstr()
    {
    return _candidateBranchingGenericConstr;
    }

  /// To hold CS branching constraints
  virtual ColClassesVector & treeOfColClasses()
    {
      return _treeOfColClasses;
    }
  
  void sortTreeOfColClasses();

  virtual const std::vector< ColGenSpConf * > & colGenSubProbConfPts() const
    {
      return _colGenSubProbConfPts;
    }
  const std::list<Problem *> & problemList() const
    {
      return _problemList;
    }

  virtual ColGenSpConf * colGenSubProbConf(const int & index) const;

  virtual MastColumn * recordSubproblemSolution(Solution * spSolPtr,
                                                bool inPhaseOne,
                                                const int & insertionLevel = 2,
                                                Solution * masterSolPtr = NULL,
                                                bool changeEnumeratedFlag = false)
  {return NULL;}

  virtual const MultiIndexNames multiIndexNames()  const {return _multiIndexNames;}
  virtual void  multiIndexNames(const MultiIndexNames & mi)  {_multiIndexNames = mi;}
  virtual void upperBoundPtr(Double * ubPtr){}
  virtual void lowerBoundPtr(Double * lbPtr){}
  virtual Double * upperBoundPtr() const {return NULL;}
  virtual Double * lowerBoundPtr() const {return NULL;}
  virtual const IndexCell & id() const;

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;
  const ControlParameters& param() const;
  ControlParameters& param();

  const ProgStatus& progStatus() const;
  ProgStatus& progStatus();

  const std::list< BcFormulation > & colGenSubProbFormList() const {return  _colGenSubProbFormList;}

  virtual bool isTypeOf(const PcId::PcIdentifier & pcIdentifier) const
  {
    return false;
  }

  void setNetworkFlowPtr(NetworkFlow * networkFlowPtr)
  {
    _networkFlowPtr = networkFlowPtr;
  }
  
  NetworkFlow * networkFlowPtr()
  {
    return _networkFlowPtr;
  }

#ifdef BCP_RCSP_IS_FOUND
  const bcp_rcsp::GraphData * rcspGraphPtr() const
  {
      return _rcspGraphPtr;
  }
#endif

   const std::unordered_map<InstanciatedVar *, int> & instVarToIdMap() const
   {
      return _instVarToIdMap;
   }

   const std::vector<InstanciatedVar *> & instVarPts() const
   {
      return _instVarPts;
   }

    std::vector<double> & instVarRedCosts()
    {
        return _instVarRedCosts;
    }

    const std::unordered_map<int, InstanciatedVar *> & resIdToVarIdMap() const
   {
      return _resIdToVarIdMap;
   }


    const NetworkFlow * networkFlowPtr() const
  {
    return _networkFlowPtr;
  }
    
  void setMipRequiredStatus(const SolutionStatus & newStatus);

#ifdef BCP_RCSP_IS_FOUND
  bool fillRCSPGraph();
#endif //BCP_RCSP_IS_FOUND
};

inline std::ostream& operator<<(std::ostream& os, const ProbConfig & that)
{
  return that.print(os);
}

#endif // ProbConfigClasses_h
