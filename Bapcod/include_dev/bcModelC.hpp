/**
 *
 * This file bcModelC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BcModelClasses_h
#define BcModelClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBapcodInit.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcNodeC.hpp"
#include "bcProbConfigC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcMasterConfC.hpp"
#include "bcOvfConfC.hpp"
#include "bcProblemC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcTimeC.hpp"
#include "bcVarConstrC.hpp"

class BcSolutionFoundCallback;
class BcModel;

class Model
{
  friend class BcObjectiveArray;
  friend class BcObjective;
  friend class BcColGenSpArray;
  friend class BcMasterArray;
  friend class BcModel;
  BapcodInit * _bapcodInitPtr;
  BcSolutionFoundCallback * _solutionFoundCallbackPtr;

  BcModel * _bcModelPtr;

  int _modelMasterCnt;
  int _modelOvfCnt;
  int _modelColGenSpCnt;
  int _modelConstrCnt;
  int _modelVarCnt;
  int _modelMastColCnt;
  int _modelBrConstrGenCnt;
  int _modelProblemCnt;
  int _modelGenVarConstrCnt;

  std::vector<ProbConfig *> _sideProbConfPts;
  MasterConf * _masterConfPtr;
  ColGenSpConf * _curColGenSpConfPtr;
  MapColGenSpConfPtrByNameAndMultiIndex _mapOfColGenSpConf;

  BcObjStatus::MinMaxIntFloat _objectiveSense;

  bool _modelIsSetup;

  PlainVarConstrPointerSet _garbageCollector;

 protected:

  Bound _initialPrimalBound;
  Bound _initialDualBound;

  std::string _modelName;
  Double _artVarCost;

  /**
   * Problem Specific Functions reqired for generic subroutines but requiring user input
   */
 public:
    int modelConstrCnt() const {return _modelConstrCnt;}
    int modelVarCnt() const {return _modelVarCnt;}
    int modelMastColCnt() const {return _modelMastColCnt;}
    int modelBrConstrGenCnt() const {return _modelBrConstrGenCnt;}
    int modelProblemCnt() const {return _modelProblemCnt;}
    int modelGenVarConstrCnt() const {return _modelGenVarConstrCnt;}
    void increaseModelConstrCnt();
    void increaseModelVarCnt();
    void increaseModelMastColCnt();
    void increaseModelBrConstrGenCnt();
    void increaseModelProblemCnt();
    void increaseModelGenVarConstrCnt();

    /**
   * Create Generic Variables and Constraints, Master, SubProblems, OVF or plain probConfig
   */
  void setup();
  virtual Solution * initPrimalHeur()
  {
    return NULL;
  }
  Model(BapcodInit *  bapcodInitPtr,
	    const std::string & modelName = "Model",
	    const BcObjStatus::MinMaxIntFloat & objectiveSense = BcObjStatus::minFloat);
  virtual ~Model();
  const BcObjStatus::MinMaxIntFloat & objectiveSense() const;
  void objectiveSense(const BcObjStatus::MinMaxIntFloat & newObjectiveSense);

  void setSolutionFoundCallback(BcSolutionFoundCallback * solutionFoundCallbackPtr)
  {
    _solutionFoundCallbackPtr = solutionFoundCallbackPtr;
  }

  bool checkIfSolutionIsFeasibleUsingCallback(Solution * solPtr) const;

  Solution * enumerateAllColumns(int & nbEnumColumns);
  Solution * solve();
  const std::string & name() const
  {
    return _modelName;
  }

  void setName(const std::string & name)
  {
    _modelName = name;
  }

  void mastConfPtr(MasterConf * MCptr);
  MasterConf * master() const;
  OvfConf * ovf() const;
  MasterConf * masterConfPtr() const {return master();}
  OvfConf * ovfConfPtr() const {return ovf();}

  bool ovfConfPtr(OvfConf * ptr)
  {
    if (_masterConfPtr != NULL)
	{
		_masterConfPtr->ovfConfPtr(ptr);
		return true;
	}
    return false;
  }

  const Double & artVarCost() const
  {
     return _artVarCost;
  }
  std::ostream & print(std::ostream& os = std::cout) const;
  std::ostream& printSol(std::ostream& os) const;
  const Bound & bestPrimalBound() const {return _masterConfPtr->bestPrimalBound();} //(Issam)
  const Bound & bestDualBound() const {return _masterConfPtr->bestDualBound();} //(Issam)
  void initialPrimalBound(const Bound & val) {_initialPrimalBound = val;} //(Issam)
  void initialDualBound(const Bound & val) {_initialDualBound = val;} //(Issam)
  const Bound & initialPrimalBound() const {return _initialPrimalBound;} //(Issam)
  const Bound & initialDualBound() const {return _initialDualBound;} //(Issam)

  BapcodInit & bapcodInit() const;
  BapcodInit * bapcodInitPtr() const;

  inline const ControlParameters& param() const
  {
    return _bapcodInitPtr->param();
  }

  inline ControlParameters& param()
  {
    return _bapcodInitPtr->param();
  }

  const ProgStatus& progStatus() const;
  ProgStatus& progStatus();
  PlainVarConstrPointerSet& garbageCollector();

  const BcModel & bcModel() const;
  BcModel & bcModel();

  void setArtCostValue(const Double & primalBound);

  OvfConf * createOvfConf(const std::string & name,
				          const IndexCell & id);

  ProbConfig * createProbConf(const std::string & name,
				              const IndexCell & id);

  MasterConf * createMasterConf(const std::string & name,
				                const IndexCell & id);

  ColGenSpConf * createColGenSubproblem(const std::string & spname = std::string("ColGenSp"),
						                const MultiIndex & id = MultiIndex(),
						                const bool & implicitlyFixCardinality = false,
						                const Double & upperBound = 1,
						                const Double & lowerBound = 0,
						                const Double & fixedCost = 0,
						                const Double & defaultDualVal = 1);

  GenericVar * createGenericVar(ProbConfig * probConfPtr,
					            const BcVarConstrType::BcVcType & vctype,
					            const std::string & name,
					            const MultiIndexNames & multiIndexNames = MultiIndexNames(),
					            const char & type = 'I',
					            const Double & cost = 0,
					            const Double & ub = BapcodInfinity,
					            const SelectionStrategy & branchingPriorityRule = SelectionStrategy::MostFractional,
					            const Double & genericBranchingOnAggregateVarPL = 10.0,
					            const Double & compBoundSetBranchingPL = 1.0,
					            const char & flag = 's', /// 's' for static, 'd' for dynamic
					            const char & sense = 'P',
                                int firstIndexMax = -1, int secondIndexMax = -1, int thirdIndexMax = -1);


    GenericConstr * createGenericConstr(ProbConfig * probConfPtr,
                                        const BcVarConstrType::BcVcType & vctype,
                                        const std::string & name,
                                        const MultiIndexNames & multiIndexNames = MultiIndexNames(),
                                        const char & sense = 'E',
                                        const Double & rhs = 0,
                                        const Double & val = 0,
                                        const bool & toBeUsedInPreprocessing = true,
                                        const char & flag = 's',
                                        const char & type = 'C',
                                        const char & kind = 'E',
                                        const SelectionStrategy & priority = SelectionStrategy::NotConsideredForSelection,
                                        const Double & priorityLevel = -1);

    /// type  'C' for core (required for the IP formulation,
    /// 'F' for facultative (only helpfull for the LP formulation),
    /// 'S' for constraints defining a subsystem in column generation for extended formulation approach
    GenericCutConstr* createGenericCut(ProbConfig * probConfigPtr,
                                       const std::string & name,
                                       const char & type,
                                       const SelectionStrategy & separationPriorityRule,
                                       const Double & nonRootPriorityLevel,
                                       const Double & rootPriorityLevel,
                                       const char & sense = 'G',
                                       const Double & rhs = 0,
                                       const bool & toBeUsedInPreprocessing = true);

  GenericCustomNonLinearCutConstr *
  createGenericCustomNonLinearCutConstr(BcCustomNonLinearCutArrayFunctor * functorPtr,
                                        ProbConfig * probConfigPtr,
                                        const std::string & name,
					                    const char & type,
					                    const SelectionStrategy & separationPriorityRule,
					                    const Double & nonRootPriorityLevel,
					                    const Double & rootPriorityLevel,
                                        const char & sense = 'G',
					                    const Double & rhs = 0);

	GenericSoftConflictsCutConstr *
	createGenericSoftConflictsCutConstr(BcSoftConflictsCutArrayFunctor * functorPtr,
									    ProbConfig * probConfigPtr,
									    const std::string & name,
									    const char & type,
									    const SelectionStrategy & separationPriorityRule,
									    const Double & nonRootPriorityLevel,
									    const Double & rootPriorityLevel,
									    const char & sense = 'L',
									    const Double & rhs = 0);

    GenericBranchingConstr *
    createGenericBranching(ProbConfig * probConfigPtr,
                           const std::string & name,
                           const char & type,
                           const SelectionStrategy & separationPriorityRule,
                           const Double & nonRootPriorityLevel = 1.0,
                           const Double & rootPriorityLevel = 1.0,
                           const char & sense = ' ',
                           const Double & rhs = 0,
                           const bool & toBeUsedInPreprocessing = true);

    Variable * recordSubproblemSolution(ColGenSpConf * cgSpConfPtr,
                                        Solution * solPtr,
                                        const int & insertionLevel = 2) const;

    InstanciatedConstr * createConstraint(ProbConfig * probConfPtr,
                                          GenericConstr * genConstrPtr,
                                          const MultiIndex & id,
                                          const Double & rhs,
                                          const char & sense,
                                          const Double & val,
                                          const std::string & name,
                                          const bool & toBeUsedInPreprocessing = true,
                                          const bool & considerAsEqualityInPreprocessing = false);

    InstanciatedConstr * createConstraint(ProbConfig * probConfPtr,
                                          GenericConstr * genConstrPtr,
                                          const MultiIndex & id,
                                          const Double & rhs,
                                          const char & sense,
                                          const Double & val);

    InstanciatedConstr * createConstraint(ProbConfig * probConfPtr,
                                          GenericConstr * genConstrPtr,
                                          const MultiIndex & id,
                                          const Double & rhs,
                                          const char & sense);

    InstanciatedConstr * createConstraint(ProbConfig * probConfPtr,
                                          GenericConstr * genConstrPtr,
                                          const MultiIndex & id,
                                          const Double & rhs);

    InstanciatedConstr * createConstraint(ProbConfig * probConfPtr,
                                          GenericConstr * genConstrPtr,
                                          const MultiIndex & id);

    InstanciatedVar * createVariable(ProbConfig * probConfPtr,
                                     GenericVar * genVarPtr,
                                     const MultiIndex & id,
                                     const Double & cost,
                                     const char & type,
                                     const std::string & name,
                                     const Double & ub = BapcodInfinity,
                                     const Double & lb = 0,
                                     const Double & priority = 1,
                                     const char & directive = 'U',
                                     const Double & globalUb = BapcodInfinity,
                                     const Double & globalLb = 0,
                                     const char & kind = 'E',
                                     const char & sense = 'P',
                                     const char & flag = 's');

    InstanciatedVar * createVariable(ProbConfig * probConfPtr,
                                     GenericVar * genVarPtr,
                                     const MultiIndex & id,
                                     const Double & cost,
                                     const char & type,
                                     const Double & ub,
                                     const Double & lb = 0,
                                     const Double & priority = 1,
                                     const char & directive = 'U',
                                     const Double & globalUb = BapcodInfinity,
                                     const Double & globalLb = 0,
                                     const char & kind = 'E',
                                     const char & sense = 'P',
                                     const char & flag = 's');

    InstanciatedVar * createVariable(ProbConfig * probConfPtr,
                                     GenericVar * genVarPtr,
                                     const MultiIndex & id,
                                     const Double & cost,
                                     const char & type);

    InstanciatedVar * createVariable(ProbConfig * probConfPtr,
                                     GenericVar * genVarPtr,
                                     const MultiIndex & id,
                                     const Double & cost);

    InstanciatedVar * createVariable(ProbConfig * probConfPtr,
                                     GenericVar * genVarPtr,
                                     const MultiIndex & id);


    bool addCoefficient(Constraint * constrPtr, Variable * varPtr, const Double & coef = 1);

    void prepareModel();
    ColGenSpConf * getColGenSubproblem(const std::string & name, const MultiIndex & id);

    void getEnumeratedSolutions(std::vector<std::tuple<double, double, BcSolution> > & enumSolutions);

};
#endif // BcModelClasses_h
