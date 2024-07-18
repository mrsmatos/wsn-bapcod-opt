/**
 *
 * This file bcModelRCSPSolver.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELRCSPSOLVER_H_
#define BCMODELRCSPSOLVER_H_

#define DSSR_MODE_NONE                          0 /// 0000
#define DSSR_MODE_FIRST_COVERGENCE_OF_A_NODE    1 /// 0001
#define DSSR_MODE_WHEN_NG_CONVERGES             2 /// 0010
#define DSSR_MODE_ON_ROLLBACK                   4 /// 0100
#define DSSR_MODE_AFTER_EACH_RELAXATION_IMPROVE 8 /// 1000

#include "bcModelFormulationC.hpp"
#include "bcModelNetworkFlow.hpp"
#include "bcModelCutConstrC.hpp"
#include "bcControlParameters.hpp"
#include <lemon/list_graph.h>

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif /* BCP_RCSP_IS_FOUND */

/*****************************************************************************/
/*****************************************************************************/
/*****                        TODO: Cuts                                 *****/
/*****************************************************************************/
/*****************************************************************************/


class BcLimMemRankOneCutConstrArray: public BcCutConstrArray
{
	GenericCutConstr* _genCutCostrPtr;

public:
    BcLimMemRankOneCutConstrArray(const BcFormulation & formulation,
                                  const double & rootPriorityLevel = 1.0,
                                  const double & nonRootPriorityLevel = 1.0,
                                  const bool & isFacultative = true,
                                  const std::string & spVarName = "");

    void runStandAloneSeparation(const std::string& filePath);

    virtual ~BcLimMemRankOneCutConstrArray();
};

class BcCapacityCutConstrArray: public BcCutConstrArray
{
public:
    BcCapacityCutConstrArray(const BcFormulation & formulation,
                             const int & maxCapacity,
                             const std::vector<int> & demands,
                             const bool & isFacultative = true,
                             const bool & equalityCase = true,
                             const int & twoPathCutsResId = -1,
                             const double & rootPriorityLevel = 1.0,
                             const double & nonRootPriorityLevel = 1.0);
    virtual ~BcCapacityCutConstrArray();
};


#ifdef BCP_RCSP_IS_FOUND

void getActiveRankOneCuts(const BcFormulation & spPtr,
                          std::vector<std::pair<const bcp_rcsp::RankOneCut *, double> > & rankOneCuts);

class BcExtendedArcCut : public bcp_rcsp::DiscreteCutInterface
{
  const NetworkFlow * _networkPtr;
  ExtendedArcCut * _cutPtr;
public:
  BcExtendedArcCut(const NetworkFlow * networkPtr, ExtendedArcCut * cutPtr);
  ~BcExtendedArcCut() {};
  double curDualVal() const;
  bool inCurProb() const;
  virtual bcp_rcsp::DiscreteCutInterface * createCopy() const;
  virtual int id() const;
  virtual void nicePrint() const;
  virtual bool customCut() const;
  virtual double getArcCoefficient(const int & tailVertId, const int & headVertId, const double * tailResCons) const;
  virtual double getArcCoefficient(const int & arcId, const double * resCons, const bool & isTailResCons) const;
  virtual double getVertRouteCoefficient(const std::vector<int> & routeVertIds,
                                         const std::vector<std::vector<double> > & routeResCons) const;
  virtual double getArcRouteCoefficient(const std::vector<int> & routeArcIds,
                                        const std::vector<std::vector<double> > & routeResCons) const;
  virtual void reserveCut() const; /// this prevents the cut from being deleted by BaPCod
  virtual void releaseCut() const ; /// this releases the cut (it can now be deleted by BaPCod)
  inline bool operator<(const BcExtendedArcCut & otherCut) const
  {
    return (id() < otherCut.id());
  }
  inline bool operator==(const BcExtendedArcCut & otherCut) const
  {
    return (id() == otherCut.id());
  }
};
#endif /* BCP_RCSP_IS_FOUND */

struct BcCustomExtendedArcCutInfo
{
  BcCustomExtendedArcCutInfo()
  {
  }
  virtual ~BcCustomExtendedArcCutInfo()
  {
  }
  virtual double getArcCoefficient(const BcArcInfo * arcInfoPtr, const double * resCons,
                                   const bool & isTailResCons) const = 0;
  virtual double getRouteCoefficient(const std::vector<const BcArcInfo *> & arcInfoPts,
                                     const std::vector<std::vector<double> > & resConsumption) const = 0;

  virtual void nicePrint(std::ostream& os) const = 0;
};

class GenericExtendedArcCutConstr;

class BcCustomExtendedArcCutSeparationFunctor
{
  friend class BcCustomExtendedArcCutArray;
  GenericExtendedArcCutConstr * _genExtArcCutCostrPtr;
protected:
  BcConstr createNewCut(BcCustomExtendedArcCutInfo * cutInfoPtr, const char & sense, const double & rhs);
public:
  BcCustomExtendedArcCutSeparationFunctor();
  virtual ~BcCustomExtendedArcCutSeparationFunctor();
  virtual int cutSeparationRoutine(BcFormulation formPtr,
				     	           BcSolution & projectedSol,
                                   std::list<std::pair<double, BcSolution> > & columnsInSol,
					               const double & violationTolerance,
					               std::list<BcConstr> & cutList);
};

class BcCustomExtendedArcCutArray
{
protected:
  GenericExtendedArcCutConstr * _genExtArcCutCostrPtr;
public:
  BcCustomExtendedArcCutArray(const BcFormulation & formulation,
                              const std::string & name,
                              const char & type = 'F',
                              const SelectionStrategy & priorityRule = SelectionStrategy::MostViolated,
                              const double & rootPriorityLevel = 1.0,
                              const double & nonRootPriorityLevel = 1.0);
  virtual ~BcCustomExtendedArcCutArray();
  const BcCustomExtendedArcCutArray & attach(BcCustomExtendedArcCutSeparationFunctor * separationFunctorPtr);
  inline const BcCustomExtendedArcCutArray & operator<<(BcCustomExtendedArcCutSeparationFunctor * separationFunctorPtr)
  {return attach(separationFunctorPtr);}
};


void getActiveDiscreteCuts(const BcFormulation & spPtr, std::vector<const BcExtendedArcCut *> & cutPts);

class BcResConsumptionKnapsackCutConstrArray: public BcCutConstrArray
{
    GenericCutConstr * _genCutCostrPtr;

public:
    BcResConsumptionKnapsackCutConstrArray(const BcFormulation & formulation,
                                           const double & rootPriorityLevel = 1.0,
                                           const double & nonRootPriorityLevel = 1.0);

    virtual ~BcResConsumptionKnapsackCutConstrArray();
};

#ifdef BCP_RCSP_IS_FOUND

void getActiveRouteLoadKnapsackCuts(const BcFormulation & spPtr,
                                    std::vector<std::pair<const bcp_rcsp::RouteLoadKnapsackCut *, double> > & rlkCuts);

#endif

/*****************************************************************************/
/*****************************************************************************/
/*****                   TODO: Branching constraints                     *****/
/*****************************************************************************/
/*****************************************************************************/

/*
** Paths per network branching
*/

class GenPathsPerNetworkBranchingConstr;

class BcPathsPerNetworkBranching
{
  GenPathsPerNetworkBranchingConstr * _genPathsPerNetworkBranchingConstrPtr;

public:
  BcPathsPerNetworkBranching(const BcFormulation & formulation,
		                     const double & priorityLevel = 1.0,
		                     const bool & toBeUsedInPreprocessing = true);

  virtual ~BcPathsPerNetworkBranching();
};

/*
** Elementarity sets assignment branching
*/

class GenPackSetAssignBranchingConstr;

class BcPackSetAssignBranching
{
  GenPackSetAssignBranchingConstr * _genPackSetAssignBranchingConstrPtr;

public:
  BcPackSetAssignBranching(const BcFormulation & formulation,
		                   const double & priorityLevel = 1.0,
		                   const bool & toBeUsedInPreprocessing = true);

  virtual ~BcPackSetAssignBranching();
};

/*
** Elementarity sets resource consumption branching
*/

class PackSetResConsGenBranchConstr;

class BcPackSetResConsumptionBranching
{
#ifdef BCP_RCSP_IS_FOUND
    PackSetResConsGenBranchConstr * _PackSetResConsGenBranchConstrPtr;
#endif
public:
    BcPackSetResConsumptionBranching(const BcFormulation & formulation,
                             const double & priorityLevel = 1.0);

    virtual ~BcPackSetResConsumptionBranching();
};


/*
** Packing sets Ryan&Foster branching
*/

class PackSetRyanFosterGenBranchConstr;

class BcPackSetRyanFosterBranching
{
#ifdef BCP_RCSP_IS_FOUND
    PackSetRyanFosterGenBranchConstr * _packSetRyanFosterGenBranchConstrPtr;
#endif
public:
    BcPackSetRyanFosterBranching(const BcFormulation & formulation,
                                 const double & priorityLevel = 1.0,
                                 const bool & usePackingSets = true);

    virtual ~BcPackSetRyanFosterBranching();
};

#ifdef BCP_RCSP_IS_FOUND

void getPackSetResConsActiveBranchConstrList(const BcFormulation & spPtr,
                                             std::vector<const bcp_rcsp::AccumResConsBranchConstraint *> & constrPts);

void getPackSetRyanFosterActiveBranchConstrList(const BcFormulation & spPtr,
                                                std::vector<const bcp_rcsp::RyanFosterBranchConstraint *> & constrPts);

/*****************************************************************************/
/*****************************************************************************/
/*****                          TODO: RCSP solver                        *****/
/*****************************************************************************/
/*****************************************************************************/

class RCSPOracleInfo : public BcSolverOracleInfo
{
    const bcp_rcsp::SolverRecord * _solverRecordPtr;

public:
    RCSPOracleInfo(const bcp_rcsp::SolverRecord * solverStatePtr) :
            _solverRecordPtr(solverStatePtr)
    {}

    const bcp_rcsp::SolverRecord * getSolverState() const
    {
        return _solverRecordPtr;
    }

    virtual ~RCSPOracleInfo()
    {
        delete _solverRecordPtr;
    }
};

class RCSPLabelExtensionCostFunctor : public bcp_rcsp::LabelExtensionCostFunctor
{
    NetworkFlow * _networkPtr;
    InstanciatedVar * _costInstVarPtr;

    virtual double getResConsDependentCost(const int arcId, const std::vector<double> & resConsumption,
                                           const bool isTailResCons) const override;
public:
    RCSPLabelExtensionCostFunctor(const BcFormulation & bcForm, const BcVar & costVar);

    InstanciatedVar * costInstVarPtr() {return _costInstVarPtr;}
    virtual double operator()(const BcArcInfo * arcInfoPtr,
                              const double * resConsumption, const bool isTailResCons) const = 0;
};

class RCSPFeasibilityCheckFunctor : public bcp_rcsp::FeasibilityCheckFunctor
{
    NetworkFlow * _networkPtr;

    virtual bool isFeasible(const bcp_rcsp::Solution & rcspSol) const override;
public:
    RCSPFeasibilityCheckFunctor(const BcFormulation & bcForm);

    virtual bool operator()(const std::vector<int> & vertIds, const std::vector<double> & waitingTimes,
            const std::vector<std::vector<double> > & resourceConsumptions) const = 0;

};


class BcRCSPFunctor : public BcSolverOracleFunctor
{
    ProbConfig * _probConfPtr;
    int _graphId;
    bool _useMetaSolver;
    bcp_rcsp::SolverParameters _solverParams;
    bcp_rcsp::SolverInterface * _solverInterface;
    std::vector<bcp_rcsp::ColGenPhaseConfig> _solverColGenPhaseConfig;
    BcRCSPFunctor * _verificationFunctorPtr;
    RCSPFeasibilityCheckFunctor * _feasibilityCheckFunctorPtr;
    RCSPLabelExtensionCostFunctor * _labelExtensionCostFunctorPtr;

    int _messageIdToCutGeneration;
    bool _paramFormulationIsBuilt;
    InstanciatedVar * _pureCostInstVarPtr;

    std::vector<const BcExtendedArcCut *> _discreteCutPts;
#ifdef CLIQUE_SEP_IS_FOUND
    std::vector<const BcCliqueCut *> _cliqueCutPts;
#endif

    void clearCutPts();
    void fillRCSPSolverParameters(const ControlParameters & bapcodParams);
    bool fillRCSPInput(BcFormulation spPtr, int colGenPhase, const std::vector<InstanciatedVar *> & instVarPts,
                       bcp_rcsp::SolverInput & input);
    void addToSolution(const NetworkFlow * networkPtr, const bcp_rcsp::Solution * rcspSolPtr,
                       const std::unordered_map<int, InstanciatedVar *> & resIdToVarPtrMap,
                       BcSolution & solution) const;

public:

    BcRCSPFunctor(const BcFormulation & spForm, const ControlParameters & params);

    explicit BcRCSPFunctor(const BcFormulation & spForm);

    ~BcRCSPFunctor() override;

    bool runAsStandaloneRCSPsolver(const std::string & fileName, int colGenPhase = 0);

    bool prepareSolver() override;

    bool setupNode(BcFormulation spPtr, const BcSolverOracleInfo * infoPtr) override;
    BcSolverOracleInfo * recordSolverOracleInfo(const BcFormulation spPtr) override;

    bool operator() (IN_ BcFormulation spPtr,
                     IN_ int phaseOfStageApproach,
                     OUT_ double & objVal,
                     OUT_ double & dualBound,
                     OUT_ BcSolution & primalSol) override;

    void reducedCostFixingAndEnumeration(IN_ BcFormulation spPtr,
                                         IN_ const int & enumerationMode,
                                         IN_ const double & threshold) override;

    int getMessageIdToCutGeneration() const override;

    bool getEnumeratedStatus() const override;

    int getNumberOfEnumeratedSolutions() const override;

    void checkEnumeratedSolutions(IN_ BcFormulation spPtr,
                                  IN_ const std::vector<Solution *> & solPts,
                                  OUT_ std::vector<bool> & solIsEnumerated) override;

    void getEnumeratedSolutions(IN_ BcFormulation spPtr,
                                IN_ const int & maxNumberOfSolutions,
                                OUT_ BcSolution & enumeratedSol,
                                OUT_ std::vector<double> & reducedCosts) override;

    bool getDebugSolution(BcFormulation spPtr, BcSolution & primalSol) override;

    bool setDebugSolution(const std::vector<std::vector<int> > & ids, bool vertexBased) override;

    bool isProperSolution(IN_ BcFormulation spPtr, IN_ BcSolution & bcSol) override;

    bool solSatisfiesCurrentSpRelaxation(IN_ BcFormulation spPtr, IN_ const BcSolution & solution) override;

    bool improveCurrentSpRelaxation(IN_ BcFormulation spPtr,
                                    IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                    IN_ const bool & masterConverged) override;

    bool drawPrimalSolutionToDotFile(IN_ BcFormulation spPtr,
                                     IN_ const std::vector<std::pair<BcSolution, double> > & colsInMasterSolution,
                                     IN_ const std::string & filename) override;

    bool lightenCurrentSpRelaxation(IN_ BcFormulation spPtr, const int & masterConverged,
                                    const int & callMode) override;

    void columnGenerationTerminated(IN_ BcFormulation spPtr, bool afterRedCostFixing,  int nodeOrder, int nodeDepth,
                                    int cutSepRound, double dualBound, double elapsedTime,
                                    bool masterConverged) override;

    void runWithDuals(BcFormulation spPtr, std::vector<double> & duals);


    void setDiscreteCase(const bool & value);
    void setVerificationFunctor(BcRCSPFunctor * functorPtr);
    void setLabelExtensionCostFunctor(RCSPLabelExtensionCostFunctor * functorPtr);
    void setPureCostBcVar(BcVar bcVar);
    void setFeasibilityCheckFunctor(RCSPFeasibilityCheckFunctor * functorPtr);
    void setParamPhasesConfig(const std::vector<bcp_rcsp::ColGenPhaseConfig> & colGenPhaseConfigVector);
    void setParamUseMoreTimers(const bool & value);
    void imposeSameResConsumptionInBucketCase();
    void setParamFormulationIsBuilt(const bool & value);
    void setParamPrintLevel(const int & value);
    void saveToStandaloneRCSPfile(const std::string & fileName);
    void enumeratedCoreFileName(const std::string & fileName, BcVar outOfEnumCoreVar);

    void addToSolution(bcp_rcsp::Solution * rcspSolPtr, BcSolution & bcSol) const;
};

bool runRCSPoracleAsStandaloneSolver(const BcModel & bcModel, const std::string & fileName);
#endif /* BCP_RCSP_IS_FOUND */

#endif //BCMODELGLOBALCUSTOMSOLVER_H_

