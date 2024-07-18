/**
 *
 * This file bcUserControlParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCUSERCONTROLPARAMETERS_HPP_
#define BCUSERCONTROLPARAMETERS_HPP_

#include "bcAbstractControlParameters.hpp"
#include "bcModelParameterC.hpp"

class UserControlParameters : public AbstractControlParameters
{
public:
    UserControlParameters();
    virtual ~UserControlParameters();

    /// Main parameters
    /// ---------------

    ApplicationParameter<long> GlobalTimeLimitInTick;
    ApplicationParameter<long> GlobalTimeLimit;
    /// default = 0: Wallclock time, 1: CPU time
    ApplicationParameter<int> ClockType;
    ApplicationParameter<int> MaxNbOfBBtreeNodeTreated;
    ApplicationParameter<double> optimalityGapTolerance;
    ApplicationParameter<double> relOptimalityGapTolerance;
    ApplicationParameter<bool> ApplyPreprocessing;
    ApplicationParameter<bool> PreprocessVariablesLocalBounds;
    ApplicationParameter<int> treeSearchStrategy;
    ApplicationParameter<int> OpenNodesLimit;
    ApplicationParameter<int> DEFAULTPRINTLEVEL;

    /// MIP solver parameters
    /// ---------------------

    /// for the moment, only CPLEX_SOLVER is available for solverName
    ApplicationParameter<std::string> solverName;
    ApplicationParameter<int> MipSolverMaxBBNodes;
    ApplicationParameter<int> MipSolverMaxTime;
    ApplicationParameter<int> MipSolverMultiThread;

    /// Column generation parameters
    /// ----------------------------

    ApplicationParameter<char> SolverSelectForMast;
    ApplicationAdvancedParameter<SolutionMethod, int> colGenSubProbSolMode;
    ApplicationAdvancedParameter<MasterInitMode, int> mastInitMode;
    /// artificial variables cost update factor for Phase 1, (default is 1.2)
    ApplicationParameter<float> ArtVarPenaltyUpdateFactor;
    /// maximum number of updates of the penalties of artificial variables cost updates before going to pure Phase 1
    ApplicationParameter<int> ArtVarMaxNbOfPenaltyUpdates;
    ApplicationParameter<int> MaxNbOfStagesInColGenProcedure;
    ApplicationParameter<bool> GenerateProperColumns;
    /// default = true # if true all generated columns are condidate for insertion in the master formulation;
    /// if false, only the best reduced cost column is candidate, the others go to the pool
    ApplicationParameter<bool> InsertAllGeneratedColumnsInFormRatherThanInPool;
    /// default = false # if true all candidate entering columns (as defined above)  are inserted in
    /// the master formulation wheiter they have negative reduced cost or not; if false, only the negative
    /// reduced cost columns among the candidates are inserted, the others go to the pool
    ApplicationParameter<bool> InsertNewNonNegColumnsDirectlyInFormRatherThanInPool;
    /// default = 10000, if the number of columns is at least this number, the columns clean-up will be called
    /// if ColumnCleanupThreshold is <= 0, then no column clean is performed
    ApplicationParameter<float> ColumnCleanupThreshold;
    /// default = 0.66, when cleaning up the master at least this ratio of columns will be left
    ApplicationParameter<float> ColumnCleanupRatio;
    /// default = 0 # we call reduced cost fixing if the current gap is below this ratio from the gap when
    /// the reduced cost fixing was called the last time (if 0 then never called, if >= 1, always called)
    ApplicationParameter<float> ReducedCostFixingThreshold;

    /// Cut generation parameters
    /// -------------------------

    ApplicationParameter<int> MasterCuttingPlanesDepthLimit;
    ApplicationParameter<int> MaxNbOfCutGeneratedAtEachIter;
    /// default = 0, if the gap improvement after the latest cut round is smaller than this value
    /// (in per cent from the latest gap) the cut tailing off counter will be increased
    ApplicationParameter<float> CutTailingOffThreshold;
    /// default = 0, if cut tailing off counter exceeds this value, cut generation will be stopped at this node
    ApplicationParameter<int> CutTailingOffCounterThreshold;
    /// default = 1000, if the number of master cuts is at least this number, the cuts clean-up will be called
    ApplicationParameter<float> CutCleanupThreshold;
    /// default = 0.66, when cleaning up the master this ratio of columns will be left
    ApplicationParameter<float> CutCleanupRatio;
    /// default = 2.0, priority of sp. relaxation improvement (in comparison to cut generators)
    ApplicationParameter<float> ColGenSpRelaxationImprovementPriority;

    /// Stabilization parameters
    /// ------------------------

    /// default = 1 ; 0 - no smoothing; 1 - automatic smoothing # alpha in [0,1):
    /// \pi^smooth = alpha {\pi^best, \pi^old, \pi^active} + (1- alpha) \pi^kelley
    ApplicationParameter<float> colGenDualPriceSmoothingAlphaFactor;
    /// default = 0 ; 0 - no smoothing; 1 - automatic smoothing # beta is the smoothing factor in computing
    /// the direction of the sep poApplicationParameter<int>
    ApplicationParameter<float> colGenDualPriceSmoothingBetaFactor;
    /// default = 0 - curvature based behaviour, 1 - half interval and slope based behaviour, 2 - multi-point mode
    ApplicationAdvancedParameter<StabilizationFunctionType, int> colGenStabilizationFunctionType;
    /// default = 0 - curvature based behaviour, 1 - half interval and slope based behaviour, 2 - multi-point mode
    ApplicationAdvancedParameter<ColGenProximalStabilizationMode, int> colGenProximalStabilizationRule;
    /// parameter for automatic setting of initial curvature of half intervals (default is 1), 0 if manual setting is
    /// used
    ApplicationParameter<float> StabilFuncKappa;
    /// maximum B&P tree depth when the smoothing is applied
    ApplicationParameter<long> colGenStabilizationMaxTreeDepth;
    /// default = 0, minimum phase in which stabilization is applied
    ApplicationParameter<int> StabilizationMinPhaseOfStage;

    /// Primal heuristics parameters
    /// ----------------------------

    ApplicationParameter<bool> UseInitialPrimalHeur;

    ApplicationParameter<int> MaxTimeForRestrictedMasterIpHeur;
    ApplicationParameter<int> CallFrequencyOfRestrictedMasterIpHeur;
    ApplicationParameter<int> MIPemphasisInRestrictedMasterIpHeur;
    ApplicationParameter<double> PolishingAfterTimeInRestrictedMasterIpHeur;

    ApplicationParameter<int> DivingHeurUseDepthLimit;
    ApplicationParameter<int> CallFrequencyOfDivingHeur;
    ApplicationAdvancedParameter<MultitokenSelectionStrategyVector, SelectionStrategy> RoundingColSelectionCriteria;
    ApplicationParameter<bool> FixIntValBeforeRoundingHeur;
    ApplicationParameter<int> MaxNbOfCgIteDuringRh;
    ApplicationParameter<int> MaxLDSbreadth;
    ApplicationParameter<int> MaxLDSdepth;
    ApplicationParameter<bool> DivingHeurStopsWithFirstFeasSol;
    /// default false, if true, we call the preprocessing before choosing a candidate to verify whether it causes
    /// infeasibility
    ApplicationParameter<bool> DivingHeurPreprocessBeforeChoosingVar;
    ApplicationParameter<int> StrongDivingCandidatesNumber;
    ApplicationAdvancedParameter<StrongBranchingPhaseParameter, std::string> EvalAlgParamsInDiving;

    ApplicationParameter<int> LocalSearchHeurUseDepthLimit;
    ApplicationAdvancedParameter<MultitokenSelectionStrategyVector, SelectionStrategy> LocalSearchColSelectionCriteria;
    ApplicationParameter<double> MaxFactorOfColFixedByLocalSearchHeur;
    /// to make control the number of improvement trials of the localSearch heuristic
    ApplicationParameter<int> MaxLocalSearchIterationCounter;

    /// Strong branching parameters
    /// ---------------------------

    ApplicationAdvancedParameter<StrongBranchingPhaseParameter, std::string> StrongBranchingPhaseOne;
    ApplicationAdvancedParameter<StrongBranchingPhaseParameter, std::string> StrongBranchingPhaseTwo;
    ApplicationAdvancedParameter<StrongBranchingPhaseParameter, std::string> StrongBranchingPhaseThree;
    ApplicationAdvancedParameter<StrongBranchingPhaseParameter, std::string> StrongBranchingPhaseFour;
    ApplicationParameter<bool> SimplifiedStrongBranchingParameterisation;
    ApplicationParameter<int> StrongBranchingPhaseOneCandidatesNumber;
    ApplicationParameter<float> StrongBranchingPhaseOneTreeSizeEstimRatio;
    ApplicationParameter<int> StrongBranchingPhaseTwoCandidatesNumber;
    ApplicationParameter<float> StrongBranchingPhaseTwoTreeSizeEstimRatio;

    /// Debug output
    /// ------------

    ApplicationParameter<int> printMasterPrimalSols;

    /// VRPSolver extension parameters
    /// ------------------------------

    /// time threshold in RCSP pricing after which non-robust cut generation is stopped
    ApplicationParameter<double> RCSPstopCutGenTimeThresholdInPricing;

    /// time threshold in RCSP pricing after which rollback is performed
    ApplicationParameter<double> RCSPhardTimeThresholdInPricing;

    /// time threshold in RCSP for the bucket arc elimination
    ApplicationParameter<double> RCSPredCostFixingTimeThreshold;

    /// # max. of buckets per vertex, should be >= 1 (used if step size provided by the user is non-positive)
    ApplicationParameter<int> RCSPnumberOfBucketsPerVertex;

    /// 0 - static bucket steps (given by the user)
    /// 1 - dynamic mode, "aggregated" adjustment of bucket steps (same for all vertices)
    /// 2 - dynamic mode, vertex specific adjustment of bucket steps
    ApplicationParameter<int> RCSPdynamicBucketSteps;

    /// 0 - monodirectional search in the RCSP pricing solver
    /// 1 - bidirectional search in the RCSP solver and in the RCSP enumeration
    /// 2 - bidirectional search only in exact RCSP solver and in the RCSP enumeration
    /// 3 - bidirectional search in RCSP solver and monodirectional RCSP enumeration
    /// 4 - bidirectional search only in exact RCSP solver and monodirectional RCSP enumeration
    ApplicationParameter<int> RCSPuseBidirectionalSearch;

    /// 0 - no bucket arc elimination in the RCSP solver
    /// 1 - normal (exaustive) reduced cost fixing
    /// 2 - light reduced cost fixing (check only completion bounds)
    /// 3 - normal reduced cost fixing of original arcs only
    /// 4 - light reduced cost fixing of original arcs only
    ApplicationParameter<int> RCSPapplyReducedCostFixing;

    /// maximum number of generated columns from each heuristic RCSP pricing subproblem
    ApplicationParameter<int> RCSPmaxNumOfColsPerIteration;

    /// maximum number of generated columns from each exact RCSP pricing subproblem
    ApplicationParameter<int> RCSPmaxNumOfColsPerExactIteration;

    /// maximum number of the same sense labels in the enumeration procedure for the RCSP pricing solver
    ApplicationParameter<int> RCSPmaxNumOfLabelsInEnumeration;

    /// maximum number of enumerated solutions in the RCSP pricing solver
    ApplicationParameter<int> RCSPmaxNumOfEnumeratedSolutions;

    /// maximum total number of enumerated solutions to close the node by MIP solver (in the middle of the node)
    ApplicationParameter<int> RCSPmaxNumOfEnumSolutionsForMIP;

    /// maximum total number of enumerated solutions to to close the node by MIP solver (in the end of the node)
    ApplicationParameter<int> RCSPmaxNumOfEnumSolsForEndOfNodeMIP;

    /// initial neighbourhood size for elementarity sets for the case the distance matrix is given
    ApplicationParameter<int> RCSPinitNGneighbourhoodSize;

    /// maximum neighbourhood size for elementarity sets in the case of dynamic NG
    ApplicationParameter<int> RCSPmaxNGneighbourhoodSize;

    /// maximum number of rows in the rank-1 cuts (if equals to 0, rank-1 cuts are not generated)
    ApplicationParameter<int> RCSPrankOneCutsMaxNumRows;

    /// maximum number of rank-1 cuts with the same number of rows generated per round
    ApplicationParameter<int> RCSPrankOneCutsMaxNumPerRound;

    /// memory type for rank-1 cuts :
    /// 0 - automatic limited memory, both arc memory and node memory are tried;
    /// 1 - arc memory;
    /// 2 - node memory;
    /// 3 - full memory
    ApplicationParameter<int> RCSPrankOneCutsMemoryType;

    /// 0 - use exhaustive separation
    /// < 0 - use exhaustive separation of 4 and 5-rows cuts and local search with this number of iterations for >5 rows
    /// > 0 - use local search with this number of iterations
    ApplicationParameter<int> RCSPrankOneCutsLSnumIterations;

    /// whether the columns covering the same set of vertice may be generated at the same time by the RCSP solver
    ApplicationParameter<bool> RCSPallowRoutesWithSameVerticesSet;

    /// if absolute value <= 1.0, it has not effect; if absolute value > 1.0, the reduced cost fixing and enumeration
    /// are performed with the primal-dual gap divided by this value, thus making the whole approach heuristic,
    /// if negative, false gap is applied only when rank-1 cuts are present
    ApplicationParameter<double> RCSPredCostFixingFalseGap;

    ApplicationParameter<long long int> SafeDualBoundScaleFactor;

    /// maximum number of the same sense labels in the enumeration procedure for the RCSP pricing solver
    /// launched inside the diving heuristic
    ApplicationParameter<int> RCSPmaxNumOfLabelsInHeurEnumeration;

    ApplicationParameter<int> MaxNumEnumSolsInRestrictedMasterIpHeur;

    virtual void addParameters(ParameterManager& parameterManager);
    virtual std::ostream& printParameters(std::ostream& os) const;
    virtual std::ostream& printVRPSolverParameters(std::ostream& os) const;
    virtual void postTreatment();
};

#ifndef SOLVER_NAMES
#define SOLVER_NAMES
/**
 * Solver's name constant.
 */
static const std::string CPLEX_SOLVER = "CPLEX_SOLVER";
static const std::string GLPK_SOLVER = "GLPK_SOLVER";
static const std::string XPRESSMP_SOLVER = "XPRESSMP_SOLVER";
static const std::string XPRESS_SOLVER = "XPRESS_SOLVER";
static const std::string CLP_SOLVER = "CLP_SOLVER";
static const std::string GUROBI_SOLVER = "GUROBI_SOLVER";
#endif

#endif /* BCUSERCONTROLPARAMETERS_HPP_ */
