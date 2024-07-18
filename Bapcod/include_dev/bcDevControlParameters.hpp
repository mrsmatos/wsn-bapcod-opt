/**
 *
 * This file bcDevControlParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCDEVCONTROLPARAMETERS_HPP_
#define BCDEVCONTROLPARAMETERS_HPP_
#include "bcAbstractControlParameters.hpp"
#include "bcModelParameterC.hpp"

class DevControlParameters : public AbstractControlParameters
{
public:
    DevControlParameters();
    virtual ~DevControlParameters();

    ApplicationParameter<int> DEFAULTTESTLEVEL;

    ApplicationParameter<bool> VRPSEimposeCapacityResourceByCuts;
    /// -1 : automatic detection
    /// 0 : time resource is critical
    /// 1 : capacity resource is critical
    ApplicationParameter<int> VRPSEcriticalResource;

    ApplicationAdvancedParameter<SolutionMethod, int> ovfSolMode;
    ApplicationAdvancedParameter<SolutionMethod, int> masterSolMode;
    ApplicationParameter<int> MaxNbOfCgIterations;
    ApplicationParameter<bool> Search4NegRedCostColInInactivePool;

    ApplicationParameter<bool> MipSolverRecordNamesInFormulation;
    ApplicationAdvancedParameter<SolutionStatus, int> RequiredSolStatForOvf;
    ApplicationParameter<bool> PreprocessorOnForOvf;
    ApplicationParameter<bool> ProbingOnForOvf;
    ApplicationParameter<bool> AutomaticCuttingPlanesOnForOvf;
    ApplicationParameter<char> SolverSelectForOvf;
    ApplicationAdvancedParameter<SolutionStatus, int> RequiredSolStatForMast;
    ApplicationParameter<bool> PreprocessorOnForMast;
    ApplicationParameter<bool> ProbingOnForMast;
    ApplicationParameter<bool> AutomaticCuttingPlanesOnForMast;
    ApplicationParameter<char> SolverSelectForMastAfterBarrierConvergence;
    ApplicationAdvancedParameter<SolutionStatus, int> RequiredSolStatForColGenSp;
    ApplicationParameter<bool> PreprocessorOnForColGenSp;
    ApplicationParameter<bool> ProbingOnForColGenSp;
    ApplicationParameter<bool> AutomaticCuttingPlanesOnForColGenSp;
    ApplicationParameter<char> SolverSelectForColGenSp;
    ApplicationParameter<double> MasterMipSolverRightHandSideZeroTol;
    ApplicationParameter<double> MasterMipSolverReducedCostTolerance;
    ApplicationParameter<double> MasterMipSolverBarrierConvergenceTolerance;
    ApplicationParameter<double> ColGenSpMipSolverRightHAndSideZeroTol;
    ApplicationParameter<double> ColGenSpMipSolverReducedCostTolerance;

    ApplicationParameter<int> MaxDepthInBBtree;
    ApplicationParameter<bool> StopWhenBranchingFails;
    ApplicationParameter<double> BapCodReducedCostTolerance;
    ApplicationParameter<double> BapCodCutViolationTolerance;
    ApplicationParameter<double> BapCodIntegralityTolerance;

    ApplicationParameter<bool> TestAggregateMasterSol4Integrality;
    ApplicationParameter<bool> VerifyColsIntegralityInTestSolForIntegrality;
    ///  default = true, Hierarchic Branching, first consider Pure Master Var, then generated columns and the associated subproblem variables
    ApplicationParameter<bool> BranchFirstOnPureMasterVariables;
    /// default = false, Hierarchic Branching, first generated columns and the associated subproblem variables, then consider Pure Master Var
    ApplicationParameter<bool> BranchFirstOnSubProblemVariables;
    ApplicationParameter<bool> ApplyStrongBranchingEvaluation; /// for backwards compatibility
    ApplicationParameter<float> StrongBranchingParameter4NotPromissingConservativeness;
    ApplicationParameter<bool> StrongBranchingExactPhaseInterrupt;
    ApplicationParameter<bool> StrongBranchingUseHistorySubtreeSize;
    ApplicationParameter<bool> StrongBranchingUseHistory;
    ApplicationParameter<bool> TerminateCgWhenRoundedDbCannotImprove;
    ApplicationParameter<bool> UseObjScalingFact;

    ApplicationParameter<bool> LocArtVarInConvexityConstr;

    /// curvature or angles are divided by this factor when stabilization variables are in the solution (default is 10)
    ApplicationParameter<float> StabilFuncArtVarInSolUpdateFactor;
    /// initial curvature of the penalty function (default is 1)
    ApplicationParameter<float> StabilFuncCurvature;
    /// advance rate along the curvature (default is 1)
    ApplicationParameter<float> StabilFuncCurvatureAdvanceRate;
    /// initial length of the outer flat half-interval of the stabilization function (>0, default=100)
    ApplicationParameter<float> StabilFuncOuterHalfInterval;
    /// initial length of the inner flat half-interval of the stabilization function
    /// (>=0, <=StabilFuncOuterFlatInterval, default=StabilFuncOuterHalfInterval/2)
    ApplicationParameter<float> StabilFuncInnerHalfInterval;
    /// multiplication factor of inner half interval when passing to a child node
    ApplicationParameter<float> StabilFuncHalfIntervalChildNodeFactor;
    /// initial angle of the outer piece (additionally to inner one) of the stabilization function (>=0, default is 0.9)
    ApplicationParameter<float> StabilFuncOuterAngle;
    /// initial angle multiplied by rhs of the inner piece of the stabilization function (>=0, default is 0.1)
    ApplicationParameter<float> StabilFuncInnerAngle;
    /// default = 0 (off) # stop smoothing when optimality gap is less than colGenDualPriceSmoothingMinGap
    ApplicationParameter<float> colGenDualPriceSmoothingMinGap;
    /// default = 0 (off) # start smoothing when optimality gap is less than colGenDualPriceSmoothingMaxGap
    ApplicationParameter<float> colGenDualPriceSmoothingMaxGap;
    ApplicationParameter<int> colGenDualPriceSmoothingMaxNbOfUpdate;
    /// (restricting \pi to be in convex hull of trust \pi's).

    ApplicationParameter<bool> SplitColIntoDissagregateSpVar;
    ApplicationParameter<bool> SplitColIntoDissagregateSpVarInHeadIn;
    ApplicationParameter<bool> PriceAllSubproblemsForBestDualBound;
    ApplicationParameter<bool> CyclicSpScanning;

    /// default = true, if true, each time a new column is generated, we verify whether it exists in a pool,
    /// if yes, we do not regenerate it
    ApplicationParameter<bool> UseColumnsPool;
    /// if difference between upper bound and lower bound is less than this value, cut generation is stopped
    ApplicationParameter<double> CutTailingOffAbsTolerance;
    ApplicationParameter<int> MinNumOfCutRoundsBeforeStopBySp; /// default = 1,
    /// default = true, the minimum level of priority can be increased
    ApplicationParameter<bool> CutDynamicPriorityLevel;
    /// default = false, if true, column generation is not terminated early because of the small primal-dual gap
    ApplicationParameter<bool> runColGenUntilFullConvergence;

    /// Primal Heuristic Control Parameters
    ApplicationParameter<int> UseGreedyHeur;
    ApplicationParameter<bool> UseCustomFracSolBasedHeur;
    ApplicationParameter<bool> runColGenAfterFixingPureMastVarInDiving; /// default = true
    /// The false gap is multiplied by this parameter when trying to enumerate columns for the restricted master heur.
    ApplicationParameter<double> EnumHeuristicFalseGapMultiplier;
    /// Number of false gap multiplications by this parameter when trying enumeration for the restr. master heur.
    ApplicationParameter<int> EnumHeuristicNumberOfTries;
    /// default = false. When the master includes both pure original variables and generated colmns,
    /// if true this forces to dive by fixing only master columns
    ApplicationParameter<bool> UseDivingHeurOnMasterColOnly;
    /// default = false. When the master includes both pure original variables and generated colmns,
    /// if true this forces to dive by fixing only Pure master variables
    ApplicationParameter<bool> UseDivingHeurOnPureMastVarOnly;
    ApplicationParameter<int> MaxNbOfPenaltUpdatesDuringRH;
    ApplicationParameter<bool> IgnoreIntValWhenSelectingRoundedVarInRH;
    ApplicationParameter<bool> ActivateAllColumnsForRestrictedMasterIpHeur;
    ApplicationParameter<int> MaxNumEnumSolsInUserHeuristicFunctor;

    ApplicationParameter<int> PricingStrategy;
    ApplicationParameter<int> MaxNbPromisingSpFound;
    ApplicationParameter<int> MaxNbUnpromisingSpFound;

    ApplicationParameter<std::string> NameOfOutputSummaryFile;

    /// experimental parameter : whether a new generation meta-solver should be used instead of the "traditional" one
    /// 0 - do not use meta-solver
    /// 1 - use meta-solver alone
    /// 2 - use meta-solver as a verification solver for the "traditional" one
    ApplicationParameter<int> RCSPuseMetaSolver;

    ApplicationParameter<int> RCSPprintLevel;
    /// whether to check dominance between labels in
    /// different buckets in the RCSP pricing solver
    ApplicationParameter<bool> RCSPcheckDominInOtherBuckets;
    /// 0 - do not use completion bounds in the RCSP solver
    /// 1 - use standard competion bounds in the RCSP solver
    /// 2 - use exact (exaustive) completion bounds
    ApplicationParameter<int> RCSPuseCompletionBoundsInPricing;
    /// number of dominance checks threshold in RCSP pricing after which non-robust cut generation is stopped
    ApplicationParameter<double> RCSPdomChecksThresholdInPricing;
    /// 0 (default) : 1 (phase 2) - inf with simple dom. (phase 1) - exact (phase 0)
    /// 1 : 1 (phase 2) - 4 (phase 1) - exact (phase 0)
    /// 2 : 1 (phase 2) - 4 (phase 1) - 8 (phase 0)
    /// 3 : 1 (phase 2) - 4 (phase 1) - 16 (phase 0)
    ApplicationParameter<int> RCSPheurLabelingStrategy;
    /// experimental parameter : the split strategy to apply to labels
    /// 0 - no split (default)
    /// 1 - split by linear pieces
    /// 2 - split by buckets
    ApplicationParameter<int> RCSPlabelSplitStrategy;

    /// advanced parameter (can be modified only by an expert)
    ApplicationParameter<double> RCSPdynBuckStepAdjustRatioThreshold;
    /// advanced parameter (can be modified only by an expert)
    ApplicationParameter<double> RCSPdynBuckStepAdjustNumBuckArcsThreshold;
    /// 0 - do not use completion bounds
    /// 1 - use standard competion bounds
    /// 2 - use exact (exaustive) completion bounds
    ApplicationParameter<int> RCSPuseComplBoundsInRedCostFixing;

    /// DSSR serves to decrease NG-memory
    ///  first bit 001 - call after the first convergence in a node
    /// second bit 010 - call after the NG convergence in the root
    ///  third bit 100 - call on the rollback
    ApplicationParameter<int> RCSPuseDSSRInMode;
    /// maximum number of columns with forbidden cycles in one round of DSSR
    ApplicationParameter<int> RCSPmaxNumOfColumnsInDSSR;
    /// all cycles with up to this size are forbidden in DSSR
    ApplicationParameter<int> RCSPmaxCycleSizeInDSSR;
    /// 0 - no dynamic NG
    /// +/-1 - dynamic NG with increase on vertices
    /// +/-2 - dynamic NG with increase on arcs
    /// negative - increase NG only when there are no rank-1 cuts
    ApplicationParameter<int> RCSPdynamicNGmode;
    /// maximum number of columns with forbidden cycles in one round of NG neghbourhood increase
    ApplicationParameter<int> RCSPmaxNumOfColumnsInDynamicNG;
    /// all cycles with up to this size are forbidden in NG neighbourhood increas
    ApplicationParameter<int> RCSPmaxCycleSizeInDynamicNG;
    /// maximum average neighbourhood size for elementarity sets in the case of dynamic NG
    ApplicationParameter<int> RCSPmaxNGaverNeighbourhoodSize;

    /// generate capacity cuts (if defined by the model)
    ApplicationParameter<bool> RCSPuseCapacityCuts;
    /// 0 - CVRPSEP separator (need to put CVRPSEP in Tools/separation_libs/CVRPSEP and rerun cmake)
    /// 1 - BCP_RCSP library separator
    ApplicationParameter<int> RCSPcapacityCutsSeparator;
    /// maximum number of capacity cuts generated per round
    ApplicationParameter<int> RCSPcapCutsMaxNumPerRound;
    /// generate resource consumption knapsack cuts (if defined by the model)
    /// -1 - do not generate RLCK cuts
    /// 0 - generate RLCK cuts only by rounding
    /// 6\8\10 - generate RLCK cuts by rounding and 1/6\1/8\1/10-inequalities by enumeration
    /// -6\-8\-10 - generate only 1/6\1/8\1/10-inequalities by enumeration
    ApplicationParameter<int> RCSPresConsKnapsackCutsMode;
    /// maximum number of capacity cuts generated per round
    ApplicationParameter<int> RCSPresConsKnapsackCutsMaxNumPerRound;
    /// maximum number of clique cuts generated per round
    /// (need to put CliqueSep in Tools/separation_libs/CliqueSep and rerun cmake)
    ApplicationParameter<int> RCSPcliqueCutsMaxNumPerRound;
    /// maximum separation time for clique cuts
    ApplicationParameter<double> RCSPcliqueCutsMaxSepTime;
    /// 0 - smallest arc memory + loop arcs
    /// 1 - smallest arc memory + loop arcs + arcs connecting row sets
    /// 2 - smallest arc memory + loop arcs + all arcs in the row set
    ApplicationParameter<int> RCSPrankOneCutsArcMemoryType;
    /// whether subproblem dependent memory is used for rank-1 cuts
    ApplicationParameter<bool> RCSPrankOneCutsSpDependentMemory;
    /// packing sets neighbourhood size for 4 and 5-rows cut separation by local search
    /// (for the case when the distance matrix for packing sets is defined)
    ApplicationParameter<int> RCSPrankOneCutsNeighbourhoodSize;
    ///maximum number of triplets in the tuple based separation
    ApplicationParameter<int> RCSPrankOneCutsMaxNumTriplets;
    /// which type of cuts to separate
    /// 0 - all cuts (packing and covering)
    /// 1 - only packing cuts
    /// 2 - only covering cuts
    ApplicationParameter<int>  RCSPrankOneCutsTypeToSeparate;
    ApplicationParameter<bool> RCSPuseCovSetsForSepOneRowPackingCuts;
    /// whether to use exact (exhaustive) completion bounds in enumeration
    ApplicationParameter<bool> RCSPuseExactComplBoundsInEnumeration;
    /// maximum possible gap to run enumeration if gap is larger the enumeration will not be called
    ApplicationParameter<double> RCSPmaxGapToRunEnumeration;
    /// for advanced use (used in formulation for split-delivery)
    ApplicationParameter<bool> RCSPuseVertexIdsInEnumeration;
    /// advanced parameter (can be modified only by an expert)
    ApplicationParameter<double> RCSPstopRatioForConcatenationInEnum;

    /// Code OUTPUT CONTROL
    ApplicationParameter<int> ColGenLogFrequency;
    ApplicationAdvancedParameter<SetParameter<std::string>, std::string> timeKeySet;
    ApplicationAdvancedParameter<VectorParameter<std::string>, std::string> outputKeyVector;

    // batch mode options
    ApplicationParameter<std::string> configFileName;
    ApplicationParameter<std::string> instance_file;
    ApplicationParameter<std::string> statistics_file;
    ApplicationParameter<std::string> baPTreeDot_file;
    ApplicationParameter<std::string> colgeninfo_file;
    ApplicationParameter<std::string> debugsolution_file;
    ApplicationParameter<std::string> fractionalSolutionDotFile;

    ApplicationParameter<bool> bestDualBoundIsConstantForTest;
    ApplicationParameter<bool> bestIncumbentIsConstantForTest;

    // advanced parameters
    ApplicationParameter<bool> useSPVarMembershipCache;

    ApplicationParameter<bool> CheckSpOracleFeasibility;
    ApplicationParameter<bool> CheckOracleOptimality;

    /***************************************************
     * END of parameters definitions
     ****************************************************/
    virtual void addParameters(ParameterManager& parameterManager);
    virtual std::ostream& printParameters(std::ostream& os) const;

    virtual void postTreatment();
};

#ifndef SOLVER_NAMES
#define SOLVER_NAMES
/**
 * Solver's name constant.
 */
static const std::string CPLEX_SOLVER = "CPLEX_SOLVER";
static const std::string GLPK_SOLVER = "GLPK_SOLVER";
static const std::string KNAPSACK_SOLVER = "KNAPSACK_SOLVER";
static const std::string XPRESSMP_SOLVER = "XPRESSMP_SOLVER";
static const std::string XPRESS_SOLVER = "XPRESS_SOLVER";
static const std::string CLP_SOLVER = "CLP_SOLVER";
static const std::string GUROBI_SOLVER = "GUROBI_SOLVER";
#endif


#endif /* BCDEVCONTROLPARAMETERS_HPP_ */
