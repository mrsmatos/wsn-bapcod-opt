/**
 *
 * This file bcUserControlParameter.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUserControlParameters.hpp"

#include "bcParameterManager.hpp"
#include <iostream>
#include "bcModelParameterC.hpp"

using namespace std;

UserControlParameters::UserControlParameters() :
    GlobalTimeLimitInTick("GlobalTimeLimitInTick", 2147483645, "Undocumented"),
    GlobalTimeLimit("GlobalTimeLimit", 0, "Undocumented"),
    ClockType("ClockType", 0, "Undocumented"),
    solverName("solverName", CPLEX_SOLVER, "Undocumented"),
    MipSolverMaxBBNodes("MipSolverMaxBBNodes", 2000000, "Undocumented"),
    MipSolverMaxTime("MipSolverMaxTime", 360000, "Undocumented"),
    colGenSubProbSolMode("colGenSubProbSolMode", 2, "Undocumented"),
    MaxNbOfBBtreeNodeTreated("MaxNbOfBBtreeNodeTreated", 100000, "Undocumented"),
    optimalityGapTolerance("optimalityGapTolerance", 0.000001, "Undocumented"),
    relOptimalityGapTolerance("relOptimalityGapTolerance", 0.000000001, "Undocumented"),
    mastInitMode("mastInitMode", MasterInitMode::localArtCol, "Undocumented"),
    DEFAULTPRINTLEVEL("DEFAULTPRINTLEVEL", -1, "Undocumented"),
    printMasterPrimalSols("printMasterPrimalSols", 0, "Undocumented"),
    ApplyPreprocessing("ApplyPreprocessing", true, "Undocumented"),
    MaxNbOfStagesInColGenProcedure("MaxNbOfStagesInColGenProcedure", 1, "Undocumented"),
    GenerateProperColumns("GenerateProperColumns", false, "Undocumented"),
    PreprocessVariablesLocalBounds("PreprocessVariablesLocalBounds", true, "Undocumented"),
    SolverSelectForMast("SolverSelectForMast", 'a', "Undocumented"),
    SafeDualBoundScaleFactor("SafeDualBoundScaleFactor", -1, "Undocumented"),
    MipSolverMultiThread("MipSolverMultiThread", 0, "Undocumented"),
    treeSearchStrategy("treeSearchStrategy", 1, "Undocumented"),
    StrongBranchingPhaseOne("StrongBranchingPhaseOne", "", "Undocumented"),
    StrongBranchingPhaseTwo("StrongBranchingPhaseTwo", "", "Undocumented"),
    StrongBranchingPhaseThree("StrongBranchingPhaseThree", "", "Undocumented"),
    StrongBranchingPhaseFour("StrongBranchingPhaseFour", "", "Undocumented"),
    SimplifiedStrongBranchingParameterisation("SimplifiedStrongBranchingParameterisation", false, "Undocumented"),
    StrongBranchingPhaseOneCandidatesNumber("StrongBranchingPhaseOneCandidatesNumber", 100, "Undocumented"),
    StrongBranchingPhaseOneTreeSizeEstimRatio("StrongBranchingPhaseOneTreeSizeEstimRatio", 0.3, "Undocumented"),
    StrongBranchingPhaseTwoCandidatesNumber("StrongBranchingPhaseTwoCandidatesNumber", 3, "Undocumented"),
    StrongBranchingPhaseTwoTreeSizeEstimRatio("StrongBranchingPhaseTwoTreeSizeEstimRatio", 0.1, "Undocumented"),
    MasterCuttingPlanesDepthLimit("MasterCuttingPlanesDepthLimit", 1000, "Undocumented"),
    MaxNbOfCutGeneratedAtEachIter("MaxNbOfCutGeneratedAtEachIter", 1000, "Undocumented"),
    OpenNodesLimit("OpenNodesLimit", 1000, "Undocumented"),
    ArtVarPenaltyUpdateFactor("ArtVarPenaltyUpdateFactor", 2, "Undocumented"),
    ArtVarMaxNbOfPenaltyUpdates("ArtVarMaxNbOfPenaltyUpdates", 5, "Undocumented"),
    colGenStabilizationFunctionType("colGenStabilizationFunctionType", 0, "Undocumented"),
    colGenProximalStabilizationRule("colGenProximalStabilizationRule", 1, "Undocumented"),
    StabilFuncKappa("StabilFuncKappa", 1, "Undocumented"),
    colGenDualPriceSmoothingAlphaFactor("colGenDualPriceSmoothingAlphaFactor", 1, "Undocumented"),
    colGenDualPriceSmoothingBetaFactor("colGenDualPriceSmoothingBetaFactor", 0, "Undocumented"),
    colGenStabilizationMaxTreeDepth("colGenStabilizationMaxTreeDepth", 10000, "Undocumented"),
    StabilizationMinPhaseOfStage("StabilizationMinPhaseOfStage", 0, "Undocumented"),
    InsertAllGeneratedColumnsInFormRatherThanInPool("InsertAllGeneratedColumnsInFormRatherThanInPool", true, "Undocumented"),
    InsertNewNonNegColumnsDirectlyInFormRatherThanInPool("InsertNewNonNegColumnsDirectlyInFormRatherThanInPool", true, "Undocumented"),
    ColumnCleanupThreshold("ColumnCleanupThreshold", 10000, "Undocumented"),
    ColumnCleanupRatio("ColumnCleanupRatio", 0.66, "Undocumented"),
    CutCleanupThreshold("CutCleanupThreshold", 1, "Undocumented"),
    CutTailingOffThreshold("CutTailingOffThreshold", 0.02, "Undocumented"),
    CutTailingOffCounterThreshold("CutTailingOffCounterThreshold", 3, "Undocumented"),
    CutCleanupRatio("CutCleanupRatio", 0.66, "Undocumented"),
    ColGenSpRelaxationImprovementPriority("ColGenSpRelaxationImprovementPriority", 0.0, "Undocumented"),
    ReducedCostFixingThreshold("ReducedCostFixingThreshold", 0.9, "Undocumented"),
    UseInitialPrimalHeur("UseInitialPrimalHeur", false, "Undocumented"),
    MaxLocalSearchIterationCounter("MaxLocalSearchIterationCounter", 3, "Undocumented"),
    EvalAlgParamsInDiving("EvalAlgParamsInDiving", "", "Undocumented"),
    MaxFactorOfColFixedByLocalSearchHeur("MaxFactorOfColFixedByLocalSearchHeur", 0.8, "Undocumented"),
    StrongDivingCandidatesNumber("StrongDivingCandidatesNumber", 1, "Undocumented"),
    MaxNbOfCgIteDuringRh("MaxNbOfCgIteDuringRh", 5000, "Undocumented"),
    MaxLDSbreadth("MaxLDSbreadth", 0, "Undocumented"),
    MaxLDSdepth("MaxLDSdepth", 0, "Undocumented"),
    DivingHeurUseDepthLimit("DivingHeurUseDepthLimit", -1, "Undocumented"),
    LocalSearchHeurUseDepthLimit("LocalSearchHeurUseDepthLimit", -1, "Undocumented"),
    DivingHeurStopsWithFirstFeasSol("DivingHeurStopsWithFirstFeasSol", false, "Undocumented"),
    DivingHeurPreprocessBeforeChoosingVar("DivingHeurPreprocessBeforeChoosingVar", false, "Undocumented"),
    FixIntValBeforeRoundingHeur("FixIntValBeforeRoundingHeur", true, "Undocumented"),
    RoundingColSelectionCriteria("RoundingColSelectionCriteria", "Undocumented"),
    LocalSearchColSelectionCriteria("LocalSearchColSelectionCriteria", "Undocumented"),
    MaxTimeForRestrictedMasterIpHeur("MaxTimeForRestrictedMasterIpHeur", -1, "Undocumented"),
    PolishingAfterTimeInRestrictedMasterIpHeur("PolishingAfterTimeInRestrictedMasterIpHeur", -1, "Undocumented"),
    MaxNumEnumSolsInRestrictedMasterIpHeur("MaxNumEnumSolsInRestrictedMasterIpHeur", 5000, "Undocumented"),
    MIPemphasisInRestrictedMasterIpHeur("MIPemphasisInRestrictedMasterIpHeur", 1, "Undocumented"),
    CallFrequencyOfRestrictedMasterIpHeur("CallFrequencyOfRestrictedMasterIpHeur", 0, "Undocumented"),
    CallFrequencyOfDivingHeur("CallFrequencyOfDivingHeur", 0, "Undocumented"),

    RCSPuseBidirectionalSearch("RCSPuseBidirectionalSearch", 2, "Undocumented"),
    RCSPstopCutGenTimeThresholdInPricing("RCSPstopCutGenTimeThresholdInPricing", 10, "Undocumented"),
    RCSPhardTimeThresholdInPricing("RCSPhardTimeThresholdInPricing", 20, "Undocumented"),
    RCSPredCostFixingTimeThreshold("RCSPredCostFixingTimeThreshold", 100, "Undocumented"),
    RCSPmaxNumOfColsPerIteration("RCSPmaxNumOfColsPerIteration", 30, "Undocumented"),
    RCSPmaxNumOfColsPerExactIteration("RCSPmaxNumOfColsPerExactIteration", 150, "Undocumented"),
    RCSPallowRoutesWithSameVerticesSet("RCSPallowRoutesWithSameVerticesSet", true, "Undocumented"),
    RCSPnumberOfBucketsPerVertex("RCSPnumberOfBucketsPerVertex", 25, "Undocumented"),
    RCSPdynamicBucketSteps("RCSPdynamicBucketSteps", 1, "Undocumented"),
    RCSPapplyReducedCostFixing("RCSPapplyReducedCostFixing", 1, "Undocumented"),
    RCSPredCostFixingFalseGap("RCSPredCostFixingFalseGap", 0.0, "Undocumented"),
    RCSPinitNGneighbourhoodSize("RCSPinitNGneighbourhoodSize", 8, "Undocumented"),
    RCSPmaxNGneighbourhoodSize("RCSPmaxNGneighbourhoodSize", 0, "Undocumented"),
    RCSPrankOneCutsMaxNumRows("RCSPrankOneCutsMaxNumRows", 5, "Undocumented"),
    RCSPrankOneCutsMaxNumPerRound("RCSPrankOneCutsMaxNumPerRound", 100, "Undocumented"),
    RCSPrankOneCutsMemoryType("RCSPrankOneCutsMemoryType", 2, "Undocumented"),
    RCSPrankOneCutsLSnumIterations("RCSPrankOneCutsLSnumIterations", 1000, "Undocumented"),
    RCSPmaxNumOfLabelsInEnumeration("RCSPmaxNumOfLabelsInEnumeration", 1000000, "Undocumented"),
    RCSPmaxNumOfLabelsInHeurEnumeration("RCSPmaxNumOfLabelsInHeurEnumeration", 0, "Undocumented"),
    RCSPmaxNumOfEnumeratedSolutions("RCSPmaxNumOfEnumeratedSolutions", 1000000, "Undocumented"),
    RCSPmaxNumOfEnumSolutionsForMIP("RCSPmaxNumOfEnumSolutionsForMIP", 10000, "Undocumented"),
    RCSPmaxNumOfEnumSolsForEndOfNodeMIP("RCSPmaxNumOfEnumSolsForEndOfNodeMIP", 0, "Undocumented")
{
    int f[4] = {SelectionStrategy::LeastFractional, SelectionStrategy::Closest2RoundUp,
                SelectionStrategy::LeastGreedyCost, SelectionStrategy::LeastCost};
    for (int i = 0; i < 4; i++)
        RoundingColSelectionCriteria().push_back(f[i]);
    for (int i = 0; i < 4; i++)
        LocalSearchColSelectionCriteria().push_back(f[i]);
}

UserControlParameters::~UserControlParameters()
{
}

void UserControlParameters::addParameters(ParameterManager& parameterManager)
{
    parameterManager.addParameter(GlobalTimeLimitInTick);
    parameterManager.addParameter(GlobalTimeLimit);
    parameterManager.addParameter(ClockType);
    parameterManager.addParameter(colGenSubProbSolMode, false);
    parameterManager.addParameter(MipSolverMaxBBNodes);
    parameterManager.addParameter(MipSolverMaxTime);
    parameterManager.addParameter(MaxNbOfStagesInColGenProcedure);
    parameterManager.addParameter(MaxNbOfBBtreeNodeTreated);
    parameterManager.addParameter(ApplyPreprocessing);
    parameterManager.addParameter(mastInitMode, false);
    parameterManager.addParameter(GenerateProperColumns);
    parameterManager.addParameter(DEFAULTPRINTLEVEL);
    parameterManager.addParameter(printMasterPrimalSols);
    parameterManager.addParameter(solverName);
    parameterManager.addParameter(optimalityGapTolerance);
    parameterManager.addParameter(relOptimalityGapTolerance);
    parameterManager.addParameter(SimplifiedStrongBranchingParameterisation);
    parameterManager.addParameter(PreprocessVariablesLocalBounds);
    parameterManager.addParameter(SolverSelectForMast);
    parameterManager.addParameter(SafeDualBoundScaleFactor);
    parameterManager.addParameter(MipSolverMultiThread);
    parameterManager.addParameter(treeSearchStrategy);
    parameterManager.addParameter(StrongBranchingPhaseOne, false);
    parameterManager.addParameter(StrongBranchingPhaseTwo, false);
    parameterManager.addParameter(StrongBranchingPhaseThree, false);
    parameterManager.addParameter(StrongBranchingPhaseFour, false);
    parameterManager.addParameter(StrongBranchingPhaseOneCandidatesNumber);
    parameterManager.addParameter(StrongBranchingPhaseOneTreeSizeEstimRatio);
    parameterManager.addParameter(StrongBranchingPhaseTwoCandidatesNumber);
    parameterManager.addParameter(StrongBranchingPhaseTwoTreeSizeEstimRatio);
    parameterManager.addParameter(MasterCuttingPlanesDepthLimit);
    parameterManager.addParameter(MaxNbOfCutGeneratedAtEachIter);
    parameterManager.addParameter(OpenNodesLimit);
    parameterManager.addParameter(ArtVarPenaltyUpdateFactor);
    parameterManager.addParameter(ArtVarMaxNbOfPenaltyUpdates);
    parameterManager.addParameter(colGenStabilizationFunctionType, false);
    parameterManager.addParameter(colGenProximalStabilizationRule, false);
    parameterManager.addParameter(StabilFuncKappa);
    parameterManager.addParameter(colGenDualPriceSmoothingAlphaFactor);
    parameterManager.addParameter(colGenDualPriceSmoothingBetaFactor);
    parameterManager.addParameter(colGenStabilizationMaxTreeDepth);
    parameterManager.addParameter(StabilizationMinPhaseOfStage);
    parameterManager.addParameter(InsertAllGeneratedColumnsInFormRatherThanInPool);
    parameterManager.addParameter(InsertNewNonNegColumnsDirectlyInFormRatherThanInPool);
    parameterManager.addParameter(ColumnCleanupThreshold);
    parameterManager.addParameter(ColumnCleanupRatio);
    parameterManager.addParameter(CutCleanupThreshold);
    parameterManager.addParameter(CutTailingOffThreshold);
    parameterManager.addParameter(CutTailingOffCounterThreshold);
    parameterManager.addParameter(CutCleanupRatio);
    parameterManager.addParameter(ColGenSpRelaxationImprovementPriority);
    parameterManager.addParameter(ReducedCostFixingThreshold);
    parameterManager.addParameter(MaxLocalSearchIterationCounter);
    parameterManager.addParameter(EvalAlgParamsInDiving, false);
    parameterManager.addParameter(MaxFactorOfColFixedByLocalSearchHeur);
    parameterManager.addParameter(StrongDivingCandidatesNumber);
    parameterManager.addParameter(MaxNbOfCgIteDuringRh);
    parameterManager.addParameter(MaxLDSbreadth);
    parameterManager.addParameter(MaxLDSdepth);
    parameterManager.addParameter(UseInitialPrimalHeur);
    parameterManager.addParameter(DivingHeurUseDepthLimit);
    parameterManager.addParameter(LocalSearchHeurUseDepthLimit);
    parameterManager.addParameter(DivingHeurStopsWithFirstFeasSol);
    parameterManager.addParameter(DivingHeurPreprocessBeforeChoosingVar);
    parameterManager.addParameter(FixIntValBeforeRoundingHeur);
    parameterManager.addParameter(RoundingColSelectionCriteria, false);
    parameterManager.addParameter(LocalSearchColSelectionCriteria, false);
    parameterManager.addParameter(CallFrequencyOfRestrictedMasterIpHeur);
    parameterManager.addParameter(CallFrequencyOfDivingHeur);
    parameterManager.addParameter(MaxTimeForRestrictedMasterIpHeur);
    parameterManager.addParameter(PolishingAfterTimeInRestrictedMasterIpHeur);
    parameterManager.addParameter(MaxNumEnumSolsInRestrictedMasterIpHeur);
    parameterManager.addParameter(MIPemphasisInRestrictedMasterIpHeur);
    parameterManager.addParameter(RCSPuseBidirectionalSearch);
    parameterManager.addParameter(RCSPstopCutGenTimeThresholdInPricing);
    parameterManager.addParameter(RCSPhardTimeThresholdInPricing);
    parameterManager.addParameter(RCSPredCostFixingTimeThreshold);
    parameterManager.addParameter(RCSPmaxNumOfColsPerIteration);
    parameterManager.addParameter(RCSPmaxNumOfColsPerExactIteration);
    parameterManager.addParameter(RCSPallowRoutesWithSameVerticesSet);
    parameterManager.addParameter(RCSPnumberOfBucketsPerVertex);
    parameterManager.addParameter(RCSPdynamicBucketSteps);
    parameterManager.addParameter(RCSPapplyReducedCostFixing);
    parameterManager.addParameter(RCSPredCostFixingFalseGap);
    parameterManager.addParameter(RCSPinitNGneighbourhoodSize);
    parameterManager.addParameter(RCSPmaxNGneighbourhoodSize);
    parameterManager.addParameter(RCSPrankOneCutsMaxNumRows);
    parameterManager.addParameter(RCSPrankOneCutsMaxNumPerRound);
    parameterManager.addParameter(RCSPrankOneCutsMemoryType);
    parameterManager.addParameter(RCSPrankOneCutsLSnumIterations);
    parameterManager.addParameter(RCSPmaxNumOfLabelsInEnumeration);
    parameterManager.addParameter(RCSPmaxNumOfLabelsInHeurEnumeration);
    parameterManager.addParameter(RCSPmaxNumOfEnumeratedSolutions);
    parameterManager.addParameter(RCSPmaxNumOfEnumSolutionsForMIP);
    parameterManager.addParameter(RCSPmaxNumOfEnumSolsForEndOfNodeMIP);
}

std::ostream& UserControlParameters::printVRPSolverParameters(std::ostream& os) const
{
    os << "--- VRPSOLVER PARAMETERS ---" << endl;
    os << GlobalTimeLimit << endl;
    os << MaxNbOfBBtreeNodeTreated << endl;
    os << treeSearchStrategy << endl;

    os << RCSPstopCutGenTimeThresholdInPricing << endl;
    os << RCSPhardTimeThresholdInPricing << endl;
    os << RCSPredCostFixingTimeThreshold << endl;

    os << RCSPnumberOfBucketsPerVertex << endl;
    os << RCSPdynamicBucketSteps << endl;

    os << RCSPuseBidirectionalSearch << endl;
    os << RCSPapplyReducedCostFixing << endl;
    os << RCSPredCostFixingFalseGap << endl;

    os << RCSPmaxNumOfColsPerIteration << endl;
    os << RCSPmaxNumOfColsPerExactIteration << endl;

    os << StabilizationMinPhaseOfStage << endl;

    os << RCSPmaxNumOfLabelsInEnumeration << endl;
    os << RCSPmaxNumOfEnumeratedSolutions << endl;
    os << RCSPmaxNumOfEnumSolutionsForMIP << endl;

    os << RCSPinitNGneighbourhoodSize << endl;
    os << RCSPmaxNGneighbourhoodSize << endl;

    os << RCSPrankOneCutsMaxNumPerRound << endl;
    os << RCSPrankOneCutsMaxNumRows << endl;
    os << RCSPrankOneCutsMemoryType << endl;

    os << CutTailingOffThreshold << endl;
    os << CutTailingOffCounterThreshold << endl;

    os << SafeDualBoundScaleFactor << endl;

    os << StrongBranchingPhaseOneCandidatesNumber << endl;
    os << StrongBranchingPhaseOneTreeSizeEstimRatio << endl;
    os << StrongBranchingPhaseTwoCandidatesNumber << endl;
    os << StrongBranchingPhaseTwoTreeSizeEstimRatio << endl;

    os << MaxTimeForRestrictedMasterIpHeur << endl;
    os << DivingHeurUseDepthLimit << endl;
    os << MaxLDSbreadth << endl;
    os << MaxLDSdepth << endl;
  return os;
}

std::ostream& UserControlParameters::printParameters(std::ostream& os) const
{
    os << "--- MAIN PARAMETERS ---" << endl;
    os << GlobalTimeLimitInTick << endl;
    os << GlobalTimeLimit << endl;
    os << ClockType << endl;
    os << MaxNbOfBBtreeNodeTreated << endl;
    os << optimalityGapTolerance << endl;
    os << relOptimalityGapTolerance << endl;
    os << ApplyPreprocessing << endl;
    os << PreprocessVariablesLocalBounds << endl;
    os << treeSearchStrategy << endl;
    os << OpenNodesLimit << endl;
    os << DEFAULTPRINTLEVEL << endl;
    os << "--- MIP SOLVER PARAMETERS ---" << endl;
    os << solverName << endl;
    os << MipSolverMaxBBNodes  << endl;
    os << MipSolverMaxTime << endl;
    os << MipSolverMultiThread << endl;
    os << "--- COLUMN GENERATION PARAMETERS ---" << endl;
    os << SolverSelectForMast << endl;
    os << "  " << colGenSubProbSolMode << endl;
    os << mastInitMode << endl;
    os << ArtVarPenaltyUpdateFactor << endl;
    os << ArtVarMaxNbOfPenaltyUpdates << endl;
    os << MaxNbOfStagesInColGenProcedure << endl;
    os << GenerateProperColumns << endl;
    os << InsertAllGeneratedColumnsInFormRatherThanInPool << endl;
    os << InsertNewNonNegColumnsDirectlyInFormRatherThanInPool << endl;
    os << ColumnCleanupThreshold << endl;
    os << ColumnCleanupRatio << endl;
    os << ReducedCostFixingThreshold << endl;
    os << "--- CUT GENERATION PARAMETERS ---" << endl;
    os << MasterCuttingPlanesDepthLimit << endl;
    os << MaxNbOfCutGeneratedAtEachIter << endl;
    os << CutTailingOffThreshold << endl;
    os << CutTailingOffCounterThreshold << endl;
    os << CutCleanupThreshold << endl;
    os << CutCleanupRatio << endl;
    os << ColGenSpRelaxationImprovementPriority << endl;
    os << "--- STABILIZATION PARAMETERS ---" << endl;
    os << colGenDualPriceSmoothingAlphaFactor << endl;
    os << colGenDualPriceSmoothingBetaFactor << endl;
    os << colGenStabilizationFunctionType << endl;
    os << colGenProximalStabilizationRule << endl;
    os << StabilFuncKappa << endl;
    os << colGenStabilizationMaxTreeDepth << endl;
    os << StabilizationMinPhaseOfStage << endl;
    os << "--- PRIMAL HEURISTIC PARAMETERS ---" << endl;
    os << UseInitialPrimalHeur << endl;
    os << MaxTimeForRestrictedMasterIpHeur << endl;
    os << CallFrequencyOfRestrictedMasterIpHeur << endl;
    os << MIPemphasisInRestrictedMasterIpHeur << endl;
    os << PolishingAfterTimeInRestrictedMasterIpHeur << endl;
    os << DivingHeurUseDepthLimit << endl;
    os << CallFrequencyOfDivingHeur << endl;
    for (MultitokenSelectionStrategyVector::const_iterator it = RoundingColSelectionCriteria().begin();
         it != RoundingColSelectionCriteria().end(); ++it)
        os << *it << endl;
    os << FixIntValBeforeRoundingHeur << endl;
    os << MaxNbOfCgIteDuringRh << endl;
    os << MaxLDSbreadth << endl;
    os << MaxLDSdepth << endl;
    os << DivingHeurStopsWithFirstFeasSol << endl;
    os << DivingHeurPreprocessBeforeChoosingVar << endl;
    os << StrongDivingCandidatesNumber << endl;
    if (EvalAlgParamsInDiving().active())
        os << "EvalAlgParamsInDiving =" << EvalAlgParamsInDiving() << std::endl;
    os << LocalSearchHeurUseDepthLimit << endl;
    for (MultitokenSelectionStrategyVector::const_iterator it = LocalSearchColSelectionCriteria().begin();
         it != LocalSearchColSelectionCriteria().end(); ++it)
        os << *it << endl;
    os << MaxFactorOfColFixedByLocalSearchHeur << endl;
    os << MaxLocalSearchIterationCounter << endl;
    os << "--- STRONG BRANCHING PARAMETERS ---" << endl;
    os << SimplifiedStrongBranchingParameterisation << endl;
    if (SimplifiedStrongBranchingParameterisation())
    {
        os << StrongBranchingPhaseOneCandidatesNumber << std::endl;
        os << StrongBranchingPhaseOneTreeSizeEstimRatio << std::endl;
        os << StrongBranchingPhaseTwoCandidatesNumber << std::endl;
        os << StrongBranchingPhaseTwoTreeSizeEstimRatio << std::endl;
    } else
    {
        if (StrongBranchingPhaseOne().active())
            os << "StrongBranchingPhaseOne =" << StrongBranchingPhaseOne() << std::endl;
        if (StrongBranchingPhaseTwo().active())
            os << "StrongBranchingPhaseTwo =" << StrongBranchingPhaseTwo() << std::endl;
        if (StrongBranchingPhaseThree().active())
            os << "StrongBranchingPhaseThree =" << StrongBranchingPhaseThree() << std::endl;
        if (StrongBranchingPhaseFour().active())
            os << "StrongBranchingPhaseFour =" << StrongBranchingPhaseFour() << std::endl;
    }
    os << SafeDualBoundScaleFactor << endl;
    os << "StrongBranchingPhaseOne =" << StrongBranchingPhaseOne() << std::endl;
    os << "StrongBranchingPhaseTwo =" << StrongBranchingPhaseTwo() << std::endl;
    os << "StrongBranchingPhaseThree =" << StrongBranchingPhaseThree() << std::endl;
    os << "StrongBranchingPhaseFour =" << StrongBranchingPhaseFour() << std::endl;
    os << "--- DEBUG OUTPUT PARAMETERS ---" << endl;
    os << printMasterPrimalSols << endl;
    os << "--- VRPSOLVER PARAMETERS ---" << endl;
    os << RCSPstopCutGenTimeThresholdInPricing << endl;
    os << RCSPhardTimeThresholdInPricing << endl;
    os << RCSPredCostFixingTimeThreshold << endl;
    os << RCSPnumberOfBucketsPerVertex << endl;
    os << RCSPdynamicBucketSteps << endl;
    os << RCSPuseBidirectionalSearch << endl;
    os << RCSPapplyReducedCostFixing << endl;
    os << RCSPmaxNumOfColsPerIteration << endl;
    os << RCSPmaxNumOfColsPerExactIteration << endl;
    os << RCSPmaxNumOfLabelsInEnumeration << endl;
    os << RCSPmaxNumOfLabelsInHeurEnumeration << endl;
    os << RCSPmaxNumOfEnumeratedSolutions << endl;
    os << RCSPmaxNumOfEnumSolutionsForMIP << endl;
    os << RCSPmaxNumOfEnumSolsForEndOfNodeMIP << endl;
    os << RCSPinitNGneighbourhoodSize << endl;
    os << RCSPmaxNGneighbourhoodSize << endl;
    os << RCSPrankOneCutsMaxNumRows << endl;
    os << RCSPrankOneCutsMaxNumPerRound << endl;
    os << RCSPrankOneCutsMemoryType << endl;
    os << RCSPrankOneCutsLSnumIterations << endl;
    os << RCSPallowRoutesWithSameVerticesSet << endl;
    os << RCSPredCostFixingFalseGap << endl;
    os << SafeDualBoundScaleFactor << endl;
    os << RCSPmaxNumOfLabelsInHeurEnumeration << endl;
    os << MaxNumEnumSolsInRestrictedMasterIpHeur << endl;

    return os;
}

void UserControlParameters::postTreatment()
{
}
