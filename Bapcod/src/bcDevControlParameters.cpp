/**
 *
 * This file bcDevControlParameters.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcDevControlParameters.hpp"
#include "bcParameterManager.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include "bcModelParameterC.hpp"

using namespace std;

DevControlParameters::DevControlParameters():
  DEFAULTTESTLEVEL("DEFAULTTESTLEVEL", 1, "Undocumented"),
  MipSolverRecordNamesInFormulation("MipSolverRecordNamesInFormulation", true, "Undocumented"),
  ovfSolMode("ovfSolMode", SolutionMethod::none, "Undocumented"),
  masterSolMode("masterSolMode", SolutionMethod::mipSolver, "Undocumented"),
  VRPSEimposeCapacityResourceByCuts("VRPSEimposeCapacityResourceByCuts", false, "Undocumented"),
  VRPSEcriticalResource("VRPSEcriticalResource", -1, "Undocumented"),
  MaxNbOfCgIterations("MaxNbOfCgIterations", 100000, "Undocumented"),
  Search4NegRedCostColInInactivePool("Search4NegRedCostColInInactivePool", false, "Undocumented"),
  RequiredSolStatForOvf("RequiredSolStatForOvf", "Undocumented"),
  PreprocessorOnForOvf("PreprocessorOnForOvf", true, "Undocumented"),
  ProbingOnForOvf("ProbingOnForOvf", true, "Undocumented"),
  AutomaticCuttingPlanesOnForOvf("AutomaticCuttingPlanesOnForOvf", true, "Undocumented"),
  SolverSelectForOvf("SolverSelectForOvf", 'd', "Undocumented"),
  RequiredSolStatForMast("RequiredSolStatForMast", "Undocumented"),
  PreprocessorOnForMast("PreprocessorOnForMast", true, "Undocumented"),
  ProbingOnForMast("ProbingOnForMast", true, "Undocumented"),
  AutomaticCuttingPlanesOnForMast("AutomaticCuttingPlanesOnForMast", true, "Undocumented"),
  SolverSelectForMastAfterBarrierConvergence("SolverSelectForMastAfterBarrierConvergence", 'd',
                                             "Undocumented"),
  RequiredSolStatForColGenSp("RequiredSolStatForColGenSp", "Undocumented"),
  PreprocessorOnForColGenSp("PreprocessorOnForColGenSp", true, "Undocumented"),
  ProbingOnForColGenSp("ProbingOnForColGenSp", true, "Undocumented"),
  AutomaticCuttingPlanesOnForColGenSp("AutomaticCuttingPlanesOnForColGenSp", true, "Undocumented"),
  SolverSelectForColGenSp("SolverSelectForColGenSp", 'd', "Undocumented"),
  MasterMipSolverRightHandSideZeroTol("MasterMipSolverRightHAndSideZeroTol", 1e-7, "Undocumented"),
  MasterMipSolverReducedCostTolerance("MasterMipSolverReducedCostTolerance", 1e-7, "Undocumented"),
  MasterMipSolverBarrierConvergenceTolerance("MasterMipSolverBarrierConvergenceTolerance", 1e-2,
                                             "Undocumented"),
  ColGenSpMipSolverRightHAndSideZeroTol("ColGenMipSolverRightHAndSideZeroTol", 1e-8, "Undocumented"),
  ColGenSpMipSolverReducedCostTolerance("ColGenMipSolverReducedCostTolerance", 1e-6, "Undocumented"),
  MaxDepthInBBtree("MaxDepthInBBtree", 5000, "Undocumented"),
  StopWhenBranchingFails("StopWhenBranchingFails", true, "Undocumented"),
  BapCodReducedCostTolerance("BapCodReducedCostTolerance", 1e-6, "Undocumented"),
  BapCodCutViolationTolerance("BapCodCutViolationTolerance", 0.02, "Undocumented"),
  BapCodIntegralityTolerance("BapCodIntegralityTolerance", 1e-8, "Undocumented"),
  TestAggregateMasterSol4Integrality("TestAggregateMasterSol4Integrality", false, "Undocumented"),
  VerifyColsIntegralityInTestSolForIntegrality("VerifyColsIntegralityInTestSolForIntegrality", false,
                                               "Undocumented"),
  BranchFirstOnPureMasterVariables("BranchFirstOnPureMasterVariables", true, "Undocumented"),
  BranchFirstOnSubProblemVariables("BranchFirstOnSubProblemVariables", false,
                                   "Undocumented"),
  ApplyStrongBranchingEvaluation("ApplyStrongBranchingEvaluation", false, "Undocumented"),
  StrongBranchingParameter4NotPromissingConservativeness("StrongBranchingParameter4NotPromissingConservativeness",
                                                         1, "Undocumented"),
  StrongBranchingExactPhaseInterrupt("StrongBranchingExactPhaseInterrupt", false,
                                     "Undocumented"),
  StrongBranchingUseHistorySubtreeSize("StrongBranchingUseHistorySubtreeSize", false,
                                       "Undocumented"),
  StrongBranchingUseHistory("StrongBranchingUseHistory", true, "Undocumented"),
  TerminateCgWhenRoundedDbCannotImprove("TerminateCgWhenRoundedDbCannotImprove", false,
                                        "Undocumented"),
  UseObjScalingFact("UseObjScalingFact", false, "Undocumented"),
  LocArtVarInConvexityConstr("LocArtVarInConvexityConstr", true, "Undocumented"),
  StabilFuncArtVarInSolUpdateFactor("StabilFuncArtVarInSolUpdateFactor", 10,
                                    "Undocumented"),
  StabilFuncCurvature("StabilFuncCurvature", 1, "Undocumented"),
  StabilFuncCurvatureAdvanceRate("StabilFuncCurvatureAdvanceRate", 1, "Undocumented"),
  StabilFuncOuterHalfInterval("StabilFuncOuterHalfInterval", 100, "Undocumented"),
  StabilFuncInnerHalfInterval("StabilFuncInnerHalfInterval", 10, "Undocumented"),
  StabilFuncHalfIntervalChildNodeFactor("StabilFuncHalfIntervalChildNodeFactor", 10,
                                        "Undocumented"),
  StabilFuncOuterAngle("StabilFuncOuterAngle", 0.9, "Undocumented"),
  StabilFuncInnerAngle("StabilFuncInnerAngle", 0.9, "Undocumented"),
  colGenDualPriceSmoothingMinGap("colGenDualPriceSmoothingMinGap", 0, "Undocumented"),
  colGenDualPriceSmoothingMaxGap("colGenDualPriceSmoothingMaxGap", 0, "Undocumented"),
  colGenDualPriceSmoothingMaxNbOfUpdate("colGenDualPriceSmoothingMaxNbOfUpdate", 10,
                                        "Undocumented"),
  SplitColIntoDissagregateSpVar("SplitColIntoDissagregateSpVar", false, "Undocumented"),
  SplitColIntoDissagregateSpVarInHeadIn("SplitColIntoDissagregateSpVarInHeadIn", true,
                                        "Undocumented"),
  PriceAllSubproblemsForBestDualBound("PriceAllSubproblemsForBestDualBound", true,
                                      "Undocumented"),
  CyclicSpScanning("CyclicSpScanning", false, "Undocumented"),
  UseColumnsPool("UseColumnsPool", false, "Undocumented"),
  CutTailingOffAbsTolerance("CutTailingOffAbsTolerance", -1),
  MinNumOfCutRoundsBeforeStopBySp("MinNumOfCutRoundsBeforeStopBySp", 1, "Undocumented"),
  CutDynamicPriorityLevel("CutDynamicPriorityLevel", false, "Undocumented"),
  runColGenUntilFullConvergence("runColGenUntilFullConvergence", false, "Undocumented"),
  UseGreedyHeur("UseGreedyHeur", 0, "Undocumented"),
  UseCustomFracSolBasedHeur("UseCustomFracSolBasedHeur", false, "Undocumented"),
  runColGenAfterFixingPureMastVarInDiving("runColGenAfterFixingPureMastVarInDiving", true,
                                          "Undocumented"),
  EnumHeuristicFalseGapMultiplier("EnumHeuristicFalseGapMultiplier", 0.5, "Undocumented"),
  EnumHeuristicNumberOfTries("EnumHeuristicNumberOfTries", 7, "Undocumented"),
  UseDivingHeurOnMasterColOnly("UseDivingHeurOnMasterColOnly", true, "Undocumented"),
  UseDivingHeurOnPureMastVarOnly("UseDivingHeurOnPureMastVarOnly", false, "Undocumented"),
  MaxNbOfPenaltUpdatesDuringRH("MaxNbOfPenaltUpdatesDuringRH", 10, "Undocumented"),
  IgnoreIntValWhenSelectingRoundedVarInRH("IgnoreIntValWhenSelectingRoundedVarInRH", true,
                                          "Undocumented"),
  ActivateAllColumnsForRestrictedMasterIpHeur("ActivateAllColumnsForRestrictedMasterIpHeur", false,
                                              "Undocumented"),
  MaxNumEnumSolsInUserHeuristicFunctor("MaxNumEnumSolsInUserHeuristicFunctor", 10000,
                                       "Undocumented"),
  PricingStrategy("PricingStrategy", 0, "Undocumented"),
  MaxNbPromisingSpFound("MaxNbPromisingSpFound", 1000000, "Undocumented"),
  MaxNbUnpromisingSpFound("MaxNbUnpromisingSpFound", 1, "Undocumented"),
  NameOfOutputSummaryFile("NameOfOutputSummaryFile", "",
                          "Print the solutions to the filename (if non-empty)."),
  RCSPuseMetaSolver("RCSPuseMetaSolver", 0, "Undocumented"),
  RCSPlabelSplitStrategy("RCSPlabelSplitStrategy", 0, "Undocumented"),
  RCSPprintLevel("RCSPprintLevel", -2, "Undocumented"),
  RCSPcheckDominInOtherBuckets("RCSPcheckDominInOtherBuckets", true, "Undocumented"),
  RCSPuseCompletionBoundsInPricing("RCSPuseCompletionBoundsInPricing", 0, "Undocumented"),
  RCSPdomChecksThresholdInPricing("RCSPdomChecksThresholdInPricing", 1e15, "Undocumented"),
  RCSPheurLabelingStrategy("RCSPheurLabelingStrategy", 0, "Undocumented"),
  RCSPdynBuckStepAdjustRatioThreshold("RCSPdynBuckStepAdjustRatioThreshold", 500,
                                      "Undocumented"),
  RCSPdynBuckStepAdjustNumBuckArcsThreshold("RCSPdynBuckStepAdjustNumBuckArcsThreshold", 10000,
                                            "Undocumented"),
  RCSPuseComplBoundsInRedCostFixing("RCSPuseComplBoundsInRedCostFixing", 2, "Undocumented"),
  RCSPuseDSSRInMode("RCSPuseDSSRInMode", 0, "Undocumented"),
  RCSPmaxNumOfColumnsInDSSR("RCSPmaxNumOfColumnsInDSSR", 100, "Undocumented"),
  RCSPmaxCycleSizeInDSSR("RCSPmaxCycleSizeInDSSR", 5, "Undocumented"),
  RCSPdynamicNGmode("RCSPdynamicNGmode", -1, "Undocumented"),
  RCSPmaxNumOfColumnsInDynamicNG("RCSPmaxNumOfColumnsInDynamicNG", 100, "Undocumented"),
  RCSPmaxCycleSizeInDynamicNG("RCSPmaxCycleSizeInDynamicNG", 5, "Undocumented"),
  RCSPuseCapacityCuts("RCSPuseCapacityCuts", true, "Undocumented"),
  RCSPcapacityCutsSeparator("RCSPcapacityCutsSeparator", 1, "Undocumented"),
  RCSPcapCutsMaxNumPerRound("RCSPcapCutsMaxNumPerRound", 100, "Undocumented"),
  RCSPresConsKnapsackCutsMode("RCSPresConsKnapsackCutsMode", -1, "Undocumented"),
  RCSPresConsKnapsackCutsMaxNumPerRound("RCSPresConsKnapsackCutsMaxNumPerRound", 100, "Undocumented"),
  RCSPcliqueCutsMaxNumPerRound("RCSPcliqueCutsMaxNumPerRound", 0, "Undocumented"),
  RCSPcliqueCutsMaxSepTime("RCSPcliqueCutsMaxSepTime", 5.0, "Undocumented"),
  RCSPrankOneCutsSpDependentMemory("RCSPrankOneCutsSpDependentMemory", false,
                                   "Undocumented"),
  RCSPmaxNGaverNeighbourhoodSize("RCSPmaxNGaverNeighbourhoodSize", 0, "Undocumented"),
  RCSPrankOneCutsArcMemoryType("RCSPrankOneCutsArcMemoryType", 2, "Undocumented"),
  RCSPrankOneCutsNeighbourhoodSize("RCSPrankOneCutsNeighbourhoodSize", 16, "Undocumented"),
  RCSPrankOneCutsMaxNumTriplets("RCSPrankOneCutsMaxNumTriplets", 10000, "Undocumented"),
  RCSPrankOneCutsTypeToSeparate("RCSPrankOneCutsTypeToSeparate", 0, "Undocumented"),
  RCSPuseCovSetsForSepOneRowPackingCuts("RCSPuseCovSetsForSepOneRowPackingCuts", false, "Undocumented"),
  RCSPuseExactComplBoundsInEnumeration("RCSPuseExactComplBoundsInEnumeration", true,
                                       "Undocumented"),
  RCSPmaxGapToRunEnumeration("RCSPmaxGapToRunEnumeration", 1.0, "Undocumented"),
  RCSPuseVertexIdsInEnumeration("RCSPuseVertexIdsInEnumeration", false, "Undocumented"),
  RCSPstopRatioForConcatenationInEnum("RCSPstopRatioForConcatenationInEnum", 30,
                                      "Undocumented"),
  ColGenLogFrequency("ColGenLogFrequency", 1, "Undocumented"),
  timeKeySet("timeKeySet", "Undocumented"),
  outputKeyVector("outputKeyVector", "Undocumented"),
  configFileName("bcParameters,b", "NOT_SPECIFIED",
                 "Path to the Bapcod's parameters file."),
  instance_file("instance,i", "", "Path to the data instance file."),
  statistics_file("statistics,s", "", "Path to the file to output the statistics."),
  baPTreeDot_file("BaPTree,t", "BaPTree.dot",
                  "Path to the file to output the BaP Tree in dot format."),
  colgeninfo_file("colgeninfo,f", "",
                  "Path to the file to output the detailed information about col. gen. iterations"),
  debugsolution_file("debugsolution,d", "", "Path to the file which contains the debug solution"),
  fractionalSolutionDotFile("fractionalSolutionDotFile", "", "Undocumented"),
  bestDualBoundIsConstantForTest("bestDualBoundIsConstantForTest", true,
                                 "Set bestDualBound parameter to be tested or not."),
  bestIncumbentIsConstantForTest("bestIncumbentIsConstantForTest", true,
                                 "Set bestIncumbent parameter to be tested or not."),
  useSPVarMembershipCache("useSPVarMembershipCache", false,
                          "Allows fast reduced cost computation. "
                          "Uses more memory, and cannot be used with dynamic formulations."), /// to remove
  CheckSpOracleFeasibility("CheckSpOracleFeasibility", false, "Undocumented"),
  CheckOracleOptimality("CheckOracleOptimality", false, "Undocumented")
{
  int a[2] = {SolutionStatus::Optimum,SolutionStatus::PrimalFeasSolFound};
  RequiredSolStatForOvf =  SolutionStatus(a, a+2);

  int b[3] = {SolutionStatus::Optimum, SolutionStatus::DualFeasSolFound, SolutionStatus::OptimumUnscalInfeas};
  RequiredSolStatForMast =  SolutionStatus(b, b+3);

  int c[2] = {SolutionStatus::Optimum,SolutionStatus::PrimalFeasSolFound};
  RequiredSolStatForColGenSp =  SolutionStatus(c, c+2);
}

DevControlParameters::~DevControlParameters()
{
}

void DevControlParameters::addParameters(ParameterManager& parameterManager)
{
    parameterManager.addParameter(configFileName, false, false, true);
    parameterManager.addParameter(instance_file, false, false, true);
    parameterManager.addParameter(statistics_file, false, false, true);
    parameterManager.addParameter(baPTreeDot_file, false, false, true);
    parameterManager.addParameter(colgeninfo_file, false, false, true);

    parameterManager.addParameter(debugsolution_file);

    parameterManager.addParameter(fractionalSolutionDotFile);

    parameterManager.addParameter(DEFAULTTESTLEVEL);

    parameterManager.addParameter(VRPSEimposeCapacityResourceByCuts);
    parameterManager.addParameter(VRPSEcriticalResource);

    parameterManager.addParameter(ovfSolMode, false);
    parameterManager.addParameter(masterSolMode, false);
    parameterManager.addParameter(MipSolverRecordNamesInFormulation);
    parameterManager.addParameter(MaxNbOfCgIterations);
    parameterManager.addParameter(Search4NegRedCostColInInactivePool);

    parameterManager.addParameter(RequiredSolStatForOvf, true);
    parameterManager.addParameter(PreprocessorOnForOvf);
    parameterManager.addParameter(ProbingOnForOvf);
    parameterManager.addParameter(AutomaticCuttingPlanesOnForOvf);
    parameterManager.addParameter(SolverSelectForOvf);
    parameterManager.addParameter(RequiredSolStatForMast, true);
    parameterManager.addParameter(PreprocessorOnForMast);
    parameterManager.addParameter(ProbingOnForMast);
    parameterManager.addParameter(AutomaticCuttingPlanesOnForMast);
    parameterManager.addParameter(SolverSelectForMastAfterBarrierConvergence);
    parameterManager.addParameter(RequiredSolStatForColGenSp, true);
    parameterManager.addParameter(PreprocessorOnForColGenSp);
    parameterManager.addParameter(ProbingOnForColGenSp);
    parameterManager.addParameter(AutomaticCuttingPlanesOnForColGenSp);
    parameterManager.addParameter(SolverSelectForColGenSp);
    parameterManager.addParameter(MasterMipSolverRightHandSideZeroTol);
    parameterManager.addParameter(MasterMipSolverReducedCostTolerance);
    parameterManager.addParameter(MasterMipSolverBarrierConvergenceTolerance);
    parameterManager.addParameter(ColGenSpMipSolverRightHAndSideZeroTol);
    parameterManager.addParameter(ColGenSpMipSolverReducedCostTolerance);
    parameterManager.addParameter(MaxDepthInBBtree);
    parameterManager.addParameter(StopWhenBranchingFails);
    parameterManager.addParameter(BapCodReducedCostTolerance);
    parameterManager.addParameter(BapCodCutViolationTolerance);
    parameterManager.addParameter(BapCodIntegralityTolerance);
    parameterManager.addParameter(TestAggregateMasterSol4Integrality);
    parameterManager.addParameter(VerifyColsIntegralityInTestSolForIntegrality);
    parameterManager.addParameter(BranchFirstOnPureMasterVariables);
    parameterManager.addParameter(ApplyStrongBranchingEvaluation);
    parameterManager.addParameter(StrongBranchingParameter4NotPromissingConservativeness);
    parameterManager.addParameter(StrongBranchingExactPhaseInterrupt);
    parameterManager.addParameter(StrongBranchingUseHistorySubtreeSize);
    parameterManager.addParameter(StrongBranchingUseHistory);
    parameterManager.addParameter(TerminateCgWhenRoundedDbCannotImprove);
    parameterManager.addParameter(UseObjScalingFact);
    parameterManager.addParameter(LocArtVarInConvexityConstr);
    parameterManager.addParameter(StabilFuncArtVarInSolUpdateFactor);
    parameterManager.addParameter(StabilFuncCurvature);
    parameterManager.addParameter(StabilFuncCurvatureAdvanceRate);
    parameterManager.addParameter(StabilFuncOuterHalfInterval);
    parameterManager.addParameter(StabilFuncInnerHalfInterval);
    parameterManager.addParameter(StabilFuncHalfIntervalChildNodeFactor);
    parameterManager.addParameter(StabilFuncOuterAngle);
    parameterManager.addParameter(StabilFuncInnerAngle);
    parameterManager.addParameter(colGenDualPriceSmoothingMinGap);
    parameterManager.addParameter(colGenDualPriceSmoothingMaxGap);
    parameterManager.addParameter(colGenDualPriceSmoothingMaxNbOfUpdate);
    parameterManager.addParameter(SplitColIntoDissagregateSpVar);
    parameterManager.addParameter(SplitColIntoDissagregateSpVarInHeadIn);
    parameterManager.addParameter(PriceAllSubproblemsForBestDualBound);
    parameterManager.addParameter(CyclicSpScanning);
    parameterManager.addParameter(UseColumnsPool);
    parameterManager.addParameter(CutTailingOffAbsTolerance);
    parameterManager.addParameter(MinNumOfCutRoundsBeforeStopBySp);
    parameterManager.addParameter(CutDynamicPriorityLevel);
    parameterManager.addParameter(runColGenUntilFullConvergence);
    parameterManager.addParameter(UseGreedyHeur);
    parameterManager.addParameter(UseCustomFracSolBasedHeur);
    parameterManager.addParameter(runColGenAfterFixingPureMastVarInDiving);
    parameterManager.addParameter(EnumHeuristicFalseGapMultiplier);
    parameterManager.addParameter(EnumHeuristicNumberOfTries);
    parameterManager.addParameter(UseDivingHeurOnMasterColOnly);
    parameterManager.addParameter(UseDivingHeurOnPureMastVarOnly);
    parameterManager.addParameter(MaxNbOfPenaltUpdatesDuringRH);
    parameterManager.addParameter(IgnoreIntValWhenSelectingRoundedVarInRH);
    parameterManager.addParameter(ActivateAllColumnsForRestrictedMasterIpHeur);
    parameterManager.addParameter(MaxNumEnumSolsInUserHeuristicFunctor);

    parameterManager.addParameter(PricingStrategy);
    parameterManager.addParameter(MaxNbPromisingSpFound);
    parameterManager.addParameter(MaxNbUnpromisingSpFound);

    parameterManager.addParameter(RCSPuseMetaSolver);
    parameterManager.addParameter(RCSPlabelSplitStrategy);
    parameterManager.addParameter(RCSPprintLevel);
    parameterManager.addParameter(RCSPcheckDominInOtherBuckets);
    parameterManager.addParameter(RCSPuseCompletionBoundsInPricing);
    parameterManager.addParameter(RCSPdomChecksThresholdInPricing);
    parameterManager.addParameter(RCSPheurLabelingStrategy);
    parameterManager.addParameter(RCSPlabelSplitStrategy);
    parameterManager.addParameter(RCSPdynBuckStepAdjustRatioThreshold);
    parameterManager.addParameter(RCSPdynBuckStepAdjustNumBuckArcsThreshold);
    parameterManager.addParameter(RCSPuseComplBoundsInRedCostFixing);
    parameterManager.addParameter(RCSPuseDSSRInMode);
    parameterManager.addParameter(RCSPmaxNumOfColumnsInDSSR);
    parameterManager.addParameter(RCSPmaxCycleSizeInDSSR);
    parameterManager.addParameter(RCSPdynamicNGmode);
    parameterManager.addParameter(RCSPmaxNumOfColumnsInDynamicNG);
    parameterManager.addParameter(RCSPmaxCycleSizeInDynamicNG);
    parameterManager.addParameter(RCSPmaxNGaverNeighbourhoodSize);
    parameterManager.addParameter(RCSPuseCapacityCuts);
    parameterManager.addParameter(RCSPcapacityCutsSeparator);
    parameterManager.addParameter(RCSPcapCutsMaxNumPerRound);
    parameterManager.addParameter(RCSPresConsKnapsackCutsMode);
    parameterManager.addParameter(RCSPresConsKnapsackCutsMaxNumPerRound);
    parameterManager.addParameter(RCSPcliqueCutsMaxNumPerRound);
    parameterManager.addParameter(RCSPcliqueCutsMaxSepTime);
    parameterManager.addParameter(RCSPrankOneCutsArcMemoryType);
    parameterManager.addParameter(RCSPrankOneCutsSpDependentMemory);
    parameterManager.addParameter(RCSPrankOneCutsNeighbourhoodSize);
    parameterManager.addParameter(RCSPrankOneCutsMaxNumTriplets);
    parameterManager.addParameter(RCSPrankOneCutsTypeToSeparate);
    parameterManager.addParameter(RCSPuseCovSetsForSepOneRowPackingCuts);
    parameterManager.addParameter(RCSPuseExactComplBoundsInEnumeration);
    parameterManager.addParameter(RCSPmaxGapToRunEnumeration);
    parameterManager.addParameter(RCSPuseVertexIdsInEnumeration);
    parameterManager.addParameter(RCSPstopRatioForConcatenationInEnum);

    parameterManager.addParameter(ColGenLogFrequency);
    parameterManager.addParameter(timeKeySet, false);
    parameterManager.addParameter(outputKeyVector, true);

    parameterManager.addParameter(NameOfOutputSummaryFile);

    parameterManager.addParameter(bestDualBoundIsConstantForTest);
    parameterManager.addParameter(bestIncumbentIsConstantForTest);
    parameterManager.addParameter(useSPVarMembershipCache);

    parameterManager.addParameter(CheckOracleOptimality);
    parameterManager.addParameter(CheckSpOracleFeasibility);
}

std::ostream& DevControlParameters::printParameters(ostream& os) const
{
    os << "!-- ADVANCED PARAMETERS --!" << endl;
    os << "  " << ovfSolMode << endl;
    os << "  " << masterSolMode << endl;
    os << VRPSEimposeCapacityResourceByCuts << endl;
    os << VRPSEcriticalResource << endl;
    os << MipSolverRecordNamesInFormulation << endl;
    os << MaxNbOfCgIterations << endl;
    os << Search4NegRedCostColInInactivePool << endl;
    os << RequiredSolStatForOvf << endl;
    os << PreprocessorOnForOvf << endl;
    os << ProbingOnForOvf << endl;
    os << AutomaticCuttingPlanesOnForOvf << endl;
    os << SolverSelectForOvf << endl;
    os << RequiredSolStatForMast << endl;
    os << PreprocessorOnForMast << endl;
    os << ProbingOnForMast << endl;
    os << AutomaticCuttingPlanesOnForMast << endl;
    os << SolverSelectForMastAfterBarrierConvergence << endl;
    os << RequiredSolStatForColGenSp << endl;
    os << PreprocessorOnForColGenSp << endl;
    os << ProbingOnForColGenSp << endl;
    os << AutomaticCuttingPlanesOnForColGenSp << endl;
    os << MasterMipSolverRightHandSideZeroTol << endl;
    os << MasterMipSolverReducedCostTolerance << endl;
    os << ColGenSpMipSolverRightHAndSideZeroTol << endl;
    os << ColGenSpMipSolverReducedCostTolerance << endl;
    os << MaxDepthInBBtree << endl;
    os << StopWhenBranchingFails << endl;
    os << BapCodReducedCostTolerance << endl;
    os << BapCodCutViolationTolerance << endl;
    os << BapCodIntegralityTolerance << endl;
    os << TestAggregateMasterSol4Integrality << endl;
    os << VerifyColsIntegralityInTestSolForIntegrality << endl;
    os << BranchFirstOnPureMasterVariables << endl;
    os << BranchFirstOnSubProblemVariables << endl;
    os << StrongBranchingParameter4NotPromissingConservativeness << endl;
    os << StrongBranchingExactPhaseInterrupt << endl;
    os << StrongBranchingUseHistorySubtreeSize << endl;
    os << StrongBranchingUseHistory << endl;
    os << TerminateCgWhenRoundedDbCannotImprove << endl;
    os << UseObjScalingFact << endl;
    os << LocArtVarInConvexityConstr << endl;
    os << StabilFuncArtVarInSolUpdateFactor << endl;
    os << StabilFuncCurvature << endl;
    os << StabilFuncCurvatureAdvanceRate << endl;
    os << StabilFuncOuterHalfInterval << endl;
    os << StabilFuncInnerHalfInterval << endl;
    os << StabilFuncHalfIntervalChildNodeFactor << endl;
    os << StabilFuncOuterAngle << endl;
    os << StabilFuncInnerAngle << endl;
    os << colGenDualPriceSmoothingMinGap << endl;
    os << colGenDualPriceSmoothingMaxGap << endl;
    os << colGenDualPriceSmoothingMaxNbOfUpdate << endl;
    os << SplitColIntoDissagregateSpVar << endl;
    os << SplitColIntoDissagregateSpVarInHeadIn << endl;
    os << PriceAllSubproblemsForBestDualBound << endl;
    os << CyclicSpScanning << endl;
    os << UseColumnsPool << endl;
    os << CutTailingOffAbsTolerance << endl;
    os << MinNumOfCutRoundsBeforeStopBySp << endl;
    os << CutDynamicPriorityLevel << endl;
    os << runColGenUntilFullConvergence << endl;
    os << UseGreedyHeur << endl;
    os << UseCustomFracSolBasedHeur << endl;
    os << ActivateAllColumnsForRestrictedMasterIpHeur << endl;
    os << MaxNumEnumSolsInUserHeuristicFunctor << endl;
    os << UseDivingHeurOnMasterColOnly << endl;
    os << UseDivingHeurOnPureMastVarOnly << endl;
    os << MaxNbOfPenaltUpdatesDuringRH << endl;
    os << IgnoreIntValWhenSelectingRoundedVarInRH << endl;
    os << useSPVarMembershipCache << endl;
    os << CheckSpOracleFeasibility << endl;
    os << CheckOracleOptimality << endl;
    os << ColGenLogFrequency << endl;
    os << RCSPuseMetaSolver << endl;
    os << RCSPlabelSplitStrategy << endl;
    os << RCSPprintLevel << endl;
    os << RCSPcheckDominInOtherBuckets << endl;
    os << RCSPuseCompletionBoundsInPricing << endl;
    os << RCSPdomChecksThresholdInPricing << endl;
    os << RCSPheurLabelingStrategy << endl;
    os << RCSPlabelSplitStrategy << endl;
    os << RCSPdynBuckStepAdjustRatioThreshold << endl;
    os << RCSPdynBuckStepAdjustNumBuckArcsThreshold << endl;
    os << RCSPuseComplBoundsInRedCostFixing << endl;
    os << RCSPuseDSSRInMode << endl;
    os << RCSPmaxNumOfColumnsInDSSR << endl;
    os << RCSPmaxCycleSizeInDSSR << endl;
    os << RCSPdynamicNGmode << endl;
    os << RCSPmaxNumOfColumnsInDynamicNG << endl;
    os << RCSPmaxCycleSizeInDynamicNG << endl;
    os << RCSPmaxNGaverNeighbourhoodSize << endl;
    os << RCSPuseCapacityCuts << endl;
    os << RCSPcapacityCutsSeparator << endl;
    os << RCSPresConsKnapsackCutsMode << endl;
    os << RCSPresConsKnapsackCutsMaxNumPerRound << endl;
    os << RCSPcapCutsMaxNumPerRound << endl;
    os << RCSPcliqueCutsMaxNumPerRound << endl;
    os << RCSPcliqueCutsMaxSepTime << endl;
    os << RCSPrankOneCutsArcMemoryType << endl;
    os << RCSPrankOneCutsSpDependentMemory << endl;
    os << RCSPrankOneCutsNeighbourhoodSize << endl;
    os << RCSPrankOneCutsMaxNumTriplets << endl;
    os << RCSPrankOneCutsTypeToSeparate << endl;
    os << RCSPuseCovSetsForSepOneRowPackingCuts << endl;
    os << RCSPuseExactComplBoundsInEnumeration << endl;
    os << RCSPmaxGapToRunEnumeration << endl;
    os << RCSPuseVertexIdsInEnumeration << endl;
    os << RCSPstopRatioForConcatenationInEnum << endl;

    for (set< string>::const_iterator it = timeKeySet().begin(); it != timeKeySet().end(); ++it)
        os <<  "time key = " << *it << endl;
    for (vector< string>::const_iterator it = outputKeyVector().begin(); it != outputKeyVector().end(); ++it)
        os <<  "output key = " << *it << endl;

    return os;
}

void DevControlParameters::postTreatment()
{
  //postTreatment is now done in bcControlParameters, since it involves both user and dev params.
}
