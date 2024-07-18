/**
 *
 * This file bcControlParameters.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcControlParameters.hpp"
#include "bcDevControlParameters.hpp"

#include "bcBoundLevC.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include "bcModelParameterC.hpp"

ControlParameters::ControlParameters() : UserControlParameters(),
                                         DevControlParameters()//1e-6, 1e-8)  // precision, highprecision
{
}

ControlParameters::~ControlParameters()
{
}

void ControlParameters::addParameters(ParameterManager& parameterManager)
{
  UserControlParameters::addParameters(parameterManager);
  DevControlParameters::addParameters(parameterManager);
}

std::ostream& ControlParameters::printDevParameters(std::ostream& os) const
{
  UserControlParameters::printParameters(os);
  DevControlParameters::printParameters(os);
  return os;
}

std::ostream& ControlParameters::printUserParameters(std::ostream& os) const
{
    UserControlParameters::printParameters(os);
    return os;
}

std::ostream& ControlParameters::printVRPSolverParameters(std::ostream& os) const
{
    UserControlParameters::printVRPSolverParameters(os);
    return os;
}

void ControlParameters::postTreatment()
{
    UserControlParameters::postTreatment();
    DevControlParameters::postTreatment();

  /*********************************************
             Parameters Consistency
  *********************************************/

//  if (printL(0))
//  {
//    std::cout << "Info about modified params should be ignored for params which are not in the config file."
//              << std::endl;
//    std::cout << "Such warnings may have been triggered because the default value of the param was modified."
//              << std::endl;
//  }

    if (solverName() == "CLP_SOLVER") /// CLP is an LP solver only
    {
        if (masterSolMode().status() != SolutionMethod::lpSolver || RCSPmaxNumOfEnumSolutionsForMIP() > 0
            || MaxTimeForRestrictedMasterIpHeur() > 0)
        {
            if (printL(0))
            {
                std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
                std::cout << "     solverName == CLP_SOLVER " << std::endl;
                std::cout << "     ===> " << std::endl;
                std::cout << "     masterSolMode = 1 && RCSPmaxNumOfEnumSolutionsForMIP = 0 "
                          << "&& MaxTimeForRestrictedMasterIpHeur == 0" << std::endl;
            }
            masterSolMode().set(SolutionMethod::lpSolver);
            RCSPmaxNumOfEnumSolutionsForMIP(0);
            MaxTimeForRestrictedMasterIpHeur(0);
        }

    }

    if ((MaxTimeForRestrictedMasterIpHeur() > 0) || (RCSPmaxNumOfEnumSolutionsForMIP() > 0))
    {
        if (masterSolMode().status() != SolutionMethod::mipSolver)
        {
            if (printL(0))
            {
                std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
                std::cout << "     MaxTimeForRestrictedMasterIpHeur > 0 || RCSPmaxNumOfEnumSolutionsForMIP > 0 "
                          << std::endl;
                std::cout << "     ===> " << std::endl;
                std::cout << "     masterSolMode = 2" << std::endl;
            }
            masterSolMode().set(SolutionMethod::mipSolver);
        }
    }

    /// if the master is initialized by an incumbent solution, we need to launch initial primal heurisitic
    /// to obtain this incumbent solution
    switch (mastInitMode().status())
    {
        case MasterInitMode::incSolCol:
        case MasterInitMode::incSolColAndGac:
        case MasterInitMode::incSolColAndLac:
            if (UseInitialPrimalHeur() != true)
            {
                if (printL(0))
                {
                    std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
                    std::cout << "     mastInitMode in {incSolCol, incSolColAndGac, incSolColAndLac} " << std::endl;
                    std::cout << "     ===> " << std::endl;
                    std::cout << "     UseInitialPrimalHeur = true " << std::endl;
                }
                UseInitialPrimalHeur(true);
            }
            break;
        default:
            break;
    }


    /// stabilization function parameter coherence
    switch (colGenStabilizationFunctionType().status())
    {
        case StabilizationFunctionType::none :
            if (StabilFuncInnerAngle() != 0 || StabilFuncOuterAngle() != 0)
            {
                StabilFuncInnerAngle(0);
                StabilFuncOuterAngle(0);
            }
            break;
        case StabilizationFunctionType::boxStep :
            if (StabilFuncInnerAngle() != 1.01 || StabilFuncOuterAngle() != 0)
            {
                StabilFuncInnerAngle(1.01);
                StabilFuncOuterAngle(0);
            }
            break;
        case StabilizationFunctionType::threePiece :
            if (StabilFuncOuterAngle() != 0)
            {
                StabilFuncOuterAngle(0);
            }
            break;
        case StabilizationFunctionType::undefined :
            if (colGenStabilizationFunctionType().status() != StabilizationFunctionType::none)
            {
                colGenStabilizationFunctionType().set(StabilizationFunctionType::none);
            }
            break;
        default:
            break;
    }

    if ((colGenStabilizationFunctionType().status() == StabilizationFunctionType::none)
        && (colGenProximalStabilizationRule().status() != ColGenProximalStabilizationMode::undefined))
        colGenProximalStabilizationRule().set(ColGenProximalStabilizationMode::undefined);
    if (StabilFuncOuterHalfInterval() < 0)
        StabilFuncOuterHalfInterval(100);
    if (StabilFuncInnerHalfInterval() < 0)
        StabilFuncOuterHalfInterval(StabilFuncOuterHalfInterval()/10);
    if (StabilFuncInnerHalfInterval() > StabilFuncOuterHalfInterval())
        StabilFuncInnerHalfInterval(StabilFuncOuterHalfInterval);
    if (StabilFuncOuterAngle() < 0)
        StabilFuncOuterAngle(0.9);
    if (StabilFuncInnerAngle() < 0)
        StabilFuncInnerAngle(0.1);
    if (StabilFuncCurvature <= 0)
        StabilFuncCurvature(1);
    if ((StabilFuncCurvatureAdvanceRate <= 0) || (StabilFuncCurvatureAdvanceRate() > 1))
        StabilFuncCurvatureAdvanceRate(1);
    if (ArtVarPenaltyUpdateFactor <= 1)
        ArtVarPenaltyUpdateFactor(1.2);
    if (ArtVarMaxNbOfPenaltyUpdates() < 0)
        ArtVarMaxNbOfPenaltyUpdates(20);
    if (StabilFuncArtVarInSolUpdateFactor <= 1)
        StabilFuncArtVarInSolUpdateFactor(10);

    if (debugsolution_file() != "")
        /// if debug solution is used then inactive columns which form this dual solution will be generated
        /// UseColumnPool should be then set to true so that the coefficients of these columns in the master constraints
        /// are calculated
        UseColumnsPool(true);

    if ( (colGenDualPriceSmoothingAlphaFactor() < 0) || (colGenDualPriceSmoothingAlphaFactor() > 1) )
        colGenDualPriceSmoothingAlphaFactor(0);

    if (colGenDualPriceSmoothingAlphaFactor() == 0)
        colGenDualPriceSmoothingBetaFactor(0);

    if ( (colGenDualPriceSmoothingBetaFactor() < 0) || (colGenDualPriceSmoothingBetaFactor() > 1) )
        colGenDualPriceSmoothingBetaFactor(0);

    if (colGenStabilizationMaxTreeDepth() < 0)
        colGenStabilizationMaxTreeDepth(0);

    if ((colGenStabilizationFunctionType().status() == StabilizationFunctionType::none) &&
        (colGenDualPriceSmoothingAlphaFactor() == 0))
        colgeninfo_file("");

  /// if strong diving or strong branching or reduced cost fixing is used, then the column generation
  /// should converge completely
  if ((StrongDivingCandidatesNumber() > 1)
      || (ReducedCostFixingThreshold() > 0) || (StrongBranchingPhaseOne().active()))
  {
    if (TerminateCgWhenRoundedDbCannotImprove())
    {
        if (printL(0))
          std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl
                    << "     StrongDivingCandidatesNumber > 1 || ReducedCostFixingThreshold > 0 "
                    << " || StrongBranchingPhaseOne.active" << std::endl
                    << "     ===> " << std::endl
                    << "     TerminateCgWhenRoundedDbCannotImprove = false" << std::endl;
      TerminateCgWhenRoundedDbCannotImprove(false);
    }
  }


  /// preprocessing is switched off for the col gen for ext form methods
  /// (we need the full formulation for a correct preprocessing)

  if (SplitColIntoDissagregateSpVar())
    {
      if (ApplyPreprocessing() || (DivingHeurUseDepthLimit() >= 0))
      {
          if (printL(0))
          {
              std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
              std::cout << "     SplitColIntoDissagregateSpVar = true" << std::endl;
              std::cout << "     ===> " << std::endl;
              std::cout << "     ApplyPreprocessing = false" << std::endl;
              std::cout << "     DivingHeurUseDepthLimit = -1" << std::endl;
          }
        ApplyPreprocessing(false);
        DivingHeurUseDepthLimit(-1);
      }
    }

  /// if we deal with proper columns only, we need to apply preprocessing in sub-problems
  /// however the opposite is not true (we can have GenerateProperColumns = false and PreprocessVariablesLocalBounds = true)
  /// the last case appears for example in the "non-proper" diving
  if (GenerateProperColumns())
  {
    if(!PreprocessVariablesLocalBounds())
    {
        if (printL(0))
        {
            std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
            std::cout << "     GenerateProperColumns = true" << std::endl;
            std::cout << "     ===> " << std::endl;
            std::cout << "     PreprocessVariablesLocalBounds = true" << std::endl;
        }
      PreprocessVariablesLocalBounds(true);
    }
  }

  if (DivingHeurUseDepthLimit() >= 0)
  {
    if(!ApplyPreprocessing())
    {
        if (printL(0))
        {
            std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
            std::cout << "     DivingHeurUseDepthLimit >= 0" << std::endl;
            std::cout << "     ===> " << std::endl;
            std::cout << "     ApplyPreprocessing = true" << std::endl;
        }
        ApplyPreprocessing(true);
    }
  }

  if (StrongBranchingPhaseOne().active()
      && (masterSolMode().status() != SolutionMethod::lpSolver)
      && (masterSolMode().status() != SolutionMethod::mipSolver))
    {
      std::cerr << "BaPCod parameters compatibility error : cannot use strong branching by phases "
                << "when master problem is not solved by column generation " << std::endl;
      exit(1);
    }

  if (VerifyColsIntegralityInTestSolForIntegrality())
  {
    if (TestAggregateMasterSol4Integrality())
    {
        if (printL(0))
        {
            std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
            std::cout << "     VerifyColsIntegralityInTestSolForIntegrality = true" << std::endl;
            std::cout << "     ===> " << std::endl;
            std::cout << "     TestAggregateMasterSol4Integrality = false" << std::endl;
        }
      TestAggregateMasterSol4Integrality(false);
    }
  }


  if ((colGenSubProbSolMode().status() == SolutionMethod::custom2mipSolver)
      && (MaxNbOfStagesInColGenProcedure() <= 1))
  {
    if (MaxNbOfStagesInColGenProcedure() != 2)
    {
        if (printL(0))
        {
            std::cout << " BaPCod info - PARAM MODIFIED : " << std::endl;
            std::cout << "     colGenSubProbSolMode() == SolutionMethod::custom2mipSolver "
                         "&& MaxNbOfStagesInColGenProcedure <= 1 " << std::endl;
            std::cout << "     ===> " << std::endl;
            std::cout << "     MaxNbOfStagesInColGenProcedure = 2" << std::endl;
        }
      MaxNbOfStagesInColGenProcedure(2);
    }
  }

  if (RCSPprintLevel() < DEFAULTPRINTLEVEL())
      RCSPprintLevel(DEFAULTPRINTLEVEL());

  if (RCSPmaxNumOfEnumSolsForEndOfNodeMIP() < RCSPmaxNumOfEnumSolutionsForMIP())
      RCSPmaxNumOfEnumSolsForEndOfNodeMIP(RCSPmaxNumOfEnumSolutionsForMIP());

  if (RCSPmaxNGaverNeighbourhoodSize() < RCSPinitNGneighbourhoodSize())
      RCSPmaxNGaverNeighbourhoodSize(RCSPinitNGneighbourhoodSize());

  if (RCSPmaxNGneighbourhoodSize() < RCSPinitNGneighbourhoodSize())
      RCSPmaxNGneighbourhoodSize(RCSPinitNGneighbourhoodSize());

  if (RCSPmaxNGneighbourhoodSize() == RCSPinitNGneighbourhoodSize())
      RCSPdynamicNGmode(0);

  if (RCSPdynamicNGmode() != 0 && ColGenSpRelaxationImprovementPriority() <= 0.0)
      ColGenSpRelaxationImprovementPriority(2.0);

  if (RCSPresConsKnapsackCutsMode() < -10)
      RCSPresConsKnapsackCutsMode(-10);

  if (RCSPresConsKnapsackCutsMode() > 10)
      RCSPresConsKnapsackCutsMode(10);

  /// ApplyStrongBranchingEvaluation parameter is kept for backward compatibility
  /// ApplyStrongBranchingEvaluation = true means that BapCod has been called from the VRPSolver
  if (ApplyStrongBranchingEvaluation())
  {
      SimplifiedStrongBranchingParameterisation(true);
      CallFrequencyOfDivingHeur(1);
      if (CallFrequencyOfRestrictedMasterIpHeur() <= 0)
          CallFrequencyOfRestrictedMasterIpHeur(5);
      PreprocessVariablesLocalBounds(false);
      EvalAlgParamsInDiving().setNonExact(1, 10000, 1,
                                          0, 0, 0, 0,
                                          1);
  }

  if (SimplifiedStrongBranchingParameterisation() && !StrongBranchingPhaseOne().active()) {
      StrongBranchingPhaseOne().setNonExact(StrongBranchingPhaseOneCandidatesNumber(),
                                            0, 0, 0,
                                            0, false, 0,
                                            StrongBranchingPhaseOneTreeSizeEstimRatio());
      StrongBranchingPhaseTwo().setNonExact(StrongBranchingPhaseTwoCandidatesNumber(),
                                            10000, 1, 0,
                                            0, false, 0,
                                            StrongBranchingPhaseTwoTreeSizeEstimRatio());
      StrongBranchingPhaseThree().setExact();
      StrongBranchingPhaseFour().setNonActive();
  }

  if (colGenSubProbSolMode().status() == SolutionMethod::mipSolver)
      GenerateProperColumns(true);

  if (GlobalTimeLimit() > 0)
  {
      GlobalTimeLimitInTick(GlobalTimeLimit() * 100);
  }
  else
  {
      GlobalTimeLimit(GlobalTimeLimitInTick() / 100);
  }
}
