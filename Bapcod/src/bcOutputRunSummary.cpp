/**
 *
 * This file bcOutputRunSummary.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcOutputRunSummary.hpp"

OutputRunSummary::OutputRunSummary(BapcodInit* bapcodInitPtr, const std::string& instanceName)
{
  if (bapcodInitPtr != NULL)
    {
      if (!bapcodInitPtr->param().NameOfOutputSummaryFile().empty())
        {
          std::ofstream ofs(bapcodInitPtr->param().NameOfOutputSummaryFile().c_str(), std::ofstream::out | std::ofstream::app);

          if (!ofs.is_open())
            {
              throw GlobalException(std::string("Cannot print the solution into the ") +
                  "file: " + bapcodInitPtr->param().NameOfOutputSummaryFile());
            }
          ofs << "Model: " << instanceName << std::endl;
          ofs << (bapcodInitPtr->param().bestDualBoundIsConstantForTest() ? "CONSTANT " : "") << "bestDualBound: " << bapcodInitPtr->statistics().getValue("bcRecBestDb")  << std::endl;
          ofs << (bapcodInitPtr->param().bestIncumbentIsConstantForTest() ? "CONSTANT " : "") << "bestIncumbent: " << bapcodInitPtr->statistics().getValue("bcRecBestInc") << std::endl;

          ofs << "bcCountMastSol: " << bapcodInitPtr->statistics().getCounter("bcCountMastSol") << std::endl;
          ofs << "bcCountSpSol: " << bapcodInitPtr->statistics().getCounter("bcCountSpSol") << std::endl;
          ofs << "bcCountNodeProc: " << bapcodInitPtr->statistics().getCounter("bcCountNodeProc") << std::endl;
          ofs << "bcTimeMain: " << bapcodInitPtr->statistics().getRecord("bcTimeMain")._time << std::endl;
          ofs << "bcTimeMastMPsol: " << bapcodInitPtr->statistics().getRecord("bcTimeMastMPsol")._time << std::endl;
          ofs << "bcTimeSpSol: " <<  bapcodInitPtr->statistics().getRecord("bcTimeSpSol")._time << std::endl;
          ofs << std::endl;
        }
    }
  else
    {
      throw GlobalException(std::string("Cannot print the solution: BapcodInit or the ") +
          "bapcodInitPtr is NULL");
    }
}


OutputRunSummary::~OutputRunSummary()
{
}

