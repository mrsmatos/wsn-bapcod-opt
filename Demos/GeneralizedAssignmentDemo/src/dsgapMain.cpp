/**
 *
 * This file dsgapMain.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <iostream>
#include <cstdlib>

#include "bcModelingLanguageC.hpp"
#include "bcApplicationException.hpp"
#include "dsgapDataC.hpp"
#include "dsgapModelC.hpp"

using namespace std;

void solveInstance(bool solveMIP, BcInitialisation & bapcodInit, GapData & dataStruct, std::string pbName)
{
    BcModel gapModel(bapcodInit, pbName);

    if (solveMIP)
    {
        BcFormulation mpForm = buildMpModel(dataStruct, gapModel);
        BcSolution solution = mpForm.solve();
        std::cout << solution << std::endl;
    }
    else
    {
        buildColGenModel(dataStruct, gapModel);
        BcSolution solution = gapModel.solve();
        std::cout << solution << std::endl;

        BcMaster master(gapModel);
        BcConstrArray covConstr(master, "COV");
        for (int jobIndex = 0; jobIndex < dataStruct.jobPts().size(); ++jobIndex)
            std::cout << "Dual value for constraint " << jobIndex << " is "
                      << ((BcConstr)covConstr[jobIndex]).curDualVal() << std::endl;
    }
    bapcodInit.outputBaPCodStatistics(pbName);
}

void testSimpleLP()
{
    BcInitialisation bapcodInit;
    bapcodInit.param().masterSolMode(SolutionMethod::lpSolver);

    BcModel model(bapcodInit, "simpleLP");
    BcFormulation form(model);
    BcObjective objective(form);
    BcVarArray vars(form, "X");
    BcConstrArray constrs(form, "Cstr");
    constrs(0) >= 25;
    constrs[0] += 3 * vars(0) + 4 * vars(1);
    constrs(1) >= 40;
    constrs[1] += 4 * vars(0) + 3 * vars(1);
    objective += 1 * vars(0) + 2 * vars(1);

    BcSolution solution = form.solve();

    std::cout << "Solution value is " << solution.cost() << std::endl;
    std::cout << "Dual value of constraint 0 is " << ((BcConstr)constrs[0]).curDualVal() << std::endl;
    std::cout << "Dual value of constraint 1 is " << ((BcConstr)constrs[1]).curDualVal() << std::endl;
    return;
}

int main(int argc, char *argv[])
{
  //Initialize a BaPCod context and read the file bcParameters.cfg
  BcInitialisation bapcodInit(argc, argv);

  ApplicationSpecificParam param("config/dsgapParameters.cfg");

  try
    {
      param.parse(argc, argv); //parse all the parameters.
    }
  catch (const ApplicationException& e) //if the file dsgapParameters.cfg is not found, an exception is thrown
    {
      //Print the message of the exception e.
      cerr << e.message() << endl;
    }

    if (param.testSimpleLp())
    {
        testSimpleLP();
        return 0;
    }

    param.print(); //print all the parameters and their values.

  if (!param.readProb())
    {
      /// random instance generation
      int nbProblems(param.nbProblems());
      for (int currentProblem = param.firstInstance();
           currentProblem < nbProblems; ++currentProblem)
        {
          GapData dataStruct(&param);
          dataStruct.generateAtRandom(currentProblem);
            solveInstance(param.justSolveMIP(), bapcodInit, dataStruct, "RandomData");
        }
    }
  else
    {
      /// reading data from a file
      if (bapcodInit.instanceFile() != "")
        {
          GapData dataStruct(&param);

          dataStruct.readData(bapcodInit.instanceFile());

            solveInstance(param.justSolveMIP(), bapcodInit, dataStruct, bapcodInit.instanceFile());
        }
      else
        {
          ifstream ifs(param.inputFile().c_str(), ios::in);
          //if the file (by default data/inputList) is not found.
          if (!ifs)
            {
              throw ApplicationException("Could not find " + param.inputFile());
            }

          while (!ifs.eof())
            {
              std::string probName;
              ifs >> probName;

              if (!probName.compare(""))
                {
                  break;
                }

              /// Data reading
              GapData dataStruct(&param);

              dataStruct.readData(probName);

              solveInstance(param.justSolveMIP(), bapcodInit, dataStruct, probName);

            }
        }
    }
  return 0;
}

