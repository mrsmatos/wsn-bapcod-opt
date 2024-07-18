/**
 *
 * This file dcspMain.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#include "dcspModelC.hpp"
#include "dcspParameters.hpp"
#include "bcModelingLanguageC.hpp"
#include "bcBapcodUtilities.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  /// The file containing bapcod parameters is parsed and the initialisation is done.
  BcInitialisation bapcodInit(argc, argv);
  
  /// The file containing the application parameters is parsed.
  ApplicationSpecificParam applicationParam("config/dcspParameters.cfg", argc, argv);

  string inputFileName = "data/inputList";

  ifstream is(inputFileName.c_str(), std::ios::in);

  if (applicationParam.readProb())
    {
      if (!is)
        {
          cerr << "Could not find " << inputFileName << endl;
          exit (EXIT_FAILURE);
        }
    }

  int nbProblems = 0;
  if (applicationParam.readProb())
    {
      if (!(is >> nbProblems))
	{
	  cerr << "Could not read nbData." << endl;
	  exit(EXIT_FAILURE);
	}
      applicationParam.firstInstance(0);      
    }
  else
    {
      nbProblems = applicationParam.nbInstances();
    }

  //For all problems
  for (int currentProblem = applicationParam.firstInstance(); currentProblem < nbProblems; ++currentProblem)
    {
      DcspData data(&applicationParam);

      string pbName = "";
      if (!applicationParam.readProb())
        {
          pbName = string("RandomData") + currentProblem + "-" + nbProblems;
          data.generateAtRandom(currentProblem); 
          data.modelName(pbName);
        }
      else
        {
	  if (!(is >> pbName))
	    {
	      cerr << "Could not read problem name: " << currentProblem << endl;
	      exit(EXIT_FAILURE);
	    }
	  data.readData(pbName);
        }

      cout << " NEXT PROBLEM n0 " << currentProblem << " IS " << pbName << endl;


      BcModel dcspModel(bapcodInit, "cspProb");
      buildModel(data, dcspModel);
      BcSolution solution = dcspModel.solve();

      cout << solution << endl;

      bapcodInit.outputBaPCodStatistics(pbName);
    }

  return EXIT_SUCCESS;
}

