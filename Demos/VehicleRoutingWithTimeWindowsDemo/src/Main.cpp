/**
 *
 * This file Main.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#include <bcModelPointerC.hpp>
#include <bcModelMasterC.hpp>

#include "Parameters.h"
#include "Loader.h"

#include "Model.h"
#include "SolutionChecker.h"

using namespace std;

int main(int argc, char** argv)
{
	BcInitialisation bapcodInit(argc, argv);

	vrpstw::Loader loader;
	if (!loader.loadParameters(bapcodInit.configFile(), argc, argv)
	    || !loader.loadData(bapcodInit.instanceFile()))
		return -1;

    vrpstw::SolutionChecker * sol_checker = new vrpstw::SolutionChecker;

    vrpstw::Model model(bapcodInit);
	model.attach(sol_checker);

	BcSolution solution(model.master());
	solution = model.solve();
	sol_checker->isFeasible(solution, true, true);

	if (bapcodInit.statFile() != "none")
	{
		std::ofstream fs(bapcodInit.statFile().c_str(), std::ios::app);
        bapcodInit.outputBaPCodStatistics(bapcodInit.instanceFile(), std::cout, fs,
                                          bapcodInit.statFile() + std::string(".out"));
		fs.close();
	}
	else
	{
        bapcodInit.outputBaPCodStatistics(bapcodInit.instanceFile());
    }

	return 0;
}
