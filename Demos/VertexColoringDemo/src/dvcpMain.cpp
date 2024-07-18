/**
 *
 * This file dvcpMain.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcModelingLanguageC.hpp"
#include "dvcpParameters.hpp"
#include "dvcpData.hpp"
#include "dvcpHeuristics.hpp"
#include "dvcpModel.hpp"

void printSolution(BcSolution solPtr)
{
    std::cout << "------------------------------------------ " << std::endl;
    if (!solPtr.defined())
    {
        std::cout << "Solution is not defined!" << std::endl;
        std::cout << "------------------------------------------ " << std::endl;
        return;
    }
    std::cout << "The best found solution of value " << solPtr.cost() << " : " << std::endl;
    
    solPtr = solPtr.next(); /// skip master solution, go straight to the subproblem solutions
    
    int stabSetIndex = 0;
    while (solPtr.defined())
    {
        std::cout << "Color " << ++stabSetIndex << " : ";
        std::set< BcVar > varSet;
        solPtr.getVar(varSet);
        
        for (std::set< BcVar >::iterator varIt = varSet.begin(); varIt != varSet.end(); ++varIt)
        {
            if (varIt->genericName() == "X")
              std::cout << varIt->id().first() + 1 << " ";
        }
        std::cout << std::endl;
        
        solPtr = solPtr.next();
    }
    std::cout << "------------------------------------------ " << std::endl;
}


int main(int argc, char *argv[])
{
    BcInitialisation bapcodInit(argc, argv);
    VCPdemoParameters applSpecParam(bapcodInit.configFile(), argc, argv);
    applSpecParam.print();

    DvcpData data;
    data.readData(bapcodInit.instanceFile());

    VCPsolution initialSolution;
    switch (applSpecParam.initialHeuristic())
    {
        case DSATUR_INIT_HEURISTIC :
#ifndef _MSC_VER
            runDsaturGreedyHeuristic(data, applSpecParam.printLevelDVCP(), initialSolution);
            break;
#endif
        case SIMPLE_INIT_HEURISTIC:
            runSimpleGreedyHeuristic(data, applSpecParam.printLevelDVCP(), initialSolution);
            break;
        case NONE_INIT_HEURISTIC :
        default:
            break;
    }
    BcModel model(bapcodInit, bapcodInit.instanceFile());
    buildVCPModel(data, initialSolution, model);
    BcSolution solPtr = model.solve();
    printSolution(solPtr);

    if (applSpecParam.printLevelDVCP() >= 1)
    {
        std::cout << "Main solution statistics : " << std::endl;
        std::cout << "Total time = " << bapcodInit.getStatisticTime("bcTimeMain") / 100 << " sec." << std::endl;
        std::cout << "Root time = " << bapcodInit.getStatisticTime("bcTimeRootEval") / 100 << " sec." << std::endl;
        std::cout << "Restricted master solution time = " << bapcodInit.getStatisticTime("bcTimeMastMPsol") / 100
                  << " sec." << std::endl;
        std::cout << "Pricing problem solution time = " << bapcodInit.getStatisticTime("bcTimeCgSpOracle")  / 100
                  << " sec." << std::endl;
        std::cout << "Number of nodes = " << bapcodInit.getStatisticCounter("bcCountNodeProc") << std::endl;
        std::cout << "Number of times the restricted master is solved = "
                  << bapcodInit.getStatisticCounter("bcCountMastSol") << std::endl;
        std::cout << "Number of times the pricing problem is solved = " << bapcodInit.getStatisticCounter("bcCountCg")
                  << std::endl;
        std::cout << "Number of pricing subproblems solved = " << bapcodInit.getStatisticCounter("bcCountSpSol")
                  << std::endl;
        std::cout << "Number of generated columns = " << bapcodInit.getStatisticCounter("bcCountCol") << std::endl;
        std::cout << "Best primal bound value = " << bapcodInit.getStatisticValue("bcRecBestInc") << std::endl;
        std::cout << "Best dual bound value = " << bapcodInit.getStatisticValue("bcRecBestDb") << std::endl;
        std::cout << "Dual bound value obtained at the root = " << bapcodInit.getStatisticValue("bcRecRootDb")
                  << std::endl;
    }

    /// this is needed for correct execution of tests
    bapcodInit.outputResultsForTests();

    return EXIT_SUCCESS;
}

