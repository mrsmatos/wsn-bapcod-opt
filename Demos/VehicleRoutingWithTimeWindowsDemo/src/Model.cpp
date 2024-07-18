/**
 *
 * This file Model.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <Parameters.h>

#include "Model.h"
#include "bcModelingLanguageC.hpp"
#include "Data.h"
#include "RCSPSolver.h"

vrpstw::Model::Model(const BcInitialisation& bc_init) : BcModel(bc_init, bc_init.instanceFile())
{
    BcObjective objective(*this);
    if (data.roundType == Data::NO_ROUND || data.roundType == Data::ROUND_ONE_DECIMAL)
        objective.setStatus(BcObjStatus::minFloat);
    else
	    objective.setStatus(BcObjStatus::minInt);

    double upperBound = params.cutOffValue();
	if (upperBound != std::numeric_limits<double>::infinity())
	{
        if (data.roundType == Data::NO_ROUND || data.roundType == Data::ROUND_ONE_DECIMAL)
            objective <= upperBound + 0.1;
        else
            objective <= upperBound + 1;
	}
    if (std::abs(upperBound) < 1e4)
        objective.setArtCostValue(std::abs(upperBound));
    else
        objective.setArtCostValue(1e4);

	std::vector<int> customerIds;
	for (int custId = 1; custId <= data.nbCustomers; ++custId)
	    customerIds.push_back(custId);

	std::vector<int> depotIds;
    for (int depotId = 1; depotId <= data.nbDepots; ++depotId)
        depotIds.push_back(depotId);

	BcMaster master(*this);

    BcConstrArray degreeConstr(master, "DEG");
    for (auto custId : customerIds)
    {
        degreeConstr(custId) == 2;
    }

    BcBranchingConstrArray vehNumberBranching(master, "VNB", SelectionStrategy::MostFractional, 1.0);
    if (data.nbDepots > 1)
        vehNumberBranching(0); /// this is for the total number of vehicles
    for (auto depotId : depotIds)
        vehNumberBranching(depotId);

    BcBranchingConstrArray edgeBranching(master, "EDGE", SelectionStrategy::MostFractional, 1.0);
    for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
            edgeBranching(firstCustId, secondCustId);

    BcColGenSpArray depotCGSp(*this);
    for (auto depotId : depotIds)
    {
        BcFormulation spForm = depotCGSp(depotId);
        spForm <= data.depots[depotId].veh_number;

        BcVarArray xVar(spForm, "X");
        xVar.type('I');
        xVar.priorityForMasterBranching(-1);
        xVar.priorityForSubproblemBranching(-1);
        xVar.defineIndexNames(MultiIndexNames('d', 'i', 'j'));

        /// customer number 0 is the depot
        for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
            for (int secondCustId = firstCustId + 1; secondCustId <= data.nbCustomers; ++secondCustId)
            {
                BcVar bcVar = xVar(depotId, firstCustId, secondCustId);
                if (firstCustId == 0)
                {
                    if (data.nbDepots > 1)
                        vehNumberBranching[0] += 0.5 * bcVar;
                    vehNumberBranching[depotId] += 0.5 * bcVar;
                    if (data.depots[depotId].veh_cost > 0.0)
                        objective += (0.5 * data.depots[depotId].veh_cost
                                        + data.getDepotToCustDistance(depotId, secondCustId)) * bcVar;
                    else
                        objective += data.getDepotToCustDistance(depotId, secondCustId) * bcVar;
                }
                else
                {
                    degreeConstr[firstCustId] += bcVar;
                    objective += data.getCustToCustDistance(firstCustId, secondCustId) * bcVar;
                }
                degreeConstr[secondCustId] += bcVar;
                edgeBranching[firstCustId][secondCustId] += bcVar;
            }

        RCSPSolver solver(spForm, depotId);
#ifdef BCP_RCSP_IS_FOUND
        spForm.attach(solver.getOracle());
#endif
    }

    std::vector<int> demands(data.nbCustomers + 1, 0);
    for (auto custId : customerIds)
        demands[custId] = data.customers[custId].demand;
    int max_capacity = 0;
    for (auto depotId : depotIds)
        max_capacity = (std::max)(max_capacity, data.depots[depotId].veh_capacity);

    BcCapacityCutConstrArray capacityCuts(master, max_capacity, demands,  true, true,
                                          -1, 3.0, 1.0);

    BcLimMemRankOneCutConstrArray limMemRank1Cuts(master);
}
