/**
 *
 * This file RCSPSolver.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "RCSPSolver.h"
#include "Data.h"
#include "Parameters.h"
#include "SolutionChecker.h"

vrpstw::RCSPSolver::RCSPSolver(BcFormulation spForm, int depotId) : spForm(std::move(spForm)), depotId(depotId)
#ifdef BCP_RCSP_IS_FOUND
    , oracle(nullptr)
#endif
{
	BcNetwork network(spForm, data.nbCustomers + 1, data.nbCustomers + 1);

    BcNetworkResource time_res(network, 0);
    time_res.setAsMainResource();

    BcNetworkResource cap_res(network, 1);
	cap_res.setAsMainResource();

    buildVertices(network, time_res, cap_res, depotId);
	buildArcs(network, time_res, cap_res, depotId);
    buildElemSetDistanceMatrix(network);

#ifdef BCP_RCSP_IS_FOUND
    oracle = new BcRCSPFunctor(spForm);
#else
    std::cerr << "Cannot use BcRCSPFunctor (VRPSolver extension) as the BCP_RCSP library is not found." << std::endl;
    std::cerr << "Please rerun CMake if you have already put this library in BapcodFramework/Tools/rcsp/." << std::endl;
    std::cerr << "Othewise, please contact ruslan.sadykov@inria.fr to obtain this library in compiled form." << std::endl;
    std::cerr << "The BCP_RCSP library is available only for Mac OS and Linux." << std::endl;
    exit(1);
#endif
}

void vrpstw::RCSPSolver::buildVertices(BcNetwork & network, BcNetworkResource & time_res, BcNetworkResource & cap_res,
                                       int depot_id)
{
	for (int custId = 0; custId <= data.nbCustomers + 1; ++custId)
	{
		BcVertex vertex = network.createVertex();

    	if ((custId >= 1) && (custId <= data.nbCustomers))
    	{
            time_res.setVertexConsumptionLB(vertex, data.customers[custId].tw_start);
            time_res.setVertexConsumptionUB(vertex, data.customers[custId].tw_end);
        } else
        {
            time_res.setVertexConsumptionLB(vertex, data.depots[depot_id].tw_start);
            time_res.setVertexConsumptionUB(vertex, data.depots[depot_id].tw_end);
        }

        cap_res.setVertexConsumptionLB(vertex, 0);
		cap_res.setVertexConsumptionUB(vertex, data.depots[depot_id].veh_capacity);

		if (custId == 0)
        {
            network.setPathSource(vertex);
        }
		else if (custId <= data.nbCustomers)
        {
		    vertex.setElementaritySet(custId);
		    vertex.setPackingSet(custId);
        }
		else
        {
		    network.setPathSink(vertex);
        }
	}
}

void vrpstw::RCSPSolver::buildArcs(BcNetwork & network, BcNetworkResource & time_res, BcNetworkResource & cap_res,
                                   int depot_id)
{
	BcVarArray xVar(spForm, "X");

    for (int firstCustId = 0; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = 1; secondCustId <= data.nbCustomers + 1; ++secondCustId)
        {
            if (firstCustId == secondCustId)
                continue;

            int minCustId = (std::min)(firstCustId, (secondCustId <= data.nbCustomers) ? secondCustId : 0);
            int maxCustId = (std::max)(firstCustId, (secondCustId <= data.nbCustomers) ? secondCustId : 0);

            if (minCustId == maxCustId)
                continue;

            BcArc arc = network.createArc(firstCustId, secondCustId, 0.0);
            arc.arcVar((BcVar)xVar[depotId][minCustId][maxCustId]);

            cap_res.setArcConsumption(arc, secondCustId <= data.nbCustomers ? data.customers[secondCustId].demand : 0);
            if (minCustId == 0)
                time_res.setArcConsumption(arc, data.getDepotToCustDistance(depot_id, maxCustId)
                                                + ((firstCustId > 0) ? data.customers[firstCustId].service_time : 0.0));
            else
                time_res.setArcConsumption(arc, data.getCustToCustDistance(firstCustId, secondCustId)
                                                + data.customers[firstCustId].service_time);
        }
}

void vrpstw::RCSPSolver::buildElemSetDistanceMatrix(BcNetwork & network)
{
    int nbElemSets = data.nbCustomers + 1;
    std::vector<std::vector<double> > distanceMatrix(nbElemSets, std::vector<double>(nbElemSets, 1e12));

    for (int firstCustId = 1; firstCustId <= data.nbCustomers; ++firstCustId )
        for (int secondCustId = 1; secondCustId <= data.nbCustomers; ++secondCustId)
            distanceMatrix[firstCustId][secondCustId] = data.getCustToCustDistance(firstCustId, secondCustId);

    network.setElemSetsDistanceMatrix(distanceMatrix);
}


