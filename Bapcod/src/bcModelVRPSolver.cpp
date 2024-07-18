/**
 *
 * This file bcModelVRPSolver.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved.
 * See License.pdf file for license terms.
 *
 * This file is available only for academic use.
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#if(JSON_IS_FOUND && BCP_RCSP_IS_FOUND)

#include <bcModelVRPSolver.hpp>
#include "bcModelingLanguageC.hpp"
#include <cmath>
#include <cfloat>

namespace API_VRP
{
    /**
     * @brief print values of map
    **/
    template<typename S>
    void print_values(S &map) noexcept
    {
        for (auto it: map)
        {
            std::cout << it.first << " : ";
            for (auto &values: it.second)
            {
                std::cout << std::get<0>(values) << " " << std::get<1>(values) << " ";
            }
            std::cout << "\n";
        }
    }

    /**
      * @brief set map with a priority queue
    **/
    template<typename S, typename D, typename Q>
    void set_priorityQueue(S &map, D key, Q value) noexcept
    {
        auto it = map.find(key);
        if (it != map.end())
            it->second.push(Q(value));
        else
        {
            std::priority_queue<Q, std::vector<Q>, Comparator<Q> > priorityQueue;
            priorityQueue.push(Q(value));
            map[key] = priorityQueue;
        }
    }


    /**
   * @brief set map with a vector
   **/
    template<typename S, typename D, typename Q>
    void set_vector(S &map, D key, Q value) noexcept
    {
        auto it = map.find(key);
        if (it != map.end())
            it->second.push_back(value);
        else
        {
            std::vector<Q> vect;
            vect.push_back(value);
            map[key] = vect;
        }
    }

    size_t get_nbExtensions(size_t startPointIndex,
                            std::map<size_t, std::vector<std::tuple<double, size_t, size_t>>> &successorsPoint) noexcept
    {
        auto it = successorsPoint.find(startPointIndex);
        if (it != successorsPoint.end())
        {
            return it->second.size();
        }
        return 0;
    }

}

API_VRP::Point::Point(size_t iID, size_t idCustomer, double iPenalty, double iServiceTime, double iTimeWindowsBegin,
                      double iTimeWindowsEnd, size_t demand, bool iArtificial)
{
	id = iID;
	idCust = idCustomer;
	costORPenalty = iPenalty;
	serviceTime = iServiceTime;
    timeWindowBegin = iTimeWindowsBegin;
    timeWindowEnd = iTimeWindowsEnd;
	demandOrCapacity = demand;
	isDepot = (idCustomer > 0) || iArtificial; /// artificial point is always a depot
    isArtificial = iArtificial;
}

API_VRP::Parameters::Parameters(double iTimeLimit, double iUB, bool iHeuristicUsed)
{
	timeLimit = iTimeLimit;
	initialUB = iUB;
	heuristicUsed = iHeuristicUsed;
}

API_VRP::Parameters::Parameters(rapidjson::GenericObject<false, rapidjson::Value>& iParam)
{
	auto it = iParam.FindMember(timeLimitStr);
	if (it != iParam.MemberEnd())
		timeLimit = it->value.GetDouble();
	
	it= iParam.FindMember(initialUBStr);
	if (it != iParam.MemberEnd())
		initialUB = it->value.GetDouble();
	
	it = iParam.FindMember(heuristicUsedStr);
	if (it != iParam.MemberEnd())
		heuristicUsed = it->value.GetBool();
	
	it = iParam.FindMember(timeLimitHeuristicStr);
	if (it != iParam.MemberEnd())
		timeLimitHeuristic = it->value.GetDouble();

	it = iParam.FindMember(configFileStr);
	if (it != iParam.MemberEnd())
		configFile = it->value.GetString();

	it = iParam.FindMember(solverNameStr);
	if (it != iParam.MemberEnd())
		solverName = it->value.GetString();

	it = iParam.FindMember(actionStr);
	std::string action{ "solve" };
	if (it != iParam.MemberEnd())
	{
		action = it->value.GetString();
	}
	enumerate = (action == "enumAllFeasibleRoutes");

	it = iParam.FindMember(printlevelStr);
	if (it != iParam.MemberEnd())
		printlevel = it->value.GetInt();
}


API_VRP::Point::Point(rapidjson::GenericArray<false, rapidjson::Value>& iPoint, size_t index)
{
	auto it = iPoint[index].FindMember(nameStr);
	if (it != iPoint[index].MemberEnd())
		name = it->value.GetString();

	id = iPoint[index][idStr].GetInt();

	it = iPoint[index].FindMember(idCustomerStr);
	if (it != iPoint[index].MemberEnd())
		idCust = it->value.GetInt();

	it = iPoint[index].FindMember(penaltyOrCostStr);
	if (it != iPoint[index].MemberEnd())
		costORPenalty = it->value.GetDouble();

	it = iPoint[index].FindMember(serviceTimeStr);
	if (it != iPoint[index].MemberEnd())
		serviceTime = it->value.GetDouble();
	
	it = iPoint[index].FindMember(twBeginStr);
	if (it != iPoint[index].MemberEnd())
        timeWindowBegin = it->value.GetDouble();

	it = iPoint[index].FindMember(twEndStr);
	if (it != iPoint[index].MemberEnd())
        timeWindowEnd = it->value.GetDouble();

	it = iPoint[index].FindMember(demandOrCapacityStr);
	if (it != iPoint[index].MemberEnd())
		demandOrCapacity = it->value.GetInt();
	
	if (idCust == 0)
        isDepot = true;

    it = iPoint[index].FindMember(incompatibleVehiclesStr);
    if (it != iPoint[index].MemberEnd())
    {
        auto array = iPoint[index][incompatibleVehiclesStr].GetArray();
        set_incompatibility(array);
    }

}

API_VRP::VehicleType::VehicleType(rapidjson::GenericArray<false, rapidjson::Value> & iVehicle, size_t index)
{
	auto it = iVehicle[index].FindMember(startPointIdStr);
	if (it != iVehicle[index].MemberEnd())
        startPointId = it->value.GetInt();

    it = iVehicle[index].FindMember(endPointIdStr);
    if (it != iVehicle[index].MemberEnd())
        endPointId = it->value.GetInt();

	it = iVehicle[index].FindMember(nameStr);
	if (it != iVehicle[index].MemberEnd())
		name = it->value.GetString();

	id = iVehicle[index][idStr].GetInt();
	index = index;

	it = iVehicle[index].FindMember(capacityStr);
	if (it != iVehicle[index].MemberEnd())
		capacity = it->value.GetInt();

	it = iVehicle[index].FindMember(varCostDistStr);
	if (it != iVehicle[index].MemberEnd())
		varCostDist = it->value.GetDouble();

	it = iVehicle[index].FindMember(varCostTimeStr);
	if (it != iVehicle[index].MemberEnd())
		varCostTime = it->value.GetDouble();

	it = iVehicle[index].FindMember(maxNumberStr);
	if (it != iVehicle[index].MemberEnd())
        maxNbVehicles = it->value.GetInt();
	
	it = iVehicle[index].FindMember(twBeginStr);
	if (it != iVehicle[index].MemberEnd())
        timeWindowBegin = it->value.GetDouble();

	it = iVehicle[index].FindMember(twEndStr);
	if (it != iVehicle[index].MemberEnd())
        timeWindowEnd = it->value.GetDouble();

	it = iVehicle[index].FindMember(fixCostStr);
	if (it != iVehicle[index].MemberEnd())
		fixedCost = it->value.GetDouble();
}

API_VRP::Link::Link(rapidjson::GenericArray<false, rapidjson::Value>& iLink, size_t index)
{
	//Assert that start and endPointId are known
	auto it = iLink[index].FindMember(startPointIdStr);
    if (it != iLink[index].MemberEnd())
        startPointId = it->value.GetInt();

    it = iLink[index].FindMember(endPointIdStr);
    if (it != iLink[index].MemberEnd())
        endPointId = it->value.GetInt();

	it = iLink[index].FindMember(nameStr);
	if (it != iLink[index].MemberEnd())
		name = it->value.GetString();

	id = index;

	it = iLink[index].FindMember(isDirectedStr);
	if (it != iLink[index].MemberEnd())
		isDirected = it->value.GetBool();

	it = iLink[index].FindMember(distanceStr);
	if (it != iLink[index].MemberEnd())
		distance = it->value.GetDouble();

	it = iLink[index].FindMember(timeStr);
	if (it != iLink[index].MemberEnd())
		time = it->value.GetDouble();
	
	it = iLink[index].FindMember(fixCostStr);
	if (it != iLink[index].MemberEnd())
		fixedCost = it->value.GetDouble();
}

API_VRP::Link::Link(size_t iId, bool iDirected, size_t iStartPointId, size_t iEndPointId, double iDistance,
                    double iTime, double iFixedCost, const std::string & iName, bool iArtificial)
{
    id = iId;
    isDirected = iDirected;
    startPointId = iStartPointId;
    endPointId = iEndPointId;
    name = iName;
    distance = iDistance;
    time = iTime;
    fixedCost = iFixedCost;
    isArtificial = iArtificial;
}


API_VRP::ExceptionVRPPython::ExceptionVRPPython(ExceptionType error_type, const std::string& iMessage) : std::exception(), error_type(error_type)
{
	std::ostringstream os;
	os << "{ \"Error Type\": " << std::to_string(static_cast<int>(error_type));
	if (!iMessage.empty())
		os << R"(, "Description": ")" << iMessage << '"';
	os << "}";

	message = new char[os.str().size() + 1];
	strcpy(message, os.str().data());
	message[os.str().size()] = '\0';

	errorStr = new char[iMessage.size()];
	strcpy(errorStr, iMessage.data());
}

void API_VRP::Data::setParametersFromBapcod(BcInitialisation & bcInit)
{
    imposeCapacityResourceByCuts = bcInit.param().VRPSEimposeCapacityResourceByCuts();
    if (bcInit.param().VRPSEcriticalResource() >= 0)
    {
        for (auto & vehType : vehicleTypes)
        {
            if (bcInit.param().VRPSEcriticalResource() == 0 && vehType.get_capacityResourceId() == 0
                && vehType.get_timeResourceId() == 1)
            {
                vehType.set_timeResourceId(0);
                vehType.set_capacityResourceId(1);
            }
            if (bcInit.param().VRPSEcriticalResource() == 1 && vehType.get_capacityResourceId() == 1
                && vehType.get_timeResourceId() == 0)
            {
                vehType.set_timeResourceId(1);
                vehType.set_capacityResourceId(0);
            }
        }
    }
}

std::vector<double> API_VRP::Data
                    ::Dijkstra(const VehicleType & vehicleType,
                               const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId) const
{
    std::vector<double> shortestTimeToEndPoint(maxPointId + 1, INT_MAX);
    typedef std::pair<double, size_t> Pair;
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> queue;
    shortestTimeToEndPoint[vehicleType.get_endPointId()] = 0;
    auto * endPointPtr = pointPtrFromId[vehicleType.get_endPointId()];
    queue.push(std::make_pair(0.0 + endPointPtr->get_serviceTime(), vehicleType.get_endPointId()));
    while (!queue.empty())
    {
        auto currPointId = queue.top().second;
        queue.pop();
        /// get all predecessors of i
        for (auto & predPair : predecessorsFromPointId[currPointId])
        {
            auto & link = links[predPair.first];
            if (!link.get_isCompatibleWithVehType(vehicleType.get_id()))
                continue;

            auto prevPointId = predPair.second ? link.get_startPointId() : link.get_endPointId();
            auto * prevPointPtr = pointPtrFromId[prevPointId];
            auto time = link.get_time() + prevPointPtr->get_serviceTime();
            if (shortestTimeToEndPoint[prevPointId] > shortestTimeToEndPoint[currPointId] + time)
            {
                shortestTimeToEndPoint[prevPointId] = shortestTimeToEndPoint[currPointId] + time;
                queue.push(std::make_pair(shortestTimeToEndPoint[prevPointId], prevPointId));
            }
        }
    }
    return shortestTimeToEndPoint;
}

bool API_VRP::Data
     ::determineSolverParameterisation(bool timeResourceIsUsed, bool capacityResourceIsUsed,
                                       bool oneOptionalCustomerExists,
                                       const std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                                       const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId)
{
    if ((!capacityResourceIsUsed) && (!timeResourceIsUsed))
    {
        const char* message = "At least one resource must be used in model";
        std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
        setErrorCodeAndMsg(ExceptionType::INCONSISTENCY, message);
        return false;
    }

    if (capacityResourceIsUsed && !timeResourceIsUsed)
    {
        for (auto & vehicleType : vehicleTypes)
        {
            vehicleType.set_capacityResourceId(0);
            vehicleType.set_timeResourceId(-1);
        }
        size_t totalCustomerDemand = 0;
        for (auto & point : points)
            if (point.get_isCustomer())
                totalCustomerDemand += point.get_demandOrCapacity();
        double averDemandPerCustomer = (double) totalCustomerDemand / (double)nbCustomers;
        size_t maxCapacity = 0;
        for (auto & vehicleType : vehicleTypes)
            if (maxCapacity < vehicleType.get_capacity())
                maxCapacity = vehicleType.get_capacity();
        useArcMemoryForRankOneCuts = (double)maxCapacity / averDemandPerCustomer > rankOneCutsArcMemoryThreshold;
        if (parameters.get_printlevel() == 0)
            std::cout << ((useArcMemoryForRankOneCuts) ? "Arc " : "Node ") << "memory for rank-1 cuts will be used ("
                      << (double) maxCapacity / averDemandPerCustomer << ")" << std::endl;
        return true;
    }

    /// TO DO : we need to consider the multi-trip case individually (time would always be the critical resource)
    ///         and the time should always be used for multi-trip instances

    /// we do not initialize generator on purpose to keep the code deterministic
    std::default_random_engine generator;
    std::geometric_distribution<int> distribution(0.5);
    size_t nbGeneratedRandomPaths = 0;
    double sumPathDemandRatios = 0.0;
    bool timeResourceIsAlwaysCritical = true;
    int capacitySample = -1;
    bool allVehicleCapacitiesAreTheSame = true;
    std::vector<double> averPathLengths;
    for (auto & vehicleType : vehicleTypes)
    {
        size_t sumPathLengths = 0;
        if (vehicleType.get_maxNbVehicles() == 0)
            continue;

        if (capacityResourceIsUsed)
        {
            if (capacitySample == -1)
                capacitySample = (int) vehicleType.get_capacity();
            else if (capacitySample != vehicleType.get_capacity())
                allVehicleCapacitiesAreTheSame = false;
        }

        auto timeToEndPointFromPointId = Dijkstra(vehicleType, predecessorsFromPointId);
        auto * startPointPtr = pointPtrFromId[vehicleType.get_startPointId()];
        auto * endPointPtr = pointPtrFromId[vehicleType.get_endPointId()];
        double startTime = (std::max)(vehicleType.get_timeWindowBegin(), startPointPtr->get_timeWindowBegin());
        double endTime = (std::min)(vehicleType.get_timeWindowEnd(), endPointPtr->get_timeWindowEnd());

        int nbTimeCriticalPaths = 0;
        int nbCapacityCriticalPaths = 0;
        int nbGeneratedPathsForThisVehicleType = 0;
        for (int randomPathIndex = 0; randomPathIndex < nbRandomPathsInAutoParameterisation; ++randomPathIndex)
        {
            /// we follow a randomized greedy path from the start node in the following manner:
            /// - we collect feasible extensions from the current point to non-end points
            ///   (extension is feasible if vehicle capacity is respected and it is enough time to go to the end point)
            /// - if the majority of extensions cannot be done because of the vehicle capacity, then
            ///   the current path is capacity critical
            /// - if no time-feasible extension, then the current path is time-critical
            /// - we choose randomly the extension, the higher is the current time in the extension,
            ///   the lower is probability to choose it (use use geometric distribution with 0.5 parameter for that)

            size_t currDemand = 0;
            double currTime = startTime;
            size_t currPointId = startPointPtr->get_id();
            size_t pathLength = 0;
            std::vector<bool> visitedCustomers(maxCustId + 1, false);
            std::vector<size_t> visitedPointsSequence(1, currPointId);
            do
            {
                int nbInfeasibleExtensionsDueToCapacity = 0;
                std::vector<std::tuple<double, double, size_t>> feasibleExtensions;
                auto * currPointPtr = pointPtrFromId[currPointId];
                for (auto & succPair : successorsFromPointId[currPointId])
                {
                    auto & link = links[succPair.first];
                    if (!link.get_isCompatibleWithVehType(vehicleType.get_id()))
                        continue;

                    size_t nextPointId = (succPair.second) ? link.get_startPointId() : link.get_endPointId();
                    auto * nextPointPtr = pointPtrFromId[nextPointId];

                    if (nextPointId == endPointPtr->get_id())
                        continue;

                    if (nextPointPtr->get_isCustomer() && visitedCustomers[nextPointPtr->get_idCust()])
                        continue;

                    double extTime = (std::max)(currTime + currPointPtr->get_serviceTime() + link.get_time(),
                                                nextPointPtr->get_timeWindowBegin());

                    if (extTime + timeToEndPointFromPointId[nextPointId] > endTime)
                        continue; /// extension infeasible due to the time
                    if (capacityResourceIsUsed &&
                        currDemand + nextPointPtr->get_demandOrCapacity() > vehicleType.get_capacity())
                    {
                        nbInfeasibleExtensionsDueToCapacity += 1;
                        continue;
                    }
                    double extCost = link.get_fixedCost() + link.get_time() * vehicleType.get_varCostTime()
                                     + link.get_distance() * vehicleType.get_varCostDist();
                    feasibleExtensions.emplace_back(extCost, extTime, nextPointId);
                }
                if (nbInfeasibleExtensionsDueToCapacity > successorsFromPointId[currPointId].size() / 2)
                {
                    nbCapacityCriticalPaths += 1;
                    break;
                }
                if (feasibleExtensions.empty())
                {
                    nbTimeCriticalPaths += 1;
                    break;
                }
                std::stable_sort(feasibleExtensions.begin(), feasibleExtensions.end());
                int chosenIndex = distribution(generator);
                if (chosenIndex >= feasibleExtensions.size())
                    chosenIndex = (int)feasibleExtensions.size() - 1;

                currPointId = std::get<2>(feasibleExtensions[chosenIndex]);
                visitedCustomers[pointPtrFromId[currPointId]->get_idCust()] = true;
                visitedPointsSequence.push_back(currPointId);
                if (capacityResourceIsUsed)
                    currDemand += pointPtrFromId[currPointId]->get_demandOrCapacity();
                currTime = std::get<1>(feasibleExtensions[chosenIndex]);
                pathLength += 1;
            }
            while (true);

            if (currTime + timeToEndPointFromPointId[currPointId] <= endTime)
            {
                /// path is feasible
                nbGeneratedRandomPaths += 1;
                nbGeneratedPathsForThisVehicleType += 1;
            }
            if (capacityResourceIsUsed)
                sumPathDemandRatios += (double) currDemand / (double) vehicleType.get_capacity();
            sumPathLengths += pathLength;
        }

        if (nbGeneratedPathsForThisVehicleType == 0)
        {
            vehicleType.set_maxNbVehicles(0);
            continue;
        }
        if (nbTimeCriticalPaths > nbCapacityCriticalPaths)
        {
            vehicleType.set_timeResourceId(0);
            vehicleType.set_capacityResourceId(1);
            if (parameters.get_printlevel() == 0)
                std::cout << "Time resource is critical for vehicle type " << vehicleType.get_id() << " ("
                          << nbTimeCriticalPaths << "," << nbCapacityCriticalPaths << ")" << std::endl;
        }
        else
        {
            vehicleType.set_capacityResourceId(0);
            vehicleType.set_timeResourceId(1);
            timeResourceIsAlwaysCritical = false;
            if (parameters.get_printlevel() == 0)
                std::cout << "Capacity resource is critical for vehicle type " << vehicleType.get_id() << " ("
                          << nbTimeCriticalPaths << "," << nbCapacityCriticalPaths << ")" << std::endl;
        }
        averPathLengths.push_back( (double)sumPathLengths / (double)nbGeneratedPathsForThisVehicleType);
    }

    if (nbGeneratedRandomPaths == 0)
    {
        /// this my happen because of time windows of vehicles types
        const char* message = "No feasible path exists. Check time windows of vehicle types.";
        std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
        setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
        return false;
    }

    if (capacityResourceIsUsed && allVehicleCapacitiesAreTheSame && timeResourceIsAlwaysCritical
        && !oneOptionalCustomerExists)
    {
        /// we now determine whether the capacity resource should be imposed by rounded capacity cuts and not explicitly
        double averagePathDemandRatio = sumPathDemandRatios / (double)nbGeneratedRandomPaths;
        if (averagePathDemandRatio < pathDemandRatioThresholdToImposeCapacityByCuts)
        {
            imposeCapacityResourceByCuts = true;
            if (parameters.get_printlevel() >= 0)
                std::cout << "Capacity resource will be imposed by rounded capacity cuts ("
                          << averagePathDemandRatio << ")" << std::endl;
            for (auto & vehicleType : vehicleTypes)
                vehicleType.set_capacityResourceId(-1);
        }
        else if (parameters.get_printlevel() >= 0)
        {
            std::cout << "Capacity resource will be imposed in pricing (" << averagePathDemandRatio << ")" << std::endl;
        }
    }

    std::stable_sort(averPathLengths.rbegin(), averPathLengths.rend());
    useArcMemoryForRankOneCuts = averPathLengths.front() > rankOneCutsArcMemoryThreshold;
    if (parameters.get_printlevel() == 0)
        std::cout << ((useArcMemoryForRankOneCuts) ? "Arc " : "Node ") << "memory for rank-1 cuts will be used ("
                  << averPathLengths.front() << ")" << std::endl;
    return true;
}

bool API_VRP::Data::determineIfSymmetricCase(bool timeResourceIsUsed)
{
    if (timeResourceIsUsed)
    {
        auto timeWindowBeginSample = points.front().get_timeWindowBegin();
        auto timeWindowEndSample = points.front().get_timeWindowEnd();

        for (auto &point: points)
            if (point.get_timeWindowBegin() != timeWindowBeginSample
                || point.get_timeWindowEnd() != timeWindowEndSample)
                return false;
        for (auto & vehicleType : vehicleTypes)
            if (vehicleType.get_timeWindowBegin() != timeWindowBeginSample
                || vehicleType.get_timeWindowEnd() != timeWindowEndSample)
                return false;
    }

    for (auto & link : links)
        if (link.get_isDirected())
            return false;

    return true;
}

bool isNotInteger(double x)
{
    return trunc(x) != x;
}

bool API_VRP::Data::determineIfIntegerObjective()
{
    for (auto & vehicleType : vehicleTypes)
        if (isNotInteger(vehicleType.get_fixedCost())
            || isNotInteger(vehicleType.get_varCostTime())
            || isNotInteger(vehicleType.get_varCostDist()))
            return false;

    for (auto & point : points)
        if (isNotInteger(point.get_costOrPenalty()))
            return false;

    for (auto & link : links)
        if (isNotInteger(link.get_fixedCost()) || isNotInteger(link.get_time()) || isNotInteger(link.get_distance()))
            return false;

    return true;
}

void API_VRP::Data::getDistanceMatrixBetweenCustomers(std::vector<std::vector<double>> & matrix,
                                                      const VehicleType & vehType) const
{
    std::vector<std::vector<std::pair<int, double>>>
            linksBetweenData(maxCustId + 1, std::vector<std::pair<int, double>>(maxCustId + 1, std::make_pair(0, 0.0)));

    for (auto & link : links)
    {
        if (!link.get_isFeasible() || !link.get_isCompatibleWithVehType(vehType.get_id()))
            continue;

        double linkCost = link.get_fixedCost();
        linkCost += vehType.get_varCostDist() * link.get_distance();
        linkCost += vehType.get_varCostTime() * link.get_time();

        auto tailCustId = getTailPoint(link).get_idCust();
        auto headCustId = getHeadPoint(link).get_idCust();
        if (tailCustId > 0 && headCustId > 0)
        {
            linksBetweenData[tailCustId][headCustId].first += 1;
            linksBetweenData[tailCustId][headCustId].second += linkCost;
            if (!link.get_isDirected())
            {
                linksBetweenData[headCustId][tailCustId].first += 1;
                linksBetweenData[headCustId][tailCustId].second += linkCost;
            }
        }
    }

    matrix = std::vector<std::vector<double>>(maxCustId + 1, std::vector<double>(maxCustId + 1, 1e12));
    for (auto firstCustId : customerIds)
        for (auto secondCustId : customerIds)
        {
            auto numLinksBetween = linksBetweenData[firstCustId][secondCustId].first;
            auto totalDistanceBetween = linksBetweenData[firstCustId][secondCustId].second;
            if (numLinksBetween > 0)
                matrix[firstCustId][secondCustId] = totalDistanceBetween / numLinksBetween;
        }
}

void API_VRP::Data::preprocessLinks(std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                                    std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId)
{
	/// Initialisation of links
	for (size_t linkId = 0; linkId < links.size(); linkId++)
    {
        auto & link = links[linkId];
        auto startPointId = link.get_startPointId();
        auto startPointPtr = pointPtrFromId[startPointId];
        auto endPointId = link.get_endPointId();
		auto endPointPtr = pointPtrFromId[endPointId];

        if (startPointPtr->get_timeWindowBegin() + startPointPtr->get_serviceTime() + link.get_time()
            > endPointPtr->get_timeWindowEnd() - endPointPtr->get_serviceTime())
        {
            link.set_dirIsFeasible(false);
        }
        else
        {
            successorsFromPointId[startPointId].push_back(std::make_pair(linkId, false));
            predecessorsFromPointId[endPointId].push_back(std::make_pair(linkId, true));
        }

        if (!link.get_isDirected())
        {
            if (endPointPtr->get_timeWindowBegin() + endPointPtr->get_serviceTime() + link.get_time()
                > startPointPtr->get_timeWindowEnd() - startPointPtr->get_serviceTime())
            {
                link.set_oppIsFeasible(false);
            }
            else
            {
                successorsFromPointId[endPointId].emplace_back(linkId, true);
                predecessorsFromPointId[startPointId].emplace_back(linkId, false);
            }
        }
	}
}

bool API_VRP::Data
     ::preprocessVehicleTypes(bool timeResourceIsUsed,
                              const std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                              const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId)
{
    nbVehicleTypesAvailable = 0;
    totalNbVehicles = 0;
    totalVehicleCapacity = 0;
    bool multiTripInstance = false;
    for (auto & vehicleType : vehicleTypes)
    {
        if (vehicleType.get_maxNbVehicles() == 0)
            continue;

        if (timeResourceIsUsed
            && (vehicleType.get_timeWindowEnd() - vehicleType.get_timeWindowBegin() <= API_VRP::epsilonParam))
        {
            vehicleType.set_minNbVehicles(0);
            if (parameters.get_printlevel() >= -1)
            {
                std::cout << "VRPSolverEasy WARNING : vehicle type with id " << vehicleType.get_id()
                          << " cannot be used as its time window is empty. " << std::endl;
            }
            continue;
        }

        std::list<size_t> queueOfPointIds;

        /// first in the forward sense
        std::vector<bool> forwVisitedPointIds(maxPointId + 1, false);
        queueOfPointIds.push_back(vehicleType.get_startPointId());
        while (!queueOfPointIds.empty())
        {
            auto pointId = queueOfPointIds.front();
            queueOfPointIds.pop_front();
            for (auto & succPair : successorsFromPointId[pointId])
            {
                auto & link = links[succPair.first];
                auto nextPointId = (succPair.second) ? link.get_startPointId() : link.get_endPointId();
                if (link.get_isCompatibleWithVehType(vehicleType.get_id()) && !forwVisitedPointIds[nextPointId])
                {
                    forwVisitedPointIds[nextPointId] = true;
                    if (nextPointId != vehicleType.get_endPointId())
                        queueOfPointIds.push_back(nextPointId);
                }
            }
        }

        /// then, in the backward sense
        std::vector<bool> backVisitedPointIds(maxPointId + 1, false);
        queueOfPointIds.clear();
        queueOfPointIds.push_back(vehicleType.get_endPointId());
        while (!queueOfPointIds.empty())
        {
            auto pointId = queueOfPointIds.front();
            queueOfPointIds.pop_front();
            for (auto & predPair : predecessorsFromPointId[pointId])
            {
                auto & link = links[predPair.first];
                auto prevPointId = (predPair.second) ? link.get_startPointId() : link.get_endPointId();
                if (link.get_isCompatibleWithVehType(vehicleType.get_id()) && !backVisitedPointIds[prevPointId])
                {
                    backVisitedPointIds[prevPointId] = true;
                    if (prevPointId != vehicleType.get_startPointId())
                        queueOfPointIds.push_back(prevPointId);
                }
            }
        }

        for (auto & point : points)
        {
            size_t id = point.get_id();
            if (id == vehicleType.get_startPointId() || id == vehicleType.get_endPointId())
                continue;

            if (!forwVisitedPointIds[id] && !backVisitedPointIds[id])
                point.set_incompatibility(vehicleType.get_id());
            else if (point.get_isDepot())
                multiTripInstance = true;
        }

		if (!forwVisitedPointIds[vehicleType.get_endPointId()])
        {
            /// end node cannot be reached -> this vehicle type cannot be used
            vehicleType.set_maxNbVehicles(0);
            if (parameters.get_printlevel() >= -1)
            {
                std::cout << "VRPSolverEasy WARNING : vehicle type with id " << vehicleType.get_id()
                          << " cannot be used as there is no feasible path between its start and end points. "
                          << std::endl;
            }
        }

        if (forwVisitedPointIds[vehicleType.get_endPointId()])
        {
            nbVehicleTypesAvailable += 1;
            totalNbVehicles += vehicleType.get_maxNbVehicles();
            totalVehicleCapacity += vehicleType.get_maxNbVehicles() * vehicleType.get_capacity();
        }
	}

    if (multiTripInstance)
        return true;

    int totalDemandOfCompulsoryCustomers = 0;
    for (auto & point : points)
    {
        if (point.get_isCustomer() && !point.get_isOptional())
            totalDemandOfCompulsoryCustomers += (int)point.get_demandOrCapacity();
    }

    if (totalDemandOfCompulsoryCustomers > totalVehicleCapacity)
    {
        const char* message = "After preprocessing, the total capacity of vehicles is not enough to cover "
            "the total demand of non-optional customers.";
        std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
        setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
        return false;
    }

    for (auto & vehicleType : vehicleTypes)
    {
        int obligatoryDemand = totalDemandOfCompulsoryCustomers;
        if (vehicleTypes.size() > 1)
            obligatoryDemand -= (int)(totalVehicleCapacity - vehicleType.get_capacity());
        if (obligatoryDemand > 0)
            vehicleType.set_minNbVehicles((obligatoryDemand - 1)/vehicleType.get_capacity() + 1);
    }
    return true;
}

inline bool API_VRP::Data::checkAndSetIncompatibleVehicles()
{
	/// We check for each point the set of incompatible vehicle types
	for (auto & point : points)
        for (auto vehicleTypeId : point.get_incompatibleVehicleTypeIds())
            if (vehicleTypeId > maxVehicleTypeId || vehicleTypePtrFromId[vehicleTypeId] == nullptr)
            {
                std::string message = "Incompatible vehicle type " + std::to_string(vehicleTypeId) + " does not exist";
                std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
                setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
                return false;
            }


    for (auto & link : links)
    {
        auto & startPointIncompats = pointPtrFromId[link.get_startPointId()]->get_incompatibleVehicleTypeIds();
        auto & endPointIncompats = pointPtrFromId[link.get_endPointId()]->get_incompatibleVehicleTypeIds();
        link.set_incompatibilities(startPointIncompats, endPointIncompats);
    }
    return true;
}

API_VRP::Data::Data(rapidjson::Document *json)
{
	auto pointsArray = (*json)[POINTS].GetArray();
	auto vehicleTypesArray = (*json)[VEHICLE_TYPES].GetArray();
	auto linksArray = (*json)[LINKS].GetArray();
	auto param = (*json)[PARAMETERS].GetObj();
    auto it = (*json).FindMember(MaxTotalVehiclesNumberStr);
    if (it != (*json).MemberEnd())
        maxTotalVehiclesNumber = it->value.GetInt();

	size_t nbDepots = 0;
    nbVehicleTypes = vehicleTypesArray.Size();
	nbLinks = linksArray.Size();
	nbPoints = pointsArray.Size();

    /// Initialisation of customers and depots
	maxCustId = 0;
    maxPointId = 0;
    points.reserve(nbPoints + 2); /// possibly two more artificial depots
	for (size_t pointIndex = 0; pointIndex < nbPoints; pointIndex++)
    {
        points.emplace_back(pointsArray, pointIndex);
        auto & point = points.back();
        checkPoint(point);
        maxCustId = (std::max)(maxCustId, point.get_idCust());
        maxPointId = (std::max)(maxPointId, point.get_id());
    }

    pointPtrFromId.resize(maxPointId + 1, nullptr);
    pointPtsFromCustomerId.resize(maxCustId + 1);
    nbCustomers = 0;

    bool timeResourceIsUsed = false;
    bool capacityResourceIsUsed = false;
    bool oneOptionalCustomerExists = false;
    double minTime = DBL_MAX;
    double maxTime = DBL_MIN;
    for (auto & point : points)
    {
        auto id = point.get_id();
        if (pointPtrFromId[id] != nullptr)
        {
            std::string message = "Each point mush have an unique id.";
            std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
            setErrorCodeAndMsg(ExceptionType::POINTS_ERROR, message);
            return;
        }
        pointPtrFromId[id] = &point;
        if (point.get_timeWindowBegin() < point.get_timeWindowEnd())
            timeResourceIsUsed = true;
        if (minTime > point.get_timeWindowBegin())
            minTime = point.get_timeWindowBegin();
        if (maxTime < point.get_timeWindowEnd())
            maxTime = point.get_timeWindowEnd();
        if (point.get_isCustomer())
        {
            auto custId = point.get_idCust();
            if (pointPtsFromCustomerId[custId].empty())
            {
                nbCustomers += 1;
                customerIds.push_back(custId);
                if (point.get_costOrPenalty() > 0)
                    oneOptionalCustomerExists = true;
            }
            pointPtsFromCustomerId[custId].push_back(&point);
            if (point.get_demandOrCapacity() > 0)
                capacityResourceIsUsed = true;
        }
	}

    std::stable_sort(customerIds.begin(), customerIds.end());

	/// Initialisation of vehicleTypes at one for model
    vehicleTypes.reserve(nbVehicleTypes);
    maxVehicleTypeId = 0;
	for (size_t vehicleTypeIndex = 0; vehicleTypeIndex < nbVehicleTypes; vehicleTypeIndex++)
    {
        vehicleTypes.emplace_back(vehicleTypesArray, vehicleTypeIndex);
        maxVehicleTypeId = (std::max)(maxVehicleTypeId, vehicleTypes.back().get_id());
	}

    /// we create artificial depot if needed
    int artificialDepotId = -1;
    for (auto & vehicleType : vehicleTypes)
        if ((vehicleType.get_startPointId() < 0 || vehicleType.get_endPointId() < 0) && artificialDepotId < 0)
        {
            nbPoints += 1;
            points.emplace_back(++maxPointId, 0, 0, 0, minTime, maxTime, 0, true);
            auto & point = points.back();
            pointPtrFromId.push_back(&point);
            artificialDepotId = (int)point.get_id();
        }

    vehicleTypePtrFromId.resize(maxVehicleTypeId + 1, nullptr);
    for (auto & vehicleType : vehicleTypes)
    {
        auto vehTypeId = vehicleType.get_id();
        vehicleTypePtrFromId[vehTypeId] = &vehicleType;

        if (artificialDepotId >= 0 && vehicleType.get_startPointId() >= 0 && vehicleType.get_endPointId() >= 0)
            pointPtrFromId[artificialDepotId]->set_incompatibility(vehTypeId);
        if (vehicleType.get_startPointId() < 0)
            vehicleType.set_startPointId(artificialDepotId);
        if (vehicleType.get_endPointId() < 0)
            vehicleType.set_endPointId(artificialDepotId);
        if (vehicleType.get_maxNbVehicles() > maxTotalVehiclesNumber)
            vehicleType.set_maxNbVehicles(maxTotalVehiclesNumber);

        auto startPointId = vehicleType.get_startPointId();
        auto endPointId = vehicleType.get_endPointId();
        if (startPointId > maxPointId || pointPtrFromId[startPointId] == nullptr)
        {
            std::string message = "Point with vehTypeId : " + std::to_string(startPointId) + " is not found.";
            std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
            setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
            return;
        }
        if (endPointId > maxPointId || pointPtrFromId[endPointId] == nullptr)
        {
            std::string message = "Point with vehTypeId : " + std::to_string(endPointId) + " is not found.";
            std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
            setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
            return;
        }
    }

    /// Initialisation of links
    links.reserve(nbLinks);
	for (size_t linkIndex = 0; linkIndex < nbLinks; linkIndex++)
    {
		links.emplace_back(linksArray, linkIndex);
        auto & link = links.back();
        auto startPointid = link.get_startPointId();
        auto endPointId = link.get_endPointId();
        if (startPointid > maxPointId || pointPtrFromId[startPointid] == nullptr)
        {
            std::string message = "Point with id : " + std::to_string(startPointid) + " is not found.";
            std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
            setErrorCodeAndMsg(ExceptionType::LINKS_ERROR, message);
            return;
        }
        if (endPointId > maxPointId || pointPtrFromId[endPointId] == nullptr)
        {
            std::string message = "Point with id : " + std::to_string(endPointId) + " is not found.";
            std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
            setErrorCodeAndMsg(ExceptionType::LINKS_ERROR, message);
            return;
        }
	}

    /// we create links to the artificial depot if needed
    if (artificialDepotId >= 0)
        for (auto & point : points)
            if (!point.get_isArtificial())
                links.emplace_back(nbLinks++, false, artificialDepotId, point.get_id(), 0.0, 0.0, 0.0, "", true);

    /// get parameters from json
    API_VRP::Parameters param_(param);
    parameters = param_;

    /// every successor (predecessor) is a pair
    /// (link id, <next/prev point is the start (true) of end (false) point of link>)
    std::vector<std::vector<std::pair<size_t, bool>>> successorsFromPointId(maxPointId + 1);
    std::vector<std::vector<std::pair<size_t, bool>>> predecessorsFromPointId(maxPointId + 1);

	/// automatic detections
    if (!checkSolverNbVehicleTypes())
        return;
    if (!checkIntermediateDepots())
        return;
    if (!checkDemandsAndPenaltiesOfCustomers())
        return;
    symmetricCase = determineIfSymmetricCase(timeResourceIsUsed);
    objFuncValueIsInteger = determineIfIntegerObjective();
    preprocessLinks(successorsFromPointId, predecessorsFromPointId);
    if (!checkAndSetIncompatibleVehicles())
        return;
    /// should be done after checkAndSetIncompatibleVehicles
    if (!preprocessVehicleTypes(timeResourceIsUsed, successorsFromPointId, predecessorsFromPointId))
        return;

    if (!determineSolverParameterisation(timeResourceIsUsed, capacityResourceIsUsed, oneOptionalCustomerExists,
                                         successorsFromPointId, predecessorsFromPointId))
        return;
}

std::string API_VRP::runModel(rapidjson::Document & document)
{
    /// get an instance and test if it conforms to the json schema
    try
    {
        rapidjson::Document document_shema;
        rapidjson::SchemaDocument schema(document_shema.Parse(API_VRP::shema));
        rapidjson::SchemaValidator validator(schema);

        if (!document.Accept(validator))
        {
            API_VRP::Serializer json(static_cast<int>(API_VRP::ExceptionType::JSON_ERROR),
                                     "The JSON input does not satisfy the scheme, "
                                     "please save the input to a JSON file and send it to solver developers.");
            return json.get_output();
        }
        else
        {
            API_VRP::Data data(&document);
            if (data.getErrorCode() < 0)
            {
                API_VRP::Serializer json(data.getErrorCode(),
                    data.getMsgError().c_str());
                return json.get_output();
            }
            std::string parametersFileName = data.getParameters().get_configFile();
            if (parametersFileName.empty())
                parametersFileName = "NOT_SPECIFIED";
            bool printHeader = (data.getParameters().get_printlevel() >= -1);
            BcInitialisation bcInit(parametersFileName, false, printHeader, printHeader);
            if (parametersFileName == "NOT_SPECIFIED")
                API_VRP::parameterizeBapcod(data, bcInit);
            else
                data.setParametersFromBapcod(bcInit);
            BcModel model(bcInit);
            API_VRP::createBapcodModel(data, bcInit, model);
            if (data.getParameters().get_enumerate())
            {
                int nbEnumColumns = 0;
                BcSolution enumSol = model.enumerateAllColumns(nbEnumColumns);
                if (nbEnumColumns == 0)
                {
                    API_VRP::Serializer json(static_cast<int>(API_VRP::EnumerationStatus::ENUMERATION_INFEASIBLE),
                                             "No feasible route was found.");
                    return json.get_output();
                }
                else if (nbEnumColumns == -1)
                {
                    API_VRP::Serializer json(static_cast<int>(API_VRP::EnumerationStatus::ENUMERATION_NOT_SUCCEEDED),
                                             "Enumeration cannot be performed on this instance");
                    return json.get_output();
                }
                else
                {
                    API_VRP::Solution solution;
                    if (!solution.obtainFromBcSolution(bcInit, enumSol, data, false))
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::INFEASIBLE_SOLUTION_RETURNED),
                                                 "Solution checker detected infeasibility.");
                        return json.get_output();
                    }
                    else
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::EnumerationStatus::ENUMERATION_SUCCEEDED),
                                                 "Feasible routes was found.", solution,
                                                 data.getParameters().get_enumerate());
                        return json.get_output();
                    }
                }
            }
            else
            {
                BcSolution bcSol = model.solve();
                API_VRP::Solution solution;
                if (bcSol.defined())
                {
                    if (!solution.obtainFromBcSolution(bcInit, bcSol, data))
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::INFEASIBLE_SOLUTION_RETURNED),
                                                 "Solution checker detected infeasibility.");
                        return json.get_output();
                    }

                    if (bcInit.getStatisticCounter("bcFailToSolveModel") == 0)
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::OPTIMAL_SOL_FOUND),
                                                 "The solution found is optimal.", solution, false);
                        return json.get_output();
                    }

                    API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::BETTER_SOL_FOUND),
                                             "The solution is feasible and better than the cut-off value.",
                                             solution, false);
                    return json.get_output();
                }
                else
                {
                    solution.obtainFromBcSolution(bcInit, bcSol, data); /// just to get statistics
                    if (bcInit.getStatisticCounter("bcFailToSolveModel") == 0)
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::BETTER_SOL_DOES_NOT_EXISTS),
                                                 "Better solution doesn't exist.", solution, false);
                        return json.get_output();
                    }
                    else if (solution.get_statistics().solutionTime >= data.getParameters().get_timeLimit())
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::BETTER_SOL_NOT_FOUND),
                                                 "Better solution was not found within the time limit.",
                                                 solution, false);
                        return json.get_output();
                    }
                    else
                    {
                        API_VRP::Serializer json(static_cast<int>(API_VRP::SolutionStatus::INTERRUPTED_BY_ERROR),
                                                 "Solver finished with an error.", solution, false);
                        return json.get_output();
                    }
                }
            }
        }
    }
    catch (API_VRP::ExceptionVRPPython& error)
    {
        API_VRP::Serializer json(static_cast<int>(error.get_type()), error.get_errorStr());
        return json.get_output();
    }
}


 char * solveModel(const char * jsonModel)
{
    rapidjson::Document document;
    document.Parse(jsonModel);
    auto output = API_VRP::runModel(document);
    char* charOutput = new char[INT_MAX];
    strcpy(charOutput, output.c_str());
    charOutput[output.size()] = '\0';
    return charOutput;
}

void freeMemory(char* ptr)
{
    if (ptr)
        free(ptr);
}
void API_VRP::buildArcs(BcFormulation & spForm, const API_VRP::Data & data, const VehicleType & vehType,
                        BcNetwork & network, BcNetworkResource & timeResource, BcNetworkResource & capResource,
                        bool oppositeSense, bool & atLeastOneArcIsUndirected)
{
    int vehTypeId = (int)vehType.get_id();
    int nbVertices = (int)data.getMaxPointId() + 1;
    bool startAndEndPointsAreTheSame = vehType.get_endPointId() == vehType.get_startPointId();
    if (startAndEndPointsAreTheSame)
        nbVertices += 1;
    BcVarArray xVar(spForm, "X");

    atLeastOneArcIsUndirected = false;
    for (int linkId = 0; linkId < (int)data.getLinks().size(); ++linkId)
    {
        auto & link = data.getLinks()[linkId];

        if (!link.get_isFeasibleForSense(oppositeSense) || !link.get_isCompatibleWithVehType(vehTypeId))
        {
            /// we create fake arc just to increase current arc id so that we have the correspondence between link ids
            /// and arc ids in the solver
            network.createArc(nbVertices, nbVertices, 0.0, true);
            continue;
        }

        auto & tailPoint = (oppositeSense) ? data.getHeadPoint(link) : data.getTailPoint(link);
        auto & headPoint = (oppositeSense) ? data.getTailPoint(link) : data.getHeadPoint(link);

        if (!startAndEndPointsAreTheSame
            && (headPoint.get_id() == vehType.get_startPointId() || tailPoint.get_id() == vehType.get_endPointId()))
        {
            network.createArc(nbVertices, nbVertices, 0.0, true);
            continue;
        }

        if (!link.get_isDirected())
            atLeastOneArcIsUndirected = true;

        auto tailCustId = tailPoint.get_idCust();
        auto headCustId = headPoint.get_idCust();

        BcArc arc;
        if (headPoint.get_id() == vehType.get_startPointId())
            arc = network.createArc((int)tailPoint.get_id(), nbVertices - 1, 0.0);
        else
            arc = network.createArc((int)tailPoint.get_id(), (int)headPoint.get_id(), 0.0);

        arc.setName(link.get_name());

        arc.arcVar((BcVar)xVar[vehTypeId][linkId]);

        /// set arc consumption for the capacity resource
        if (capResource.isDefined())
        {
            double capConsumption = 0.0;
            if (data.getIfSymmetricCase())
            {
                if (tailCustId == 0)
                {
                    if (tailPoint.get_id() != vehType.get_startPointId())
                        capConsumption -= (int)vehType.get_capacity() * 0.5;
                }
                else
                {
                    capConsumption += (int)tailPoint.get_demandOrCapacity() * 0.5;
                }
                if (headCustId == 0)
                {
                    if (headPoint.get_id() != vehType.get_endPointId())
                        capConsumption -= (int) vehType.get_capacity() * 0.5;
                }
                else
                {
                    capConsumption += (int)headPoint.get_demandOrCapacity() * 0.5;
                }
            }
            else /// non-symmetric case
            {
                if (headCustId == 0)
                {
                    if (headPoint.get_id() != vehType.get_endPointId())
                        capConsumption -= (int) vehType.get_capacity();
                }
                else
                {
                    capConsumption += (int)headPoint.get_demandOrCapacity();
                }
            }
            capResource.setArcConsumption(arc, capConsumption);
        }

        /// set arc consumption for the time resource
        if (timeResource.isDefined())
        {
            double timeConsumption = 0.0;
            if (data.getIfSymmetricCase())
            {
                timeConsumption += link.get_time();
                if (tailPoint.get_id() == vehType.get_startPointId())
                    timeConsumption += tailPoint.get_serviceTime();
                else
                    timeConsumption += tailPoint.get_serviceTime() * 0.5;
                if (headPoint.get_id() == vehType.get_endPointId())
                    timeConsumption += headPoint.get_serviceTime();
                else
                    timeConsumption += headPoint.get_serviceTime() * 0.5;
            }
            else /// non-symmetric case
            {
                timeConsumption += link.get_time() + tailPoint.get_serviceTime();
            }
            timeResource.setArcConsumption(arc, timeConsumption);
        }
    }
}


BcRCSPFunctor* API_VRP::createRCSPOracle(BcFormulation& spForm, const API_VRP::Data & data, const VehicleType & vehType)
{
    int vehTypeId = (int)vehType.get_id();
    int nbSets = (int)data.getMaxCustId() + 1;

	BcNetwork network(spForm, nbSets, nbSets);

	BcNetworkResource timeResource;
	BcNetworkResource capResource;

	if (vehType.get_timeResourceId() >= 0)
	{
		timeResource = BcNetworkResource(network, vehType.get_timeResourceId());
		timeResource.setAsMainResource();
	}

    if (vehType.get_capacityResourceId() >= 0)
	{
		capResource = BcNetworkResource(network, vehType.get_capacityResourceId() );
		capResource.setAsMainResource();
	}

	int nbVertices = (int)data.getMaxPointId() + 1;
	bool startAndEndPointsAreTheSame = vehType.get_endPointId() == vehType.get_startPointId();

	for (size_t pointId = 0; pointId < nbVertices; ++pointId)
	{
		auto * pointPtr = data.getPointPtr(pointId);

        if (pointPtr == nullptr || !pointPtr->isCompatible(vehTypeId))
        {
            network.createVertex(true); /// we create fake vertex just to increase current vertex id
            continue;
        }

        BcVertex vertex = network.createVertex();
        vertex.setName(pointPtr->get_name());

        if (capResource.isDefined())
        {
            if (pointPtr->get_isDepot() && pointPtr->get_id() != vehType.get_startPointId()
                && pointPtr->get_id() != vehType.get_endPointId())
            {
                /// if an intermediate depot
                if (data.getIfSymmetricCase())
                {
                    capResource.setVertexConsumptionLB(vertex, (int) vehType.get_capacity() * 0.5);
                    capResource.setVertexConsumptionUB(vertex, (int) vehType.get_capacity() * 0.5);
                }
                else
                {
                    capResource.setVertexConsumptionLB(vertex, 0);
                    capResource.setVertexConsumptionUB(vertex, 0);
                }
            }
            else
            {
                capResource.setVertexConsumptionLB(vertex, 0);
                capResource.setVertexConsumptionUB(vertex, (int) vehType.get_capacity());
            }
        }

		if (timeResource.isDefined())
		{
            double lb = pointPtr->get_timeWindowBegin();
            double ub = pointPtr->get_timeWindowEnd();
            if (pointId == vehType.get_startPointId())
                lb = (std::max)(lb, vehType.get_timeWindowBegin());
            else if (pointId == vehType.get_endPointId())
                ub = (std::min)(ub, vehType.get_timeWindowEnd());
			timeResource.setVertexConsumptionLB(vertex, lb);
            timeResource.setVertexConsumptionUB(vertex, ub - pointPtr->get_serviceTime() + API_VRP::epsilonParam);
		}

		if (pointId == vehType.get_startPointId())
		{
			network.setPathSource(vertex);
		}
		else if (pointId == vehType.get_endPointId())
		{
			if (!startAndEndPointsAreTheSame)
				network.setPathSink(vertex);
		}
		else if (!pointPtr->get_isDepot())
		{
			vertex.setElementaritySet((int)pointPtr->get_idCust());
			vertex.setPackingSet((int)pointPtr->get_idCust());
		}
	}

    /// we add an additional vertex for the sink if startAndEndPointsAreTheSame == true
    if (startAndEndPointsAreTheSame)
    {
        BcVertex vertex = network.createVertex();

        if (capResource.isDefined())
        {
            capResource.setVertexConsumptionLB(vertex, 0);
            capResource.setVertexConsumptionUB(vertex, (int) vehType.get_capacity());
        }

        if (timeResource.isDefined())
        {
            auto * pointPtr = data.getPointPtr(vehType.get_endPointId());
            double ub = (std::min)(pointPtr->get_timeWindowEnd(), vehType.get_timeWindowEnd());

            timeResource.setVertexConsumptionLB(vertex, pointPtr->get_timeWindowBegin());
            timeResource.setVertexConsumptionUB(vertex, ub - pointPtr->get_serviceTime() + API_VRP::epsilonParam);
        }

        network.setPathSink(vertex);
        nbVertices += 1;
    }

    network.createVertex(true); /// we always create fake vertex with id = nbVertices

    bool atLeastOneArcIsUndirected;
    buildArcs(spForm, data, vehType, network, timeResource, capResource, false, atLeastOneArcIsUndirected);

    /// if there are undirected links compatible with this vehicle type, we should also add arcs of opposite sense
    /// for every undirected link, id of the opposite sense arc for a link is 'links.size() + link.id'
    if (atLeastOneArcIsUndirected)
        buildArcs(spForm, data, vehType, network, timeResource, capResource, true, atLeastOneArcIsUndirected);

    std::vector<std::vector<double> > distanceMatrix;
    data.getDistanceMatrixBetweenCustomers(distanceMatrix, vehType);
    network.setElemSetsDistanceMatrix(distanceMatrix);

    return new BcRCSPFunctor(spForm);
}

void API_VRP::parameterizeBapcod(const Data & data, BcInitialisation & bcInit)
{
    /// Parameterisation
    bcInit.param().MaxNbOfStagesInColGenProcedure(3);
    bcInit.param().colGenSubProbSolMode().set(SolutionMethod::customSolver);
    bcInit.param().MipSolverMultiThread(1);
    bcInit.param().ApplyStrongBranchingEvaluation(true);

    auto & params = data.getParameters();

    if (data.getIfUseArcMemoryForRankOneCuts())
        bcInit.param().RCSPrankOneCutsMemoryType(1);
    else
        bcInit.param().RCSPrankOneCutsMemoryType(2);

    /// this is for heuristics
    if (params.get_heuristicUsed())
    {
        bcInit.param().RCSPmaxNumOfLabelsInHeurEnumeration(100000);
        bcInit.param().MaxNumEnumSolsInRestrictedMasterIpHeur(5000);
        bcInit.param().MaxTimeForRestrictedMasterIpHeur((int)params.get_timeLimitHeuristic());
        bcInit.param().CallFrequencyOfRestrictedMasterIpHeur(5);
    }

    if (params.get_enumerate())
        bcInit.param().RCSPmaxNumOfLabelsInHeurEnumeration(100000);

    bcInit.param().DEFAULTPRINTLEVEL(params.get_printlevel());

    bcInit.param().GlobalTimeLimitInTick((long)(params.get_timeLimit()*100));

    bcInit.param().solverName(params.get_solverName());
    if (params.get_solverName() == "CLP_SOLVER")
    {
        bcInit.param().RCSPmaxNumOfEnumSolutionsForMIP(0);
        bcInit.param().MaxTimeForRestrictedMasterIpHeur(0);
        bcInit.param().masterSolMode().set(SolutionMethod::lpSolver);
        bcInit.param().relOptimalityGapTolerance(1e-7);
    }

    bcInit.param().RCSPmaxNumOfLabelsInEnumeration(200000);

    /// this is to analyse parameters and do their post-treatment
    bcInit.bcReset();
}

void API_VRP::createBapcodModel(const API_VRP::Data & data, BcInitialisation & bcInit, BcModel & model)
{

    BcObjective objective(model);
    if (data.getIfObjFunctionValueIsInteger())
        objective.setStatus(BcObjStatus::minInt);
    else
        objective.setStatus(BcObjStatus::minFloat);

    if (data.getParameters().get_initialUB() < 10000)
        objective.setArtCostValue(data.getParameters().get_initialUB());
    else if (data.getParameters().get_initialUB() / (double)data.getNbCustomers() * 2.0 > 10000)
        objective.setArtCostValue(data.getParameters().get_initialUB() / (double)data.getNbCustomers() * 2.0);
    else
        objective.setArtCostValue(10000);

    objective <= data.getParameters().get_initialUB();

    BcMaster master(model);

    BcVarArray yVar(master, "Y");
    yVar.type('B');
    yVar.priorityForMasterBranching(1);
    yVar.priorityForSubproblemBranching(-1);

    BcConstrArray degreeConstr(master, "DEG");
    for (auto custId : data.getCustomerIds())
    {
        degreeConstr((int)custId) == 2;
        if (data.getIfOptionalCustomer(custId))
        {
            degreeConstr[(int)custId] += 2 * yVar((int)custId);
            objective += data.getPenalty(custId) * yVar[(int)custId];
        }
    }

    BcConstrArray totalVehNumberConstr(master, "TVN");
    size_t sumVehicleNumbers = 0;
    for (auto & vehType : data.getVehicleTypes())
        sumVehicleNumbers += vehType.get_maxNbVehicles();
    bool imposeTotalVehNumberConstraint = sumVehicleNumbers > data.getMaxTotalVehiclesNumber();
    if (imposeTotalVehNumberConstraint)
        totalVehNumberConstr(0) <= (int)data.getMaxTotalVehiclesNumber();

    BcBranchingConstrArray vehNumberBranching(master, "VNB", SelectionStrategy::MostFractional, 1.0);

    bool totalVehNumBranchDefined = !data.onlyOneVehicleAvailable();
    bool severalVehicleTypes = data.getNbVehicleTypesAvailable() > 1;
    if (totalVehNumBranchDefined)
        vehNumberBranching(0); /// branching on the total number of vehicles type
	if (severalVehicleTypes)
		for (auto &vehType : data.getVehicleTypes())
			vehNumberBranching((int)vehType.get_id());

    BcBranchingConstrArray assignBranching(master, "ASB", SelectionStrategy::MostFractional, 1.0);
    if (severalVehicleTypes)
		for (auto & vehType : data.getVehicleTypes())
		    if (vehType.get_isAvailable())
                for (auto custId : data.getCustomerIds())
                    assignBranching((int)vehType.get_id(), (int)custId);

    BcBranchingConstrArray edgeBranching(master, "EDGE", SelectionStrategy::MostFractional, 1.0);
	for (auto & link : data.getLinks())
    {
        if (link.get_isFeasible())
        {
            if (link.get_startPointId() < link.get_endPointId())
                edgeBranching((int)link.get_startPointId(), (int)link.get_endPointId());
            else
                edgeBranching((int)link.get_endPointId(), (int)link.get_startPointId());
        }
    }

    BcColGenSpArray vehTypeCGSp(model);
    bool resourceCapacityIsUsedInGraphs = false;
	for (auto & vehType : data.getVehicleTypes())
    {
        if (!vehType.get_isAvailable())
            continue;

		int vehTypeId = (int)vehType.get_id();
        BcFormulation spForm = vehTypeCGSp(vehTypeId);

		spForm >= (int) vehType.get_minNbVehicles();
        spForm <= (int) vehType.get_maxNbVehicles();

        BcVarArray xVar(spForm, "X");
        xVar.type('I');
        xVar.priorityForMasterBranching(0.5);
        xVar.priorityForSubproblemBranching(-1);
        xVar.defineIndexNames(MultiIndexNames('k', 'i'));

        for (int linkId = 0; linkId < (int)data.getLinks().size(); ++linkId)
        {
			auto & link = data.getLinks()[linkId];

            if (!link.get_isFeasible() || !link.get_isCompatibleWithVehType(vehTypeId))
                continue;

            auto & tailPoint = data.getTailPoint(link);
            auto & headPoint = data.getHeadPoint(link);
            int tailCustId = (int)tailPoint.get_idCust();
            int headCustId = (int)headPoint.get_idCust();

            BcVar bcVar = xVar(vehTypeId, linkId);
            if (!tailPoint.get_isDepot())
            {
                degreeConstr[tailCustId] += bcVar;
                if (severalVehicleTypes)
                    assignBranching[vehTypeId][tailCustId] += 0.5 * bcVar;
            }
            if (!headPoint.get_isDepot())
            {
                degreeConstr[headCustId] += bcVar;
                if (severalVehicleTypes)
                    assignBranching[vehTypeId][headCustId] += 0.5 * bcVar;
            }

            if (tailCustId <= headCustId)
                edgeBranching[tailCustId][headCustId] += bcVar;
            else
                edgeBranching[headCustId][tailCustId] += bcVar;

            if (tailPoint.get_isDepot() || headPoint.get_isDepot())
            {
                if (totalVehNumBranchDefined)
                    vehNumberBranching[0] += 0.5 * bcVar;
                if (severalVehicleTypes)
                    vehNumberBranching[vehTypeId] += 0.5 * bcVar;
                if (imposeTotalVehNumberConstraint)
                    totalVehNumberConstr[0] += 0.5 * bcVar;
            }

            double arcCost = link.get_fixedCost();
            arcCost += vehType.get_varCostDist() * link.get_distance();
            arcCost += vehType.get_varCostTime() * link.get_time();
            if (tailPoint.get_isDepot() || headPoint.get_isDepot())
                arcCost += vehType.get_fixedCost() * 0.5;
            objective += arcCost * bcVar;
        }

        spForm.attach(createRCSPOracle(spForm, data, vehType));
        if (vehType.get_capacityResourceId() >= 0)
            resourceCapacityIsUsedInGraphs = true;
    }

    if (resourceCapacityIsUsedInGraphs || data.getIfImposeCapacityResourceByCuts())
    {
        int maxCapacity = 0;
        for (auto & vehicleType : data.getVehicleTypes())
        {
            if (maxCapacity < vehicleType.get_capacity())
                maxCapacity = (int)vehicleType.get_capacity(); 
        }

        std::vector<int> demands(data.getMaxCustId() + 1, 0);
        for (auto custId : data.getCustomerIds())
        {
            if (data.getIfOptionalCustomer(custId))
            {
                /// we set demand to zero for optional customers, as they cannot participate in the covering constraints
                demands[custId] = 0;
                continue;
            }
            demands[custId] = (int)data.getDemand(custId);
        }

        bool capacityCutsAreFaculatative = !data.getIfImposeCapacityResourceByCuts();
        BcCapacityCutConstrArray capacityCuts(master, maxCapacity, demands, capacityCutsAreFaculatative,
                                              true, -1, 3.0, 1.0);
    }

    BcLimMemRankOneCutConstrArray limMemRank1Cuts(master);
}

API_VRP::Route::Route(size_t id_) :
    id(id_), vehicleTypeId(0), cost(0.0)
{
}

bool API_VRP::Route::obtainFromBcSolution(const BcSolution & solution, const API_VRP::Data & data)
{
    vehicleTypeId = solution.formulation().id().first();
    cost = solution.cost();
	const BcNetwork network(solution.formulation().network());
	pointIds.reserve(solution.orderedIds().size() + 1);

	if (solution.orderedIds().empty())
		return false;

	auto * vehTypePtr = data.getVehicleTypePtr(vehicleTypeId);
    if (vehTypePtr == nullptr)
        return false;

	auto firstArcId = solution.orderedIds().front();
	auto firstLinkId = firstArcId < data.getLinks().size() ? firstArcId : firstArcId - data.getLinks().size() - 1;
    if (firstLinkId < 0 || firstLinkId >= data.getLinks().size())
        return false;
	auto startPointId = network.getArc(firstArcId).tail().ref();
    auto startPointPtr = data.getPointPtr(startPointId);
    if (startPointPtr == nullptr)
        return false;

    bool timeResourceIsUsed = vehTypePtr->get_timeResourceId() >= 0;
    bool capResourceIsUsed = data.getIfImposeCapacityResourceByCuts() || vehTypePtr->get_capacityResourceId() >= 0;
	double currTime = timeResourceIsUsed ? (std::max)(startPointPtr->get_timeWindowBegin(),
                                                      vehTypePtr->get_timeWindowBegin())
                                           + startPointPtr->get_serviceTime() : 0.0;
	int currDemand = 0;

    if (!startPointPtr->get_isArtificial())
    {
        pointIds.push_back(startPointId);
        pointNames.push_back(startPointPtr->get_name());
        timeConsumption.push_back(currTime);
        capConsumption.push_back(currDemand);
        incomingArcNames.emplace_back("");
    }

	for (auto arcId : solution.orderedIds())
	{
        size_t linkId = arcId < data.getLinks().size() ? arcId : arcId - data.getLinks().size();
        if (linkId < 0 || linkId >= data.getLinks().size())
            return false;
        auto & link = data.getLinks()[linkId];
        int currPointId = network.getArc(arcId).head().ref();
        if (currPointId == data.getMaxPointId() + 1 && vehTypePtr->get_startPointId() == vehTypePtr->get_endPointId())
            currPointId = (int) vehTypePtr->get_startPointId();
        auto currPointPtr = data.getPointPtr(currPointId);
        if (currPointPtr == nullptr)
            return false;

        if (capResourceIsUsed)
        {
            currDemand += (int) currPointPtr->get_demandOrCapacity();
            if (currDemand > vehTypePtr->get_capacity())
                return false;
        }
        if (timeResourceIsUsed)
        {
            currTime = (std::max)(currTime + link.get_time(), currPointPtr->get_timeWindowBegin())
                       + currPointPtr->get_serviceTime();
            if (currTime > currPointPtr->get_timeWindowEnd())
                return false;
        }

        if (!currPointPtr->get_isArtificial())
        {
            pointIds.push_back(currPointId);
            pointNames.push_back(currPointPtr->get_name());
            capConsumption.push_back(currDemand);
            timeConsumption.push_back(currTime);
            incomingArcNames.push_back(link.get_name());
        }
	}

    if (timeResourceIsUsed && currTime > vehTypePtr->get_timeWindowEnd())
        return false;

    /// TO DO : to check that Bapcod solution cost equals to the computed cost

    return true;
}

API_VRP::Statistics::Statistics():
  solutionValue(0.0), solutionTime(0.0), bestLB(0.0), rootLB(0.0), rootTime(0.0), nbBranchAndBoundNodes(0)
{
}

void API_VRP::Statistics::getFromBapcod(BcInitialisation & bcInit, BcSolution & bcSol)
{
    if (bcSol.defined())
        solutionValue = bcSol.cost();
    solutionTime = bcInit.getStatisticTime("bcTimeMain") / 100;
    bestLB = bcInit.getStatisticValue("bcRecBestDb") ;
    rootLB = bcInit.getStatisticValue("bcRecRootDb");
    rootTime = bcInit.getStatisticTime("bcTimeRootEval") / 100;
    nbBranchAndBoundNodes = (int)bcInit.getStatisticCounter("bcCountNodeProc");
}

bool API_VRP::Solution::obtainFromBcSolution(BcInitialisation & bcInit, BcSolution & bcSolution, API_VRP::Data & data,
                                             bool checkForGlobalFeasibility)
{

    statistics.getFromBapcod(bcInit, bcSolution);

    if (!bcSolution.defined())
        return false;

	BcSolution currSol = bcSolution.next();

	size_t routeId = 0;
	while (currSol.defined())
	{
		routes.emplace_back(routeId++);
        if (!routes.back().obtainFromBcSolution(currSol, data))
            return false;
        currSol = currSol.next();
	}

    if (!checkForGlobalFeasibility)
        return true;

	std::vector<size_t> nbUsedVehiclesForVehTypeId(data.getMaxVehicleTypeId() + 1, 0);
    std::vector<size_t> nbCustomerAppearences(data.getMaxCustId() + 1, 0);
	for (auto & route : routes)
	{
        /// existence of this vehicle type is already verified while obtaining the route
        auto * vehTypePtr = data.getVehicleTypePtr(route.vehicleTypeId);

        /// check the availabilities of vehicle type
        nbUsedVehiclesForVehTypeId[vehTypePtr->get_id()] += 1;

        for (auto pointId : route.pointIds)
        {
            auto * pointPtr = data.getPointPtr(pointId);
            if (pointPtr->get_isCustomer())
                nbCustomerAppearences[pointPtr->get_idCust()] += 1;
        }

        if (nbUsedVehiclesForVehTypeId[vehTypePtr->get_id()] > vehTypePtr->get_maxNbVehicles())
            return false;

	}

    for (auto & vehType : data.getVehicleTypes())
        if (nbUsedVehiclesForVehTypeId[vehType.get_id()] > vehType.get_maxNbVehicles())
            return false;

    for (auto & custId : data.getCustomerIds())
    {
        if (data.getIfOptionalCustomer(custId) && nbCustomerAppearences[custId] > 1)
            return false;
        if (!data.getIfOptionalCustomer(custId) && nbCustomerAppearences[custId] != 1)
            return false;
    }

    /// TO DO : to check that Bapcod solution cost equals to the computed cost

    return true;
}

API_VRP::Serializer::Serializer(int code, const char* iMessage)
{
	writer.Reset(buffer);
	writer.StartObject();
	writer.String("Status");
	writer.StartObject();
	add_properties(writer, "code", code);
	add_properties(writer, "message", iMessage);
	writer.EndObject();
	writer.EndObject();

	output = buffer.GetString();
}

API_VRP::Serializer::Serializer(int code, const char* message, const API_VRP::Solution& iSolution, bool enumeration)
{
	// Write Status
	//Serializer(code, message);
	writer.Reset(buffer);
	writer.StartObject();
		writer.String("Status");
		writer.StartObject();
			add_properties(writer, "code", code);
			add_properties(writer, "message", message);
		writer.EndObject();

        auto solutionValue = 0.0;
		// Write Statistics
		if (!enumeration)
		{
			auto statis = iSolution.get_statistics();
			writer.String("Statistics");
			writer.StartObject();
			add_properties(writer, "bestLB", statis.bestLB);
			add_properties(writer, "nbBranchAndBoundNodes", statis.nbBranchAndBoundNodes);
			add_properties(writer, "rootLB", statis.rootLB);
			add_properties(writer, "rootTime", statis.rootTime);
			add_properties(writer, "solutionTime", statis.solutionTime);		
			writer.EndObject();
            solutionValue = statis.solutionValue;

		}
	

	// Write Solution
	auto routes = iSolution.get_routes();
	writer.String("Solution");
    
	writer.StartObject();
    add_properties(writer, "bestSolutionValue", solutionValue);
    

    writer.Key("Routes");
    writer.StartArray();
	for(const auto&route : routes)
	{ 
		writer.StartObject();
		add_properties(writer, "vehicleTypeId", (int)route.vehicleTypeId);
		add_properties(writer, "routeCost", route.cost);
		writer.String("visitedPoints");
		writer.StartArray();
			for (size_t i=0;i<route.pointIds.size();i++)
			{
				writer.StartObject();
					add_properties(writer, "incomingArcName", route.incomingArcNames[i]);
					add_properties(writer, "pointId", route.pointIds[i]);
					add_properties(writer, "pointName", route.pointNames[i]); 
					add_properties(writer, "endTime", route.timeConsumption[i]); 
					add_properties(writer, "load", route.capConsumption[i]); 
				writer.EndObject();
			}
		writer.EndArray();
		writer.EndObject();
	}

	writer.EndArray();
    writer.EndObject();
	writer.EndObject();


	output = buffer.GetString();
}

#endif
