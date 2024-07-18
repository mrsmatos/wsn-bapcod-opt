#ifndef VRPSTW_DATA_H
#define VRPSTW_DATA_H

#include "Singleton.h"

#include <vector>
#include <string>
#include <set>
#include <limits>
#include <cmath>

namespace vrpstw
{

    class Customer
    {
    public:
        Customer(int id = -1, int demand = 0, double x = 0.0, double y = 0.0,
                 double tw_start = -1e12, double tw_end = 1e12, double service_time = 0.0) :
                id(id), demand(demand), x(x), y(y), tw_start(tw_start), tw_end(tw_end), service_time(service_time)
        {}

        int id;
        int demand;
        double x;
        double y;
        double tw_start;
        double tw_end;
        double service_time;
    };

    class Depot
    {
    public:
        Depot(int id = -1, int veh_capacity = 0, int veh_number = 10000, double veh_cost = 0.0,
              double x = 0.0, double y = 0.0, double tw_start = 0.0, double tw_end = 1e12) :
                id(id), veh_capacity(veh_capacity), veh_number(veh_number), veh_cost(veh_cost),  x(x), y(y),
                tw_start(tw_start), tw_end(tw_end)
        {}

        int id;
        int veh_capacity;
        int veh_number;
        double veh_cost;
        double x;
        double y;
        double tw_start;
        double tw_end;
    };

    class Data : public Singleton<Data>
	{
		friend class Singleton<Data>;
	public:

		std::string name;

        /// customers and depots are indexed starting from 1, customers[0] and depots[0] are fictive
		int nbCustomers;
		int nbDepots;
		std::vector<Customer> customers;
        std::vector<Depot> depots;

        enum RoundType
        {
            NO_ROUND,
            ROUND_CLOSEST,
            ROUND_UP,
            ROUND_ONE_DECIMAL
        };

        RoundType roundType;

        double getCustToCustDistance(int firstCustId, int secondCustId) const
        {
            return getDistance(customers[firstCustId].x, customers[firstCustId].y,
                               customers[secondCustId].x, customers[secondCustId].y);
        }

        double getDepotToCustDistance(int depotId, int custId) const
        {
            return getDistance(depots[depotId].x, depots[depotId].y, customers[custId].x, customers[custId].y);
        }

	private:
		Data() :
		    name(), nbCustomers(0), nbDepots(0), customers(1, Customer()),
		    depots(1, Depot()), roundType(RoundType::NO_ROUND)
		{}

		double getDistance(double x1, double y1, double x2, double y2) const
        {
		    double distance = sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) );
		    if (roundType == ROUND_CLOSEST)
		        return round(distance);
		    else if (roundType == ROUND_UP)
		        return ceil(distance);
		    else if (roundType == ROUND_ONE_DECIMAL)
                return floor(distance * 10) / 10;
		    return distance;
        }
	};
}

#endif
