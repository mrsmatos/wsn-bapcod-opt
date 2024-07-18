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
#ifndef BCMODELVRPSOLVER_H_
#define BCMODELVRPSOLVER_H_

#include <exception>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <random>

#include <bcModelPointerC.hpp>
#include <bcModelRCSPSolver.hpp>
#include "bcParameterParserC.hpp"


#pragma GCC visibility push(default)

#ifdef _MSC_VER
#define EXPORTED  __declspec( dllexport )
#pragma warning( push )
#pragma warning( disable : 4190 )
#else
#define EXPORTED
#endif

#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type-c-linkage"
#endif

#if(JSON_IS_FOUND && BCP_RCSP_IS_FOUND )

#include "rcsp_interface.hpp"
#include <rapidjson.h>
#include <document.h>
#include <schema.h>
#include <writer.h>
#include <prettywriter.h>


namespace API_VRP {




    const char* const shema = R"__(
{
    "$schema": "https://json-schema.org/draft/2019-09/schema",
    "$id": "http://example.com/example.json",
    "type": "object",
    "default": {},
    "required": ["Points","VehicleTypes","Links","Parameters"],
    "properties": {
		"MaxTotalVehiclesNumber": {
            "type": "integer","default": 10000,"minimum": 1
        },
        "Points": {
            "type": "array",
            "default": [],
            "title": "The Point Schema",
            "minItems": 1,
            "maxItems": 1022,
            "items": {
                "type": "object",
				"required": ["id"],
                "properties": {
                    "name": {
                        "type": "string","default": "" 
                    },
                    "id": {
                        "type": "integer","default": 0,"minimum": 0
                    },
                    "idCustomer": {
                        "type": "integer","default": 0,"minimum": 0,"maximum":1022
                    },
                    "penaltyOrCost": {
                        "type": "number", "default": 0, "minimum" : 0
                    },
                    "serviceTime": {
                        "type": "number","default": 0, "minimum" : 0
                    },
                    "twBegin": {
                        "type": "number","default": 0                       
                    },
                    "twEnd": {
                        "type": "number","default": 0                      
                    },
                    "demandOrCapacity": {
						"type": "integer","default": 0,"minimum": 0
                    },
                    "incompatibleVehicles": {
						"type": "array","default": [], "items": {"type": "integer"}
                    }
                },
				"additionalProperties": false
            }
        },
        "VehicleTypes": {
            "type": "array",
            "default": [],
            "title": "The Vehicles Schema",
            "minItems": 1,
            "maxItems": 50,
            "items": {
                "type": "object",
                "default": {},
                "title": "A Schema",
                "required": ["id","startPointId","endPointId"],
                "properties": {
                    "name": {
                        "type": "string","default": ""   
                    },
                    "id": {
                        "type": "integer", "default": 1, "minimum" : 1
                    },
                    "capacity": {
                        "type": "integer","default": 0, "minimum" : 0
                    },
                    "fixedCost": {
                        "type": "number", "default": 0, "minimum" : 0
                    },
                    "varCostDist": {
                        "type": "number", "default": 0, "minimum" : 0
                    },
                    "varCostTime": {
                        "type": "number", "default": 0, "minimum" : 0
                    },
                    "maxNumber": {
                        "type": "integer", "default": 1000, "minimum" : 1
                    },
                    "startPointId": {
                        "type": "integer", "default": 0, "minimum" : -1
                    },
                    "endPointId": {
                        "type": "integer", "default": 0, "minimum" : -1
                    },
                    "twBegin": {
                        "type": "number","default": 0
                    },
                    "twEnd": {
                        "type": "number","default": 0
                    }
                },
				"additionalProperties": false
            }
        },
        "Links": {
            "type": "array",
            "default": [],
            "title": "The Link Schema",
            "minItems": 1,
            "items": {
                "type": "object",
                "properties": {
					"name": {
                        "type": "string","default": ""   
                    },
                    "isDirected": {
                        "type": "boolean", "default": false
                    },
                    "startPointId": {
                        "type": "integer", "minimum": 0
                    },
                    "endPointId": {
                        "type": "integer", "minimum": 0
                    },
                    "distance": {
                        "type": "number", "default": 0, "minimum": 0
                    },
                    "time": {
                        "type": "number", "default": 0, "minimum": 0
                    },
                    "fixedCost": {
                        "type": "number", "default": 0, "minimum": 0
                    }
                },
				"additionalProperties": false
            }
        },
        "Parameters": {
            "type": "object",
            "default": {},
            "properties": {
                "timeLimit": {
                    "type": "number",
                    "default": 300 
                },
                "upperBound": {
                    "type": "number",
                    "default": 1000000 
                },
                "heuristicUsed": {
                    "type": "boolean",
                    "default": false
                },
				"timeLimitHeuristic": {
                    "type": "number",
                    "default": 20 
                },
                "configFile": {
                    "type": "string",
                    "default": ""
                },
                "solverName": {
                    "type": "string",
                    "default": "CLP",
					"enum": ["CLP","CPLEX"]
                },
                "printLevel": {
                    "type": "integer",
                    "default": -1,
					"enum": [-2,-1,0,1,2]
                },
                "action": {
                    "type": "string",
                    "default": false,
					"enum": ["solve","enumAllFeasibleRoutes"]
                }
            },
			"additionalProperties": false
        }
    }
}
)__";

const char* const shema_output = R"__(
{
    "$schema": "https://json-schema.org/draft/2019-10/schema",
    "$id": "http://example.com/example.json",
    "type": "object",
    "default": {},
    "required": ["Status","Statistics","Solution"],
    "properties": {
		"Status": {
				"type": "object",
				"required":["code","message"],
				"properties": {
					"code": {
						"type": "integer","default": -999999 
					},
					"message": {
						"type": "string","default": ""
					}
				}
        },
        "Statistics": {
            "type": "object",
            "default": {},
            "title": "A Schema",
            "required": ["solutionTime","bestSolutionValue","bestLB","rootLB","rootTime","numberBranchAndBoundNodes"],
            "properties": {
                "solutionTime": {
                    "type": "number","default": 0.0   
                },
                "bestLB": {
                    "type": "number","default": 0
                },
                "rootLB": {
                    "type": "number","default": 0
                },
                "rootTime": {
                    "type": "number","default": 0
                },
                "numberBranchAndBoundNodes": {
                    "type": "integer","default": 0
                }
            }
        },
        "Solution": {
            "type": "object",
            "default": [],
            "title": "The Solution Schema",
			"properties": {
			    "Routes": {
					"items": {
						"type": "object",
						"required": ["vehicleTypeId","routeCost","visitedPoints"],
						"properties": {
							"vehicleTypeId": {
								"type": "integer"  
							},
							"routeCost": {
								"type": "number"
							},
							"visitedPoints": {
								"type": "array",
								"items": {
									"type": "object",
									"required": ["incomingArcName","pointId","pointName","time","load"],
									"properties": {
										"incomingArcName": {
											"type": "string"  
										},
										"pointId": {
											"type": "number"
										},
										"pointName": {
											"type": "string"
										},
										"time": {
											"type": "number"
										},
										"load": {
											"type": "number"
										}
									}
								}
							}
						}
					}
				},
				"bestSolutionValue": {
                    "type": "number", "default": 0
                }
			}
        }
    }
}
)__";

    /// Declaration of strings for json
    const char* const POINTS = "Points";
    const char* const VEHICLE_TYPES = "VehicleTypes";
    const char* const LINKS = "Links";
    const char* const PARAMETERS = "Parameters";
	const char* const MaxTotalVehiclesNumberStr = "MaxTotalVehiclesNumber";
    const char* const nameStr = "name";
    const char* const idStr = "id";
    const char* const idCustomerStr = "idCustomer";
    const char* const penaltyOrCostStr = "penaltyOrCost";
    const char* const serviceTimeStr = "serviceTime";
    const char* const twBeginStr = "twBegin";
    const char* const twEndStr = "twEnd";
    const char* const fixCostStr = "fixedCost";
    const char* const capacityStr = "capacity";
    const char* const demandOrCapacityStr = "demandOrCapacity";
    const char* const isDirectedStr = "isDirected";
    const char* const startPointIdStr = "startPointId";
    const char* const endPointIdStr = "endPointId";
    const char* const distanceStr = "distance";
    const char* const incompatibleVehiclesStr = "incompatibleVehicles";
    const char* const timeStr = "time";
    const char* const varCostDistStr = "varCostDist";
    const char* const varCostTimeStr = "varCostTime";
    const char* const maxNumberStr = "maxNumber";
    const char* const timeLimitStr = "timeLimit";
    const char* const initialUBStr = "upperBound";
    const char* const heuristicUsedStr = "heuristicUsed";
    const char* const timeLimitHeuristicStr = "timeLimitHeuristic";
    const char* const configFileStr = "configFile";
    const char* const solverNameStr = "solverName";
    const char* const actionStr = "action";
    const char* const printlevelStr = "printLevel";

    /// Declaration of json default values
    const int defaultUB = 999999999;

    /// Parameters for automatic detection
    const double rankOneCutsArcMemoryThreshold = 12.0;
    const int nbRandomPathsInAutoParameterisation = 21;
    const double pathDemandRatioThresholdToImposeCapacityByCuts = 0.35;
    const double epsilonParam = 1e-6;

    /// take minimum in tuple using the first element
    template<typename S>
    struct Comparator {
        bool operator()(S& left, S& right) {
            return std::get<0>(left) > std::get<0>(right);
        }
    };


    enum class ExceptionType
	{
		JSON_ERROR = -3,
		VALIDATOR_JSON_FAILED = -4,
		INCONSISTENCY = -5,
		CUSTOMERS_ERROR = -6,
		DEPOTS_ERROR = -7,
		VEHICLES_ERROR = -8,
		LINKS_ERROR = -9,
        POINTS_ERROR = -10
	};

	enum class SolutionStatus
	{
        INFEASIBLE_SOLUTION_RETURNED = -2,
        INTERRUPTED_BY_ERROR = -1,
        OPTIMAL_SOL_FOUND = 0,
		BETTER_SOL_FOUND = 1,
        BETTER_SOL_DOES_NOT_EXISTS = 2,
        BETTER_SOL_NOT_FOUND = 3
	};

	enum class EnumerationStatus
	{
        ENUMERATION_SUCCEEDED = 0,
        ENUMERATION_NOT_SUCCEEDED = 1,
        ENUMERATION_INFEASIBLE = 2
	};

	struct Statistics
	{
		double solutionTime;
		double solutionValue;
		double bestLB;
		double rootLB;
		double rootTime;
		int nbBranchAndBoundNodes;

        Statistics();
		void getFromBapcod(BcInitialisation & bcInit, BcSolution & bcSol);
	};

	enum class typeHeuristic
	{
		NoHeuristic=0,
		Greedy,
		Enumerate
	};

	enum class typePoint
	{
		Customer = 0,
		Depot,
	};

	enum class typeResource
	{
		noResource = 0,
		demand,
		time,
	};

	/**
	* @brief defines all the characteristics of a point
	**/
	class Point
    {
	private:
		std::string name;
		size_t id = 0;
        size_t idCust = 0;
		bool isDepot = false;
		double costORPenalty = 0;
		double serviceTime = 0;
		double timeWindowBegin = 0;
		double timeWindowEnd = 0;
		size_t demandOrCapacity = 0;
        bool isArtificial = false;
		std::set<size_t> incompatibleVehicleTypeIds;
		Point() = default;
	public:
		/**
		* @brief defines a point from customer components
		**/
		Point(size_t iID, size_t idCustomer, double ipenalty = 0, double iServiceTime = 0,
              double iTimeWindowsBegin = 0.0, double iTimeWindowsEnd = 0.0, size_t demand = 0,
              bool iArtificial = false);
		/**
		* @brief defines a point from depot components
		**/
		//Point(size_t iID, double icost=0, double iServiceTime=0, double iTimeWindowsBegin = INT_MIN, double iTimeWindowsEnd = INT_MAX, double capacity= INT_MAX);
		Point(rapidjson::GenericArray<false, rapidjson::Value> & iPoint, size_t index);
		virtual ~Point() = default;

		size_t get_id() const noexcept;
		std::string get_name() const noexcept;
		size_t get_idCust() const noexcept;
		double get_timeWindowBegin() const noexcept;
		double get_timeWindowEnd() const noexcept;
		size_t get_demandOrCapacity() const noexcept;
		double get_serviceTime() const noexcept;
		bool get_isDepot() const noexcept;
        bool get_isCustomer() const noexcept;
		double get_costOrPenalty() const noexcept;
        bool get_isOptional() const noexcept;
        bool get_isArtificial() const noexcept;
		const std::set<size_t> &  get_incompatibleVehicleTypeIds() const noexcept;
		
		void set_incompatibility(rapidjson::GenericArray<false, rapidjson::Value>& iVector);
        void set_incompatibility(size_t vehicleTypeId);
		bool isCompatible(size_t vehicleTypeId) const noexcept;
	};




	/**
	* @brief defines all the characteristics of a vehicleTypes
	**/
	class VehicleType
    {

	private:
		std::string name;
		size_t id =  0;
		size_t index = 0;
		size_t capacity = 0;
		double fixedCost = 0;
		double varCostDist = 0;
		double varCostTime = 0;
		size_t maxNbVehicles = 1;
        size_t minNbVehicles = 0;
		int startPointId = 0;
		int endPointId = 0;
		double timeWindowBegin = 0;
		double timeWindowEnd = 0;
        int timeResourceId = -1;
        int capacityResourceId = -1;

	public:
		VehicleType() = default;
		VehicleType(rapidjson::GenericArray<false, rapidjson::Value> & iVehicle, size_t index);
		virtual ~VehicleType() = default;

		size_t get_id() const noexcept;
		size_t get_Index() const noexcept;
		size_t get_capacity() const noexcept;
		size_t get_maxNbVehicles() const noexcept;
        size_t get_minNbVehicles() const noexcept;
		int get_startPointId() const noexcept;
		int get_endPointId() const noexcept;
		double get_timeWindowBegin() const noexcept;
		double get_timeWindowEnd() const noexcept;

		double get_fixedCost() const noexcept;
		double get_varCostDist() const noexcept;
		double get_varCostTime() const noexcept;
        int get_timeResourceId() const noexcept;
        int get_capacityResourceId() const noexcept;
        bool get_isAvailable() const noexcept;

		void set_maxNbVehicles(size_t value);
        void set_minNbVehicles(size_t value);
        void set_timeResourceId(int resId);
        void set_capacityResourceId(int resId);
        void set_startPointId(int pointId);
        void set_endPointId(int pointId);
    };



	/**
	* @brief defines all the characteristics of an edge
	**/
	class Link
    {

	private:
		std::string name;
		size_t id = 0;
		bool isDirected = false;
		size_t startPointId = 0;
		size_t endPointId = 0;
		double distance = 0;
		double time = 0;
		double fixedCost = 0;
		bool dirIsFeasible = true;
        bool oppIsFeasible = true;
        bool isArtificial = false;
		std::set<size_t> incompatibleVehicleTypeIds;

	public:
		Link() = default;
		//Link(size_t id, size_t startPointId, size_t endPointId, double distance = 0, double time = 0); // On renvoie une erreur si on a pas de startPointId
		Link(rapidjson::GenericArray<false, rapidjson::Value>& iPoint, size_t index);
        Link(size_t iId, bool iDirected, size_t iStartPointId, size_t iEndPointId, double iDistance, double iTime, double iFixedCost,
             const std::string & iName = "", bool iArtificial = false);
		size_t get_id() const noexcept;
		std::string get_name() const noexcept;
		bool get_isDirected() const noexcept;
		size_t get_startPointId() const noexcept;
		size_t get_endPointId() const noexcept;
		double get_distance() const noexcept;
		double get_time() const noexcept;
		double get_fixedCost() const noexcept;
		bool get_isCompatibleWithVehType(size_t vehicleTypeId) const noexcept;
		void set_incompatibilities(const std::set<size_t> & first_point_incompatibilities,
                                   const std::set<size_t> & second_point_incompatibilities);

		/**
		* @brief return true if the link is not feasible
		**/
        bool get_isFeasible() const noexcept;
        bool get_isFeasibleForSense(bool oppositeSense) const noexcept;
        bool get_isArtificial() const noexcept;
		void set_dirIsFeasible(bool result);
        void set_oppIsFeasible(bool result);
		virtual ~Link() = default;

	};

	/**
	* @brief defines all the parameters of model
	**/
	class Parameters
    {
    private:
        double timeLimit = 300;
        double initialUB = 1000000;
        bool heuristicUsed = false;
        double timeLimitHeuristic = 20;
        std::string configFile;
        std::string solverName = "CLP";
        bool enumerate = false;
        int printlevel = -1;

	public:
		Parameters() = default;
		Parameters(double iTimeLimit, double iUB, bool iHeuristicUsed) ;
		explicit Parameters(rapidjson::GenericObject<false, rapidjson::Value> & iParam);
		double get_initialUB() const noexcept;
		double get_timeLimit() const noexcept;
		int get_printlevel() const noexcept;
		bool get_enumerate() const noexcept;
		std::string get_configFile() const noexcept;
		bool get_heuristicUsed() const noexcept;
		double get_timeLimitHeuristic() const noexcept;
		std::string get_solverName() const noexcept;
		virtual ~Parameters() = default;
	};


	/**
	* @brief defines all the data and produce automatic detections
	**/
	class Data
    {
	private:
		std::vector<Link> links;
		std::vector<Point> points;
		std::vector<VehicleType> vehicleTypes;
		Parameters parameters;
		size_t nbPoints = 0;
		size_t nbCustomers = 0;
		size_t nbVehicleTypes = 0;
		size_t nbVehicleTypesAvailable = 0;
		size_t nbLinks = 0;
		size_t maxCustId = 0;
		size_t maxPointId = 0;
        size_t maxVehicleTypeId = 0;
		size_t maxTotalVehiclesNumber = 10000;
		int errorCode = 0;
		std::string errorMsg;

        std::vector<size_t> customerIds; /// all customer ids
		std::vector<std::vector<Point *>> pointPtsFromCustomerId; /// for each customer id, vector of its point ids
		std::vector<Point *> pointPtrFromId; /// for each point id, pointer to its data
                                             /// (nullptr if does not exist)
		std::vector<VehicleType *> vehicleTypePtrFromId;  /// for each vehicle type id, pointer to its data
                                                          /// (nullptr if does not exist)

		/**
		* @brief indicates if objective function is integer or not
		**/
		bool objFuncValueIsInteger = false;


		/**
		* @brief indicates if the graph is symetric or not
		* @details We check 2 informations:
		*			-All times windows are the same
		*			-All arcs are non directed
		**/
		bool symmetricCase = true;

		/**
		* @brief indicates the demands are feasible
		* @details We check 2 informations:
		*			-Check if user didn't indicate the penalty
		*			-if the first check passed, we check if the sum of demand is higher than total capacity
		**/
		size_t totalVehicleCapacity = 0;
		size_t totalNbVehicles = 0;

        bool useArcMemoryForRankOneCuts = false;
        bool imposeCapacityResourceByCuts = false;

		/**
		* @brief will check if the links are feasible and fill vectors and maps for set the critical resources
		**/
		void preprocessLinks(std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                             std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId);

		/**
		* @brief will check if the endPointId is reachable from the start point for each vehicle type
		**/
		bool preprocessVehicleTypes(bool timeResourceIsUsed,
                                    const std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                                    const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId);

		/**
		* @brief will check if the objective function is integer or not
		**/
		bool determineIfIntegerObjective();

		/**
		* @brief will check if the graph is symmetric
		**/
		bool determineIfSymmetricCase(bool timeResourceIsUsed);

        /**
         * @brief find shortest path between all points and the end point of given vehicle type
        **/
        std::vector<double> Dijkstra(const VehicleType & vehicleType,
                                     const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId) const;

        /**
         *  @brief determine 1) which memory type should be used for rank-1 cuts, 2) which resource is critical, and
         *                   3) whether capacity resource should be imposed by rounded capacity cuts
         **/
        bool determineSolverParameterisation(bool timeResourceIsUsed, bool capacityResourceIsUsed,
                                             bool oneOptionalCustomerExists,
                                             const std::vector<std::vector<std::pair<size_t, bool>>> & successorsFromPointId,
                                             const std::vector<std::vector<std::pair<size_t, bool>>> & predecessorsFromPointId);

    public:
		explicit Data(rapidjson::Document *json);

        const Parameters & getParameters() const {return parameters;}

        const std::vector<Link> & getLinks() const {return links;}

        const std::vector<Point> & getPoints() const {return points;}

        const std::vector<VehicleType> & getVehicleTypes() const {return vehicleTypes;}

        bool getIfSymmetricCase() const {return symmetricCase;}

        bool getIfObjFunctionValueIsInteger() const {return objFuncValueIsInteger;}

        size_t getNbCustomers() const {return nbCustomers;}

        size_t getNbVehicleTypesAvailable() const {return nbVehicleTypesAvailable;}

        bool getIfUseArcMemoryForRankOneCuts() const {return useArcMemoryForRankOneCuts;}
        bool getIfImposeCapacityResourceByCuts() const {return imposeCapacityResourceByCuts;}

        const std::vector<size_t> & getCustomerIds() const {return customerIds;}

        const API_VRP::Point & getTailPoint(const Link & iLink) const;

		const API_VRP::Point & getHeadPoint(const Link & iLink) const;

		const API_VRP::Point * getPointPtr(size_t pointId) const;

        const API_VRP::VehicleType * getVehicleTypePtr(size_t vehTypeId) const;

        void setParametersFromBapcod(BcInitialisation & bcInit);

        /**
		* @brief return true if there is only one vehicule
		**/
		bool onlyOneVehicleAvailable() const noexcept;

		/**
		* @brief return true if the maxmimum id from customers
		**/
		size_t getMaxCustId() const;

		/**
		* @brief return demand from a custId
		**/
		size_t getDemand(size_t custId) const;

		/**
		* @brief return penalty from a custId
		**/
		double getPenalty(size_t custId) const;

        /**
          * @brief return if customer is optional from a custId
        **/
        bool getIfOptionalCustomer(size_t custId) const;

		/**
		  * @brief return the maximum total vehicles number
		**/
		size_t getMaxTotalVehiclesNumber() const;

        /**
		* @brief for each point we check a vector of incompatibility and
		 *       set incompatible vehicles for links
		**/
		 bool checkAndSetIncompatibleVehicles();

		 /**
		 * @brief return the greatest id from points
		 **/
		 size_t getMaxPointId() const;

         size_t getMaxVehicleTypeId() const;

		 int getErrorCode() const;

		 const std::string& getMsgError() const;

         void getDistanceMatrixBetweenCustomers(std::vector<std::vector<double>> & matrix,
                                                const VehicleType & vehType) const;
		 
		 /**
		 * @brief checks that for each group of customers, all have the same request and the same penalty
		 **/
		 bool checkDemandsAndPenaltiesOfCustomers();

		 /**
		 * @brief checks the number of vehicle types depending on the solver.
		 * @details If the solver is CPLEX we allow many vehicle type but 
		 * if the solver is CLP we must have only one vehicle type.
		 **/
		 bool checkSolverNbVehicleTypes();

         bool checkIntermediateDepots();

        /**
		 * @brief set error code and message if an is thrown
		 **/
		 void setErrorCodeAndMsg(ExceptionType iException,std::string iMsg) ;

		 /**
		 * @brief check if a point is valid
		 * @details -the time window begin must be greater than time window end
		 * -the addition of time window begin and service time must be less than time window end
		 **/
		 void checkPoint(Point& iPoint);


		 /**
		 * @brief check if a vehicle type is valid
		 * @details the time window begin must be is greater than time window end
		 **/
		 void checkVehicleType(VehicleType& iVehicleType);

		 virtual ~Data() = default;;

	};

	struct Route
	{
		size_t id;
		size_t vehicleTypeId;
		double cost;
		std::vector<int> pointIds;
		std::vector<std::string> pointNames;
		std::vector<double> capConsumption;
		std::vector<double> timeConsumption;
		std::vector<std::string> incomingArcNames;

		explicit Route(size_t id_);
        bool obtainFromBcSolution(const BcSolution & solution, const API_VRP::Data & data);
	};

	/**
	* @brief defines a solution for all vehicules type
	**/
	class Solution {

	private:
        double cost;
		std::vector<Route> routes;
		Statistics statistics;
	public:
        Solution() : cost(0.0) {}
        /**
        * @brief return false if solution is infeasible
        **/
		bool obtainFromBcSolution(BcInitialisation & bcInit, BcSolution& bcSolution, API_VRP::Data& data,
                                  bool checkForGlobalFeasibility = true);
		const Statistics & get_statistics() const noexcept;
		const std::vector<Route> & get_routes() const noexcept;
        const double get_cost() const noexcept;
		virtual ~Solution() = default;
	};

	/**
	* @brief construct output json
	**/
	class Serializer {

	private:
		rapidjson::Document document;
		rapidjson::StringBuffer buffer;
		rapidjson::PrettyWriter<rapidjson::StringBuffer> writer;
		const char* output;

	public:
		Serializer() : output("") {}
		Serializer(int code, const char* message); /// first constructor is called if there are no solutions
		Serializer(int code, const char* message,const API_VRP::Solution& iSolution,bool enumeration); //is called for action solve and enumeration if it's feasible
		std::string get_output() const;
		void add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, int value);
		void add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, const std::string& value);
		void add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, const char* value);
		void add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, double value);

		virtual ~Serializer() = default;
	};

	/**
	* @brief defines all the exceptions of model
	**/
	class ExceptionVRPPython : public std::exception
    {
	public:
		ExceptionVRPPython(ExceptionType error_type, const std::string& message);
		char* get_message() const noexcept;
		const char* get_errorStr() const noexcept;
		ExceptionType get_type() const noexcept;
	private:
		ExceptionType error_type;
		char* message = nullptr;
		char* errorStr = nullptr;
	};


    void buildArcs(BcFormulation & spForm, const API_VRP::Data & data, const VehicleType & vehType,
                   BcNetwork & network, BcNetworkResource & timeResource, BcNetworkResource & capResource,
                   bool oppositeSense, bool & atLeastOneArcIsUndirected);
    BcRCSPFunctor * createRCSPOracle(BcFormulation& spForm, const Data & data, const VehicleType & vehType);
    void parameterizeBapcod(const Data & data, BcInitialisation & bcInit);
    void createBapcodModel(const Data & data, BcInitialisation & bcInit, BcModel & model);
}
inline const std::vector<API_VRP::Route>& API_VRP::Solution::get_routes() const noexcept
{
	return routes;
}

inline const double API_VRP::Solution::get_cost() const noexcept
{
	return cost;
}

inline const API_VRP::Statistics& API_VRP::Solution::get_statistics() const noexcept
{
	return statistics;
}

inline std::string API_VRP::Serializer::get_output() const
{
//	std::string out_ = output;
//	char* output_ = new char[out_.size()];
//	strcpy(output_, out_.c_str());
	return {output};
}

inline size_t API_VRP::Data::getMaxPointId() const
{
	return maxPointId;
}

inline size_t API_VRP::Data::getMaxVehicleTypeId() const
{
    return maxVehicleTypeId;
}

inline int API_VRP::Data::getErrorCode() const
{
	return errorCode;
}

inline const std::string& API_VRP::Data::getMsgError() const
{
	return errorMsg;
}

inline bool API_VRP::Data::checkDemandsAndPenaltiesOfCustomers()
{
	for (auto custId : customerIds)
        if (!pointPtsFromCustomerId[custId].empty())
        {
            auto * firstPointPtr = pointPtsFromCustomerId[custId].front();
            auto demand = firstPointPtr->get_demandOrCapacity();
            auto penalty = firstPointPtr->get_costOrPenalty();
            for (auto * pointPtr : pointPtsFromCustomerId[custId])
                if (pointPtr->get_demandOrCapacity() != demand || pointPtr->get_costOrPenalty() != penalty)
                {
					std::string message = "All points with the same customer id must have the same penalty "
                                          "and the same demand";
					std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
					setErrorCodeAndMsg(ExceptionType::CUSTOMERS_ERROR, message);
                    return false;
                }
        }
    return true;
}

inline bool API_VRP::Data::checkSolverNbVehicleTypes()
{
    /// deactivated : more than one vehicle type is allowed even in the free version
//#ifndef _CPLEX_FOUND
//    if (nbVehicleTypes > 1)
//		setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR,
//                           "More than one type of vehicle is allowed only in the academic version.");
//#endif
    return true;
}

inline bool API_VRP::Data::checkIntermediateDepots()
{
    for (auto & point : points)
        if (point.get_isDepot())
            for (auto & vehType : vehicleTypes)
            {
                if (vehType.get_startPointId() != point.get_id() && vehType.get_endPointId() != point.get_id()
                    && !point.get_incompatibleVehicleTypeIds().count(vehType.get_id()))
                {
                    std::string message = "Intermediate depots are not allowed for the moment.";
                    std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
                    setErrorCodeAndMsg(ExceptionType::DEPOTS_ERROR, message);
                    return false;
                }
            }
    return true;
}

inline void API_VRP::Data::setErrorCodeAndMsg(ExceptionType iException, std::string iMsg)
{
	errorCode = static_cast<int>(iException);
	errorMsg = std::move(iMsg);
}

inline void API_VRP::Data::checkPoint(Point& iPoint)
{
	if (iPoint.get_timeWindowBegin() > iPoint.get_timeWindowEnd())
	{
		std::string message = "Time windows begin is greater than time windows end.";
		std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
		setErrorCodeAndMsg(iPoint.get_idCust() == 0 ? ExceptionType::DEPOTS_ERROR : ExceptionType::CUSTOMERS_ERROR,
                           message);
	}


	if (iPoint.get_timeWindowBegin() + iPoint.get_serviceTime() > iPoint.get_timeWindowEnd())
	{
		std::string message = "There is insufficient time to serve a node";
		std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
		setErrorCodeAndMsg(iPoint.get_idCust() == 0 ? ExceptionType::DEPOTS_ERROR : ExceptionType::CUSTOMERS_ERROR,
                           message);
	}
}

inline void API_VRP::Data::checkVehicleType(VehicleType& iVehicleType)
{
	if (iVehicleType.get_timeWindowBegin() > iVehicleType.get_timeWindowEnd())
	{
		std::string message = "Time windows begin is greater than time windows end.";
		std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
		setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
	}
	if (iVehicleType.get_capacity() < 0)
	{
		std::string message = "The capacity must be greater than 0.";
		std::cout << "VRPSolverEasy ERROR : " << message << std::endl;
		setErrorCodeAndMsg(ExceptionType::VEHICLES_ERROR, message);
	}
}

inline void API_VRP::Serializer::add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, double value)
{
	writer.String(nameProperty);
	writer.Double(value);
}
inline void API_VRP::Serializer::add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, const char* value)
{
	writer.String(nameProperty);
	writer.String(value);
}
inline void API_VRP::Serializer::add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty, int value)
{
	writer.String(nameProperty);
	writer.Int(value);
}

inline void API_VRP::Serializer::add_properties(rapidjson::PrettyWriter<rapidjson::StringBuffer>& writer, const char* nameProperty,const std::string& value)
{
	writer.String(nameProperty);
	writer.String(value.c_str());
}

inline const API_VRP::Point* API_VRP::Data::getPointPtr(size_t pointId) const
{
    return (pointId <= maxPointId) ? pointPtrFromId[pointId] : nullptr;
}

inline const API_VRP::VehicleType * API_VRP::Data::getVehicleTypePtr(size_t vehTypeId) const
{
    return (vehTypeId <= maxVehicleTypeId) ? vehicleTypePtrFromId[vehTypeId] : nullptr;
}

inline size_t API_VRP::Data::getDemand(size_t custId) const
{
    if (custId > maxCustId || pointPtsFromCustomerId[custId].empty())
        return 0;

    return pointPtsFromCustomerId[custId].front()->get_demandOrCapacity();
}


inline double API_VRP::Data::getPenalty(size_t custId) const
{
    if (custId > maxCustId || pointPtsFromCustomerId[custId].empty())
        return 0;

    return pointPtsFromCustomerId[custId].front()->get_costOrPenalty();
}

inline bool API_VRP::Data::getIfOptionalCustomer(size_t custId) const
{
    if (custId > maxCustId || pointPtsFromCustomerId[custId].empty())
        return true;

    return (pointPtsFromCustomerId[custId].front()->get_costOrPenalty() > 0.0);
}

inline size_t API_VRP::Data::getMaxTotalVehiclesNumber() const
{
	return maxTotalVehiclesNumber;
}

inline size_t API_VRP::Data::getMaxCustId() const
{
	return maxCustId;
}

inline bool API_VRP::Data::onlyOneVehicleAvailable() const noexcept
{
	return totalNbVehicles == 1;
}

inline const API_VRP::Point& API_VRP::Data::getTailPoint(const Link & iLink) const
{
	return *pointPtrFromId[iLink.get_startPointId()];
}

inline const API_VRP::Point& API_VRP::Data::getHeadPoint(const Link & iLink) const
{
	return *pointPtrFromId[iLink.get_endPointId()];
}

inline std::string API_VRP::Parameters::get_configFile() const noexcept
{
	return configFile;
}
inline bool API_VRP::Parameters::get_heuristicUsed() const noexcept
{
	return heuristicUsed;
}

inline double API_VRP::Parameters::get_timeLimitHeuristic() const noexcept
{
	return timeLimitHeuristic;
}

inline std::string API_VRP::Parameters::get_solverName() const noexcept
{
	return solverName + "_SOLVER";
}

inline bool API_VRP::Parameters::get_enumerate() const noexcept
{
	return enumerate;
}

inline double API_VRP::Parameters::get_initialUB() const noexcept
{
	return initialUB;
}

inline double API_VRP::Parameters::get_timeLimit() const noexcept
{
	return timeLimit;
}

inline int API_VRP::Parameters::get_printlevel() const noexcept
{
	return printlevel;
}

inline const char* API_VRP::ExceptionVRPPython::get_errorStr() const noexcept
{
	return errorStr;
}
inline API_VRP::ExceptionType API_VRP::ExceptionVRPPython::get_type() const noexcept
{
	return error_type;
}

inline char* API_VRP::ExceptionVRPPython::get_message() const noexcept
{
	return message;
}

inline size_t API_VRP::Link::get_id() const noexcept
{
	return id;
}

inline double API_VRP::Link::get_fixedCost() const noexcept
{
	return fixedCost;
}
inline double API_VRP::Link::get_distance() const noexcept
{
	return distance;
}
inline double API_VRP::Link::get_time() const noexcept
{
	return time;
}
inline size_t API_VRP::Link::get_endPointId() const noexcept
{
	return endPointId;
}

inline size_t API_VRP::Link::get_startPointId() const noexcept
{
	return startPointId;
}

inline bool API_VRP::Link::get_isDirected() const noexcept
{
	return isDirected;
}

inline bool API_VRP::Link::get_isArtificial() const noexcept
{
    return isArtificial;
}

inline bool API_VRP::Link::get_isFeasible() const noexcept
{
    return (isDirected && dirIsFeasible) || (!isDirected && (dirIsFeasible || oppIsFeasible));
}

inline bool API_VRP::Link::get_isFeasibleForSense(bool oppositeSense) const noexcept
{
    if (oppositeSense)
        return oppIsFeasible;
    return dirIsFeasible;
}

inline bool API_VRP::Link::get_isCompatibleWithVehType(size_t vehicleTypeId) const noexcept
{
    return !incompatibleVehicleTypeIds.count(vehicleTypeId);
}

inline void API_VRP::Link::set_dirIsFeasible(bool result)
{
    dirIsFeasible = result;
}

inline void API_VRP::Link::set_oppIsFeasible(bool result)
{
    oppIsFeasible = result;
}

inline std::string API_VRP::Link::get_name() const noexcept
{
	return name;
}

inline void API_VRP::Link::set_incompatibilities(const std::set<size_t> & first_point_incompatibilities,
                                                 const std::set<size_t> & second_point_incompatibilities)
{
	for (auto vehicleTypeId : first_point_incompatibilities)
        incompatibleVehicleTypeIds.insert(vehicleTypeId);
    for (auto vehicleTypeId : second_point_incompatibilities)
        incompatibleVehicleTypeIds.insert(vehicleTypeId);
}

inline void API_VRP::VehicleType::set_maxNbVehicles(size_t value)
{
    maxNbVehicles = value;
}

inline void API_VRP::VehicleType::set_minNbVehicles(size_t value)
{
    minNbVehicles = value;
}

inline void API_VRP::VehicleType::set_timeResourceId(int resId)
{
    timeResourceId = resId;
}

inline void API_VRP::VehicleType::set_capacityResourceId(int resId)
{
    capacityResourceId = resId;
}

inline void API_VRP::VehicleType::set_startPointId(int pointId)
{
    startPointId = pointId;
}

inline void API_VRP::VehicleType::set_endPointId(int pointId)
{
    endPointId = pointId;
}

inline size_t API_VRP::VehicleType::get_id() const noexcept
{
	return id;
}

inline size_t API_VRP::VehicleType::get_Index() const noexcept
{
	return index;
}

inline double API_VRP::VehicleType::get_timeWindowEnd() const noexcept
{
	return timeWindowEnd;
}

inline double API_VRP::VehicleType::get_timeWindowBegin() const noexcept
{
	return timeWindowBegin;
}

inline int API_VRP::VehicleType::get_endPointId() const noexcept
{
	return endPointId;
}

inline int API_VRP::VehicleType::get_startPointId() const noexcept
{
	return startPointId;
}

inline size_t API_VRP::VehicleType::get_maxNbVehicles() const noexcept
{
	return maxNbVehicles;
}

inline size_t API_VRP::VehicleType::get_minNbVehicles() const noexcept
{
    return minNbVehicles;
}

inline size_t API_VRP::VehicleType::get_capacity() const noexcept
{
	return capacity;
}

inline  double API_VRP::VehicleType::get_fixedCost() const noexcept
{
	return fixedCost; // fixedCost contains int, not float: will throw
}

inline double API_VRP::VehicleType::get_varCostDist() const noexcept
{
	return varCostDist;
}

inline double API_VRP::VehicleType::get_varCostTime() const noexcept
{
	return varCostTime;
}

inline int API_VRP::VehicleType::get_timeResourceId() const noexcept
{
    return timeResourceId;
}

inline int API_VRP::VehicleType::get_capacityResourceId() const noexcept
{
    return capacityResourceId;
}

inline bool API_VRP::VehicleType::get_isAvailable() const noexcept
{
    return maxNbVehicles > 0;
}

inline void API_VRP::Point::set_incompatibility(rapidjson::GenericArray<false, rapidjson::Value>& iVector)
{
	for (size_t IndexVeh = 0; IndexVeh < iVector.Size(); IndexVeh++)
        incompatibleVehicleTypeIds.insert(iVector[IndexVeh].GetInt());
}

inline void API_VRP::Point::set_incompatibility(size_t vehicleTypeId)
{
    incompatibleVehicleTypeIds.insert(vehicleTypeId);
}

inline bool API_VRP::Point::isCompatible(size_t vehicleTypeId) const noexcept
{
    return !incompatibleVehicleTypeIds.count(vehicleTypeId);
}


//inline bool API_VRP::Data::capacityResourceIsImposedByCuts(const VehicleType& vehType) const
//{
//	if (oneOptionalCustomerExists)
//		return false;
//	else
//		return std::get<0>(vehicleTypesCR[vehType.get_Index()]);
//}

inline double API_VRP::Point::get_costOrPenalty() const noexcept
{
	return costORPenalty;
}

inline bool API_VRP::Point::get_isOptional() const noexcept
{
    return costORPenalty > 0.0;
}

inline bool API_VRP::Point::get_isArtificial() const noexcept
{
    return isArtificial;
}

inline const std::set<size_t> & API_VRP::Point::get_incompatibleVehicleTypeIds() const noexcept
{
	return incompatibleVehicleTypeIds;
}


inline bool API_VRP::Point::get_isDepot() const noexcept
{
	return isDepot;
}

inline bool API_VRP::Point::get_isCustomer() const noexcept
{
    return !isDepot;
}

inline double API_VRP::Point::get_serviceTime() const noexcept
{
	return serviceTime;
}

inline double API_VRP::Point::get_timeWindowBegin() const noexcept
{
	return timeWindowBegin;
}

inline double API_VRP::Point::get_timeWindowEnd() const noexcept
{
	return timeWindowEnd;
}

inline size_t API_VRP::Point::get_demandOrCapacity() const noexcept
{
	return demandOrCapacity;
}

inline size_t API_VRP::Point::get_id() const noexcept
{
	return id;
}

inline std::string API_VRP::Point::get_name() const noexcept
{
	return name;
}

inline size_t API_VRP::Point::get_idCust() const noexcept
{
	return idCust;
}

namespace API_VRP
{
    std::string runModel(rapidjson::Document &document);
}

extern "C" {
	EXPORTED  char* solveModel(const char * jsonModel);
	EXPORTED  void freeMemory(char* ptr);

}

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#ifdef __clang__
#pragma GCC diagnostic pop
#endif

#pragma GCC visibility pop

#endif
#endif
