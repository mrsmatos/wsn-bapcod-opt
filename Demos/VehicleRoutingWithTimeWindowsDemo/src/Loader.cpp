/**
 *
 * This file Loader.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "Loader.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <locale>
#include <cmath>
#include <cstring>

#include "Data.h"
#include "Parameters.h"

vrpstw::Loader::Loader() : data(Data::getInstance()), parameters(Parameters::getInstance()) {}

bool vrpstw::Loader::loadParameters(const std::string & file_name, int argc, char* argv[])
{
    return parameters.loadParameters(file_name, argc, argv);
}

bool vrpstw::Loader::loadData(const std::string & file_name)
{
    std::ifstream ifs(file_name.c_str(), std::ios::in);
    if (!ifs)
    {
        std::cerr << "Instance reader error : cannot open file " << file_name << std::endl;
        return false;
    }

    if (ifs.eof())
    {
        std::cout << "Instance reader error : empty input file " << file_name << std::endl;
        ifs.close();
        return false;
    }

    std::size_t found = file_name.find("VRPTW/");

    data.name = file_name;

    bool success = false;
    if (found != std::string::npos)
        success = loadVRPTWFile(ifs);
    else
        std::cerr << "Instance reader error : cannot determine the instance type " << file_name << std::endl;

    ifs.close();
    return success;
}


bool vrpstw::Loader::loadVRPTWFile(std::ifstream & ifs)
{
    data.roundType = Data::ROUND_ONE_DECIMAL;

    std::string line;
    std::getline(ifs, line);

    switch (line.at(0)) {
        case 'r':
        case 'R':
        case 'c':
        case 'C':
            break;
        default:
            return false;
    }
    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);
    int maxNbVehicles, capacity;
    ifs >> maxNbVehicles >> capacity;

    if ((maxNbVehicles <= 0) || (capacity <= 0) || ifs.eof() || ifs.bad())
        return false;


    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);
    std::getline(ifs, line);

    while (true)
    {
        int customerId = -1, demand, readyTime, dueDate, serviceTime;
        double xCoord, yCoord;
        ifs >> customerId >> xCoord >> yCoord >> demand >> readyTime >> dueDate >> serviceTime;
        if ((customerId < 0) || ifs.eof() || ifs.bad())
            break;
        if (customerId == 0)
        {
            data.depots.push_back(Depot(1, capacity,maxNbVehicles, 0.0, xCoord, yCoord, readyTime, dueDate));
        }
        else
        {
            data.customers.push_back(Customer(customerId, demand, xCoord, yCoord, readyTime, dueDate, serviceTime));
        }
    }
    if (data.depots.empty() || data.customers.empty())
        return false;

    data.nbCustomers = data.customers.size() - 1;
    data.nbDepots = data.depots.size() - 1;

    std::cout << "VRPTW data file detected" << std::endl;

    return true;
}

