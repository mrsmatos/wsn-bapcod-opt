/**
 *
 * This file Parameters.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "Parameters.h"

vrpstw::Parameters::Parameters() :
        silent("silent", false),
        cutOffValue("cutOffValue", std::numeric_limits<double>::infinity())
{
}

bool vrpstw::Parameters::loadParameters(const std::string & parameterFileName, int argc, char* argv[])
{
    setParameterFileName(parameterFileName);
    addApplicationParameter(silent);
    addApplicationParameter(cutOffValue);
    parse(argc, argv);

    return true;
}
