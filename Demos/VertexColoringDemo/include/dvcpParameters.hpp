/**
 *
 * This file dvcpParameters.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DvcpParameters_h
#define DvcpParameters_h

#include "bcParameterParserC.hpp"

#define NONE_INIT_HEURISTIC 0
#define SIMPLE_INIT_HEURISTIC 1
#define DSATUR_INIT_HEURISTIC 2

/* PROBLEM SPECIFIC PARAMETERS */
class VCPdemoParameters : public ParameterParser
{
public:
    ApplicationParameter<int> printLevelDVCP;
    ApplicationParameter<int> initialHeuristic;

    VCPdemoParameters(std::string parameterFileName, int argc, char* argv[], std::ostream & os = std::cout);
    virtual ~VCPdemoParameters()
    {
    }

    std::ostream & print(std::ostream & os = std::cout) const;
};

#endif // DvcpParameters_h
