/**
 *
 * This file dvcpParameters.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dvcpParameters.hpp"


VCPdemoParameters::VCPdemoParameters(std::string parameterFileName, int argc, char * argv[], std::ostream & os) :
    ParameterParser(parameterFileName),
    printLevelDVCP("printLevelDVCP", 0, ""),
    initialHeuristic("initialHeuristic", 2, "")
{
  addApplicationParameter(printLevelDVCP);
  addApplicationParameter(initialHeuristic);

  parse(argc, argv);

}

std::ostream & VCPdemoParameters::print(std::ostream& os) const
{
    os << "Vertex Coloring Demo Parameters  " << std::endl;
    
    os << printLevelDVCP << std::endl;
    os << initialHeuristic << std::endl;

    
    return (os);
}
