/**
 *
 * This file dvcpHeuristics.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DvcpHeuristics_h
#define DvcpHeuristics_h

#include "dvcpData.hpp"

void runSimpleGreedyHeuristic(DvcpData & data, const int printLevel, VCPsolution & solution);

#ifndef _MSC_VER
void runDsaturGreedyHeuristic(DvcpData & data, const int printLevel, VCPsolution & solution);
#endif

#endif // DvcpHeuristics_h
