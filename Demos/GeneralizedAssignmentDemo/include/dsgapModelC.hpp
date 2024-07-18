/**
 *
 * This file dsgapModelC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef StandardCutStockProbModelClass_h
#define StandardCutStockProbModelClass_h

#include "bcUsefulHeadFil.hpp"
#include "bcModelFormulationC.hpp"

#include "dsgapDataC.hpp"


BcFormulation buildMpModel(GapData & dataStruct, BcModel & gapModel);

void buildColGenModel(GapData & dataStruct, BcModel & gapModel);

void initPrimalHeur(GapData & dataStruct, BcModel & gapModel);


#endif // StandardCutStockProbModelClass_h
