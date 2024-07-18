/**
 *
 * This file dvcpModel.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DvcpModel_h
#define DvcpModel_h

#include "dvcpData.hpp"
#include "bcModelingLanguageC.hpp"

void buildVCPModel(DvcpData & data, const VCPsolution & initSolution, BcModel & model);

#endif // DvcpModel_h
