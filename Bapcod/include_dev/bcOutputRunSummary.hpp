/**
 *
 * This file bcOutputRunSummary.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCOUTPUTRUNSUMMARY_HPP_
#define BCOUTPUTRUNSUMMARY_HPP_

#include "bcBapcodInit.hpp"
#include "bcModelC.hpp"

/**
 * Print the primal and dual bound solution inside a file.
 * The filename is the parameter NameOfOutputSummaryFile.
 *
 * TODO: Change this class to be more generic: without
 * dependencies to BapcodInit or Model.
 */
class OutputRunSummary
{
public:
  OutputRunSummary(BapcodInit* bapcodInitPtr, Model* modelPtr);
  OutputRunSummary(BapcodInit* bapcodInitPtr, const std::string& instanceName);
  virtual ~OutputRunSummary();
};

#endif /* BCOUTPUTRUNSUMMARY_HPP_ */
