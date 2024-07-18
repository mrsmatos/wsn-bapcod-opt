/**
 *
 * This file bcModelMasterC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcModelC.hpp"
#include "bcModelingLanguageC.hpp"
#include "bcPrintC.hpp"
#include "bcProbConfigC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcSolutionC.hpp"
/**
 * Generic code
 */
using namespace std;

BcMaster::BcMaster(BcModel & model, const std::string & name):
    BcFormulation()
{
    Model * modelPtr = (Model *)model;
    _probConfPtr = modelPtr->createMasterConf(name, MultiIndex());
}

BcMaster::BcMaster(BcFormulation & bcForm) :
    BcFormulation(bcForm)
{
    if (!bcForm.isMaster())
    {
        std::cerr << "BaPCod error : formulation is not master" << std::endl;
        exit(1);
    }
}

BcMasterArray::BcMasterArray(BcModel & model, const std::string & name) :
  BcFormulationArray(model, name)
{
    _modelPtr->createMasterConf(_name, MultiIndex());
}

BcFormulation & BcMasterArray::getElement(const MultiIndex & multiIndexId)
{
    if (_curForm.isDefined())
    {
        if (_curForm.id() == multiIndexId)
            return _curForm;

        _curForm = BcFormulation();
    }
    return _curForm;
}

BcFormulation & BcMasterArray::createElement(const MultiIndex & multiIndexId)
{
    if (_curForm.isDefined() && (_curForm.id() == multiIndexId))
        return _curForm;

    return _curForm = BcFormulation(_modelPtr->createMasterConf(_name, multiIndexId));
}

