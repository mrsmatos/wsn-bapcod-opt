/**
 *
 * This file bcGenMastVarConstrC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenMastVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcModelC.hpp"
#include "bcModelVarC.hpp"
#include "bcMastVarConstrC.hpp"

#include "bcPrintC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcVarConstrC.hpp"
#include "bcGlobalException.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcFormC.hpp"

/**
 * Generic code
 */
using namespace std;

/***************************************************************************
 *************   TODO: Methods for ConvexityGenConstr      *****************
 ***************************************************************************/

ConvexityGenConstr::ConvexityGenConstr(Model * modelPtr, MasterConf * masterPtr,
                                       const std::vector<ColGenSpConf *> & colGenSubProbConfPts) :
    GenericConstr(modelPtr, BcVarConstrType::local2Formulation, masterPtr, "")
{
  IndexCell id(-1);
  InstMastConvexityConstr * mccPtr(NULL);

  char implicitOrExplicitType((param().SplitColIntoDissagregateSpVar() ? 'I' : 'E'));

  for (std::vector<ColGenSpConf *>::const_iterator spcPt = colGenSubProbConfPts.begin();
       spcPt != colGenSubProbConfPts.end(); spcPt++)
  {
    if (printL(7))
      std::cout << "sp conf " << (*spcPt)->name() << std::endl;

    /// we force creation of all convexity constraints later redundant constraints will be deactivated
    if ((*spcPt)->upperBoundPtr() == NULL)
      (*spcPt)->upperBoundPtr(new Double(InstMastConvexityConstr::upperBoundWhenInactive));
    if ((*spcPt)->lowerBoundPtr() == NULL)
      (*spcPt)->lowerBoundPtr(new Double(InstMastConvexityConstr::lowerBoundWhenInactive));

    id = IndexCell((*spcPt)->ref());
    if ((*spcPt)->upperBoundPtr() != NULL)
    {
      if (printL(7))
        std::cout << "upper bound " << *(*spcPt)->upperBoundPtr() << std::endl;

        std::string name("convU");
        name = name + (*spcPt)->ref();
        (*spcPt)->id().appendRef2name(name);
        mccPtr = new InstMastConvexityConstr(id,
                                             this,
                                             masterPtr,
                                             *spcPt,
                                             name,
                                             *(*spcPt)->upperBoundPtr(),
                                             'L',
                                             'X', /// this is view as a subproblem constraint
                                             implicitOrExplicitType, //'E',
                                             's',
                                             -1,
                                             (*spcPt)->defaultDualVal4UbConstr(),
                                             BapcodInfinity,
                                             0);
        if (printL(5))
          std::cout << " new convexity constr  " << *mccPtr << std::endl;

        (*spcPt)->upperBoundMastConstrPtr(mccPtr);

        if ((param().mastInitMode().status() == MasterInitMode::localAndGlobAc)
                || (param().mastInitMode().status() == MasterInitMode::localArtCol)
                || (param().mastInitMode().status() == MasterInitMode::incSolColAndLac))
        {
          mccPtr->addLocalArtVar(modelPtr->objectiveSense());
        }
        if (masterPtr->probPtr()->posGlobalArtVarPtr() != NULL)
          mccPtr->addMember(masterPtr->probPtr()->posGlobalArtVarPtr());
        if (masterPtr->probPtr()->negGlobalArtVarPtr() != NULL)
          mccPtr->addMember(masterPtr->probPtr()->negGlobalArtVarPtr());
    }

    if ((*spcPt)->lowerBoundPtr() != NULL)
    {
        std::string name("convL");
        name = name + (*spcPt)->ref();
        (*spcPt)->id().appendRef2name(name);
        mccPtr = new InstMastConvexityConstr(id,
                this,
                masterPtr,
                *spcPt,
                name,
                *(*spcPt)->lowerBoundPtr(),
                'G',
                'X', /// this is view as a subproblem constraint
                implicitOrExplicitType, //'E',
                's',
                -1,
                (*spcPt)->defaultDualVal4LbConstr(),
                BapcodInfinity,
                0);
        if (printL(5))
          std::cout << " new lowerbound convexity constr  " << *mccPtr << std::endl;

        (*spcPt)->lowerBoundMastConstrPtr(mccPtr);

        if ((param().mastInitMode().status() == MasterInitMode::localAndGlobAc)
                || (param().mastInitMode().status() == MasterInitMode::localArtCol)
                || (param().mastInitMode().status() == MasterInitMode::incSolColAndLac))

        {
          mccPtr->addLocalArtVar(modelPtr->objectiveSense());
        }
        if (masterPtr->probPtr()->posGlobalArtVarPtr() != NULL)
          mccPtr->addMember(masterPtr->probPtr()->posGlobalArtVarPtr());
        if (masterPtr->probPtr()->negGlobalArtVarPtr() != NULL)
          mccPtr->addMember(masterPtr->probPtr()->negGlobalArtVarPtr());
    }
  }

  return;
}

ConvexityGenConstr::~ConvexityGenConstr()
{
  return;
}

bool ConvexityGenConstr::genericCount(const InstanciatedConstr * const iconstrPtr,
                                      const InstanciatedVar * const ivarPtr) const
{
  bapcodInit().check(true, "ConvexityGenConstr::genericCount should not be called",
                     ProgStatus::run);

  return false;
}

const LpCoef ConvexityGenConstr::genericCoef(const InstanciatedConstr * const iconstrPtr,
                                             const InstanciatedVar * const ivarPtr) const
{
  bapcodInit().check(true, "ConvexityGenConstr::genericCoef should not be called",
                     ProgStatus::run);

  return LpCoef::ZeroCoef;
}

std::ostream& ConvexityGenConstr::print(std::ostream& os) const
{
  os << "ConvexityGenConstr" << std::endl;

  return (GenericConstr::print(os));
}

