/**
 *
 * This file bcModelNonPublicCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

//
//  bcCapacityCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 11/09/2017.
//
//

#include "bcModelNonPublicCuts.hpp"
#include "bcNonPublicCuts.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

const int MaxOldCutsMemory = 1000;

BcCliqueCutConstrArray::BcCliqueCutConstrArray(const BcFormulation & formulation,
                                               const double & rootPriorityLevel,
                                               const double & nonRootPriorityLevel):
        BcCutConstrArray()
{
    if (printL(5))
        std::cout << " BcCliqueCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcCliqueCutConstrArray =  CLQ" << std::endl;

    const ControlParameters & params = formulation.probConfPtr()->param();
    if ((params.RCSPcliqueCutsMaxNumPerRound() <= 0) || (params.RCSPcliqueCutsMaxSepTime() <= 0))
        return;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("CLQ");

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << " BcCliqueCutConstrArray() : need to create cut  " << std::endl;

#if defined(CLIQUE_SEP_IS_FOUND) && defined(BCP_RCSP_IS_FOUND)
        _genericCutConstrPtr = new GenericCliqueCutConstr(formulation.probConfPtr()->modelPtr(),
                                                        formulation.probConfPtr(), "CLQ",
                                                        nonRootPriorityLevel, rootPriorityLevel,
                                                        params.RCSPcliqueCutsMaxNumPerRound(),
                                                        params.RCSPcliqueCutsMaxSepTime());
#else
        std::cerr << "BaPCod error : cannot use clique cuts, as CliqueSep or BCP_RCSP libraries are not found."
                  << std::endl;
        exit(1);
#endif
        _genericCutConstrPtr->defaultSense('L');
        _genericCutConstrPtr->defaultCostRhs(1.0);
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcCliqueCutConstrArray::~BcCliqueCutConstrArray()
{
}

BcStrongKPathCutConstrArray::BcStrongKPathCutConstrArray(const BcFormulation & formulation,
                                                         const int & maxCapacity,
                                                         const std::vector<int> & demands,
                                                         const bool & isFacultative,
                                                         const bool & equalityCase,
                                                         const int & twoPathCutsResId,
                                                         const double & rootPriorityLevel,
                                                         const double & nonRootPriorityLevel)
{
    if (printL(5))
        std::cout << " BcStrongKPathCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcStrongKPathCutConstrArray = SKP" << std::endl;

    const ControlParameters & params = formulation.probConfPtr()->param();
    if (!params.RCSPuseCapacityCuts())
        return;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("SKP");

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << "BcStrongKPathCutConstrArray() : need to create cut" << std::endl;

#ifdef BCP_RCSP_IS_FOUND
        _genericCutConstrPtr = new GenericLimMemStrongKPathCutConstr(formulation.probConfPtr()->modelPtr(),
                                                                     formulation.probConfPtr(), "SKP",
                                                                     nonRootPriorityLevel, rootPriorityLevel,
                                                                     isFacultative, equalityCase, maxCapacity, demands,
                                                                     twoPathCutsResId);
#else
        std::cerr << "BaPCod error : cannot use strong k-path cuts, as BCP_RCSP library is not found." << std::endl;
        exit(1);
#endif

        _genericCutConstrPtr->defaultSense('G');
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcStrongKPathCutConstrArray::~BcStrongKPathCutConstrArray()
{
}

#ifdef BCP_RCSP_IS_FOUND

void getActiveStrongKPathCuts(const BcFormulation & spPtr,
                              std::vector<std::pair<const bcp_rcsp::StrongKPathCut *, double> > & kPathCuts)
{
    if (spPtr.probConfPtr() == NULL)
    {
        std::cerr << "ERROR Model BcFormulation == NULL in getActiveStrongKPathCuts" << std::endl;
        exit(1);
    }
    kPathCuts.clear();

    MasterConf * mastConfPtr = spPtr.probConfPtr()->modelPtr()->master();
    GenericCutConstr * genKPathCutConstrPtr = mastConfPtr->getGenericCutConstr("SKP");
    if (genKPathCutConstrPtr == NULL)
        return;

    long long int scaleFactor = spPtr.probConfPtr()->param().SafeDualBoundScaleFactor();

    const IndexCell2InstancConstrPtrMap & constrPtrMap = genKPathCutConstrPtr->indexCell2InstancConstrPtrMap();
    for (IndexCell2InstancConstrPtrMap::const_iterator it = constrPtrMap.begin(); it != constrPtrMap.end(); ++it)
        if ((it->second->vcIndexStatus() == VcIndexStatus::Active)
            && it->second->isTypeOf(VcId::LimMemoryKPathCutConstrMask))
        {
            LimMemKPathCut * cutPtr = static_cast<LimMemKPathCut *>(it->second);
            double curDualValue = (scaleFactor > 0) ? (double)ceil(cutPtr->valOrSepPointVal()._val * scaleFactor)
                                                    : (double)cutPtr->valOrSepPointVal();
            kPathCuts.push_back(std::make_pair(cutPtr->rcspCutPtr(), curDualValue));
        }
}

#endif

BcHomExtCapCutConstrArray::BcHomExtCapCutConstrArray(const BcFormulation & formulation,
                                                     const std::vector<int> & demands,
                                                     const std::vector<int> & capacities,
                                                     const int & resourceId,
                                                     const double & priorityLevel,
                                                     const int & maxNumCutsPerRound):
        BcCutConstrArray()
{
    if (printL(5))
        std::cout << " BcHomExtCapCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcHomExtCapCutConstrArray = HECC" << std::endl;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("HECC");

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << "BcHomExtCapCutConstrArray() : need to create cut" << std::endl;

#if defined(RHECCSEP_IS_FOUND) && defined(BCP_RCSP_IS_FOUND)
        _genericCutConstrPtr = new GenericHomExtCapCutConstr(formulation.probConfPtr()->modelPtr(),
                                                           formulation.probConfPtr(), "HECC",
                                                           priorityLevel, priorityLevel, maxNumCutsPerRound,
                                                           demands, capacities, resourceId);
#else
        std::cerr << "BaPCod error : cannot use extended capacity cuts, as RHECC_Sep or BCP_RCSP libraries are not found."
                  << std::endl;
        exit(1);
#endif
        _genericCutConstrPtr->defaultSense('G');
        _genericCutConstrPtr->defaultCostRhs(1.0);
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);
    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcHomExtCapCutConstrArray::~BcHomExtCapCutConstrArray()
{
}

BcOverlEliminCutConstrArray::BcOverlEliminCutConstrArray(const BcFormulation & formulation,
                                                         const std::vector<int> & demands,
                                                         const std::vector<int> & capacities,
                                                         const int & resourceId,
                                                         const double & priorityLevel,
                                                         const int & maxNumCutsPerRound):
        BcCutConstrArray()
{
    if (printL(5))
        std::cout << " BcOverlEliminCutConstrArray() : ProbConfig =  " << formulation.probConfPtr()->name()
                  << " BcOverlEliminCutConstrArray = OEC" << std::endl;

    _genericCutConstrPtr = formulation.probConfPtr()->getGenericCutConstr("OEC");

    if (_genericCutConstrPtr == NULL)
    {
        if (printL(5))
            std::cout << "BcOverlEliminCutConstrArray() : need to create cut" << std::endl;

#if defined(RHECCSEP_IS_FOUND) && defined(BCP_RCSP_IS_FOUND)
        _genericCutConstrPtr = new GenericOverlEliminCutConstr(formulation.probConfPtr()->modelPtr(),
                                                             formulation.probConfPtr(), "OEC",
                                                             priorityLevel, priorityLevel, maxNumCutsPerRound,
                                                             demands, capacities, resourceId);
#else
        std::cerr << "BaPCod error : cannot use overload elimination cuts, as as RHECC_Sep or BCP_RCSP libraries are "
                  << "not found." << std::endl;
        exit(1);
#endif
        _genericCutConstrPtr->defaultSense('G');
        _genericCutConstrPtr->defaultCostRhs(1.0);
        _genericCutConstrPtr->defaultFlag('d');
        _genericCutConstrPtr->defaultVal(0);

    }
    _genericConstrPtr = _genericCutConstrPtr;
}

BcOverlEliminCutConstrArray::~BcOverlEliminCutConstrArray()
{
}

#ifdef BCP_RCSP_IS_FOUND

#ifdef CLIQUE_SEP_IS_FOUND

BcCliqueCut::BcCliqueCut(CliqueCut * cutPtr) :
        _cutPtr(cutPtr)
{
}

BcCliqueCut::BcCliqueCut(const BcCliqueCut & cut) :
    _cutPtr(cut._cutPtr)
{
}

double BcCliqueCut::curDualVal() const
{
    return _cutPtr->valOrSepPointVal();
}

const std::vector<std::vector<int> > & BcCliqueCut::setIds() const
{
  return _cutPtr->_setIds;
}

int BcCliqueCut::id() const
{
  return _cutPtr->id().first();
}


void BcCliqueCut::nicePrint() const
{
  _cutPtr->nicePrint();
}
#endif

#endif