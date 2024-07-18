/**
 *
 * This file bcModelNonPublicCuts.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELNONPUBLICCUTS_H_
#define BCMODELNONPUBLICCUTS_H_

#include "bcModelCutConstrC.hpp"

#ifdef BCP_RCSP_IS_FOUND
#include "rcsp_interface.hpp"
#endif

class BcStrongKPathCutConstrArray: public BcCutConstrArray
{
public:
    BcStrongKPathCutConstrArray(const BcFormulation & formulation,
                                const int & maxCapacity,
                                const std::vector<int> & demands,
                                const bool & isFacultative = true,
                                const bool & equalityCase = true,
                                const int & twoPathCutsResId = -1,
                                const double & rootPriorityLevel = 1.0,
                                const double & nonRootPriorityLevel = 1.0);
    virtual ~BcStrongKPathCutConstrArray();
};

#ifdef BCP_RCSP_IS_FOUND
void getActiveStrongKPathCuts(const BcFormulation & spPtr,
                              std::vector<std::pair<const bcp_rcsp::StrongKPathCut *, double> > & kPathCuts);
#endif

class BcCliqueCutConstrArray: public BcCutConstrArray
{
public:
    BcCliqueCutConstrArray(const BcFormulation & formulation,
                           const double & rootPriorityLevel = 1.0,
                           const double & nonRootPriorityLevel = 1.0);
    virtual ~BcCliqueCutConstrArray();
};

class BcHomExtCapCutConstrArray: public BcCutConstrArray
{
public:
  BcHomExtCapCutConstrArray(const BcFormulation & formulation,
                            const std::vector<int> & demands,
                            const std::vector<int> & capacities,
                            const int & resourceId,
                            const double & priorityLevel = 1.0,
                            const int & maxNumCutsPerRound = 1000000);
  virtual ~BcHomExtCapCutConstrArray();
};

class BcOverlEliminCutConstrArray: public BcCutConstrArray
{
public:
  BcOverlEliminCutConstrArray(const BcFormulation & formulation,
                              const std::vector<int> & demands,
                              const std::vector<int> & capacities,
                              const int & resourceId,
                              const double & priorityLevel = 1.0,
                              const int & maxNumCutsPerRound = 1000000);
  virtual ~BcOverlEliminCutConstrArray();
};

#ifdef BCP_RCSP_IS_FOUND
#ifdef CLIQUE_SEP_IS_FOUND

class BcCliqueCut : public bcp_rcsp::CliqueCutInterface
{
protected:
    CliqueCut * _cutPtr;
public:
    BcCliqueCut(CliqueCut * cutPtr);
    BcCliqueCut(const BcCliqueCut & cut);
    double curDualVal() const;
    virtual const std::vector<std::vector<int> > & setIds() const;
    virtual int id() const;
    virtual void nicePrint() const;
};

void getActiveCliqueCuts(const BcFormulation & spPtr, std::vector<const BcCliqueCut *> & cutPts);

#endif

#endif

#ifdef _CGL_FOUND

struct BcCutType {

    enum Cuts { // TODO faire comme l'énumération vcType
        cglClique
        ,
        cglKnapsackCoverCut
        ,
        cglOddHoleCut
        ,
        cglZeroHalf
        ,
        cglFlowCover
        ,
        cglMixedIntegerRoundingCut
        ,
        cglMixedIntegerRounding2Cut
        ,
        cglTwoMir
        ,
        cglResidualCapacity
        ,
        cglSimpleRoundingCut
        ,
        _BcCutTypeMax // end marker
    };
};

class CGLSeparationFunctor : public BcCutSeparationFunctor
{
public:

    CGLSeparationFunctor();

    virtual ~CGLSeparationFunctor() {}

    virtual int operator()(BcFormulation formPtr,
                           BcSolution & primalSol,
                           double & maxViolation,
                           std::list< BcConstr > & cutList);

    void enableCut(BcCutType::Cuts cutType);
    void disableCut(BcCutType::Cuts cutType);
    void enableCutAll();
    void disableCutAll();

private:

    int last_cut_index;

    std::vector<bool> enabledCuts;
};

#endif // _CGL_FOUND

#endif //BCMODELNONPUBLICCUTS_H_

