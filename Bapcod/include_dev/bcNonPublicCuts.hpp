/**
 *
 * This file bcNonPublicCuts.hpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcNetworkBasedCutsC.hpp
//  Project
//
//  Created by Ruslan Sadykov on 6/10/2017.
//
//

#ifndef Project_bcNonPublicCuts_hpp
#define Project_bcNonPublicCuts_hpp

#include "bcNetworkBasedCuts.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcMastVarConstrC.hpp"

#include <boost/functional/hash.hpp>
#include <unordered_set>

#ifdef _WINDOWS
#include <bitset>
#endif

class InstanciatedVar;

#ifdef CVRPSEP_IS_FOUND
#include "cnstrmgr.h"

class GenericCapacityCutConstr: public GenericCutConstr
{
    int _maxCapacity;
    std::vector<int> _demands;
    CnstrMgrPointer _myCutsCMP;
    CnstrMgrPointer _myOldCutsCMP;
    int _cutCount;
    int _cutRound;
    int _maxNumCutsPerRound;
    int _numElemSets;
    std::vector<std::vector<std::vector<InstanciatedVar *> > > _elemSetPairVarPts;
public:
    GenericCapacityCutConstr(Model * modelPtr,
                             ProbConfig * probConfPtr,
                             const std::string & name,
                             const Double & nonRootPriorityLevel,
                             const Double & rootPriorityLevel,
                             const bool & isFacultative,
                             const int & maxNumCutsPerRound,
                             const int & maxCapacity,
                             const std::vector<int> & demands);
    virtual ~GenericCapacityCutConstr();
    virtual bool prepareSeparation();
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset < InstanciatedConstr * ,
                                              CutSeparationPriorityComp > & generatedCutConstrSet);
};

#endif

#ifdef BCP_RCSP_IS_FOUND

namespace bcp_rcsp {
    struct StrongKPathCut;
    class StrongKPathCutSeparationInterface;
    struct RouteLoadKnapsackCut;
    class RouteLoadKnapsackCutSeparationInterface;
}

class GenericLimMemStrongKPathCutConstr;

class LimMemKPathCut: public InstMasterConstr, Base4NonLinearConstraint
{
    friend class GenericLimMemRankOneCutConstr;

    const bcp_rcsp::StrongKPathCut * _rcspCutPtr;

    GenericLimMemStrongKPathCutConstr * _genLimMemStrongKPathCutConstr;

public:
    LimMemKPathCut(const IndexCell& id,
                   GenericLimMemStrongKPathCutConstr * genConstrPtr,
                   ProbConfig * probConfigPtr,
                   const std::string & name,
                   const bcp_rcsp::StrongKPathCut * rcspCutPtr);
    virtual ~LimMemKPathCut();

    const bcp_rcsp::StrongKPathCut * rcspCutPtr() const {return _rcspCutPtr;}

    const std::vector<int> & setIds() const;

    virtual void setMembership();
    virtual bool isTypeOf(const VcId::VcIdentifier & vcIdentifier) const;
    void nicePrint(std::ostream& os) const;
};

class GenericLimMemStrongKPathCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
    bool _equalityCase;
    int _maxCapacity;
    std::vector<int> _demands;
    int _twoPathCutsResId;
    bcp_rcsp::StrongKPathCutSeparationInterface * _interfacePtr;
    std::vector<const ColGenSpConf *> _cgSpConfPts;
public:
    GenericLimMemStrongKPathCutConstr(Model * modelPtr,
                                      ProbConfig * probConfPtr,
                                      const std::string & name,
                                      const Double & nonRootPriorityLevel,
                                      const Double & rootPriorityLevel,
                                      const bool & isFacultative,
                                      const bool & equalityCase,
                                      const int & maxCapacity,
                                      const std::vector<int> & demands,
                                      const int & twoPathCutsResId);
    virtual ~GenericLimMemStrongKPathCutConstr();
    virtual bool prepareSeparation();
    virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
    virtual void buildMembership(InstanciatedConstr * iconstrPtr);
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset < InstanciatedConstr * ,
                                              CutSeparationPriorityComp > & generatedCutConstrSet);
    void nicePrintCut(const bcp_rcsp::StrongKPathCut * rcspCutPtr, std::ostream & os);
};


#ifdef CLIQUE_SEP_IS_FOUND

class GenericCliqueCutConstr;

class CliqueCut: public InstMasterConstr, Base4NonLinearConstraint
{
	friend class GenericCliqueCutConstr;
	friend class BcCliqueCut;
	std::vector<std::vector<int> > _setIds;
	GenericCliqueCutConstr * _genCliqueCutConstr;

public:
	CliqueCut(const IndexCell& id, GenericCliqueCutConstr * genConstrPtr, ProbConfig * probConfigPtr,
			  const std::string & name, const std::set<std::vector<int> > & setIds);
	virtual ~CliqueCut();

	virtual void setMembership();
	virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
	void nicePrint(std::ostream & os = std::cout) const;
};

class GenericCliqueCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
	int _numGeneratedCuts;
	int _maxNumCutsPerRound;
	double _separationMaxTime;
	int _nbPackingSets;

	void printColsToConstrsMatrix(int spConfNum, std::vector<std::vector<int> > & colsToConstrsMatrix);
	void minimizeColumnToSetIdMatrix(std::vector<std::vector<bool> > & columnToSetIdMatrix,
									 std::set<std::vector<int> > & cutSetIds);
public:
	GenericCliqueCutConstr(Model * modelPtr,
						   ProbConfig * probConfPtr,
						   const std::string & name,
						   const Double & nonRootPriorityLevel,
						   const Double & rootPriorityLevel,
						   const int & maxNumCutsPerRound,
	                       const double & separationMaxTime);
	virtual ~GenericCliqueCutConstr();
	virtual bool prepareSeparation();
	virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
	virtual void buildMembership(InstanciatedConstr * iconstrPtr);
	virtual void cutSeparationRoutine(const VarPtrSet & curSol,
									  std::multiset < InstanciatedConstr * ,
											  CutSeparationPriorityComp > & generatedCutConstrSet);
	const LpCoef getMastColumnCoeff(CliqueCut * cutPtr, MastColumn * colPtr) const;
};

#endif /* CLIQUE_SEP_IS_FOUND */

class CutGenerator;

namespace arturnewsep
{
	class CutGenerator2;
}


class HomogExtendedCapacityCut: public ExtendedArcCut
{
  friend class GenericHomExtCapCutConstr;
  bool _useNewSep;
  int _numerator;
  int _denominator;
  int _resourceId;
  int _sumDemands;
  std::set<int> _verticesSet;
  std::vector<bool> _verticesBitmap;

public:
  HomogExtendedCapacityCut(GenericExtendedArcCutConstr * genConstrPtr,
                           ProbConfig * probConfigPtr,
                           const std::string & name,
                           const Double & rhs,
                           const int & numerator,
                           const int & denominator,
                           const int & resourceId,
                           const std::set<int> verticesSet,
                           const int & sumDemands);
  virtual ~HomogExtendedCapacityCut();

  virtual void nicePrint(std::ostream& os = std::cout);
  virtual double getArcCoefficient(const int & tailVertId, const int & headVertId,
                                  const double * tailResCons) const;
  virtual double getRouteCoefficient(const std::vector<int> & routeVertIds,
                                    const std::vector<std::vector<double> > & routeResCons) const;
};

class GenericHomExtCapCutConstr: public GenericExtendedArcCutConstr
{
  bool _useNewSep;
  int _maxNumCutsPerRound;
  int _resourceId;
  std::vector<int> _demands; /// demands (resource consumption) of vertices
  std::vector<int> _capacities; /// capacities of different subproblems
  CutGenerator * cutGenObj;
  arturnewsep::CutGenerator2 * cutGenObjNew;

public:
  GenericHomExtCapCutConstr(Model * modelPtr,
                            ProbConfig * probConfPtr,
                            const std::string & name,
                            const Double & nonRootPriorityLevel,
                            const Double & rootPriorityLevel,
                            const int & maxNumCutsPerRound,
                            const std::vector<int> & demands,
                            const std::vector<int> & capacities,
                            const int & resourceId);

  virtual bool prepareSeparation();

  virtual ~GenericHomExtCapCutConstr();
  virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                    std::multiset < InstanciatedConstr * ,
                                                    CutSeparationPriorityComp > & generatedCutConstrSet);
};

class OverloadEliminationCut: public ExtendedArcCut
{
  friend class GenericOverlEliminCutConstr;
  int _resourceId;
  int _resConsumption;
  int _sumDemands;
  int _capacity;
  int _nbMachines;
  std::set<int> _verticesSet;
  std::vector<bool> _verticesBitmap;
    
public:
  OverloadEliminationCut(GenericExtendedArcCutConstr * genConstrPtr,
                         ProbConfig * probConfigPtr,
                         const std::string & name,
                         const Double & rhs,
                         const int & resourceId,
                         const int & resConsumption,
                         const std::set<int> verticesSet,
                         const int & sumDemands,
						 const int & capacity,
						 const int & nbMachines);
  virtual ~OverloadEliminationCut();

  virtual void nicePrint(std::ostream& os = std::cout);
  virtual double getArcCoefficient(const int & tailVertId, const int & headVertId,
                                   const double * tailResCons) const;
  virtual double getRouteCoefficient(const std::vector<int> & routeVertIds,
                                     const std::vector<std::vector<double> > & routeResCons) const;
};

class GenericOverlEliminCutConstr: public GenericExtendedArcCutConstr
{
  int _maxNumCutsPerRound;
  int _resourceId;
  std::vector<int> _demands; /// demands (resource consumption) of vertices
  std::vector<int> _capacities; /// capacities of different subproblems
  arturnewsep::CutGenerator2 * cutGenObj;

public:
  GenericOverlEliminCutConstr(Model * modelPtr,
                              ProbConfig * probConfPtr,
                              const std::string & name,
                              const Double & nonRootPriorityLevel,
                              const Double & rootPriorityLevel,
                              const int & maxNumCutsPerRound,
                              const std::vector<int> demands,
                              const std::vector<int> capacities,
                              const int & resourceId);
    
  virtual ~GenericOverlEliminCutConstr();
  virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                    std::multiset < InstanciatedConstr * ,
                                                    CutSeparationPriorityComp > & generatedCutConstrSet);
};

#endif /* BCP_RCSP_IS_FOUND */

#endif
