/**
 *
 * This file bcNetworkBasedCuts.hpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#ifndef Project_bcNetworkBasedCuts_hpp
#define Project_bcNetworkBasedCuts_hpp

#include "bcGenVarConstrC.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcMastVarConstrC.hpp"

#include <boost/functional/hash.hpp>
#include <unordered_set>

#ifdef _WINDOWS
#include <bitset>
#endif

class InstanciatedVar;

#ifdef BCP_RCSP_IS_FOUND

namespace bcp_rcsp {
    struct RoundedCapacityCut;
    class RoundCapCutSeparationInterface;
    struct RouteLoadKnapsackCut;
    class RouteLoadKnapsackCutSeparationInterface;
}

class GenericRCSPCapacityCutConstr: public GenericCutConstr
{
    int _maxCapacity;
    bool _equalityCase;
    std::vector<int> _demands;
    int _twoPathCutsResId;
    bcp_rcsp::RoundCapCutSeparationInterface * _interfacePtr;
    std::vector<const ColGenSpConf *> _cgSpConfPts;
public:
    GenericRCSPCapacityCutConstr(Model * modelPtr,
                                 ProbConfig * probConfPtr,
                                 const std::string & name,
                                 const Double & nonRootPriorityLevel,
                                 const Double & rootPriorityLevel,
                                 const bool & isFacultative,
                                 const bool & equalityCase,
                                 const int & maxCapacity,
                                 const std::vector<int> & demands,
                                 const int & twoPathCutsResId);
    virtual ~GenericRCSPCapacityCutConstr();
    virtual bool prepareSeparation();
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset < InstanciatedConstr * ,
                                              CutSeparationPriorityComp > & generatedCutConstrSet);
};


class GenericLimMemRankOneCutConstr;

namespace bcp_rcsp {
    struct RankOneCut;
    class RankOneCutSeparationInterface;
}

class LimMemRankOneCut: public InstMasterConstr, Base4NonLinearConstraint
{
  friend class GenericLimMemRankOneCutConstr;

  const bcp_rcsp::RankOneCut * _rcspCutPtr;

  GenericLimMemRankOneCutConstr * _genLimMemRankOneCutConstr;
    
public:
  LimMemRankOneCut(const IndexCell& id,
                   GenericLimMemRankOneCutConstr * genConstrPtr,
                   ProbConfig * probConfigPtr,
                   const std::string & name,
                   const bcp_rcsp::RankOneCut * rcspCutPtr);
  virtual ~LimMemRankOneCut();

  const bcp_rcsp::RankOneCut * rcspCutPtr() const {return _rcspCutPtr;}

  const std::vector<int> & setIds() const;
  const std::vector<int> & coeffs() const;
  int numRows() const;
  int denominator() const;

  virtual void setMembership();
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
  void nicePrint(std::ostream& os = std::cout) const;
};

class GenericLimMemRankOneCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
    bcp_rcsp::RankOneCutSeparationInterface * _interfacePtr;

    int _nbPackingSets;
    int _memoryType;
    int _currentPhase;
    int _maxNumberOfRows;
    std::string _standaloneFileName;

    /// needed when the subproblem is solved by MIP
    std::string _spVarName;
    std::map<ColGenSpConf *, GenericVar *> _genIndicVarPts;
    std::map<ColGenSpConf *, GenericConstr *> _genIndicConstrPts;

    void updateSubprobemsWithIndicatorVarAndConstr(const bcp_rcsp::RankOneCut * cutPtr);

public:
    GenericLimMemRankOneCutConstr(Model * modelPtr,
                                  ProbConfig * probConfPtr,
                                  const std::string & name,
                                  const Double & nonRootPriorityLevel,
                                  const Double & rootPriorityLevel,
								  const std::string & spVarName,
                                  const int & memoryType,
                                  const bool & isFacultative);
    
    virtual ~GenericLimMemRankOneCutConstr();
    
    virtual bool prepareSeparation();
    virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
    virtual void buildMembership(InstanciatedConstr * iconstrPtr);
    void cutSeparationRoutine(const std::string & fileName);
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset < InstanciatedConstr * ,
                                                      CutSeparationPriorityComp > & generatedCutConstrSet);
    const std::string & spVarName() const {return _spVarName;}

    void nicePrintCut(const bcp_rcsp::RankOneCut * rcspCutPtr, std::ostream & os);

    virtual void resetSeparationPhase();
    virtual void increaseSeparationPhase();
    virtual bool separationPhaseIsMaximum() const;
	void setArcMemory();
    void setVertexMemory();
};

class GenericExtendedArcCutConstr;
struct BcCustomExtendedArcCutInfo;

class ExtendedArcCut: public InstMasterConstr, Base4NonLinearConstraint
{
    GenericExtendedArcCutConstr * _genExtendedArcCutConstr;
    const BcCustomExtendedArcCutInfo * _cutInfoPtr; /// this is needed for custom (user) extended arc cuts
public :
    ExtendedArcCut(GenericExtendedArcCutConstr * genConstrPtr,
                   ProbConfig * probConfigPtr,
                   const std::string & name,
                   const Double & rhs,
                   const char & sense,
                   const BcCustomExtendedArcCutInfo * cutInfoPtr = NULL);
    virtual ~ExtendedArcCut();

    virtual void nicePrint(std::ostream& os = std::cout);
    virtual double getArcCoefficient(const int & tailVertId, const int & headVertId,
                                     const double * tailResCons) const;
    virtual double getRouteCoefficient(const std::vector<int> & routeVertIds,
                                       const std::vector<std::vector<double> > & routeResCons) const;
    virtual double getArcCoefficient(const NetworkFlow * netFlowPtr, const int & arcId,
                                     const double * resCons, const bool & isTailResCons) const;
    virtual double getRouteCoefficient(const NetworkFlow * netFlowPtr, const std::vector<int> & arcIds,
                                       const std::vector<std::vector<double> > & routeResCons) const;
    virtual void setMembership();
    virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
    inline const BcCustomExtendedArcCutInfo * cutInfoPtr() const {return _cutInfoPtr;}
};

class BcCustomExtendedArcCutSeparationFunctor;

namespace bcp_rcsp
{
    struct GraphData;
    struct PreprocessedGraphInfo;
}

class GenericExtendedArcCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
    friend class ExtendedArcCut;
protected:
    static int _numGeneratedCuts; /// this member serves for setting unique id for every cut derived from ExtendedArcCut
    BcCustomExtendedArcCutSeparationFunctor * _separationFunctorPtr;
    int _maxGraphId;
    std::vector<const bcp_rcsp::GraphData *> _graphPts;
    std::vector<const bcp_rcsp::PreprocessedGraphInfo *> _graphInfoPts;
public:
    GenericExtendedArcCutConstr(Model * modelPtr,
                                ProbConfig * probConfPtr,
                                const std::string & name,
                                const Double & nonRootPriorityLevel,
                                const Double & rootPriorityLevel);
    ~GenericExtendedArcCutConstr();
    virtual bool prepareSeparation();
    void setSeparationFunctor(BcCustomExtendedArcCutSeparationFunctor * separationFunctorPtr);
    virtual void buildMembership(InstanciatedConstr * iconstrPtr);
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset<InstanciatedConstr *,
                                              CutSeparationPriorityComp> & generatedCutConstrSet);
    virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
    virtual const LpCoef getMastColumnCoeff(ExtendedArcCut * cutPtr, MastColumn * colPtr) const;
};

class GenericResConsKnapsackCutConstr;

typedef std::map<std::pair<ColGenSpConf *, int>, std::map<double, double> > ResConsKnapsackMap;

class ResConsKnapsackCut: public InstMasterConstr, Base4NonLinearConstraint
{
    friend class GenericResConsKnapsackCutConstr;

    const bcp_rcsp::RouteLoadKnapsackCut * _rcspCutPtr;

    GenericResConsKnapsackCutConstr * _genResConsKnapCutConstr;

public:
    ResConsKnapsackCut(const IndexCell& id,
                       GenericResConsKnapsackCutConstr * genConstrPtr,
                       ProbConfig * probConfigPtr,
                       const std::string & name,
                       const Double & rhs,
                       const bcp_rcsp::RouteLoadKnapsackCut * rcspCutPtr);

    bool isRelatedTo(const ColGenSpConf * cgSpConfPtr);

    const bcp_rcsp::RouteLoadKnapsackCut * rcspCutPtr() const {return _rcspCutPtr;}

    virtual ~ResConsKnapsackCut();

    virtual void setMembership();
    virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
    void nicePrint(std::ostream& os = std::cout) const;
};

typedef std::map<std::pair<ColGenSpConf *, int>, std::map<double, double> > ResConsKnapsackMap;

struct ResConsConstrLhsItem {
    ColGenSpConf * colGenSpConfPtr;
    int resId;
    InstanciatedVar * varPtr;
    double coeff;

    ResConsConstrLhsItem(ColGenSpConf * colGenSpConfPtr_, int resId_, InstanciatedVar * varPtr_, double coeff_) :
            colGenSpConfPtr(colGenSpConfPtr_), resId(resId_), varPtr(varPtr_), coeff(coeff_)
    {}
};

struct ResConsConstrInfo {
    int id;
    InstMasterConstr * constrPtr;
    std::vector<ResConsConstrLhsItem> leftHandSideInfo;
    std::map<InstanciatedVar *, double> pureMastVarInfo;
    double freeCoeff;
    Model * modelPtr;
    BapcodInit * bcInitPtr;

    ResConsConstrInfo(int id_, InstMasterConstr * constrPtr_) :
            id(id_), constrPtr(constrPtr_), leftHandSideInfo(), pureMastVarInfo(), freeCoeff(0.0), bcInitPtr(nullptr),
            modelPtr(nullptr)
    {}

    ~ResConsConstrInfo();
};


class GenericResConsKnapsackCutConstr: public GenericCutConstr, public Base4NonLinearGenericConstr
{
    std::vector<ResConsConstrInfo> _constrInfos; /// each element corresponds to one master constraint which
    /// contains variables associated to resources

    bcp_rcsp::RouteLoadKnapsackCutSeparationInterface * _interfacePtr;

    /// Needed for exact Fenchel separation

    bool _useFenchelSeparation;

    bool prepareFenchelSeparation();
public:
    GenericResConsKnapsackCutConstr(Model * modelPtr,
                                    ProbConfig * probConfPtr,
                                    const std::string & name,
                                    const Double & nonRootPriorityLevel,
                                    const Double & rootPriorityLevel);

    virtual ~GenericResConsKnapsackCutConstr();

    virtual bool prepareSeparation();
    virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;
    virtual void buildMembership(InstanciatedConstr * iconstrPtr);
    virtual void cutSeparationRoutine(const VarPtrSet & curSol,
                                      std::multiset < InstanciatedConstr * ,
                                              CutSeparationPriorityComp > & generatedCutConstrSet);
};

#endif /* BCP_RCSP_IS_FOUND */


#endif
