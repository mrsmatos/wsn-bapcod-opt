/**
 *
 * This file bcCInterfaceModel.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCCINTERFACEMODEL_HPP
#define BCCINTERFACEMODEL_HPP

#include "bcModelingLanguageC.hpp"
#include "bcApplicationException.hpp"
#include "bcSolutionC.hpp"
#include "bcProbConfigC.hpp"
#include "bcNetworkFlowC.hpp"

#include "bcModelGlobalCustomSolvers.hpp"

#include <limits>
#include <string>
#include <vector>

//#define NOSTANDALONE

#pragma GCC visibility push(default)

#ifdef _MSC_VER
#define EXPORTED  __declspec( dllexport )
#pragma warning( push )
#pragma warning( disable : 4190 )
#else
#define EXPORTED
#endif

#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type-c-linkage"
#endif

extern "C" {

    struct EXPORTED CutCallback;

    struct EXPORTED InterfaceModel {
        int id; // to debug

        bool problemInit;
        bool decompositionUsed;

        BcInitialisation init;

        BcModel model;
        BcObjectiveArray objective;

        BcFormulationArray mpform;
        BcMaster master;
        BcColGenSpArray colGenSp;

        std::vector<BcVar> vars;
        std::vector<BcConstr> constrs;
        std::vector<BcBranchingConstrArray> branchingConstrs;

        std::map<int, char *> genvars; // colId to name of the generic var.

        BcSolution solution; // to keep the first ptr
        
        BcSolution lastInitialSolution;
        BcSolution initialSolution;

        std::map<std::string, std::vector<int> > varspartition; // key is subproblem name
        std::vector<CutCallback> cutCbList;
        std::vector< std::pair<int*,int> > subproblems; //int* is the id and int is the type

        InterfaceModel(char* paramFile, bool printParam, bool intObj, bool intBound, int argc, char** argv);
        ~InterfaceModel();
    };

    struct EXPORTED CutCallback {
        int id;
        int nbcuts;
        char type;
        BcCutConstrArray cuts;

        // typecut = 'C' for core (lazy) & typecut = 'F' for facultative (user)
        CutCallback(InterfaceModel* m, int id, char typecut, std::string name);
        CutCallback(const CutCallback & that);
        ~CutCallback();
    };

    struct EXPORTED bcOracleSolution {
        BcSolution _solution;
        double _objVal;
        BcSolution _lastSolution;
        BcFormulation _sp;
        int _phaseOfStageApproach;

        bcOracleSolution(BcSolution & sol, BcFormulation & sp, int phase);
        ~bcOracleSolution();
    };

}

class EXPORTED JuliaSolverOracleFunctor : public BcSolverOracleFunctor {
private:
    void* _userdata;
    int _idoracle;
    //Pointeur sur la fonction julia
    void (*_julia_fct)(const void*, const void*, const int);

public:
    JuliaSolverOracleFunctor(void (*julia_fct)(const void*, const void*, const int),
                                      void* userdata, int idoracle);

    ~JuliaSolverOracleFunctor();

#define IN_ // Indicates what is used for input
#define OUT_ // Indicates what is used for output
    virtual bool operator() (IN_ BcFormulation spPtr,
                             OUT_ double & objVal,
                             /// it represents the value of the returned solution (either primal or dual or both)
                             IN_ OUT_ double & primalBound,
                             /// on entering the procedure, primalBound represents a cutoff value (above which the resulting column cannot have negative reduced cost, for instance);
                             /// on exsiting the procedure primal bound should be reflect the value of the best feasible solution found
                             IN_ OUT_ double & dualBound,
                             /// on entering the procedure, dualBound represents an optimistic objVal evalution (for instance,the best possible value of objval ignoring the subproblem constraints);
                             /// on exsiting the procedure primal bound should be reflect the value of the best dual bound found.
                             OUT_ BcSolution & primalSol,
                             /// If (objVal == primalBound == dualbound) on exiting the procedure, it means that the problem is solved to optimality
                             /// If the returned primalBound and dualBound are different, one shall return either primalSol or dualSol but not both. In this case:
                             /// - dualSol will be considered as the relevent returned solution if objVal == dualBound
                             /// - primalSol will be considered as the relevent returned solution if objVal == primalBound
                             IN_ OUT_ int & phaseOfStageApproach);
};

class EXPORTED JuliaSeparationRoutineFunctor : public BcCutSeparationFunctor {
private:
    void* _userdata;
    int _idfunctor;
    void (*_julia_fct)(BcSolution&, const void*, std::list<BcConstr>&, int);

public:
    JuliaSeparationRoutineFunctor(void(*julia_fct)(BcSolution&, const void*, std::list<BcConstr>&, int),
                                  void* userdata, int idfunctor);
    ~JuliaSeparationRoutineFunctor();

    virtual int operator() (BcFormulation formPtr,
                            BcSolution & primalSol,
                            double & maxViolation,
                            std::list<BcConstr>& cutList);
};

typedef std::list<std::pair<BcConstr, double> > CstrsMembership;

class JuliaDefineCstrFunctor : public BcAddConstrFunctor {
private:
    BcConstrArray _cstrs;
    void* _userdata;
    void (*_julia_fct)(int*, void*, BcConstr&);

public:
    JuliaDefineCstrFunctor(BcConstrArray& cstrs, void* userdata, void(*julia_fct)(int*, void*, BcConstr&));
    virtual void operator() (const MultiIndex & indexArray);
};


class EXPORTED JuliaInformationalCallback : public BcSolutionFoundCallback {
private:
    void* _userdata;
    void (*_julia_fct)(BcSolution&, const void*);
    
public :
    JuliaInformationalCallback(void (*julia_fct)(BcSolution&, const void*), void* userdata);
    virtual void operator() (BcSolution& newSolution);
};

extern "C" EXPORTED MultiIndex arrayToMultiIndex(int* array);

extern "C" EXPORTED const BcFormulation& getProblem(InterfaceModel* m, int type, int* spId);

extern "C" EXPORTED const std::string getProblemName(const BcFormulation & f);

extern "C" EXPORTED InterfaceModel* bcInterfaceModel_new(char* paramFile, bool printParam, bool intObj, bool intBound,
                                                         int argc, char** argv);

extern "C" EXPORTED void bcInterfaceModel_delete(InterfaceModel* m);

extern "C" EXPORTED void bcInterfaceModel_initModel(InterfaceModel* m, int nbRows, int nbCols);

extern "C" EXPORTED int bcInterfaceModel_subProblemMult(InterfaceModel* m, int l, int u, int spType, int* spIndex);

extern "C" EXPORTED int bcInterfaceModel_registerSubProblem(InterfaceModel* m, int type, int* spid);

extern "C" EXPORTED int bcInterfaceModel_registerSubProblem(InterfaceModel* m, int type, int* spid);

extern "C" EXPORTED int bcInterfaceModel_registerVar(InterfaceModel* m, char* name,
                                                     int columnId, int spType, int* spId, int* varId);

extern "C" EXPORTED int bcInterfaceModel_registerGenericVar(InterfaceModel* m, char* name, int columnId);

extern "C" EXPORTED int bcInterfaceModel_addDynVar(InterfaceModel* m, int columnId, char* name, int* varId, int spType,
                                                   int* spId);

extern "C" EXPORTED int bcInterfaceModel_addMembershipToCstr(InterfaceModel* m, CstrsMembership& list, int rowId,
                                                             double coeff);

extern "C" EXPORTED void bcInterfaceModel_initVars(InterfaceModel* m, double* l, double* u, double* c);

extern "C" EXPORTED int bcInterfaceModel_getVarLb(InterfaceModel* m, double* lb, int size);

extern "C" EXPORTED int bcInterfaceModel_setVarLb(InterfaceModel* m, double* lb, int size);

extern "C" EXPORTED int bcInterfaceModel_getVarUb(InterfaceModel* m, double* ub, int size);

extern "C" EXPORTED int bcInterfaceModel_setVarUb(InterfaceModel* m, double* ub, int size);

extern "C" EXPORTED int bcInterfaceModel_getVarCosts(InterfaceModel* m, double* c, int size);

extern "C" EXPORTED int bcInterfaceModel_setVarCosts(InterfaceModel* m, double* c, int size);

extern "C" EXPORTED int bcInterfaceModel_getVarType(InterfaceModel* m, char* t, int size);

extern "C" EXPORTED int bcInterfaceModel_setVarType(InterfaceModel* m, char* t, int size);

extern "C" EXPORTED int bcInterfaceModel_setVarPriorityInMaster(InterfaceModel* m, char* varName,
                                                                int spType, int* spId, double priority);

extern "C" EXPORTED int bcInterfaceModel_setVarPriorityInMaster(InterfaceModel* m, char* varName,
                                                                int spType, int* spId, double priority);

extern "C" EXPORTED int bcInterfaceModel_setVarPriorityInSp(InterfaceModel* m, char* varName, int spType, int* spId,
                                                            double priority);

extern "C" EXPORTED int bcInterfaceModel_addAggrSubProbVarBranching(InterfaceModel* m, char* varName,
                                                                    double highestPriority, double priority,
                                                                    int ignoredIndices, bool toBeUsedInPreprocessing);

extern "C" EXPORTED int bcInterfaceModel_registerBranchingExpression(InterfaceModel* m, char* expName, double priority);

extern "C" EXPORTED int bcInterfaceModel_addBranchingExpression(InterfaceModel* m, int arrayId, int* expId, int* varIds,
                                                       double* coeffs, int length);

extern "C" EXPORTED int bcInterfaceModel_registerCstr(InterfaceModel* m, char* name,
                                             int rowId, int spType, int* spId, int* constrId);

extern "C" EXPORTED int bcInterfaceModel_registerGenericCstr(InterfaceModel* m, char* name, int rowId);

extern "C" EXPORTED int bcInterfaceModel_setConstrBcType(BcConstr* c, char type);

extern "C" EXPORTED int bcInterfaceModel_attachCstrFunctor(InterfaceModel* m, char* name,
                                                           void (*julia_fct)(int*, void *, BcConstr&), void* jcstr);

extern "C" EXPORTED int bcInterfaceModel_setConstrArrayBcType(BcConstrArray* c, char type);

extern "C" EXPORTED int bcInterfaceModel_setConstrBcType(BcConstr* c, char type);

extern "C" EXPORTED int bcInterfaceModel_setConstrArrayBcType(BcConstrArray* c, char type);

extern "C" EXPORTED int bcInterfaceModel_addDynCstr(InterfaceModel* m, int rowId,
                                           char* name, int* cstrId);

extern "C" EXPORTED BcConstr& bcInterfaceModel_getCstr(InterfaceModel* m, int rowId);

extern "C" EXPORTED int bcInterfaceModel_addCstrTerm(InterfaceModel* m, BcConstr& cstr,
                                            int varcol, double coeff);

extern "C" EXPORTED int bcInterfaceModel_addCstrTerms(InterfaceModel* m, BcConstr& cstr, int* varcols, int len,
                                                      double coeff);

extern "C" EXPORTED void bcInterfaceModel_addCstrRhs(BcConstr& cstr, char sense,
                                            double rhs);

extern "C" EXPORTED int bcInterfaceModel_initCstrs(InterfaceModel* m, int* colPtr,
                                          int* rowId, double* nonZeroVal, double* lb, double* ub);

extern "C" EXPORTED int bcInterfaceModel_cstrUsedInPreprocessing(InterfaceModel* m,
                                                        char* cstrName, int spType, int* spId, bool used);

extern "C" EXPORTED int bcInterfaceModel_cstrUsedInPreprocessing(InterfaceModel* m, char* cstrName, int spType,
                                                                 int* spId, bool used);

extern "C" EXPORTED void bcInterfaceModel_setArtCostValue(InterfaceModel* m, double artcostvalue);

extern "C" EXPORTED void bcInterfaceModel_setObjMagnitude(InterfaceModel* m, double magnitude);

extern "C" EXPORTED void bcInterfaceModel_setObjLb(InterfaceModel* m, double lb);

extern "C" EXPORTED void bcInterfaceModel_setObjUb(InterfaceModel* m, double ub);


extern "C" EXPORTED int bcInterfaceModel_getSense(InterfaceModel* m);

extern "C" EXPORTED int bcInterfaceModel_setSubproblemPriority(InterfaceModel* m, int* blockgroup, double priority);

//To add in a new file (Solver Part)

extern "C" EXPORTED void bcSolution_enumerateAllColumns(InterfaceModel* m, BcSolution* solution,
                                                        int & nbEnumeratedCols);

extern "C" EXPORTED void bcInterfaceSolve_optimize(InterfaceModel* m, BcSolution* solution);

extern "C" EXPORTED void bcInterfaceSolve_setParameter(char* param, char* value);

extern "C" EXPORTED int bcInterfaceSolve_initOracle(InterfaceModel* m,
                                                    void (*julia_func)(const void*, const void*, const int),
                                                    int sp_type, int* sp_id, void* userdata, int idoracle);
extern "C" EXPORTED int bcInterfaceSolve_getMasterDual(InterfaceModel* m, int constr, double& dual) ;

extern "C" EXPORTED int bcInterfaceSolve_getVarCurUB(InterfaceModel* m, int var, double& curUB) ;

extern "C" EXPORTED int bcInterfaceSolve_getVarCurLB(InterfaceModel* m, int var, double& curLB);

extern "C" EXPORTED int bcInterfaceSolve_getOptimalityGapTolerance(InterfaceModel* m, double& ogt);

extern "C" EXPORTED int bcInterfaceSolve_getVarCurCost(InterfaceModel* m, int var, double& curCost);

extern "C" EXPORTED int bcInterfaceSolve_getDynVarCurCost(InterfaceModel* m, char* varname, int* varid, int sptype,
                                                          int* spid, double& curCost);

extern "C" EXPORTED int bcInterfaceSolve_getDynVarCurCost(InterfaceModel* m, char* varname, int* varid, int sptype,
                                                          int* spid, double& curCost);

extern "C" EXPORTED int bcInterfaceSolve_addToOracleSol(InterfaceModel* m, bcOracleSolution* s, int varId,
                                                        double varVal);

extern "C" EXPORTED int bcInterfaceSolve_addToOracleDualSol(InterfaceModel* m, bcOracleSolution* s, BcConstr* c);

extern "C" EXPORTED int bcInterfaceSolve_newOracleSol(bcOracleSolution* s);

extern "C" EXPORTED int bcInterfaceSolve_getOraclePhaseOfStage(bcOracleSolution* s);

extern "C" EXPORTED int bcInterfaceSolve_setObjValOfOracleSol(bcOracleSolution* s, double objval);

extern "C" EXPORTED double bcInterfaceSolve_getStatisticValue(InterfaceModel* m, char* name);
extern "C" EXPORTED long bcInterfaceSolve_getStatisticTime(InterfaceModel* m, char* name);
extern "C" EXPORTED long bcInterfaceSolve_getStatisticCounter(InterfaceModel* m, char* name);

extern "C" EXPORTED int bcInterfaceSolve_initSepRoutine(InterfaceModel* m,
                                                        void (*julia_func)(BcSolution&, const void*,
                                                                           std::list<BcConstr>&, int),
                                                        void* userdata, char typecut);

extern "C" EXPORTED int bcInterfaceSolve_initInfoRoutine(InterfaceModel* m,
                                                         void (*julia_func)(BcSolution&, const void*),
                                                         void* userdata);

extern "C" EXPORTED char bcInterfaceSolve_getSepCbType(InterfaceModel* m, int id);

extern "C" EXPORTED int bcInterfaceSolve_addSepCut(InterfaceModel* m, int idcb,
                                                   std::list<BcConstr>& cutList, double* coeffs, int* colId, int size,
                                                   char sens, double rhs);

extern "C" EXPORTED void bcInterfaceSolve_addCut(std::list<BcConstr>& cutList, BcConstr& cstr);

extern "C" EXPORTED void bcInterfaceSolve_setInitialSol(InterfaceModel* m, int* var_idx, double* var_vals, int nb_var,
                                                        int spid);
extern "C" EXPORTED void bcInterfaceSolve_provideInitialSol(InterfaceModel* m);

extern "C" EXPORTED BcSolution* bcSolution_new();
extern "C" EXPORTED void bcSolution_delete(BcSolution* solution);
extern "C" EXPORTED int bcSolution_next(BcSolution* solution);
extern "C" EXPORTED int bcSolution_start(BcSolution* solution, InterfaceModel* m);
extern "C" EXPORTED int bcSolution_getMultiplicity(BcSolution* solution, int & value);
extern "C" EXPORTED int bcSolution_getCost(BcSolution* solution, double & value);
extern "C" EXPORTED int bcSolution_getTrueCost(BcSolution* solution, double & value);

extern "C" EXPORTED int bcSolution_getAggSolution(InterfaceModel* m, BcSolution* solution, double array[], int size);
extern "C" EXPORTED int bcSolution_getValues(InterfaceModel* m, BcSolution* solution, double array[], int size);
extern "C" EXPORTED int bcSolution_getValueOfVar(InterfaceModel* m, BcSolution* solution, int varcol, double& value);
extern "C" EXPORTED const char* bcSolution_getNameOfSolForm(BcSolution* solution);

extern "C" EXPORTED int bcSolution_getNbNodes(BcSolution* solution);
extern "C" EXPORTED int bcSolution_getResConsumption(BcSolution* solution, double resCons[], int nodeIds[], int size,
                                                     int resId);
extern "C" EXPORTED int bcSolution_getArcsIds(BcSolution* solution, int arcsIds[], int size);
extern "C" EXPORTED int bcSolution_getProblemFirstId(BcSolution* solution);
extern "C" EXPORTED int bcSolution_print(BcSolution* solution);
////////////////////////////////////////


extern "C" {
    struct EXPORTED RCSPNetworkInterface {
        BcNetwork network;
        std::vector<BcVertex> vertices;
        std::vector<BcArc> arcs;
        //std::vector<BcNetworkElementaritySet> elementaritySets;
        std::map<int, BcNetworkResource> resources;

        RCSPNetworkInterface(InterfaceModel* m, int spType, int* spIndex, int nbNodes, int nbElementaritySets,
                                      int nbPackingSets, int nbCoveringSets);
    };
}

extern "C" EXPORTED RCSPNetworkInterface* bcRCSP_new(InterfaceModel* m, int spType, int* spIndex, int nbNodes,
                                                     int nbElementaritySets, int nbPackingSets, int nbCoveringSets);
extern "C" EXPORTED void bcRCSP_delete(RCSPNetworkInterface* n);

extern "C" EXPORTED int bcRCSP_newResource(RCSPNetworkInterface* n, int id);
extern "C" EXPORTED void bcRCSP_setAsMainResource(RCSPNetworkInterface* n, int resId, double stepValue);
extern "C" EXPORTED void bcRCSP_setAsNonDisposableResource(RCSPNetworkInterface* n, int resId);
extern "C" EXPORTED void bcRCSP_setSpecialResourceAsNonDisposable(RCSPNetworkInterface* n, int resId);

extern "C" EXPORTED int bcRCSP_setSource(RCSPNetworkInterface* n, int nodeId);
extern "C" EXPORTED int bcRCSP_setSink(RCSPNetworkInterface* n, int nodeId);

extern "C" EXPORTED int bcRCSP_addElementaritySets(RCSPNetworkInterface* n, int nbSets);
extern "C" EXPORTED int bcRCSP_attachElementaritySetToNode(RCSPNetworkInterface* n, int nodeId, int elementaritySetId);
extern "C" EXPORTED int bcRCSP_attachElementaritySetToEdge(RCSPNetworkInterface* n, int edgeId, int elementaritySetId);

extern "C" EXPORTED int bcRCSP_newArc(RCSPNetworkInterface* n, int tailId, int headId, double cost);

extern "C" EXPORTED int bcRCSP_attachBcVarToArc(RCSPNetworkInterface* n, int arcId, InterfaceModel* m, int colId,
                                                double coeff);
extern "C" EXPORTED int bcRCSP_addBinaryResourceConsumption(RCSPNetworkInterface* n, int nodeId, int binaryResId,
                                                            int consumption, int lowerBound, int upperBound);
extern "C" EXPORTED int bcRCSP_setVertexConsumptionLB(RCSPNetworkInterface * n, int nodeId, int resId, double lb);
extern "C" EXPORTED int bcRCSP_setVertexConsumptionUB(RCSPNetworkInterface* n, int nodeId, int resId, double ub);
extern "C" EXPORTED int bcRCSP_setVertexSpecialConsumptionLB(RCSPNetworkInterface* n, int nodeId, int resId, double lb);
extern "C" EXPORTED int bcRCSP_setVertexSpecialConsumptionUB(RCSPNetworkInterface* n, int nodeId, int resId, double ub);
extern "C" EXPORTED int bcRCSP_setArcConsumptionLB(RCSPNetworkInterface* n, int arcId, int resId, double lb);
extern "C" EXPORTED int bcRCSP_setArcConsumptionUB(RCSPNetworkInterface* n, int arcId, int resId, double ub);
extern "C" EXPORTED int bcRCSP_setEdgeConsumptionValue(RCSPNetworkInterface* n, int edgeId, int resId,
                                                       double consumption);
extern "C" EXPORTED int bcRCSP_setEdgeSpecialConsumptionValue(RCSPNetworkInterface* n, int edgeId, int resId,
                                                              double consumption);
extern "C" EXPORTED int bcRCSP_addVertexToMemOfElementaritySet(RCSPNetworkInterface* n, int nodeId,
                                                               int elementaritySetId);
extern "C" EXPORTED int bcRCSP_addEdgeToMemOfElementaritySet(RCSPNetworkInterface* n, int edgeId,
                                                             int elementaritySetId);
extern "C" EXPORTED int bcRCSP_addVertexToPackingSet(RCSPNetworkInterface* n, int vertexId, int packingSetId);
extern "C" EXPORTED int bcRCSP_addEdgeToPackingSet(RCSPNetworkInterface* n, int edgeId, int packingSetId);
extern "C" EXPORTED int bcRCSP_addVertexToCoveringSet(RCSPNetworkInterface* n, int vertexId, int coveringSetId);
extern "C" EXPORTED int bcRCSP_addEdgeToCoveringSet(RCSPNetworkInterface* n, int edgeId, int coveringSetId);
extern "C" EXPORTED int bcRCSP_addPackingSetToPackingSetCutNeighbourhood(RCSPNetworkInterface* n, int elemSetId,
                                                                         int packingSetId);
extern "C" EXPORTED int bcRCSP_setElemSetsDistanceMatrix(RCSPNetworkInterface* n, double** dists, int num_packsets);
//extern "C" EXPORTED int bcRCSP_addElementaritySetToMemOfElementaritySetCut(RCSPNetworkInterface* n,
//                                                                           int elementaritySetIdToAdd,
//                                                                           int elementaritySetId);
#ifdef NOSTANDALONE
extern "C" int EXPORTED bcRCSP_createOracle(RCSPNetworkInterface* n, InterfaceModel* m, int spType, int* spId);
#else
extern "C" int EXPORTED bcRCSP_createOracle(RCSPNetworkInterface* n, InterfaceModel* m, int spType, int* spId,
                                            bool saveStandalone, char* standaloneFileName);
#endif
extern "C" int EXPORTED bcRCSP_addGenericCapacityCut(InterfaceModel* m, int maxCapacity, int* demands, int demands_size,
                                                     bool isFacultative, double rootPrioLvl, double nonRootPrioLvl,
                                                     int twoPathCutResId = -1);
extern "C" int EXPORTED bcRCSP_addGenericStrongKPathCut(InterfaceModel* m, int maxCapacity, int* demands,
                                                        int demands_size, bool isFacultative, double rootPrioLvl,
                                                        double nonRootPrioLvl);
extern "C" int EXPORTED bcRCSP_addGenericLimMemOneCut(InterfaceModel* m);
extern "C" int EXPORTED bcRCSP_addGenericCliqueCut(InterfaceModel* m);
extern "C" int EXPORTED bcRCSP_addPathsPerNetworkBranching(InterfaceModel* m, double priority);
extern "C" int EXPORTED bcRCSP_addPackSetAssignBranching(InterfaceModel* m, double priority);
extern "C" int EXPORTED bcRCSP_addElemSetResourceConsumptionBranching(InterfaceModel* m, double priority);
extern "C" int EXPORTED bcRCSP_addPackSetRyanAndFosterBranching(InterfaceModel* m, double priority);
extern "C" int EXPORTED bcRCSP_addPermanentRyanAndFosterConstraint(RCSPNetworkInterface* n, int firstPackSetId,
                                                                   int secondPackSetId, bool together);

extern "C" int EXPORTED bcRCSP_addAssociatedVarToResource(RCSPNetworkInterface* n, int resId, InterfaceModel* m,
                                                          int colId);
extern "C" int EXPORTED bcRCSP_addGenericKnapsackConsumptionCut(InterfaceModel* m, double rootPriorityLvl,
                                                                double nonRootPriorityLvl);

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#ifdef __clang__
#pragma GCC diagnostic pop
#endif

#pragma GCC visibility pop

#endif /* BCCMODELINTERFACE_HPP */
