/**
 *
 * This file bcCInterfaceModel.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcCInterfaceModel.hpp"
#include "bcModelObjectiveC.hpp"
#include "bcModelRCSPSolver.hpp"
#include "bcInstanciatedVarConstrC.hpp"

#ifdef USE_NON_PUBLIC_CUTS
#include "bcModelNonPublicCuts.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <typeinfo>

// ### Data structures
// argc & argv contains command line parameters (-t BaPCodTreeName.dot for instance)
InterfaceModel::InterfaceModel(char* paramFile, bool printParam, bool intObj, bool intBound, int argc, char** argv) :
        problemInit(false),
        decompositionUsed(false),
        init(argc, argv, paramFile, printParam, true),
        model(init, "BaPCodModel", (intObj) ? BcObjStatus::minInt : BcObjStatus::minFloat),
        objective(model),
        mpform(model),
        master(model),
        colGenSp(model),
        solution(NULL),
        initialSolution(NULL),
        lastInitialSolution(NULL)
{
    std::srand(std::time(0));
    id = std::rand();

    if (intObj)
        objective.setMinMaxStatus(BcObjStatus::minInt);
    else
        objective.setMinMaxStatus(BcObjStatus::minFloat);

#ifdef JULIADEBUG
    std::cout << "\e[1;32m@InterfaceModel\e[00m building the model." << std::endl;
    std::cout << "\e[1m interfacemodel id = \e[00m " << id << std::endl;
#endif
}

InterfaceModel::~InterfaceModel()
{}

CutCallback::CutCallback(InterfaceModel* m, int id, char typecut, std::string name):
    id(id), nbcuts(0), type(typecut),
    cuts(m->master, name, typecut, 3.0,1.0, false)
{
    cuts.defineIndexNames(MultiIndexNames('i'));
}

CutCallback::CutCallback(const CutCallback& that) = default;

CutCallback::~CutCallback()
{}

bcOracleSolution::bcOracleSolution(BcSolution& sol, BcFormulation& sp, int phase)
    : _solution(sol), _objVal(0.0), _lastSolution(sol), _sp(sp), _phaseOfStageApproach(phase)
{
}

bcOracleSolution::~bcOracleSolution()
{}

// ### Functors
// #### Oracles
JuliaSolverOracleFunctor::JuliaSolverOracleFunctor(void(*julia_fct)(const void*, const void*, const int),
                                                   void* userdata, int idoracle) :
    _julia_fct(julia_fct), _userdata(userdata), _idoracle(idoracle)
{
}

JuliaSolverOracleFunctor::~JuliaSolverOracleFunctor()
{}


bool JuliaSolverOracleFunctor::operator() (BcFormulation spPtr,
                                           double & objVal,
                                           double & primalBound,
                                           double & dualBound,
                                           BcSolution & primalSol,
                                           int & phaseOfStageApproach)
{
    bcOracleSolution * oracleSolution = new bcOracleSolution(primalSol, spPtr, phaseOfStageApproach);
    _julia_fct(oracleSolution, _userdata, _idoracle);
    primalBound = dualBound = objVal = oracleSolution->_objVal;
    delete oracleSolution;
    return true;
}

// #### Cut callback
JuliaSeparationRoutineFunctor::JuliaSeparationRoutineFunctor(void(*julia_fct)(BcSolution&, const void*,
                                                                              std::list<BcConstr>&, int),
                                                             void* userdata, int idfunctor) :
    _userdata(userdata), _idfunctor(idfunctor), _julia_fct(julia_fct)
{
}

JuliaSeparationRoutineFunctor::~JuliaSeparationRoutineFunctor()
{}

int JuliaSeparationRoutineFunctor::operator() (BcFormulation formPtr,
                                               BcSolution & primalSol,
                                               double & maxViolation,
                                               std::list<BcConstr>& cutList)
{
    _julia_fct(primalSol, _userdata, cutList, _idfunctor);
    return cutList.size();
}

JuliaDefineCstrFunctor::JuliaDefineCstrFunctor(BcConstrArray& cstrs, void* userdata,
                                               void(*julia_fct)(int*, void*, BcConstr&))
        : BcAddConstrFunctor(), _cstrs(cstrs), _userdata(userdata), _julia_fct(julia_fct)
{
//    std::cout << "&userdata = " << _userdata << std::endl;
}

void JuliaDefineCstrFunctor::operator()(const MultiIndex& indexArray)
{
    BcConstr cstr = _cstrs.getElement(indexArray);
    int index[8];
    for (int i = 0; i < 8; i++) {
        index[i] = indexArray.index(i);
    }
    _julia_fct(index, _userdata, cstr);
}

JuliaInformationalCallback::JuliaInformationalCallback(void (*julia_fct)(BcSolution&, const void*), void* userdata) :
    BcSolutionFoundCallback(), _userdata(userdata), _julia_fct(julia_fct)
{
}

void JuliaInformationalCallback::operator() (BcSolution& newSolution)
{
    std::cout << "\e[1;35m informational callback \e[00m" << std::endl;
    _julia_fct(newSolution, _userdata);
}

// ### Functions
extern MultiIndex arrayToMultiIndex(int* array)
{
    MultiIndex mi;
    int i = 0;
    while (i < 8 && array[i] >= 0) {
        mi += array[i];
        i++;
    }
    return mi;
}

extern const BcFormulation & getProblem(InterfaceModel* m, int type, int* spId)
{
#ifdef JULIADEBUGADV
    std::cout << "\e[41m getProblem( m (with id = " << m->id << ") , type = " << type << ", spId = " << arrayToMultiIndex(spId) << ") \e[00m" << std::endl;
#endif
    switch(type) {
        case 0: // MIP
            if (!m->problemInit) {
                m->problemInit = true;
                m->decompositionUsed = false;
            } else {
                if (m->decompositionUsed) {
                    std::cerr << "\e[31m BAPCOD decomposition error \e[00m" << std::endl;
                    std::cerr << "Want a MIP but decomposition is used." << std::endl;
                    exit(1);
                }
            }
            return (m->mpform)();
        case 1: // MASTER
            if (!m->problemInit) {
                m->problemInit = true;
                m->decompositionUsed = true;
            } else {
                if (!m->decompositionUsed) {
                    std::cerr << "\e[31m BAPCOD decomposition error \e[00m" << std::endl;
                    std::cerr << "Want a master but decomposition is not used." << std::endl;
                    exit(1);
                }
            }
            return m->master;
        case 2: // DW SP
            if (!m->problemInit) {
                m->problemInit = true;
                m->decompositionUsed = true;
            } else {
                if (!m->decompositionUsed) {
                    std::cerr << "\e[31m BAPCOD decomposition error \e[00m" << std::endl;
                    std::cerr << "Want a Dantzig-Wolfe subproblem but decomposition is not used." << std::endl;
                    exit(1);
                }
            }
            return  m->colGenSp.createElement(arrayToMultiIndex(spId));
        case 4: // MASTER GENERATED VAR
            if (!m->problemInit) {
                m->problemInit = true;
                m->decompositionUsed = true;
            } else {
                if (!m->decompositionUsed) {
                    std::cerr << "\e[31m BAPCOD decomposition error \e[00m" << std::endl;
                    std::cerr << "Want a master subproblem but decomposition is not used." << std::endl;
                    exit(1);
                }
            }
            return m->master;
        case 3: /// 3 is Benders, not available now
        default:
            std::cerr << "getProblem : Unrecognized type ! (type = "
            << type << ")" << std::endl;
            return (m->mpform)();
    }
}

extern const std::string getProblemName(const BcFormulation & f)
{
    std::stringstream ss;
    ss << f.name() << f.id();
    return ss.str();
}

std::string sptype(int spType)
{
    switch(spType) {
        case 0:
            return "MIP";
        case 1:
            return "MASTER";
        case 2:
            return "DW SUBPROBLEM";
        case 3:
            return "B SUBPROBLEM";
        case 4:
            return "GEN MASTER";
        default:
            return "NEITHER";
    }
}

extern "C" InterfaceModel* bcInterfaceModel_new(char* paramFile, bool printParam, bool intObj, bool intBound, int argc,
                                                char** argv)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;42m@bcInterfaceModel_new\e[00m" << std::endl;
#endif
    InterfaceModel * model = new InterfaceModel(paramFile, printParam, intObj, intBound, argc, argv);
    return model;
}

extern "C" void bcInterfaceModel_delete(InterfaceModel* m)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcm delete \e[00m" << std::endl;
    std::cout << " > delete model with id \e[1m" << m->id << "\e[0m" << std::endl;
#endif
    delete m;
}

extern "C" void bcInterfaceModel_initModel(InterfaceModel* m, int nbRows, int nbCols)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcm initModel \e[00m" << std::endl;
    std::cout << " > nbCols = " << nbCols << std::endl;
    std::cout << " > nbRows = " << nbRows << std::endl;
#endif
    m->vars.resize(nbCols, NULL);
    m->constrs.resize(nbRows, NULL);
}

extern "C" int bcInterfaceModel_subProblemMult(InterfaceModel* m, int l, int u, int spType, int* spIndex)
{
    MultiIndex spMid = arrayToMultiIndex(spIndex);
    BcFormulation sp = getProblem(m, spType, spIndex);
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcm subProblemMult : subproblem multiplicity\e[00m" << std::endl;
    std::cout << " > Subproblem with id = \e[1m" << sp.name() << "\e[00m " << std::endl
    << " > Multiplicity : lb = \e[1m" << l << "\e[00m & ub = \e[1m" << u << "\e[00m" << std::endl;
#endif
    if (!sp.isDefined())
    {
        std::cerr << "Cannot set multiplicity on subproblem" << spMid << " : not"
                  << " a subproblem." << std::endl;
        return 0;
    }
    sp >= l;
    sp <= u;
#ifdef JULIADEBUG
    std::cout << "\e[1;32m Done :\e[00m on subproblem \e[1m" << sp.name() << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_registerSubProblem(InterfaceModel* m, int type, int* id)
{
    MultiIndex mid = arrayToMultiIndex(id);
#ifdef JULIADEBUG
    std::cout << "\e[1;33@bcm registerSubProblem \e[0m" << std::endl;
    std::cout << " > registering problem id \e[1m" << mid
              << "\e[0m with type \e[1m" << sptype(type) << "\e[0m" << std::endl;
#endif
    m->subproblems.push_back(std::pair<int*,int>(id, type));
    BcFormulation f = getProblem(m, type, id);
#ifdef JULIADEBUG
    std::cout << " > \e[1m" << f.name() << "\e[0m registered." << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_registerVar(InterfaceModel* m, char* name, int columnId, int spType, int* spId,
                                            int* varId)
{
    MultiIndex varMid = arrayToMultiIndex(varId);
    MultiIndex spMid = arrayToMultiIndex(spId);
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcm registerVar \e[00m" << std::endl;
    std::cout << " > Registering variable :"
    << "    name = \e[1m" << name << "\e[00m"
    << "    id = \e[1m" << varMid << "\e[00m"
    << "    column = \e[1m" << columnId << "\e[00m" << std::endl
    << "   attached to problem :"
    << "    type = \e[1m" << sptype(spType) << "\e[00m"
    << "    id = \e[1m" << spMid << "\e[00m" << std::endl;
    if (columnId < 0 || columnId >= m->vars.size()) {
        std::cerr << "Cannot init variable with id " << columnId
        << ". Out of bounds."
        << std::endl;
        return false;
    }
#endif
    BcFormulation f = getProblem(m, spType, spId);
    if (!f.isDefined()) {
        if (spType) {
            std::cerr << "registerVar: Cannot get the subproblem with multi-index: "
            << spMid << "." << std::endl;
            return 0;
        }
        std::cerr << "registerVar : Cannot get the problem." << std::endl;
        return 0;
    }
    BcVarArray var(f, name);

    m->vars[columnId] = var.createElement(varMid);
    m->varspartition[getProblemName(f)].push_back(columnId);
#ifdef JULIADEBUG
    std::cout << "\e[1;32m Done : \e[00m variable \e[1m " << m->vars[columnId].name()
    << "\e[00m attached to problem \e[1m" <<  f.name()
    << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_registerGenericVar(InterfaceModel* m, char* name, int columnId)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm registerGenericVar\e[00m" << std::endl;
    std::cout << " > variable name \e[1m" << name << "\e[00m" << std::endl;
    std::cout << " > var array registered with key \e[1m" << columnId << "\e[00m" << std::endl;
#endif
    BcVarArray var(m->master, name);

    std::vector< std::pair<int*, int> >::const_iterator spIt = m->subproblems.begin();
    while (spIt != m->subproblems.end()) {
        if (spIt->second == 2) { // dw subproblem
            BcVarArray var(getProblem(m, 2, spIt->first), name);
        }
        spIt++;
    }
    m->vars[columnId] = BcVar(NULL);
    m->genvars[columnId] = name;
    return 1;
}

extern "C" int bcInterfaceModel_addDynVar(InterfaceModel* m, int columnId, char* name, int* varId, int spType,
                                          int* spId)
{
    MultiIndex varMid = arrayToMultiIndex(varId);
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm createVar\e[00m" << std::endl;
    std::cout << " > name of the dynamic variable \e[1m" << name << "\e[00m" << std::endl;
    std::cout << " > multi index of the dynamic variable \e[1m" << varMid << "\e[00m" << std::endl;
    std::cout << " > to be registered at \e[1m" << columnId << "\e[00m" << std::endl;
#endif
    if (columnId != m->vars.size()) {
        std::cerr << "Incorrect column index. "
        << " Must be " << m->vars.size()
        << " but it is " << columnId << "." << std::endl;
        return 0;
    }
    BcFormulation f = getProblem(m, spType, spId);
    BcVarArray var(f, name);
    if (var.isDefinedAt(varMid)) {
        std::cerr << "The variable " << var.getElement(varMid).name()
        << " has been already generated." << std::endl;
        return 0;
    }
    m->vars.push_back(var(varMid));
#ifdef JULIADEBUG
    std::cout << " > created variable is \e[1m" << m->vars[columnId].name() << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_addMembershipToCstr(InterfaceModel* m, CstrsMembership& list, int rowId, double coeff)
{
#ifdef JULIADEBUG
    std::cout << "\e[36m@bcm addMembershipToCstr\e[00m" << std::endl;
#endif
    if (rowId >= m->constrs.size()) {
        std::cerr << "Unknown constraint." << std::endl;
        return 0;
    }
    BcConstr cstr = m->constrs[rowId];
#ifdef JULIADEBUG
    std::cout << " > will add the coefficient \e[1m" << coeff << "\e[00m"
    << " to the constraint \e[1m" << cstr.name() << "\e[00m" << std::endl;
#endif
    list.push_back(std::pair<BcConstr, double>(cstr, coeff));
    return 1;
}

extern "C" void bcInterfaceModel_initVars(InterfaceModel* m, double* l, double* u, double* c)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm initVars \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    bcInterfaceModel_setVarLb(m, l, nbVars);
    bcInterfaceModel_setVarUb(m, u, nbVars);
    bcInterfaceModel_setVarCosts(m, c, nbVars);
#ifdef JULIADEBUG
    std::vector<BcVar>::const_iterator varIt;
    for (varIt = m->vars.begin(); varIt != m->vars.end(); varIt++) {
        if (varIt->isDefined()) {
            std::cout << "  initialized variable \e[1m" << std::setw(15) << varIt->name() << "\e[00m"
            << "  lb = \e[1m" << std::setw(5) << varIt->curLb() << "\e[00m"
            << "  ub = \e[1m" << std::setw(5) << varIt->curUb() << "\e[00m"
            << "  cost = \e[1m" << varIt->originalCost() << "\e[00m" << std::endl;
        }
    }
#endif
}

extern "C" int bcInterfaceModel_getVarLb(InterfaceModel* m, double* lb, int size)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;42m bcm getVarLb \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (size != nbVars) {
        std::cerr << "getVarLb : incorrect size of input array." << std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        if (m->vars[var].isDefined()) {
            lb[var] = m->vars[var].curLb();
        } else {
            lb[var] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    return 1;
}

extern "C" int bcInterfaceModel_setVarLb(InterfaceModel* m, double* lb, int size)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm setVarLb \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (size != nbVars) {
        std::cerr << "setVarLb : incorrect size of input array." << std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        if (m->vars[var].isDefined()) {
            if (!std::isinf(lb[var])) {
                m->vars[var] >= lb[var];
            } else {
                m->vars[var] >= - std::numeric_limits<double>::infinity();
            }
        }
    }
    return 1;
}

extern "C" int bcInterfaceModel_getVarUb(InterfaceModel* m, double* ub, int size)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;42m bcm getVarUb \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (size != nbVars) {
        std::cerr << "getVarUb : incorrect size of input array." << std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        if (m->vars[var].isDefined()) {
            ub[var] = m->vars[var].curUb();
        } else {
            ub[var] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    return 1;
}

extern "C" int bcInterfaceModel_setVarUb(InterfaceModel* m, double* ub, int size)
{
#ifdef JULIADEBUG

    std::cout << "\e[1;36m@bcm setVarUb \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cerr << "setVarUb : incorrect size of input array." << std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        if (!std::isinf(ub[var])) {
            m->vars[var] <= ub[var];
            m->vars[var].globalUb(ub[var]);
        }
    }
    return 1;
}

extern "C" int bcInterfaceModel_getVarCosts(InterfaceModel* m, double* c,
                                            int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm getVarCosts \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cerr << "getVarCosts : incorrect size of input array." <<std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        c[var] = m->vars[var].originalCost();
    }
    return 1;
}

extern "C" int bcInterfaceModel_setVarCosts(InterfaceModel* m, double* c,
                                            int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm setVarCosts \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cerr << "setVarCosts : incorrect size of input array." <<std::endl;
        return 0;
    }
    for (int var = 0; var < nbVars; var++) {
        (m->objective)() += c[var] * m->vars[var];
    }
    return 1;
}

extern "C" int bcInterfaceModel_getVarType(InterfaceModel* m, char* t,
                                           int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;42m bcm getVarType \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cerr << "getVarType : incorrect size of input array." << std::endl;
        return 0;
    }
    std::cout << "TODO" << std::endl;
    return 1;
}

extern "C" int bcInterfaceModel_setVarType(InterfaceModel* m, char* t,
                                           int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm setVarType \e[00m" << std::endl;
#endif
    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cerr << "setVarType : incorrect size of input array." << std::endl;
        return 0;
    }
    for (int colId = 0; colId < nbVars;  colId++) {
        char varType = t[colId];
#ifdef JULIADEBUGADV
        std::cout << "Var #\e[1m" << colId << "\e[0m";
        if (m->vars[colId].isDefined()) {
            std::cout << " defined : \e[1m" << m->vars[colId] << "\e[0m with type \e[1m" << varType << "\e[0m" << std::endl;
        } else {
            std::cout << " not defined (type was \e[1m" << varType << "\e[0m). Check in generic variables..." << std::endl;
        }
#endif
        if (m->vars[colId].isDefined()) {
            m->vars[colId].type(varType);
        } else {
            std::map<int, char*>::iterator genvarIt = m->genvars.find(colId);
            if (genvarIt != m->genvars.end()) {
                char* name = genvarIt->second;
#ifdef JULIADEBUGADV
                std::cout << "      > generic variable \e[1m" << name << "\e[0m" << std::endl;
#endif
                BcVarArray variable(m->master, name);
                variable.type(varType);
                std::vector< std::pair<int*, int> >::const_iterator spIt = m->subproblems.begin();
                while (spIt != m->subproblems.end()) {
                    if (spIt->second == 2) { // dw subproblem
                        BcVarArray variable(getProblem(m, 2, spIt->first), name);
                        variable.type(varType);
                    }
                    spIt++;
                }
            }
        }
    }
    return 1;
}

extern "C" int bcInterfaceModel_setVarPriorityInMaster(InterfaceModel* m, char* varName,
                                               int spType, int* spId, double priority) {
    BcFormulation f = getProblem(m, spType, spId);
    BcVarArray x(f, varName);
#ifdef JULIADEBUG
    std::cout << "\e[1;36m @bcm setVarPriorityInMaster \e[00m" << std::endl;
    std::cout << " > variables \e[1m" << x.genericName()
              << "\e[0m of subproblem \e[1m" << f.name()
              << "\e[0m have priority \e[1m" << priority
              << "\e[0m in master." << std::endl;
#endif
    x.priorityForMasterBranching(priority);
    return 1;
}

extern "C" int bcInterfaceModel_setVarPriorityInSp(InterfaceModel* m, char* varName,
                                                       int spType, int* spId, double priority) {
    BcFormulation f = getProblem(m, spType, spId);
    BcVarArray x(f, varName);
#ifdef JULIADEBUG
    std::cout << "\e[1;36m @bcm setVarPriorityInSp \e[00m" << std::endl;
    std::cout << " > variables \e[1m" << x.genericName()
    << "\e[0m of subproblem \e[1m" << f.name()
    << "\e[0m have priority \e[1m" << priority
    << "\e[0m in subproblems." << std::endl;
#endif
    x.priorityForSubproblemBranching(priority);
    return 1;
}

extern "C" int bcInterfaceModel_addAggrSubProbVarBranching(InterfaceModel* m, char* varName, double highestPriority,
                                                           double priority, int ignoredIndices,
                                                           bool toBeUsedInPreprocessing)
{
    BcAggrSubProbVarBranching(m->master, varName, highestPriority, priority, ignoredIndices,
                              toBeUsedInPreprocessing);
    return 1; 
}

extern "C" int bcInterfaceModel_registerBranchingExpression(InterfaceModel* m, char* expName, double priority)
{
    BcBranchingConstrArray bcBrConstr(m->master, expName, SelectionStrategy::MostFractional, priority);
    m->branchingConstrs.push_back(bcBrConstr);
    return m->branchingConstrs.size();
}

extern "C" int bcInterfaceModel_addBranchingExpression(InterfaceModel* m, int arrayId, int* expId, int* varIds,
                                                       double* coeffs, int length)
{
    MultiIndex expMid = arrayToMultiIndex(expId);
    BcBranchingConstrArray expArray = m->branchingConstrs[arrayId-1];
    expArray.createElement(expMid);
    for (int i = 0; i < length; i++) {
        expArray.getElement(expMid) +=  coeffs[i] * m->vars[varIds[i]];
    }
    return 1;
}

extern "C" int bcInterfaceModel_registerCstr(InterfaceModel* m, char* name,
                                             int rowId, int spType, int* spId, int* constrId) {
    MultiIndex spMid = arrayToMultiIndex(spId);
    MultiIndex cstrMid = arrayToMultiIndex(constrId);
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcm registerCstr \e[00m" << std::endl;
    std::cout << " > Registering constraint :"
    << "    name = \e[1m" << name << "\e[00m"
    << "    id = \e[1m" << cstrMid << "\e[00m"
    << "    row = \e[1m" << rowId << "\e[00m" << std::endl
    << "   attached to subproblem :"
    << "    type = \e[1m" << sptype(spType) << "\e[00m"
    << "    id = \e[1m" << spMid << "\e[00m" << std::endl;
    if (rowId < 0 || rowId >= m->constrs.size()) {
        std::cerr << "Cannot init constraint with id " << rowId
        << ". Out of bounds."
        << std::endl;
        return false;
    }
#endif
    BcFormulation f = getProblem(m, spType, spId);
    if (!f.isDefined())
    {
        if (spType)
        {
            std::cerr << "registerCstr : Cannot get the subproblem with "
                      << "multi-index : " << spMid << "." << std::endl;
            return 0;
        }
        std::cerr << "registerCstr : Cannot get the problem." << std::endl;
        return 0;
    }
    BcConstrArray constr(f, name);
    m->constrs[rowId] = constr.createElement(cstrMid);
#ifdef JULIADEBUG
    std::cout << "\e[1;32m Done : \e[00m constraint \e[1m " << m->constrs[rowId].name()
    << "\e[00m attached to problem \e[1m" <<  f.name()
    << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_initCstrs(InterfaceModel* m, int* colPtr,
                                          int* firstRowId, double* nonZeroVal, double* lb, double* ub)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm initCstrs \e[00m" << std::endl;
#endif
    int nbRows = m->constrs.size();
    int nbCols = m->vars.size();

    // Bounds
    for (int row = 0; row < nbRows; row++) {
        bool ub_isinf = std::isinf(ub[row]);
        bool lb_isinf = std::isinf(lb[row]);

        if (!(ub_isinf && lb_isinf) && ub[row] == lb[row]) {
            m->constrs[row] == ub[row];
            continue;
        }

        if (!ub_isinf && !lb_isinf) {
            std::cerr << "BaPCod does not support range constraints."
            << std::endl;
            return 0;
        }

        if (!ub_isinf) {
            m->constrs[row] <= ub[row];
        } else if (!lb_isinf) {
            m->constrs[row] >= lb[row];
        }
    }
    // Coefficients
    int k = 0;
    int col = 0;
    while (k != colPtr[nbCols]) {
        while (k < colPtr[col+1]) {
            int row = firstRowId[k];
            double a = nonZeroVal[k];
            m->constrs[row] += a * m->vars[col];
            k++;
        }
        col++;
    }
    return 1;
}

extern "C" int bcInterfaceModel_registerGenericCstr(InterfaceModel* m, char* name,
                                                    int rowId) {
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm registerGenericCstr \e[00m" << std::endl;
    std::cout << " > constraint name \e[1m" << name << "\e[00m" << std::endl;
    std::cout << " > cstr array registered at key \e[1m" << rowId << "\e[00m" << std::endl;
#endif
    BcConstrArray cstr(m->master, name);
    cstr.flag('d'); //constraint is dynamic
    cstr.dualVal(1); // Default dual value must be set to 1 otherwise constraints are
                                  // added to the master formulation when they are generated in the
                                  // oracle (bcProblemC.cpp: about 7070 zeroTest).
    m->constrs[rowId] = BcConstr(NULL);
    return 1;
}

extern "C" int bcInterfaceModel_setConstrBcType(BcConstr* c, char type) {
    c->type(type);
    return 1;
}

extern "C" int bcInterfaceModel_setConstrArrayBcType(BcConstrArray* c, char type) {
    c->type(type);
    return 1;
}

extern "C" int bcInterfaceModel_attachCstrFunctor(InterfaceModel* m, char* name,
                                                  void (*julia_fct)(int*, void *, BcConstr&), void* jcstr) {
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm attachCstrFunctor \e[00m" << std::endl;
    std::cout << " > to \e[1m" << name << "\e[00m" << std::endl;
#endif
    BcConstrArray cstr(m->master, name);
#ifdef JULIADEBUG
    std::cout << " > attach functor to generic cstr array \e[m" << cstr.genericName() << "\e[00m" << std::endl;
#endif
    JuliaDefineCstrFunctor * rowFunctor = new JuliaDefineCstrFunctor(cstr, jcstr, julia_fct);
    cstr.attach(rowFunctor);
    return 1;
}

extern "C" int bcInterfaceModel_addDynCstr(InterfaceModel* m, int rowId, char* name, int* cstrId)
{
    MultiIndex cstrMid = arrayToMultiIndex(cstrId);
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm addDynCstr\e[00m" << std::endl;
    std::cout << " > name of the dynamic constraint \e[1m" << name << "\e[00m" << std::endl;
    std::cout << " > multi index of the dynamic constraint \e[1m" << cstrMid << "\e[00m" << std::endl;
    std::cout << " > to be registered at \e[1m" << rowId << "\e[00m" << std::endl;
#endif
    if (rowId != m->constrs.size()) {
        std::cerr << "Incorrect row index. "
        << " Must be " << m->constrs.size()
        << " but it is " << rowId << "." << std::endl;
        return 0;
    }
    BcConstrArray cstr(m->master, name);
    m->constrs.push_back(cstr(cstrMid));
#ifdef JULIADEBUG
    std::cout << " > created constraint is \e[1m" << m->constrs[rowId].name() << "\e[00m" << std::endl;
#endif
    return 1;
}

extern BcConstr & bcInterfaceModel_getCstr(InterfaceModel* m, int rowId) {
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm getConstr\e[00m" << std::endl;
#endif
    // check inputs
    return m->constrs[rowId];
}


extern "C" int bcInterfaceModel_addCstrTerm(InterfaceModel* m, BcConstr& cstr,
                                            int varcol, double coeff) {
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm addCstrTerm \e[00m" << std::endl;
    std::cout << " > On constraint \e[1m" << cstr.name() << "\e[00m" << std::endl;
#endif
    //cstr += coeff * m->vars[varcol];
    InstanciatedConstr* icstrPtr(cstr);
    InstanciatedVar* ivarPtr(m->vars[varcol]);
    icstrPtr->includeMember(ivarPtr, coeff, true);
#ifdef JULIADEBUG
    std::cout << " > Added the term \e[1m" << coeff << " * " << m->vars[varcol].name() << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceModel_addCstrTerms(InterfaceModel* m, BcConstr& cstr, int* varcols, int nbcols, double coeff)
{
    for(int i = 0; i < nbcols; i++) {
        //cstr += coeff * m->vars[varcols[i] - 1];
        InstanciatedConstr* icstrPtr(cstr);
        InstanciatedVar* ivarPtr(m->vars[varcols[i] - 1]);
        icstrPtr->includeMember(ivarPtr, coeff, true);
    }
    return 1;
}

extern "C" void bcInterfaceModel_addCstrRhs(BcConstr& cstr, char sense,
                                            double rhs) {
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm addCstrRhs \e[00m" << std::endl;
    std::cout << " > On constraint \e[1m" << cstr.name() << "\e[00m" << std::endl;
    std::cout << " > Sense is \e[1m" << sense << "\e[00m and rhs is \e[1m" << rhs << "\e[00m" << std::endl;
#endif
    if (sense == 'L') {
        cstr <= rhs;
    } else if (sense == 'E') {
        cstr == rhs;
    } else if (sense == 'G') {
        cstr >= rhs;
    }
}

extern "C" int bcInterfaceModel_cstrUsedInPreprocessing(InterfaceModel* m, char* cstrName, int spType, int* spId,
                                                        bool used)
{
    BcFormulation f = getProblem(m, spType, spId);
    BcConstrArray cstr(f, cstrName);
    // to do error msg
#ifdef JULIADEBUG
    std::cout << "\e[31m@bcm cstrUsedInPreprocessing \e[00m" << std::endl;
    std::cout << " > Constraint \e[1m" << cstr.genericName() << "\e[00m of \e[1m" << f.name() << "\e[0m ";
    if (used) {
        std::cout << "\e[1;32m used \e[00m";
    } else {
        std::cout << "\e[1;31m not used \e[00m";
    }
    std::cout << " in preprocessing." << std::endl;
#endif
    cstr.toBeUsedInPreprocessing(used);
    return 1;
}

extern "C" void bcInterfaceModel_setArtCostValue(InterfaceModel* m,
                                                 double artcostvalue) {
    (m->objective)().setArtCostValue(artcostvalue);
}


extern "C" void bcInterfaceModel_setObjMagnitude(InterfaceModel* m,
                                                 double magnitude) {
    (m->objective)() == magnitude;
}

extern "C" void bcInterfaceModel_setObjLb(InterfaceModel* m, double lb) {
    (m->objective)() >= lb;
}

extern "C" void bcInterfaceModel_setObjUb(InterfaceModel* m, double ub) {
    (m->objective)() <= ub;
}

extern "C" int bcInterfaceModel_getSense(InterfaceModel* m) {
    return 1;
}

extern "C" int bcInterfaceModel_setSubproblemPriority(InterfaceModel* m,
                                                      int* blockgroup, double priority) {
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcm setSubproblemPriority \e[00m" << std::endl;
#endif
    MultiIndex bgMid = arrayToMultiIndex(blockgroup);
    BcFormulation f = getProblem(m, false, blockgroup);
    if (!f.isDefined()) {
        std::cerr << "setSubproblemPriority : Cannot get the subproblem with multi-index: "
                  << bgMid << "." << std::endl;
        return 0;
    }
    f.priorityLevel(priority);
    return 1;
}

extern "C" void bcSolution_enumerateAllColumns(InterfaceModel* m, BcSolution* solution, int & nbEnumeratedCols)
{
    *solution = m->model.enumerateAllColumns(nbEnumeratedCols);
    m->solution = BcSolution(*solution);
}

//To add in a new file (solver part)
extern "C" void bcInterfaceSolve_optimize(InterfaceModel* m, BcSolution* solution)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;36m@bcs optimize \e[00m" << std::endl;
    std::cout << "going to optimize \e[1m InterfaceModel with id \e[00m = " << m->id << std::endl;
#endif
    if (!m->decompositionUsed) {
        BcFormulation mpf((m->mpform)());
        *solution = mpf.solve();
        m->solution = BcSolution(*solution);
    } else {
        *solution = m->model.solve();
        m->solution = BcSolution(*solution);
    }

    if (!solution->defined()) {
    } else {
        std::cout << "\e[32mBaPCod : solution defined.\e[00m" << std::endl;
#ifdef JULIADEBUGADV
        std::cout << " solution = " << m->solution << std::endl;
#endif
    }
}

extern "C" int bcInterfaceSolve_getMasterDual(InterfaceModel* m, int constr, double& dual)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs getVarCurCost\e[00m" << std::endl;
    int nbConstrs = m->constrs.size();
    if (nbConstrs <= constr) {
        std::cerr << "getMasterDual : incorrect size of input array." << std::endl;
        return 0;
    }
#endif
    dual = - m->constrs[constr].curDualVal();
#ifdef JULIADEBUG
    std::cout << " > dual of constraint is \e[1m" << dual << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceSolve_getVarCurCost(InterfaceModel* m, int var,
                                              double& curCost) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs getVarCurCost\e[00m" << std::endl;
    int nbVars = m->vars.size();
    if (nbVars <= var) {
        std::cerr << "getCurCost : incorrect size of input array." << std::endl;
        return 0;
    }
#endif
    curCost = m->vars[var].curCost();
#ifdef JULIADEBUG
    std::cout << " > cost of the variable is \e[1m" << curCost << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceSolve_getVarCurUB(InterfaceModel* m, int var,
        double& curUB) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs getVarCurUB\e[00m" << std::endl;
    int nbVars = m->vars.size();
    if (nbVars <= var) {
        std::cerr << "getCurUB : incorrect size of input array." << std::endl;
        return 0;
    }
#endif
    curUB = m->vars[var].curUb();
#ifdef JULIADEBUG
    std::cout << " > curUB of the variable is \e[1m" << curUB << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceSolve_getVarCurLB(InterfaceModel* m, int var,
        double& curLB) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs getVarCurLB\e[00m" << std::endl;
    int nbVars = m->vars.size();
    if (nbVars <= var) {
        std::cerr << "getCurLB : incorrect size of input array." << std::endl;
        return 0;
    }
#endif
    curLB = m->vars[var].curLb();
#ifdef JULIADEBUG
    std::cout << " > curLB of the variable is \e[1m" << curLB << "\e[00m" << std::endl;
#endif
    return 1;
}

extern "C" int bcInterfaceSolve_getOptimalityGapTolerance(InterfaceModel* m, double& ogt) {
    ogt = m->init.param().optimalityGapTolerance;
    return 1;
}

extern "C" int bcInterfaceSolve_getDynVarCurCost(InterfaceModel* m, char* varname, int* varid, int sptype, int* spid,
                                                 double& curCost)
{
    BcFormulation f = getProblem(m, sptype, spid);
    BcVarArray var(f, varname);
    MultiIndex varmid = arrayToMultiIndex(varid);
    BcVarIndex vi(var, varmid);
    curCost = vi.curCost();
    return 1;
}

extern "C" int bcInterfaceSolve_initOracle(InterfaceModel* m,
                                           void (*julia_func)(const void*, const void*, const int),
                                           int sp_type, int* sp_id, void* userdata, int idoracle)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs initOracle \e[00m" << std::endl;
#endif
    JuliaSolverOracleFunctor * oraclePtr = new JuliaSolverOracleFunctor(julia_func, userdata, idoracle);
    BcFormulation f = getProblem(m, sp_type, sp_id);
    if (!f.isDefined())
    {
        MultiIndex spMid = arrayToMultiIndex(sp_id);
        std::cerr << "initOracle : Cannot get the subproblem with multi-index: "
                  << spMid << "." << std::endl;
        return 0;
    }
    f.attach(oraclePtr);
    return 1;
}

extern "C" int bcInterfaceSolve_addToOracleSol(InterfaceModel* m, bcOracleSolution* s, int varId, double varVal)
{
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs addToOracleSol \e[00m" << std::endl;
#endif
    if (NULL == s->_lastSolution) {
        std::cerr << "addToOracleSol : solution not initialized." << std::endl;
        return 0;
    }
    m->vars[varId] = varVal;
    s->_lastSolution += m->vars[varId];
    return 1;
}

extern "C" int bcInterfaceSolve_newOracleSol(bcOracleSolution* s) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs newOracleSol \e[00m" << std::endl;
#endif
    if (NULL == s->_solution) {
        std::cerr << "newOracleSol : solution not initialized." << std::endl;
        return 0;
    }
    s->_lastSolution = BcSolution(s->_sp);
    s->_solution += s->_lastSolution;
    return 1;
}

extern "C" int bcInterfaceSolve_getOraclePhaseOfStage(bcOracleSolution* s) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs getOraclePhaseOfStage" << std::endl;
#endif
    return s->_phaseOfStageApproach;
}

extern "C" int bcInterfaceSolve_setObjValOfOracleSol(bcOracleSolution* s,
                                                     double objVal) {
#ifdef JULIADEBUG
    std::cout << "\e[1;34m@bcs setObjValOfOracleSol \e[00m" << std::endl;
#endif
    if (NULL == s->_solution) {
        std::cerr << "setObjValOfOracleSol : solution not initialized."
        << std::endl;
        return 0;
    }
    s->_objVal = objVal;
    return 1;
}

// bcTimeMain bcRecBestDb
extern "C" double bcInterfaceSolve_getStatisticValue(InterfaceModel* m,
                                                     char* name) {
    double value = m->init.getStatisticValue(std::string(name));
    return value;
}

extern "C" long bcInterfaceSolve_getStatisticTime(InterfaceModel* m,
                                                  char* name) {
    long value = m->init.getStatisticTime(std::string(name));
    return value;
}

extern "C" long bcInterfaceSolve_getStatisticCounter(InterfaceModel* m,
                                                     char* name) {
    long value = m->init.getStatisticCounter(std::string(name));
    return value;
}

extern "C" int bcInterfaceSolve_initSepRoutine(InterfaceModel* m,
                                               void (*julia_func)(BcSolution&, const void*, std::list<BcConstr>&, int),
                                               void* userdata, char typecut)
{
    int id = m->cutCbList.size();
    std::stringstream ss;
    if (typecut == 'C') {
        //ss << "lazyCutCb" << id;
        ss << "UserCuts" << id;
    } else if (typecut == 'F') {
        ss << "userCutCb" << id;
    } else {
        std::cerr << "\e[41m BaPCod : Unknonwn type of cut." << std::endl;
        return 0;
    }
    std::string name = ss.str();
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcs initLazySepRoutine \e[00m" << std::endl;
    std::cout << " > name = \e[1m" << name << "\e[00m" << std::endl;
#endif
    m->cutCbList.push_back(CutCallback(m, id, typecut, name));

    JuliaSeparationRoutineFunctor* sepPtr = new JuliaSeparationRoutineFunctor(julia_func, userdata, id);
    m->cutCbList[id].cuts.attach(sepPtr);
    return 1;
}

extern "C" int bcInterfaceSolve_initInfoRoutine(InterfaceModel* m,
                                                void (*julia_func)(BcSolution&, const void*),
                                                void* userdata) {
    std::cout << "\e[1;33m@bcs init info routine \e[00m" << std::endl;
    m->model.attach(new JuliaInformationalCallback(julia_func, userdata));
    return 1;
}


extern "C" char bcInterfaceSolve_getSepCbType(InterfaceModel* m, int id) {
    return m->cutCbList[id].type;
}

//extern "C" int bcInterfaceSolve_addLazyCut(InterfaceModel* m, int idcb,
extern "C" int bcInterfaceSolve_addSepCut(InterfaceModel* m, int idcb,
                                          std::list<BcConstr>& cutList, double* coeffs, int* colId, int size,
                                          char sens, double rhs) {
    CutCallback cb = m->cutCbList[idcb];
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcs addSepCut \e[00m" << std::endl;
    std::cout << "\t > model id = \e[1m" << m->id << "\e[00m" << std::endl;
    std::cout << "\t > id callback = \e[1m" << idcb << "\e[00m" << std::endl;
    std::cout << "\t > cut number = \e[1m" << cb.nbcuts << "\e[00m" << std::endl;
    std::cout << "\t > cutlist size = \e[1m" << cutList.size() << "\e[00m" << std::endl;
#endif
    BcConstr newLazyCut = cb.cuts(MultiIndex(cb.nbcuts));
    for (int i = 0; i < size; i++) {
        newLazyCut += coeffs[i] * m->vars[colId[i]-1];
    }
    if (sens == '<') {
        newLazyCut <= rhs;
    }
    else if (sens == '>') {
        newLazyCut >= rhs;
    } else {
        newLazyCut == rhs;
    }
    cutList.push_back(newLazyCut);
    m->cutCbList[idcb].nbcuts += 1;
    return 1;
}

extern "C" void bcInterfaceSolve_addCut(std::list<BcConstr>& cutList, BcConstr& cstr) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcs addCut \e[00m" << std::endl;
    std::cout << "\t > Adding constraint \e[1m" << cstr.name() << "\e[00m" << std::endl;
#endif
    cutList.push_back(cstr);
}

extern "C" void bcInterfaceSolve_setInitialSol(InterfaceModel* m, int* var_idx, double* var_vals, int nb_var, int spid)
{
    BcMasterArray master(m->model);
    BcColGenSpArray colGenSp(m->model);
    BcFormulation spForm = colGenSp[spid];
    
    m->initialSolution = BcSolution(spForm);
    if(m->lastInitialSolution != NULL)
        m->initialSolution += m->lastInitialSolution;
    
    for(int i = 0 ; i<nb_var ; i++){
        m->vars[var_idx[i]] = var_vals[i];
        m->initialSolution += m->vars[var_idx[i]];
    }
    
    m->lastInitialSolution = m->initialSolution;
}

extern "C" void bcInterfaceSolve_provideInitialSol(InterfaceModel* m){
    std::cout<< m->initialSolution << std::endl;
    BcMasterArray master(m->model);
    master().initializeWithSolution(m->initialSolution);
}

extern "C" BcSolution* bcSolution_new() {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution new\e[00m" << std::endl;
#endif
    return new BcSolution(NULL);
}

extern "C" void bcSolution_delete(BcSolution* solution) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution delete\e[00m" << std::endl;
#endif
    delete solution;
}

extern "C" int bcSolution_next(BcSolution* solution) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution next\e[00m" << std::endl;
#endif
    if (solution == NULL) {
        return 0;
    }

    if (!solution->next().defined()) {
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > current is \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
    std::cout << " > will return \e[1m" << solution->next().defined() << "\e[00m" << std::endl;
#endif
    *solution = solution->next();
    return 1;
}

extern "C" int bcSolution_start(BcSolution* solution, InterfaceModel* m) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution reset\e[00m" << std::endl;
#endif
    if (m->solution == NULL) {
        return 0;
    }
    *solution = m->solution;
    return 1;
}


extern "C" int bcSolution_getMultiplicity(BcSolution* solution, int & value) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution getMultiplicity\e[00m" << std::endl;
#endif
    if (solution == NULL) {
        return 0;
    }

    if (!solution->defined()) {
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > solution of model \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
    std::cout << " > multiplicity is \e[1m" << solution->getMultiplicity() << "\e[00m" << std::endl;
#endif
    value = solution->getMultiplicity();
    return 1;
}

extern "C" int bcSolution_getCost(BcSolution* solution, double & value) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution getCost\e[00m" << std::endl;
#endif
    if (solution == NULL) {
        return 0;
    }

    if (!solution->defined()) {
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > solution of model \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
    std::cout << " > cost is \e[1m" << solution->cost() << "\e[00m" << std::endl;
#endif
    value = solution->cost();
    return 1;
}

extern "C" int bcSolution_getTrueCost(BcSolution* solution, double & value) {
#ifdef JULIADEBUG
    std::cout << "\e[1;33m@bcsolution getTrueCost\e[00m" << std::endl;
#endif
    if (solution == NULL) {
        return 0;
    }

    if (!solution->defined()) {
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > solution of model \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
    std::cout << " > true cost is \e[1m" << solution->resetCost() << "\e[00m" << std::endl;
#endif
    value = solution->resetCost();
    return 1;
}

extern "C" int bcSolution_getAggSolution(InterfaceModel* m, BcSolution* solution, double array[], int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;32m @bcsol getAggSolution (4args)\e[00m" << std::endl;
#endif
    for (int i = 0; i < size; i++) {
        array[i] = 0.0;
    }

    if (!solution->defined()) {
        return 0;
    }

#ifdef JULIADEBUG
    std::cout << " > compute the aggregated solution from model \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
#endif

    int nbVars = m->vars.size();

    BcSolution curSol = *solution;
    int counter = 0;
    while (curSol != NULL) {
        // Step 1 : get the name of the problem attached to the solution
        std::string key = getProblemName(curSol.formulation());
#ifdef JULIADEBUG
        std::cout << " > build solution it = \e[1m" << counter
        << "\e[00m & problem = \e[1m" << key << "\e[00m" << std::endl;
#endif
        // Step 2 : Find the vector containing columns id of the variables of this problem
        std::map<std::string, std::vector<int> >::const_iterator columnsIdPtr;
        columnsIdPtr = m->varspartition.find(key);
        if (columnsIdPtr != m->varspartition.end()) {
            // Step 3 : Interate over columns id array
            std::vector<int>::const_iterator colIdPtr;
            colIdPtr = (columnsIdPtr->second).begin();
            while (colIdPtr != (columnsIdPtr->second).end()) {
                //Step 4 : add to the solution array
                array[*colIdPtr] += curSol.getVarVal(m->vars[*colIdPtr]) * curSol.getMultiplicity();
                colIdPtr++;
            }
        }
        // next solution
        curSol = curSol.next();
        counter++;
    }
    return 1;
}

extern "C" int bcSolution_getValues(InterfaceModel* m, BcSolution* solution, double array[], int size) {
#ifdef JULIADEBUG
    std::cout << "\e[1;32m @bcsol getValues\e[00m" << std::endl;
#endif
    for (int i = 0; i < size; i++) {
        array[i] = 0.0;
    }

    if (!solution->defined()) {
        //std::cout << "Solution not defined." << std::endl;
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > solution from model \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
#endif

    int nbVars = m->vars.size();
    if (nbVars != size) {
        std::cout << "The number of variables and the size of the array are not equal." << std::endl;
        return 0;
    }

    // Step 1 : get the name of the problem attached to the solution
    std::string key = getProblemName(solution->formulation());
    // Step 2 : Find the vector containing columns id of the variables of this problem
    std::map<std::string, std::vector<int> >::const_iterator columnsIdPtr;
    columnsIdPtr = m->varspartition.find(key);
    if (columnsIdPtr != m->varspartition.end()) {
        // Step 3 : Interate over columns id array
        std::vector<int>::const_iterator colIdPtr;
        colIdPtr = (columnsIdPtr->second).begin();
        while (colIdPtr != (columnsIdPtr->second).end()) {
            //Step 4 : add to the solution array
            double value = solution->getVarVal(m->vars[*colIdPtr]);
            if (value > 10e-6) {
                array[*colIdPtr] += value;
            }
            colIdPtr++;
        }
    }
    return 1;
}

extern "C" int bcSolution_getValueOfVar(InterfaceModel* m, BcSolution* solution, int varcol, double& value) {
#ifdef JULIADEBUG
    std::cout << "\e[1;32m @bcsol getDisValueOfVar\e[00m" << std::endl;
#endif
    if (solution == NULL) {
        //std::cout << "Solution pointer is NULL" << std::endl;
        return 0;
    }

    if (!solution->defined()) {
//        std::cout << "Solution not defined." << std::endl;
        return 0;
    }
#ifdef JULIADEBUG
    std::cout << " > solution from problem \e[1m" << solution->formulation().name() << "\e[00m" << std::endl;
#endif
    BcVar variable = m->vars[varcol];
    value = solution->getVarVal(variable);
#ifdef JULIADEBUG
    std::cout << " > variable \e[1m" << variable.name() << "\e[0m value is \e[1m" << value << "\e[0m" << std::endl;
#endif
    return 1;
}

extern "C" const char* bcSolution_getNameOfSolForm(BcSolution* solution) {
#ifdef JULIADEBUG
    std::cout << "\e[1;32m@bcsol getNameOfSolForm\e[00m" << std::endl;
#endif
    if (solution == NULL) {
//        std::cout << "Solution pointer is NULL" << std::endl;
        //std::string name = "null";
        return "null";
    }
    if (!solution->defined()) {
     //   std::cout << "Solution not defined." << std::endl;
        //std::string name = "undef";
        return "undef";
    }
    return solution->formulation().name().c_str();
}

extern "C" int bcSolution_getNbNodes(BcSolution* solution) {
  return solution->orderedIds().size() + 1;
}

extern "C" int bcSolution_getResConsumption(BcSolution* solution, double resCons[], int nodeIds[], int size, int resId)
{
  Solution* sol = solution->_solutionPtr;
  const NetworkFlow * netFlowPtr = sol->probConfPtr()->networkFlowPtr();

  if (netFlowPtr == NULL) {
    std::cout << "getResConsumption : no network flow attached to the current solution." << std::endl;
    return 0;
  }

  if (size != solution->orderedIds().size() + 1) {
    std::cout << "getResConsumption : size of input arrays should be " << solution->orderedIds().size() + 1 << "; got "
              << size << std::endl;
    return 0;
  }

  std::vector<std::vector<double> >::const_iterator resConsIt = solution->resConsumption().begin();
  std::vector<int>::const_iterator arcIdIt = solution->orderedIds().begin();

  int i = 0;
  nodeIds[i] = netFlowPtr->netArcPtr(*arcIdIt)->tailVertexPtr()->id();
  resCons[i] = (*resConsIt).at(resId);
  ++resConsIt;
  ++i;
  while (arcIdIt != solution->orderedIds().end()) {
    nodeIds[i] = netFlowPtr->netArcPtr(*arcIdIt)->headVertexPtr()->id();
    resCons[i] = (*resConsIt).at(resId);
    ++i;
    ++arcIdIt;
    ++resConsIt;
  }
  return 1;
}

extern "C" int bcSolution_getArcsIds(BcSolution* solution, int arcsIds[], int size) {
  Solution* sol = solution->_solutionPtr;
  const NetworkFlow * netFlowPtr = sol->probConfPtr()->networkFlowPtr();

  if (netFlowPtr == NULL) {
    std::cout << "getArcIds : no network flow attached to the current solution." << std::endl;
    return 0;
  }

  if (size != solution->orderedIds().size()) {
    std::cout << "getArcIds : size of input arrays should be " << solution->orderedIds().size() << "; got " << size
              << std::endl;
    return 0;
  }

  std::vector<int>::const_iterator arcIdIt = solution->orderedIds().begin();

  int i = 0;
  while (arcIdIt != solution->orderedIds().end()) {
    arcsIds[i] = *arcIdIt;
    ++i;
    ++arcIdIt;
  }
  return 1;
}

extern "C" int bcSolution_getProblemFirstId(BcSolution* solution) {
  MultiIndex i = solution->formulation().id();
  return i.first();
}

extern "C" int bcSolution_print(BcSolution* solution) {
  solution->print(std::cout);
  return 1;
}

#ifdef BCP_RCSP_IS_FOUND
////////////////////// RCSP INTERFACE /////////////////////////////

RCSPNetworkInterface::RCSPNetworkInterface(InterfaceModel* m, int spType, int* spIndex, int nbNodes,
                                           int nbElementaritySets, int nbPackingSets, int nbCoveringSets) :
    network(NULL)
{
    BcFormulation f = getProblem(m, spType, spIndex);
    network = BcNetwork(f, nbElementaritySets, nbPackingSets, nbCoveringSets);
#ifdef RCSPDEBUG
    std::cout << " > creating \e[1m" << nbNodes << "\e[0m nodes. " << std::endl;
    std::cout << " > creating \e[1m" << nbElementaritySets << "\e[0m elementarity sets. " << std::endl;
    std::cout << " > creating \e[1m" << nbPackingSets << "\e[0m packing sets. " << std::endl;
    std::cout << " > creating \e[1m " << nbCoveringSets << "\e[0m covering sets. " << std::endl;
#endif
    for (int n=0; n < nbNodes; n++) {
        vertices.push_back(network.createVertex());
    }
}

extern "C" RCSPNetworkInterface* bcRCSP_new(InterfaceModel* m, int spType, int* spIndex, int nbNodes,
                                            int nbElementaritySets, int nbPackingSets, int nbCoveringSets)
{
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP new \e[0m" << std::endl;
#endif
    RCSPNetworkInterface* interface = new RCSPNetworkInterface(m, spType, spIndex, nbNodes, nbElementaritySets,
                                                               nbPackingSets, nbCoveringSets);
    return interface;
}

extern "C" void bcRCSP_delete(RCSPNetworkInterface* r) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP delete \e[0m" << std::endl;
#endif
    delete r;
}

extern "C" int bcRCSP_newResource(RCSPNetworkInterface* n, int id) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP new resource \e[0m" << std::endl;
    std::cout << " > adding resource with id \e[1m" << id << "\e[0m. " << std::endl;
#endif
    n->resources[id] = BcNetworkResource(n->network, id);
    return 1;
}

extern "C" void bcRCSP_setAsMainResource(RCSPNetworkInterface* n, int resId, double stepValue) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set as main resource \e[0m" << std::endl;
    std::cout << " > resource with id \e[1m" << resId << "\e[0m becomes the main resource with step value = \e[1m " << stepValue << "\e[0m" << std::endl;
#endif
    // to do, check if the resouce exists
    n->resources[resId].setAsMainResourceWithStep(stepValue);

}

extern "C" void bcRCSP_setAsNonDisposableResource(RCSPNetworkInterface* n, int resId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set as non-disposable \e[0m" << std::endl;
    std::cout << " > resource with id \e[1m" << resId << "\e[0m becomes non-disposable" << std::endl;
#endif
    // to do, check if the resouce exists
    n->resources[resId].setAsNonDisposableResource();
}


extern "C" void bcRCSP_setSpecialResourceAsNonDisposable(RCSPNetworkInterface* n, int resId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set special as non-disposable \e[0m" << std::endl;
    std::cout << " > special resource with id \e[1m" << resId << "\e[0m becomes non-disposable" << std::endl;
#endif
    // to do, check if the resouce exists
    n->network.setSpecialResourceNonDisposable(resId);
}

extern "C" int bcRCSP_setSource(RCSPNetworkInterface* n, int nodeId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set source \e[0m" << std::endl;
    std::cout << " > Node \e[1m" << nodeId << "\e[0m becomes the source. " << std::endl;
#endif
    // to do : check if the node exists
    BcVertex v = n->vertices[nodeId];
    n->network.setPathSource(v);
    return 1;
}

extern "C" int bcRCSP_setSink(RCSPNetworkInterface* n, int nodeId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set source \e[0m" << std::endl;
    std::cout << " > Node \e[1m" << nodeId << "\e[0m becomes the sink. " << std::endl;
#endif
    // to do : check if the node exists
    BcVertex v = n->vertices[nodeId];
    n->network.setPathSink(v);
    return 1;
}

extern "C" int bcRCSP_addElementaritySets(RCSPNetworkInterface* n, int nbSets) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP add elementarity sets \e[0m" << std::endl;
    std::cout << " > creating \e[1m" << nbSets << "\e[0m sets. " << std::endl;
#endif
    return 1;
}

extern "C" int bcRCSP_attachElementaritySetToNode(RCSPNetworkInterface* n, int nodeId, int elementaritySetId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP attach elementarity set \e[0m" << std::endl;
    std::cout << " > attaching set \e[1m" << elementaritySetId << "\e[0m  to the node \e[1m " << nodeId << "\e[0m. " << std::endl;
#endif
    // to do : check if the node exists & same elem set
    BcVertex v  = n->vertices[nodeId];
    v.setElementaritySet(elementaritySetId);
    return 1;
}

extern "C" int bcRCSP_attachElementaritySetToEdge(RCSPNetworkInterface* n, int edgeId, int elementaritySetId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP attach elementarity set \e[0m" << std::endl;
    std::cout << " > attaching set \e[1m" << elementaritySetId << "\e[0m  to the edge \e[1m " << n->arcs[edgeId] << "\e[0m. " << std::endl;
#endif
    // to do : check if the node exists & same elem set
    BcArc a = n->arcs[edgeId];
    a.setElementaritySet(elementaritySetId);
    return 1;
}

extern "C" int bcRCSP_newArc(RCSPNetworkInterface* n, int tailId, int headId, double cost) {
    int arcId = n->arcs.size();
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP new arc \e[0m" << std::endl;
    std::cout << " > creating an arc from node \e[1m" << tailId
             << "\e[0m to the node \e[1m " << headId
             << "\e[0m. with id \e[1m" << arcId
             << "\e[0m and original cost \e[1m" << cost
             << "\e[0m." <<  std::endl;
#endif
    BcArc a = n->network.createArc(tailId, headId, cost);
    n->arcs.push_back(a);
    return arcId;
}

extern "C" int bcRCSP_attachBcVarToArc(RCSPNetworkInterface* n, int arcId, InterfaceModel* m, int colId, double coeff)
{
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP attach BcVar to arc \e[0m" << std::endl;
#endif
    BcArc a = n->arcs[arcId];
    BcVar v = m->vars[colId];
#ifdef RCSPDEBUG
    std::cout << " > attaching variable \e[1m" << v.name() << "\e[0m with coeff \e[1m " << coeff << " \e[0m to the arc with id \e[1m" << arcId << "\e[0m " << std::endl;
#endif
    //a.arcVar(v);
    a.addVarAssociation(v, coeff);
    return 1;
}

extern "C" int bcRCSP_setVertexConsumptionLB(RCSPNetworkInterface* n, int nodeId, int resId, double lb) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set vertex consumption lb \e[0m" << std::endl;
    std::cout << " > node \e[1m" << nodeId << "\e[0m has consumption lb \e[1m" << lb << "\e[0m for resource \e[1m" << resId << "\e[0m" << std::endl;
#endif
    BcNetworkResource b = n->resources[resId];
    b.setVertexConsumptionLB(n->vertices[nodeId], lb);
    return 1;
}

extern "C" int bcRCSP_setArcConsumptionLB(RCSPNetworkInterface* n, int arcId, int resId, double lb) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set arc consumption lb \e[0m" << std::endl;
    std::cout << " > arc \e[1m" << arcId << "\e[0m has consumption lb \e[1m" << lb << "\e[0m for resource \e[1m" << resId << "\e[0m" << std::endl;
#endif
    BcNetworkResource b = n->resources[resId];
    b.setArcConsumptionLB(n->arcs[arcId], lb);
    return 1;
}

extern "C" int bcRCSP_setVertexSpecialConsumptionLB(RCSPNetworkInterface* n, int nodeId, int resId, double lb)
{
    BcVertex v = n->vertices[nodeId];
    v.setSpecialResourceConsumptionLB(resId, lb);
    return 1;
}


extern "C" int bcRCSP_addBinaryResourceConsumption(RCSPNetworkInterface* n, int nodeId, int binaryResId,
                                                   int consumption, int lowerBound, int upperBound)
{
    BcVertex v = n->vertices[nodeId];
    //v.addBinaryResourceConsumption(binaryResId, consumption, lowerBound, upperBound);
    return 1;
}


extern "C" int bcRCSP_setVertexConsumptionUB(RCSPNetworkInterface* n, int nodeId, int resId, double ub) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set vertex consumption ub \e[0m" << std::endl;
    std::cout << " > node \e[1m" << nodeId << "\e[0m has consumption ub \e[1m" << ub << "\e[0m for resource \e[1m" << resId << " \e[0m" << std::endl;
#endif
    BcNetworkResource b = n->resources[resId];
    b.setVertexConsumptionUB(n->vertices[nodeId], ub);
    return 1;
}

extern "C" int bcRCSP_setArcConsumptionUB(RCSPNetworkInterface* n, int arcId, int resId, double ub) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set arc consumption ub \e[0m" << std::endl;
    std::cout << " > arc \e[1m" << arcId << "\e[0m has consumption ub \e[1m" << ub << "\e[0m for resource \e[1m" << resId << "\e[0m" << std::endl;
#endif
    BcNetworkResource b = n->resources[resId];
    b.setArcConsumptionUB(n->arcs[arcId], ub);
    return 1;
}

extern "C" int bcRCSP_setVertexSpecialConsumptionUB(RCSPNetworkInterface* n, int nodeId, int resId, double ub) {
    BcVertex v = n->vertices[nodeId];
    v.setSpecialResourceConsumptionUB(resId, ub);
    return 1;
}

extern "C" int bcRCSP_setEdgeConsumptionValue(RCSPNetworkInterface* n, int edgeId, int resId, double consumption) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set edge consumption \e[0m" << std::endl;
    std::cout << " > edge \e[1m" << edgeId << "\e[0m has consumption value \e[1m" << consumption << "\e[0m "
             << " for resource id \e[1m" << resId << "\e[0m."
             << std::endl;
#endif
    BcNetworkResource b = n->resources[resId];
    b.setArcConsumption(n->arcs[edgeId], consumption);
    return 1;
}

extern "C" int bcRCSP_setEdgeSpecialConsumptionValue(RCSPNetworkInterface* n, int edgeId, int resId, double consumption)
{
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP set edge consumption \e[0m" << std::endl;
    std::cout << " > edge \e[1m" << edgeId << "\e[0m has consumption value \e[1m" << consumption << "\e[0m "
             << " for special resource id \e[1m" << resId << "\e[0m."
             << std::endl;
#endif
    n->arcs[edgeId].addBinaryResourceConsumption(resId, consumption);
    return 1;
}

extern "C" int bcRCSP_addVertexToMemOfElementaritySet(RCSPNetworkInterface* n, int nodeId, int elementaritySetId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP add vertex to memory of elem set \e[0m" << std::endl;
    std::cout << " > adding node \e[1m" << nodeId << "\e[0m  to the memory of elementarity set \e[1m " << elementaritySetId << "\e[0m. " << std::endl;
#endif
    // to do : check if the node exists & same elem set
    BcVertex v  = n->vertices[nodeId];
    v.addToMemoryOfElemSet(elementaritySetId);
    //v.setElementaritySet(n->elementaritySets[elementaritySetId]);
    return 1;
}

extern "C" int bcRCSP_addEdgeToMemOfElementaritySet(RCSPNetworkInterface* n, int edgeId, int elementaritySetId) {
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP add edge to memory of elem set \e[0m" << std::endl;
    std::cout << " > adding edge \e[1m" << edgeId << "\e[0m to the memory of elementarity set \e[1m" <<
    elementaritySetId << "\e[0m." << std::endl;
#endif
    BcArc a = n->arcs[edgeId];
    a.addToMemoryOfElemSet(elementaritySetId);
    return 1;
}

extern "C" int bcRCSP_addVertexToPackingSet(RCSPNetworkInterface* n, int vertexId, int packingSetId) {
    BcVertex v = n->vertices[vertexId];
    v.setPackingSet(packingSetId);
    return 1;
}

extern "C" int bcRCSP_addEdgeToPackingSet(RCSPNetworkInterface* n, int edgeId, int packingSetId) {
    BcArc a = n->arcs[edgeId];
    a.setPackingSet(packingSetId);
    return 1;
}

extern "C" int bcRCSP_addVertexToCoveringSet(RCSPNetworkInterface* n, int vertexId, int coveringSetId) {
    BcVertex v = n->vertices[vertexId];
    v.setCoveringSet(coveringSetId);
    return 1;
}

extern "C" int bcRCSP_addEdgeToCoveringSet(RCSPNetworkInterface* n, int edgeId, int coveringSetId) {
    BcArc a = n->arcs[edgeId];
    a.setCoveringSet(coveringSetId);
    return 1;
}

extern "C" int bcRCSP_setElemSetsDistanceMatrix(RCSPNetworkInterface* n, double** dists, int num_packsets) {
    std::vector<std::vector<double> > matrix;
    for(int i = 0; i < num_packsets; i++)
    {
	 std::vector<double> row;
	 for(int j = 0; j < num_packsets; j++)
	 {
	      row.push_back(dists[i][j]);
	 }
	 matrix.push_back(row);
    }
    n->network.setElemSetsDistanceMatrix(matrix);
    return 1;
}

extern "C" int bcRCSP_addPackingSetToPackingSetCutNeighbourhood(RCSPNetworkInterface* n, int packingSetId2add,
                                                                int packingSetId)
{
    n->network.addToPackingSetCutNeighbourhood(packingSetId, packingSetId2add);
    return 1;
}

#ifdef NOSTANDALONE
extern "C" int bcRCSP_createOracle(RCSPNetworkInterface* n, InterfaceModel* m, int spType, int* spId)
{
#else
extern "C" int bcRCSP_createOracle(RCSPNetworkInterface* n, InterfaceModel* m, int spType, int* spId,
                                   bool saveStandalone, char* standaloneFileName)
{
#endif
    BcFormulation f = getProblem(m, spType, spId);
#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP create oracle \e[0m" << std::endl;
    std::cout << " > creating oracle for subproblem \e[1m" << f.name() << "\e[0m" << std::endl;
#endif
    BcRCSPFunctor * rcspOracle = new BcRCSPFunctor(f);
#ifndef NOSTANDALONE
    if (saveStandalone) {
        rcspOracle->saveToStandaloneRCSPfile(standaloneFileName);
    }
#endif
    f.attach(rcspOracle);
    return 1;
}

extern "C" int bcRCSP_addGenericCapacityCut(InterfaceModel* m, int maxCapacity, int* demands, int demands_size,
                                            bool isFacultative, double rootPrioLvl, double nonRootPrioLvl,
                                            int twoPathCutResId)
{
    std::vector<int> v_demands;
    for (int i = 0; i < demands_size; i++) {
        v_demands.push_back(demands[i]);
    }

#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP add generic capacity cut \e[0m" << std::endl;
    std::cout << " > max capacity = \e[1m" << maxCapacity << "\e[0m & demand_length = " << demands_size << std::endl;
    std::cout << " > facultative = " << isFacultative << std::endl;
    std::cout << " > rootPrioLvl = " << rootPrioLvl << std::endl;
    std::cout << " > nonRootPrioLvl = " << nonRootPrioLvl << std::endl;
    for (int i = 0; i < demands_size; i++) {
        std::cout << "demand[" << i << "] = " << v_demands[i] << std::endl;
    }
#endif

    BcCapacityCutConstrArray capacityCuts(m->master, maxCapacity, v_demands, isFacultative, true,
                                          twoPathCutResId, rootPrioLvl, nonRootPrioLvl);
    return 1;
}

extern "C" int bcRCSP_addGenericStrongKPathCut(InterfaceModel* m, int maxCapacity, int* demands, int demands_size,
                                               bool isFacultative, double rootPrioLvl, double nonRootPrioLvl)
{
#ifdef USE_NON_PUBLIC_CUTS
    std::vector<int> v_demands;
    for (int i = 0; i < demands_size; i++) {
        v_demands.push_back(demands[i]);
    }

#ifdef RCSPDEBUG
    std::cout << "\e[1;34m @RCSP add generic capacity cut \e[0m" << std::endl;
    std::cout << " > max capacity = \e[1m" << maxCapacity << "\e[0m & demand_length = " << demands_size << std::endl;
    std::cout << " > facultative = " << isFacultative << std::endl;
    std::cout << " > rootPrioLvl = " << rootPrioLvl << std::endl;
    std::cout << " > nonRootPrioLvl = " << nonRootPrioLvl << std::endl;
    for (int i = 0; i < demands_size; i++) {
        std::cout << "demand[" << i << "] = " << v_demands[i] << std::endl;
    }
#endif

    BcStrongKPathCutConstrArray strongKpathCuts(m->master, maxCapacity, v_demands, isFacultative, true,
                                                -1, rootPrioLvl, nonRootPrioLvl);
#else
    std::cout << "BaPCod WARNING: strong K-path cuts are not supported in this library" << std::endl;
#endif
    return 1;
}

extern "C" int bcRCSP_addGenericLimMemOneCut(InterfaceModel* m) {
    BcLimMemRankOneCutConstrArray limMemRank1Cuts(m->master);
    return 1;
}

extern "C" int bcRCSP_addGenericCliqueCut(InterfaceModel* m) {
#ifdef USE_NON_PUBLIC_CUTS
    BcCliqueCutConstrArray limMemRank1Cuts(m->master);
#else
    std::cout << "BaPCod WARNING: clique cuts are not supported in this library" << std::endl;
#endif
    return 1;
}

extern "C" int bcRCSP_addPathsPerNetworkBranching(InterfaceModel* m, double priority) {
    BcPathsPerNetworkBranching(m->master, priority);
    return 1;
}


extern "C" int bcRCSP_addPackSetAssignBranching(InterfaceModel* m, double priority) {
    BcPackSetAssignBranching(m->master, priority);
    return 1;
}

extern "C" int bcRCSP_addPackSetRyanAndFosterBranching(InterfaceModel* m, double priority) {
    BcPackSetRyanFosterBranching(m->master, priority);
    return 1;
}

extern "C" int bcRCSP_addPermanentRyanAndFosterConstraint(RCSPNetworkInterface* n, int firstPackSetId,
                                                          int secondPackSetId, bool together)
{
    n->network.addPermanentRyanAndFosterConstraint(firstPackSetId, secondPackSetId, together);
    return 1;
}

extern "C" int bcRCSP_addElemSetResourceConsumptionBranching(InterfaceModel* m, double priority) {
    BcPackSetResConsumptionBranching(m->master, priority);
    return 1;
}

extern "C" int bcRCSP_addAssociatedVarToResource(RCSPNetworkInterface* n, int resId, InterfaceModel* m, int colId)
{
#ifdef RCSPDEBUG
    std::cout << "\e[1;31m associate variable " << m->vars[colId].name() << " to resource " << resId << "\e[00m" << std::endl;
#endif
    BcNetworkResource r = n->resources[resId];
    r.setAssociatedVar(m->vars[colId]);
    return 1;
}

#endif // BCP_RCSP_IS_FOUND
