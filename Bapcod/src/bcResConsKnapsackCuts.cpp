/**
 *
 * This file bcNonPublicResConsKnapsackCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifdef BCP_RCSP_IS_FOUND
#include "bcUsefulHeadFil.hpp"
#include "bcNetworkBasedCuts.hpp"

#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelingLanguageC.hpp"
#include "bcModelC.hpp"
#include "bcModelRCSPSolver.hpp"
#include "rcsp_interface.hpp"

ResConsConstrInfo::~ResConsConstrInfo()
{
    delete modelPtr;
    delete bcInitPtr;
}

/***************************************************************************
 ***********   TODO: Methods for ResConsKnapsackCut   **********************
 ***************************************************************************/

ResConsKnapsackCut::ResConsKnapsackCut(const IndexCell& id,
                                       GenericResConsKnapsackCutConstr * genConstrPtr,
                                       ProbConfig * probConfigPtr,
                                       const std::string & name,
                                       const Double & rhs,
                                       const bcp_rcsp::RouteLoadKnapsackCut * rcspCutPtr):
        InstMasterConstr(id, genConstrPtr, probConfigPtr, name, rhs, 'L',  genConstrPtr->defaultType(),
                         genConstrPtr->defaultKind(), genConstrPtr->defaultFlag()),
        Base4NonLinearConstraint(), _genResConsKnapCutConstr(genConstrPtr), _rcspCutPtr(rcspCutPtr)
{
}

bool ResConsKnapsackCut::isRelatedTo(const ColGenSpConf * cgSpConfPtr)
{
    if ((cgSpConfPtr == nullptr) || (cgSpConfPtr->rcspGraphPtr() == nullptr))
        return false;

    int graphId = cgSpConfPtr->rcspGraphPtr()->id;

    return (_rcspCutPtr->resVariableMap.find(graphId) != _rcspCutPtr->resVariableMap.end());
}


ResConsKnapsackCut::~ResConsKnapsackCut()
{
}

bool ResConsKnapsackCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::ResConsKnapsackCutConstrMask, vcIdentifier);
}

void ResConsKnapsackCut::nicePrint(std::ostream& os) const
{
  os << "Res.cons.knapsack cut " << name()  << ":" ;
  for (ConstVarConstrPtr2Double::const_iterator mapIt = member2coefMap().begin();
       mapIt != member2coefMap().end(); ++mapIt)
  {
      if (!mapIt->first->isTypeOf(VcId::MastColumnMask))
      {
          if ((mapIt != member2coefMap().begin()) && (mapIt->second > 0))
              os << "+";
          os << mapIt->second << "*" << mapIt->first->name();
      }
  }
  os << " <= " << costrhs() << "  ";

  _rcspCutPtr->nicePrint(os);
  os << std::endl;
}

void ResConsKnapsackCut::setMembership()
{
  if(!buildMembershipHasBeenPerformed())
  {
    genVarConstrPtr()->buildMembership(this);
    buildMembershipHasBeenPerformed(true);
  }

  bool cumulativeCoef(false);

  VarIndexManager::const_iterator it;
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Active, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Active, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _genResConsKnapCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _genResConsKnapCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  /// if column pool is not used, we do not generate membership of the constraint
  /// in and unsuitable column, as generated column will be active
  /// only at the node it was generated and in the subtree rooted at this node
  if (param().UseColumnsPool())
    for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
         it != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); ++it)
      if ((*it)->isTypeOf(VcId::MastColumnMask))
      {
        LpCoef lpCoeff = _genResConsKnapCutConstr->genericMastColumnCoef(this, static_cast<MastColumn *>(*it));
        if (lpCoeff.first)
          includeMember(*it, lpCoeff.second, cumulativeCoef);
      }

  Constraint::setMembership();

  return;
}

/***************************************************************************
 *******   TODO: Methods for GenericResConsKnapsackCutConstr   *************
 ***************************************************************************/

GenericResConsKnapsackCutConstr::GenericResConsKnapsackCutConstr(Model * modelPtr,
                                                                 ProbConfig * probConfPtr,
                                                                 const std::string & name,
                                                                 const Double & nonRootPriorityLevel,
                                                                 const Double & rootPriorityLevel) :
    GenericCutConstr(modelPtr, probConfPtr, name, 'F', SelectionStrategy::MostFractional,
                     nonRootPriorityLevel, rootPriorityLevel, false),
    Base4NonLinearGenericConstr(NULL), _constrInfos(), _interfacePtr(nullptr),
    _useFenchelSeparation(false)
{
    //_useFenchelSeparation = true;
}

GenericResConsKnapsackCutConstr::~GenericResConsKnapsackCutConstr()
{
    _constrInfos.clear();
    delete _interfacePtr;
}

bool GenericResConsKnapsackCutConstr::prepareSeparation()
{
    std::map<InstMasterConstr *, int> constrPtrToIdMap;

    std::set<const SubProbVariable *> resAssociatedVars;
    for (std::vector< ColGenSpConf * >::const_iterator cgspPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgspPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgspPtrIt)
    {
        const NetworkFlow * networkPtr = (*cgspPtrIt)->networkFlowPtr();
        if (networkPtr != NULL)
        {
            for (std::list<ScalableResource *>::const_iterator resPtrIt = networkPtr->sideResourceConstrList().begin();
                 resPtrIt != networkPtr->sideResourceConstrList().end(); ++resPtrIt)
            {
                InstanciatedVar * varPtr = (*resPtrIt)->associatedVar();
                if (varPtr == NULL)
                    continue;
                if (!varPtr->isTypeOf(VcId::SubProbVariableMask) || (varPtr->probConfPtr() != *cgspPtrIt))
                {
                    std::cerr << "ResConsKnapsack cuts separator error : variable associated to resource with id = "
                              << (*resPtrIt)->id() << " of the network of subproblem " << (*cgspPtrIt)->genericName()
                              << " is not from this subproblem" << std::endl;
                    return false;
                }
                const SubProbVariable * spVarPtr = static_cast<const SubProbVariable *>(varPtr);
                resAssociatedVars.insert(spVarPtr);
                for (ConstVarConstrPtr2Double::const_iterator mapIt = spVarPtr->masterConstrMember2coefMap().begin();
                     mapIt != spVarPtr->masterConstrMember2coefMap().end(); ++mapIt)
                {
                    InstMasterConstr * imConstrPtr = static_cast<InstMasterConstr *>(mapIt->first);
                    int infoId = -1;
                    auto idMapIt = constrPtrToIdMap.find(imConstrPtr);
                    if (idMapIt == constrPtrToIdMap.end())
                    {
                        infoId = (int)_constrInfos.size();
                        _constrInfos.push_back(ResConsConstrInfo(infoId, imConstrPtr));
                    }
                    else
                    {
                        infoId = idMapIt->second;
                    }
                    std::vector<ResConsConstrLhsItem> & lhsInfo = _constrInfos[infoId].leftHandSideInfo;
                    lhsInfo.push_back(ResConsConstrLhsItem(*cgspPtrIt, (*resPtrIt)->id(), varPtr, mapIt->second));
                }
            }
        }
    }

    bcp_rcsp::ResourceVarMapVector resourceVarMaps;

    /// we now verify that all constraints which contain variables associated to resources are "good" for separation
    for (auto & info : _constrInfos)
    {
       InstMasterConstr * mastConstrPtr = info.constrPtr;
       if (mastConstrPtr->sense() == 'G')
       {
           std::cerr << "ResConsKnapsack cuts separator error : master constraint " << mastConstrPtr->name()
                     << " containing a variable associated to a resource has a wrong sense "
                     << "(should be either <= or ==)" << std::endl;
           return false;
       }
       for (MapSubProbVariablePtr2Double::const_iterator spMapIt = mastConstrPtr->subProbVarMember2coefMap().begin();
            spMapIt != mastConstrPtr->subProbVarMember2coefMap().end(); ++spMapIt)
       {
           const SubProbVariable * spVarPtr = spMapIt->first;
           const Double & spVarCoeff = spMapIt->second;
           if (resAssociatedVars.find(spVarPtr) == resAssociatedVars.end())
           {
               std::cerr << "ResConsKnapsack cuts separator error : master constraint " << mastConstrPtr->name()
                         << " contains subproblem variables both associated and not to a resource" << std::endl;
               return false;
           }
           if (spVarCoeff < -Double::precision)
           {
               std::cerr << "ResConsKnapsack cuts separator error : variable " << spVarPtr->name()
                         << " associated with a resource has negative coefficient in constraint "
                         << mastConstrPtr->name() << std::endl;
               return false;
           }
       }
        for (ConstVarConstrPtr2Double::iterator mvMapIt = mastConstrPtr->member2coefMap().begin();
             mvMapIt != mastConstrPtr->member2coefMap().end(); ++mvMapIt)
        {
            if (mvMapIt->first->type() == 'C')
            {
                std::cerr << "ResConsKnapsack cuts separator error : master constraint " << mastConstrPtr->name()
                          << " having a variable associated with a resource contains pure master variable "
                          << mvMapIt->first->name() << " which is continuous" << std::endl;
                return false;
            }
            if ((mvMapIt->second < -Double::precision) && (mvMapIt->first->lb() < -Double::precision)) {
                std::cerr << "ResConsKnapsack cuts separator error : master constraint " << mastConstrPtr->name()
                          << " have a variable associated with a resource contains pure master variable "
                          << mvMapIt->first->name() << " with negative coefficient and negative lower bound"
                          << std::endl;
                return false;
            }
            InstanciatedVar * ivPtr = static_cast<InstanciatedVar *>(mvMapIt->first);
            info.pureMastVarInfo.insert(std::make_pair(ivPtr, mvMapIt->second));
        }
        if (!mastConstrPtr->costrhs().isZero())
        {
            if (mastConstrPtr->costrhs() < 0)
            {
                std::cerr << "ResConsKnapsack cuts separator error : master constraint " << mastConstrPtr->name()
                          << " containing a variable associated with a resource has negative free coefficient"
                          << std::endl;
                return false;
            }
            info.freeCoeff = mastConstrPtr->costrhs();
        }

        resourceVarMaps.push_back(std::map<const bcp_rcsp::GraphData *, std::pair<int, double> >());
        for (auto & lhsInfo : info.leftHandSideInfo)
            resourceVarMaps.back().insert(std::make_pair(lhsInfo.colGenSpConfPtr->rcspGraphPtr(),
                                                         std::make_pair(lhsInfo.resId, lhsInfo.coeff)));
    }

    bcp_rcsp::RouteLoadKnapsackCutsSeparatorParameters sepParams;
    sepParams.maxNumPerRound = param().RCSPresConsKnapsackCutsMaxNumPerRound();
    sepParams.cutViolationTolerance = param().BapCodCutViolationTolerance();
    sepParams.printLevel = param().RCSPprintLevel();
    sepParams.oneKSeparation = (param().RCSPresConsKnapsackCutsMode() > 0) ? param().RCSPresConsKnapsackCutsMode()
                                                                           : -param().RCSPresConsKnapsackCutsMode();
    sepParams.separateByRounding = (param().RCSPresConsKnapsackCutsMode() > 0);

    _interfacePtr = bcp_rcsp::createAndPrepareRouteLoadKnapsackCutSeparation(resourceVarMaps, sepParams);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "ResConsKnapsack cuts separator error : cannot create separation interface " << std::endl;
        return false;
    }

    if (_useFenchelSeparation)
      return prepareFenchelSeparation();

    return true;
}

/// experimental Fenchel separation by column generation
bool GenericResConsKnapsackCutConstr::prepareFenchelSeparation()
{
    for (auto & info : _constrInfos)
    {
        double curFreeCoeff = info.freeCoeff;
        for (auto &pair: info.pureMastVarInfo)
        {
            double curUb = pair.first->curUb();
            double varCoeff = pair.second;
            if ((varCoeff < 0.0) && (curUb < 1 + Double::precision))
                curFreeCoeff = -varCoeff;
        }
        int knapsackCapcity = (int) ceil(curFreeCoeff - Double::precision);

        info.bcInitPtr = new BapcodInit();
        info.bcInitPtr->param().DEFAULTPRINTLEVEL(-2);
        info.bcInitPtr->param().MaxNbOfBBtreeNodeTreated(1);
        info.bcInitPtr->param().masterSolMode(SolutionMethod::lpSolver);
        info.bcInitPtr->param().MipSolverMultiThread(1);
        info.modelPtr = new Model(info.bcInitPtr);
        BcModel model(info.modelPtr);
        BcObjective objective(model);
        objective.setArtCostValue(1e+6);
        BcMaster master(model);

        BcConstrArray covConstr(master, "COV");
        for (int weight = 1; weight <= knapsackCapcity; ++weight)
            covConstr(weight) >= 0;

        BcColGenSpArray knapsackCGSp(model);
        knapsackCGSp.setFixedCost(1.0);
        knapsackCGSp(0);
        BcConstrArray knapConstr(knapsackCGSp[0], "KNP");
        knapConstr(0) <= knapsackCapcity;

        BcVarArray weightVar(knapsackCGSp[0], "W");
        weightVar.type('I');
        for (int weight = 1; weight <= knapsackCapcity; ++weight)
        {
            weightVar(weight);
            covConstr(weight) += 1 * weightVar(weight);
            knapConstr(0) += weight * weightVar(weight);
        }
    }

    return true;
}

void GenericResConsKnapsackCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

void GenericResConsKnapsackCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset < InstanciatedConstr *, CutSeparationPriorityComp > & generatedCutConstrSet)
{
    bcp_rcsp::FractionalMasterSolution rcspFracSolution;
    rcspFracSolution.solPts.reserve(curSol.size());
    rcspFracSolution.values.reserve(curSol.size());

    for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    {
        if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;

        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);

        auto rcspSolPtr = colPtr->spSol()->rcspSolPtr();
        if (rcspSolPtr != nullptr)
        {
            rcspFracSolution.solPts.push_back(rcspSolPtr);
            rcspFracSolution.values.push_back(colPtr->val());
        }
    }

    std::map<int, std::pair<std::map<int, double>, int> > infosForSepInterface;
    std::map<int, InstanciatedVar *> nonFixedBinaryVarPts;

    for (auto & info : _constrInfos)
    {
        /// we now save all columns in the rcspFracSolution with non-zero coefficient in the constraint
        std::vector<std::pair<MastColumn *, double> > colsInSolMap;
        std::vector<std::vector<double> >::const_reverse_iterator colResConsIt;

        std::map<int, double> pureMasterCoeffMap; /// contribution of pure master variables to the master knapsack
        /// coefficients

        /// a valid cut can added only if
        /// - there is exactly one non-fixed pure master binary variable with negative coeff and all other pure master
        ///   integer variables with negative coeff are fixed to zero and free coefficient is zero
        /// - all pure master integer variable with negative coefficient are fixed
        int numNonFixedBinaryVars = 0;
        InstanciatedVar *nonFixedBinaryVarPtr = nullptr;
        bool someIntegerVarsAreNotFixed = false;
        double curFreeCoeff = info.freeCoeff;
        bool nonZeroFreeCoeff = (curFreeCoeff > Double::precision); /// freeCoeff cannot be negative
        for (auto &pair: info.pureMastVarInfo)
        {
            InstanciatedVar *iVarPtr = pair.first;
            double varCoeff = pair.second;
            double curUb = iVarPtr->curUb();
            double curLb = iVarPtr->curLb();
            if (varCoeff > 0)
            {
                if (!iVarPtr->val().isZero())
                {
                    int roundedDownCoeff = floor(varCoeff + Double::precision);
                    auto mapIt = pureMasterCoeffMap.find(roundedDownCoeff);
                    if (mapIt != pureMasterCoeffMap.end())
                        mapIt->second += iVarPtr->val();
                    else
                        pureMasterCoeffMap.insert(std::make_pair(roundedDownCoeff, iVarPtr->val()));
                }
            }
            else
            {
                if (curUb - curLb < Double::precision)
                {
                    if (curUb > 1 - Double::precision)
                    {
                        nonZeroFreeCoeff = true;
                        curFreeCoeff += (curUb * -varCoeff);
                    }
                }
                else
                {
                    if (curUb < 1 + Double::precision)
                    {
                        numNonFixedBinaryVars += 1;
                        nonFixedBinaryVarPtr = iVarPtr;
                        curFreeCoeff = -varCoeff;
                    }
                    else
                    {
                        someIntegerVarsAreNotFixed = true;
                    }
                }
            }
        }
        if ((numNonFixedBinaryVars > 1) || someIntegerVarsAreNotFixed
            || ((numNonFixedBinaryVars == 1) && nonZeroFreeCoeff))
            continue;

        int rightHandSideValue = (int) ceil(curFreeCoeff - Double::precision);

        infosForSepInterface.insert(std::make_pair(info.id,
                                                   std::make_pair(pureMasterCoeffMap, rightHandSideValue)));
        nonFixedBinaryVarPts.insert(std::make_pair(info.id, nonFixedBinaryVarPtr));
    }

    if (infosForSepInterface.empty())
        return;

    std::vector<const bcp_rcsp::RouteLoadKnapsackCut *> rlkCutPts;
    if (!_interfacePtr->separate(bcp_rcsp::RouteLoadKnapsackCutSeparationInput(rcspFracSolution, infosForSepInterface),
                                 rlkCutPts))
    {
        std::cerr << "BaPCod error : separation of res. cons. knapsack cuts did not succeed " << std::endl;
        exit(1);
    }

    for (auto & pair : infosForSepInterface)
    {
        int constrId = pair.first;
        ResConsConstrInfo & info = _constrInfos[constrId];

        if (info.modelPtr == nullptr)
            continue;

        std::map<int, double> coeffMap;
        for (auto & thisPair: pair.second.first)
            coeffMap[thisPair.first] = thisPair.second;

        int numRcspSols = (int) rcspFracSolution.solPts.size();
        for (int rcspSolId = 0; rcspSolId < numRcspSols; ++rcspSolId)
        {
            const auto pathPtr = rcspFracSolution.solPts[rcspSolId];
            const auto value = rcspFracSolution.values[rcspSolId];

            if (pathPtr->graphId != constrId + 1) /// this is not generic, TO DO it generic
                continue;

            int resId = 0; /// again this is not generic
            int pathAccResCons = floor(pathPtr->resConsumption.back()[resId] + Double::precision);
            auto ret = coeffMap.insert(std::pair<int, double>(pathAccResCons, value));
            if (!ret.second)
                coeffMap[pathAccResCons] += value;
        }

        BcModel model(info.modelPtr);
        BcMaster master(model);
        master.resetConstrRHS();
        BcConstrArray covConstr(master, "COV");
        std::cout << "Starting Fenchel separation, coeffs are";
        for (auto & mapPair : coeffMap)
            std::cout << " " << mapPair.first << "=>" << mapPair.second;
        int freeCoeff = (int) ceil(pair.second.second - Double::precision);
        std::cout << ", free coeff. is " << freeCoeff << std::endl;
        for (auto & mapPair : coeffMap)
            covConstr[mapPair.first] >= mapPair.second;

        PrintLevel::setPrintLevel(-2);
        info.bcInitPtr->reset();
        master.update();
        model.solve();
        PrintLevel::setPrintLevel(param().DEFAULTPRINTLEVEL());

        double masterValue = info.bcInitPtr->statistics().getValue("bcRecRootLpVal");

        if (masterValue > 1 + 1e-6)
        {
            std::cout << "Found a RLKC for constraint " << constrId << " with violation " << masterValue - 1
                      << " by Fenchel separation, coeffs are";
            double prevValue = 0.0;
            for (int weight = 1; weight <= freeCoeff; ++weight)
            {
                double value = -((BcConstr)covConstr[weight]).curDualVal();
                if (value > prevValue + Double::precision)
                {
                    prevValue = value;
                    std::cout << " " << weight << "=>" << value;
                }
            }
            std::cout << std::endl;
        }
    }

    for (auto rcspCutPtr : rlkCutPts)
    {
        MultiIndex newCutId(rcspCutPtr->id);
        std::string newCutName("RCK");
        auto mapIt = nonFixedBinaryVarPts.find(rcspCutPtr->origConstrId);
        InstanciatedVar * nonFixedBinaryVarPtr = (mapIt != nonFixedBinaryVarPts.end()) ? mapIt->second : nullptr;
        double rhsValue = (nonFixedBinaryVarPtr == nullptr) ?  rcspCutPtr->rightHandSide : 0.0;

        ResConsKnapsackCut * newCutPtr
                = new ResConsKnapsackCut(newCutId, this, probConfPtr(),
                                         newCutId.appendRef2name(newCutName, multiIndexNames()), rhsValue, rcspCutPtr);

        /// add add master variables to the cut
        for (auto & pair : _constrInfos[rcspCutPtr->origConstrId].pureMastVarInfo)
        {
            InstanciatedVar * varPtr = pair.first;
            double varCoeff = pair.second;
            if (varCoeff < 0)
            {
                if (varPtr == nonFixedBinaryVarPtr)
                    newCutPtr->includeMember(nonFixedBinaryVarPtr, -rcspCutPtr->rightHandSide, false);
            }
            else
            {
                // stopped here
                int roundedDownCoeff = (int)floor(varCoeff + Double::precision);
                auto coeffMapIt = rcspCutPtr->coeffMap.upper_bound(roundedDownCoeff);
                --coeffMapIt; /// we can always do it as there is always the key 0 in the map
                newCutPtr->includeMember(varPtr, coeffMapIt->second, false);
            }
        }
        generatedCutConstrSet.insert(newCutPtr);
    }
}

const LpCoef GenericResConsKnapsackCutConstr::genericMastColumnCoef(InstanciatedConstr * icPtr,
                                                                    MastColumn * colPtr) const
{
    if (!icPtr->isTypeOf(VcId::ResConsKnapsackCutConstrMask))
        return LpCoef(false, 0.0);

    const bcp_rcsp::RouteLoadKnapsackCut * rcspCutPtr = static_cast<ResConsKnapsackCut *>(icPtr)->rcspCutPtr();
    double coeff = _interfacePtr->coefficient(colPtr->spSol()->rcspSolPtr(), rcspCutPtr);
    if (coeff > 0.0)
        return LpCoef(true, coeff);
    return LpCoef(false, 0.0);
}


#endif //BCP_RCSP_IS_FOUND
