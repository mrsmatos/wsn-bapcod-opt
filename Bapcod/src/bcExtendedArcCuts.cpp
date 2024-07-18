/**
 *
 * This file bcExtendedArcCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcGenNonRobustCutsC.hpp
//  Project
//
//  Created by Ruslan Sadykov on 16/02/2015.
//
//

#ifdef BCP_RCSP_IS_FOUND
#include "bcUsefulHeadFil.hpp"
#include "bcNetworkBasedCuts.hpp"

#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "bcModelCutConstrC.hpp"
#include "bcModelRCSPSolver.hpp"
#include "bcModelMasterC.hpp"

#include "rcsp_data.hpp"

int GenericExtendedArcCutConstr::_numGeneratedCuts = 0;

/***************************************************************************
 *************   TODO: Methods for ExtendedArcCut   ************************
 ***************************************************************************/
 
ExtendedArcCut::ExtendedArcCut(GenericExtendedArcCutConstr * genConstrPtr,
                               ProbConfig * probConfigPtr,
                               const std::string & name,
                               const Double & rhs,
                               const char & sense,
                               const BcCustomExtendedArcCutInfo * cutInfoPtr) :
  InstMasterConstr(MultiIndex(genConstrPtr->_numGeneratedCuts), genConstrPtr,
                   probConfigPtr, std::string(name) + "_n" + (genConstrPtr->_numGeneratedCuts++),
                   rhs, sense, genConstrPtr->defaultType(), genConstrPtr->defaultKind(),
                   genConstrPtr->defaultFlag()),
  Base4NonLinearConstraint(), _genExtendedArcCutConstr(genConstrPtr), _cutInfoPtr(cutInfoPtr)
{
}

ExtendedArcCut::~ExtendedArcCut()
{
}

void ExtendedArcCut::nicePrint(std::ostream& os)
{
  if (_cutInfoPtr != NULL)
    _cutInfoPtr->nicePrint(os);
}

double ExtendedArcCut::getArcCoefficient(const int & tailVertId, const int & headVertId,
                                         const double * tailResCons) const
{
  std::cerr << "ExtendedArcCut::nicePrint() should not be called, use classes derived from ExtendedArcCut" << std::endl;
  exit(1);
  return 0.0;
}

double ExtendedArcCut::getRouteCoefficient(const std::vector<int> & arcVertIds,
                                           const std::vector<std::vector<double> > & routeResCons) const
{
  std::cerr << "ExtendedArcCut::nicePrint() should not be called, use classes derived from ExtendedArcCut" << std::endl;
  exit(1);
  return 0.0;
}

/// attention! supposes that _cutInfoPtr is not NULL and also netFlowPtr->getArcInfoPtr(arcId) is not NULL
double ExtendedArcCut::getArcCoefficient(const NetworkFlow * netFlowPtr, const int & arcId,
                                         const double * resCons, const bool & isTailResCons) const
{
  return _cutInfoPtr->getArcCoefficient(netFlowPtr->getArcInfoPtr(arcId), resCons, isTailResCons);
}

/// attention! supposes that _cutInfoPtr is not NULL and also netFlowPtr->getArcInfoPtr(arcId) is not NULL
double ExtendedArcCut::getRouteCoefficient(const NetworkFlow * netFlowPtr, const std::vector<int> & arcIds,
                                          const std::vector<std::vector<double> > & routeResCons) const
{
  std::vector<const BcArcInfo *> arcInfoPts;
  arcInfoPts.reserve(arcIds.size());
  for (std::vector<int>::const_iterator arcIdIt = arcIds.begin(); arcIdIt != arcIds.end(); ++arcIdIt)
    arcInfoPts.push_back(netFlowPtr->getArcInfoPtr(*arcIdIt));
  return _cutInfoPtr->getRouteCoefficient(arcInfoPts, routeResCons);
}

void ExtendedArcCut::setMembership()
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
        LpCoef lpCoeff = _genExtendedArcCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
        if (lpCoeff.first)
          includeMember(*it, lpCoeff.second, cumulativeCoef);
      }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
      {
        LpCoef lpCoeff = _genExtendedArcCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
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
          LpCoef lpCoeff = _genExtendedArcCutConstr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
          if (lpCoeff.first)
            includeMember(*it, lpCoeff.second, cumulativeCoef);
        }

  Constraint::setMembership();

  return;
}

bool ExtendedArcCut::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::ExtendedArcCutConstrMask, vcIdentifier);
}

/***************************************************************************
 *************   TODO: Methods for GenericExtendedArcCutConstr   ***********
 ***************************************************************************/

GenericExtendedArcCutConstr::GenericExtendedArcCutConstr(Model * modelPtr,
                                                         ProbConfig * probConfPtr,
                                                         const std::string & name,
                                                         const Double & nonRootPriorityLevel,
                                                         const Double & rootPriorityLevel):
  GenericCutConstr(modelPtr, probConfPtr, name, 'F', SelectionStrategy::MostFractional,
                   nonRootPriorityLevel, rootPriorityLevel, false),
  Base4NonLinearGenericConstr(NULL), _separationFunctorPtr(NULL), _maxGraphId(0), _graphPts()
{
}

GenericExtendedArcCutConstr::~GenericExtendedArcCutConstr()
{
  if (_separationFunctorPtr != NULL)
    delete _separationFunctorPtr;
}

bool GenericExtendedArcCutConstr::prepareSeparation()
{
    std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt;
    _maxGraphId = 0;
    for (cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() == nullptr)
            continue;
        _maxGraphId = (std::max)(_maxGraphId, (*cgSpConfPtrIt)->rcspGraphPtr()->id);
    }

    _graphPts.resize(_maxGraphId + 1, nullptr);

    for (cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() == nullptr)
            continue;

        const bcp_rcsp::GraphData * graphPtr = (*cgSpConfPtrIt)->rcspGraphPtr();

        _graphPts[graphPtr->id] = graphPtr;
    }
    return true;
}

void GenericExtendedArcCutConstr::setSeparationFunctor(BcCustomExtendedArcCutSeparationFunctor * separationFunctorPtr)
{
  _separationFunctorPtr = separationFunctorPtr;
}

void GenericExtendedArcCutConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

const LpCoef GenericExtendedArcCutConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (!icPtr->isTypeOf(VcId::ExtendedArcCutConstrMask))
    return LpCoef(false, 0.0);
    
  return getMastColumnCoeff(static_cast<ExtendedArcCut *>(icPtr), colPtr);
}

void GenericExtendedArcCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset<InstanciatedConstr *, CutSeparationPriorityComp> & generatedCutConstrSet)
{
  if ((probConfPtr() == NULL) || (_separationFunctorPtr == NULL))
    return;
  
  Solution * solPtr = probConfPtr()->getSolution(curSol);
  Solution * primalSolPtr = probConfPtr()->getAggregatedSolution(solPtr);
    
  if (printL(5))
    std::cout << "GenericCustomNonLinearCutConstr::cutSeparationRoutine: primalSol "
              << primalSolPtr << std::endl;

  std::list<BcConstr> cutList;

  BcSolution projectedSol(primalSolPtr);

  std::list<std::pair<double, BcSolution> > columnsInSol;
  for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
    if ((*varPtrIt)->isTypeOf(VcId::MastColumnMask))
      {
        MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
        columnsInSol.push_back(std::make_pair(colPtr->val(), BcSolution(colPtr->spSol())));
      }

  int counter = _separationFunctorPtr->cutSeparationRoutine(BcFormulation(modelPtr()->masterConfPtr()),
                                                            projectedSol, columnsInSol,
                                                            param().BapCodCutViolationTolerance(),
                                                            cutList);
  delete solPtr;
  primalSolPtr->deleteSolutionsChain();
  delete primalSolPtr;
    
  if (printL(5))
    std::cout << "GenericCustomNonLinearCutConstr::cutSeparationRoutine: generated CutConstraint "
              << counter << std::endl;

  if (counter >= 1)
    {
      for (std::list<BcConstr>::iterator cutIt = cutList.begin(); cutIt != cutList.end(); ++cutIt)
        {
          if (printL(5))
            std::cout << "CutConstraint " << ((InstanciatedConstr *) *cutIt)  << std::endl;
          generatedCutConstrSet.insert((InstanciatedConstr *) *cutIt);
        }
    }
}

const LpCoef GenericExtendedArcCutConstr::getMastColumnCoeff(ExtendedArcCut * cutPtr, MastColumn * colPtr) const
{
  if ((colPtr->spSol() == nullptr) || (colPtr->spSol()->rcspSolPtr() == nullptr))
    return LpCoef(false, 0.0);

  double coeff = 0.0;
  if (_separationFunctorPtr == NULL)
    {
        int graphId = colPtr->spSol()->rcspSolPtr()->graphId;
        std::vector<int> vertexIds;

        if ((_graphPts[graphId] == nullptr)
            || !_graphPts[graphId]->obtainVertexIds(colPtr->spSol()->rcspSolPtr(), vertexIds))
            return LpCoef(false, 0.0);

        const std::vector<std::vector<double> > & resCons = colPtr->spSol()->rcspSolPtr()->resConsumption;

      if (colPtr->spSol()->enumeratedFlag())
        {
          coeff = cutPtr->getRouteCoefficient(vertexIds, resCons);
        }
      else
        {
          std::vector<int>::const_iterator nextVertIt, vertIt = vertexIds.begin();
          std::vector<std::vector<double> >::const_iterator resConsIt = resCons.begin();
          if (vertIt != vertexIds.end())
            {
              nextVertIt = vertIt;
              ++nextVertIt;
              while (nextVertIt != vertexIds.end())
                {
                  coeff += cutPtr->getArcCoefficient(*vertIt, *nextVertIt, &(*resConsIt)[0]);
                  ++vertIt;
                  ++nextVertIt;
                  ++resConsIt;
                }
            }
        }
    }
  else if (cutPtr->cutInfoPtr() != NULL) /// custom extended cuts implementation
    {
      const NetworkFlow * netFlowPtr = colPtr->cgSpConfPtr()->networkFlowPtr();
      const std::vector<int> & arcIds = colPtr->spSol()->rcspSolPtr()->arcIds;
      const std::vector<std::vector<double> > & resCons = colPtr->spSol()->rcspSolPtr()->resConsumption;
      if (colPtr->spSol()->enumeratedFlag())
        {
          coeff = cutPtr->getRouteCoefficient(netFlowPtr, arcIds, resCons);
        }
      else
        {
          std::vector<std::vector<double> >::const_iterator resConsIt = resCons.begin();
          std::vector<int>::const_iterator arcIdIt = arcIds.begin();
          while (arcIdIt != arcIds.end())
            {
              coeff += cutPtr->getArcCoefficient(netFlowPtr, *arcIdIt, &((*resConsIt)[0]), true);
              ++arcIdIt;
              ++resConsIt;
            }
        }
    }
  bool nonZeroCoeff = (coeff != 0.0);
  return LpCoef(nonZeroCoeff, coeff);
}
#endif /* BCP_RCSP_IS_FOUND */
