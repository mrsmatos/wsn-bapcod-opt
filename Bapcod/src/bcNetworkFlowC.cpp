/**
 *
 * This file bcNetworkFlowC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcNetworkFlowC.hpp"
#include "bcModelNetworkFlow.hpp"

void NetworkVertex::setElementaritySet(int elemSetId)
{
    _elementaritySetPts.clear();
    _enumElementaritySetPts.clear();
    addToElementaritySet(elemSetId);
}

void NetworkVertex::addToElementaritySet(int elemSetId)
{
    if (elemSetId < _networkPtr->elementaritySetPts().size())
        _elementaritySetPts.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
}

void NetworkVertex::addToEnumerationOnlyElementaritySet(int elemSetId)
{
    if (elemSetId < _networkPtr->elementaritySetPts().size())
        _enumElementaritySetPts.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
}

void NetworkVertex::setPackingSet(int packSetId)
{
    _packingSetPts.clear();
    addToPackingSet(packSetId);
}

void NetworkVertex::addToPackingSet(int packSetId)
{
  if (packSetId < _networkPtr->packingSetPts().size())
    _packingSetPts.push_back(_networkPtr->packingSetPts()[packSetId]);
}

void NetworkVertex::setCoveringSet(int covSetId)
{
    _coveringSetPts.clear();
    addToCoveringSet(covSetId);
}

void NetworkVertex::addToCoveringSet(int covSetId)
{
  if (covSetId < _networkPtr->coveringSetPts().size())
    _coveringSetPts.push_back(_networkPtr->coveringSetPts()[covSetId]);
}

void NetworkVertex::addToMemoryOfElemSet(const int elemSetId)
{
  if ((elemSetId >= 0) && (elemSetId < _networkPtr->elementaritySetPts().size()))
    {
      _inMemoryOfElemSets.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
      _networkPtr->elementaritySetPts()[elemSetId]->addToMemory(this);
    }
}

void NetworkVertex::setSpecialResourceConsumptionLB(int binaryResId, int lowerBound)
{
  if (_specResConsumptionBounds.find(binaryResId) == _specResConsumptionBounds.end())
      _specResConsumptionBounds[binaryResId] = std::make_pair(lowerBound, 1);
  else
      _specResConsumptionBounds[binaryResId].first = lowerBound;
}

void NetworkVertex::setSpecialResourceConsumptionUB(int binaryResId, int upperBound)
{
    if (_specResConsumptionBounds.find(binaryResId) == _specResConsumptionBounds.end())
        _specResConsumptionBounds[binaryResId] = std::make_pair(0, upperBound);
    else
        _specResConsumptionBounds[binaryResId].second = upperBound;
}


NetworkVertex * NetworkArc::tailVertexPtr() const
{
  return _networkPtr->netVertexPtr(_digraph.source(_lemonArc));
}

NetworkVertex * NetworkArc::headVertexPtr() const
{
  return _networkPtr->netVertexPtr(_digraph.target(_lemonArc));
}

void NetworkArc::setElementaritySet(int elemSetId)
{
    _elementaritySetPts.clear();
    addToElementaritySet(elemSetId);
}

void NetworkArc::addToElementaritySet(int elemSetId)
{
    if (elemSetId < _networkPtr->elementaritySetPts().size())
        _elementaritySetPts.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
}

void NetworkArc::addToEnumerationOnlyElementaritySet(int elemSetId)
{
    if (elemSetId < _networkPtr->elementaritySetPts().size())
        _enumElementaritySetPts.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
}

void NetworkArc::setPackingSet(int packSetId)
{
    _packingSetPts.clear();
    addToPackingSet(packSetId);
}

void NetworkArc::addToPackingSet(int packSetId)
{
    if (packSetId < _networkPtr->packingSetPts().size())
        _packingSetPts.push_back(_networkPtr->packingSetPts()[packSetId]);
}

void NetworkArc::setCoveringSet(int covSetId)
{
    _coveringSetPts.clear();
    addToCoveringSet(covSetId);
}

void NetworkArc::addToCoveringSet(int covSetId)
{
    if (covSetId < _networkPtr->coveringSetPts().size())
        _coveringSetPts.push_back(_networkPtr->coveringSetPts()[covSetId]);
}

void NetworkArc::name(std::string name_)
{
    _name = std::move(name_);
}


void NetworkArc::addToMemoryOfElemSet(const int elemSetId)
{
  if ((elemSetId >= 0) && (elemSetId < _networkPtr->elementaritySetPts().size()))
    {
      _inMemoryOfElemSets.push_back(_networkPtr->elementaritySetPts()[elemSetId]);
      _networkPtr->elementaritySetPts()[elemSetId]->addToMemory(this);
    }
}

void NetworkArc::addBinaryResourceConsumption(int binaryResId, int consumption)
{
  _specResConsumption[binaryResId] = consumption;
}

double NetworkArc::cost() const
{
  return _networkPtr->arcCostMap()[_lemonArc];
}

ScalableResource::ScalableResource(NetworkFlow * networkFlowPtr, int id, int upperBoundOnUse) :
  _id(id), _varPtr(nullptr), _vertexConsumptionLBMap(networkFlowPtr->digraph()),
  _vertexConsumptionUBMap(networkFlowPtr->digraph()), _vertexIncompatibleValuesForSelection(networkFlowPtr->digraph()),
  _arcConsumptionMap(networkFlowPtr->digraph()), _arcConsumptionLBMap(networkFlowPtr->digraph()),
  _arcConsumptionUBMap(networkFlowPtr->digraph()), _arcDepResWaitingTimeCoeff(networkFlowPtr->digraph()),
  _arcConsumptionPWLfunctionMap(networkFlowPtr->digraph()),
  _arcConsumptionPWLlowerBoundFunctionMap(networkFlowPtr->digraph()),
  _arcConsumptionPWLupperBoundFunctionMap(networkFlowPtr->digraph()),
  _arcIncompatibleValuesForSelection(networkFlowPtr->digraph()), _step(-1), _mainResource(false), _disposable(true),
  _baseResource(false), _dependentWithBaseResId(-1), _nbValuesForSelection(0), _initAccumResConsFunction(),
  _associatedVarsForSelection()
{
}

NetworkFlow::~NetworkFlow()
{
  for (std::list<ScalableResource *>::iterator it = _sideResourceConstrList.begin();
       it != _sideResourceConstrList.end(); ++it)
    {
      delete (*it);
      *it = NULL;
    }
  _sideResourceConstrList.clear();

  for (std::vector<NetworkSet *>::iterator elemSetPtrIt = _elementaritySetPts.begin();
       elemSetPtrIt != _elementaritySetPts.end(); ++elemSetPtrIt)
    {
      delete *elemSetPtrIt;
      *elemSetPtrIt = NULL;
    }
  _elementaritySetPts.clear();

  for (std::vector<NetworkSet *>::iterator packSetPtrIt = _packingSetPts.begin();
       packSetPtrIt != _packingSetPts.end(); ++packSetPtrIt)
    {
      delete *packSetPtrIt;
      *packSetPtrIt = NULL;
    }
  _packingSetPts.clear();

  for (std::vector<NetworkSet *>::iterator covSetPtrIt = _coveringSetPts.begin();
       covSetPtrIt != _coveringSetPts.end(); ++covSetPtrIt)
    {
      delete *covSetPtrIt;
      *covSetPtrIt = NULL;
    }
  _coveringSetPts.clear();
  
  for (std::vector<BcArcInfo *>::iterator arcInfoIt = _arcIdToArcInfo.begin(); arcInfoIt != _arcIdToArcInfo.end();
       ++arcInfoIt)
    delete *arcInfoIt;

  for (lemon::ListDigraph::ArcIt lemonArc(_digraph); lemonArc != lemon::INVALID; ++lemonArc)
  {
    delete _netArcPtsMap[lemonArc];
  }

  for (lemon::ListDigraph::NodeIt lemonNode(_digraph); lemonNode != lemon::INVALID; ++lemonNode)
  {
    delete _netNodePtsMap[lemonNode];
  }

  _arcIdToArcInfo.clear();
}

void NetworkFlow::generateArcInfo()
{
  struct CompResourcesById
  {
    bool operator()(const ScalableResource * resAPtr, const ScalableResource * resBPtr) const
    {
      return (resAPtr->id() < resBPtr->id());
    }
  };
  
  std::vector<ScalableResource *> sortedByIdResPts;
  for (std::list<ScalableResource *>::iterator resPtrIt = _sideResourceConstrList.begin();
       resPtrIt != _sideResourceConstrList.end(); ++resPtrIt)
    sortedByIdResPts.push_back(*resPtrIt);
  int numResources = sortedByIdResPts.size();
  std::stable_sort(sortedByIdResPts.begin(), sortedByIdResPts.end(), CompResourcesById());

  _arcIdToArcInfo.resize(_digraph.maxArcId() + 1, NULL);
  for (lemon::ListDigraph::ArcIt lemonArc(_digraph); lemonArc != lemon::INVALID; ++lemonArc)
    {
      double * resConsumption = new double[numResources];
      double * resConsumptionUB = new double[numResources];
      for (int resOrd = 0; resOrd < numResources; ++resOrd)
      {
        resConsumption[resOrd] = sortedByIdResPts[resOrd]->arcConsumption(lemonArc);
        resConsumptionUB[resOrd] = sortedByIdResPts[resOrd]->arcConsumptionUB(lemonArc);
      }
      _arcIdToArcInfo[_digraph.id(lemonArc)] = new BcArcInfo(netArcPtr(lemonArc), resConsumption, resConsumptionUB);
    }
}
