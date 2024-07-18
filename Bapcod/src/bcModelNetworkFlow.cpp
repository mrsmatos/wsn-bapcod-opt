/**
 *
 * This file bcModelNetworkFlow.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include <utility>

#include "bcModelNetworkFlow.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcProbConfigC.hpp"
#include "bcModelFormulationC.hpp"

using namespace std;


int BcNetworkResource::id() const
{
  return _resourcePtr->id();
}

void BcNetworkResource::setAsMainResource()
{
    _resourcePtr->setAsMainResource();
}

void BcNetworkResource::setAsNonDisposableResource()
{
  _resourcePtr->disposableResource(false);
}

void BcNetworkResource::setAsMainResourceWithStep(const double & stepValue)
{
  _resourcePtr->setAsMainResource();
  _resourcePtr->step(stepValue);
}

void BcNetworkResource::setAsSelectionResourceWithNbValues(int nbValues)
{
    _resourcePtr->setAsSelectionResource(nbValues);
}

BcNetworkResource::BcNetworkResource(const BcNetwork & bcNetwork, int id, int upperBoundOnUse) :
    _resourcePtr(new ScalableResource((NetworkFlow *) bcNetwork, id, upperBoundOnUse))
{
  ((NetworkFlow *) bcNetwork)->addSideResource(_resourcePtr);
}

void BcNetworkResource::setArcConsumption(const BcArc & bcArc, const double & consumption)
{
  _resourcePtr->arcConsumption(((NetworkArc *) bcArc)->lemonArc(), consumption);
}

void BcNetworkResource::setArcConsumptionLB(const BcArc & bcArc, const double & consumptionLB)
{
  _resourcePtr->arcConsumptionLB(((NetworkArc *) bcArc)->lemonArc(), consumptionLB);
}

void BcNetworkResource::setArcConsumptionUB(const BcArc & bcArc, const double & consumptionUB)
{
  _resourcePtr->arcConsumptionUB(((NetworkArc *) bcArc)->lemonArc(), consumptionUB);
}


void BcNetworkResource::setVertexConsumptionLB(const BcVertex & bcVertex, const double & consumptionLB)
{
  _resourcePtr->vertexConsumptionLB(((NetworkVertex *) bcVertex)->lemonVertex(), consumptionLB);
}

void BcNetworkResource::setVertexConsumptionUB(const BcVertex & bcVertex, const double & consumptionUB)
{
  _resourcePtr->vertexConsumptionUB(((NetworkVertex *) bcVertex)->lemonVertex(), consumptionUB);
}

void BcNetworkResource::setVertexIncompatibleValueForSelection(const BcVertex & bcVertex, int value)
{
    _resourcePtr->addIncompatibleValue(((NetworkVertex *) bcVertex)->lemonVertex(), value);
}

void BcNetworkResource::setAssociatedVar(const BcVar & bcVar)
{
  _resourcePtr->associateVar((InstanciatedVar *) bcVar);
}

void BcNetworkResource::setAssociatedVarForSelectionValue(int value, const BcVar & bcVar)
{
    _resourcePtr->associateVarForSelectionValue(value, (InstanciatedVar *) bcVar);
}

void BcNetworkResource::setAsBaseResource()
{
    _resourcePtr->setAsBaseResource();
}

void BcNetworkResource::setAsBaseResourceWithStep(double stepValue)
{
    _resourcePtr->setAsBaseResource();
    _resourcePtr->step(stepValue);
}

void BcNetworkResource::setAsDependentResource(const BcNetworkResource & baseResource)
{
    _resourcePtr->setAsDependentResource(baseResource.id());
}

void BcNetworkResource::setInitResConsPWLfunction(std::vector<std::pair<double, double> > & function)
{
    _resourcePtr->setInitAccumResConsFunction(function);
}

void BcNetworkResource::setInitResConsLinearFunction(double initValue, double slope)
{
    std::vector<std::pair<double, double> > function(1, std::make_pair(initValue, slope));
    _resourcePtr->setInitAccumResConsFunction(function);
}

void BcNetworkResource::setWaitingTimeCoeff(const BcArc & bcArc, double coeff)
{
    _resourcePtr->arcDepResWaitingTimeCoeff(((NetworkArc *) bcArc)->lemonArc(), coeff);
}

void BcNetworkResource::setArcConsumptionPWLfunction(const BcArc & bcArc,
                                                  const std::vector<std::pair<double, double> > & function)
{
    _resourcePtr->arcConsumptionPWLfunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcConsumptionLinearFunction(const BcArc & bcArc, double initValue, double slope)
{
    std::vector<std::pair<double, double> > function(1, std::make_pair(initValue, slope));
    _resourcePtr->arcConsumptionPWLfunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcConsumptionLowerBoundPWLfunction(const BcArc & bcArc,
                                                               std::vector<std::pair<double, double> > & function)
{
    _resourcePtr->arcConsumptionPWLlowerBoundFunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcConsumptionLowerBoundLinearFunction(const BcArc & bcArc, double initValue, double slope)
{
    std::vector<std::pair<double, double> > function(1, std::make_pair(initValue, slope));
    _resourcePtr->arcConsumptionPWLlowerBoundFunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcConsumptionUpperBoundPWLfunction(const BcArc & bcArc,
                                                               std::vector<std::pair<double, double> > & function)
{
    _resourcePtr->arcConsumptionPWLupperBoundFunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcConsumptionUpperBoundLinearFunction(const BcArc & bcArc, double initValue, double slope)
{
    std::vector<std::pair<double, double> > function(1, std::make_pair(initValue, slope));
    _resourcePtr->arcConsumptionPWLupperBoundFunction(((NetworkArc *) bcArc)->lemonArc(), function);
}

void BcNetworkResource::setArcIncompatibleValueForSelection(const BcArc & bcArc, int value)
{
    _resourcePtr->addIncompatibleValue(((NetworkArc *) bcArc)->lemonArc(), value);
}


///BcVertex methods

lemon::ListDigraph::Node BcVertex::getLemonNode()
{
  return _netVertex->_lemonVertex;
}

BcVertex::BcVertex(NetworkFlow * networkPtr, bool isFake)
{
  _netVertex = networkPtr->createVertex(isFake);
}

int BcVertex::ref() const
{
  return lemon::ListDigraph::id(_netVertex->lemonVertex());
}

void BcVertex::setElementaritySet(int elemSetId)
{
  _netVertex->setElementaritySet(elemSetId);
}

void BcVertex::setPackingSet(int packSetId)
{
  _netVertex->setPackingSet(packSetId);
}

void BcVertex::setCoveringSet(int covSetId)
{
  _netVertex->setCoveringSet(covSetId);
}

void BcVertex::addToElementaritySet(int elemSetId)
{
    _netVertex->addToElementaritySet(elemSetId);
}

void BcVertex::addToEnumerationOnlyElementaritySet(int elemSetId)
{
    _netVertex->addToEnumerationOnlyElementaritySet(elemSetId);
}

void BcVertex::addToPackingSet(int packSetId)
{
    _netVertex->addToPackingSet(packSetId);
}

void BcVertex::addToCoveringSet(int covSetId)
{
    _netVertex->addToCoveringSet(covSetId);
}

void BcVertex::addToMemoryOfElemSet(int elemSetId)
{
  _netVertex->addToMemoryOfElemSet(elemSetId);
}

void BcVertex::setSpecialResourceConsumptionLB(int binaryResId, int lowerBound)
{
  _netVertex->setSpecialResourceConsumptionLB(binaryResId, lowerBound);
}

void BcVertex::setSpecialResourceConsumptionUB(int binaryResId, int upperBound)
{
  _netVertex->setSpecialResourceConsumptionUB(binaryResId, upperBound);
}

void BcVertex::setBinaryResourceConsumptionLB(int binaryResId, int lowerBound)
{
    _netVertex->setSpecialResourceConsumptionLB(binaryResId, lowerBound);
}

void BcVertex::setBinaryResourceConsumptionUB(int binaryResId, int upperBound)
{
    _netVertex->setSpecialResourceConsumptionUB(binaryResId, upperBound);
}

void BcVertex::setName(std::string name_)
{
    _netVertex->name(std::move(name_));
}

bool BcVertex::defined() 
{
  return _netVertex != nullptr;
}

/// BcArcInfo methods

BcArcInfo::BcArcInfo(const NetworkArc * const arcPtr, double * resConsumption_, double * resConsumptionUB_) :
    arcId(arcPtr->id()), arcElemSetId(-1), arcPackSetId(-1), arcCovSetId(-1), tailVertId(arcPtr->tailVertexPtr()->id()),
    tailElemSetId(-1), tailPackSetId(-1), tailCovSetId(-1), headVertId(arcPtr->headVertexPtr()->id()),
    headElemSetId(-1), headPackSetId(-1), headCovSetId(-1), resConsumption(resConsumption_),
    resConsumptionUB(resConsumptionUB_), name(arcPtr->name())
{
  if (!arcPtr->elementaritySetPts().empty())
    arcElemSetId = arcPtr->elementaritySetPts().front()->id();
  if (!arcPtr->packingSetPts().empty())
    arcPackSetId = arcPtr->packingSetPts().front()->id();
  if (!arcPtr->coveringSetPts().empty())
    arcCovSetId = arcPtr->coveringSetPts().front()->id();
  if (!arcPtr->tailVertexPtr()->elementaritySetPts().empty())
    tailElemSetId = arcPtr->tailVertexPtr()->elementaritySetPts().front()->id();
  if (!arcPtr->tailVertexPtr()->packingSetPts().empty())
    tailPackSetId = arcPtr->tailVertexPtr()->packingSetPts().front()->id();
  if (!arcPtr->tailVertexPtr()->coveringSetPts().empty())
    tailCovSetId = arcPtr->tailVertexPtr()->coveringSetPts().front()->id();
  if (!arcPtr->headVertexPtr()->elementaritySetPts().empty())
    headElemSetId = arcPtr->headVertexPtr()->elementaritySetPts().front()->id();
  if (!arcPtr->headVertexPtr()->packingSetPts().empty())
    headPackSetId = arcPtr->headVertexPtr()->packingSetPts().front()->id();
  if (!arcPtr->headVertexPtr()->coveringSetPts().empty())
    headCovSetId = arcPtr->headVertexPtr()->coveringSetPts().front()->id();
}

BcArcInfo::~BcArcInfo()
{
  delete [] resConsumption;
  delete [] resConsumptionUB;
}

/// BcArc methods

lemon::ListDigraph::Arc BcArc::getLemonArc()
{
  return _netArc->_lemonArc;
}

BcArc::BcArc(NetworkFlow * networkPtr, const BcVertex & tail, const BcVertex & head, double originalCost, bool isFake) :
    _netArc(nullptr)
{
  _netArc = networkPtr->createAndInsertArc(((NetworkVertex*) tail)->lemonVertex(),
                                           ((NetworkVertex*) head)->lemonVertex(), originalCost, isFake);
}

int BcArc::ref() const
{
    return lemon::ListDigraph::id(_netArc->lemonArc());
}

BcVertex BcArc::tail() const
{
  lemon::ListDigraph::Node lemonVertex = _netArc->digraph().source(_netArc->lemonArc());
  return BcVertex(_netArc->networkPtr()->netVertexPtr(lemonVertex));
}

BcVertex BcArc::head() const
{
  lemon::ListDigraph::Node lemonVertex = _netArc->digraph().target(_netArc->lemonArc());
  return BcVertex(_netArc->networkPtr()->netVertexPtr(lemonVertex));
}

void BcArc::arcVar(const BcVar & newvar)
{
  _netArc->var((InstanciatedVar *)newvar);
}

void BcArc::addVarAssociation(const BcVar & newvar, double coeff)
{
  _netArc->addVarAssociation((InstanciatedVar *)newvar, coeff);
}

void BcArc::addAlternativeMapping()
{
    _netArc->addAlternativeMapping();
}

bool BcArc::defined() 
{
  return _netArc != nullptr;
}

void BcArc::cost(double cost)
{
  _netArc->networkPtr()->arcCost(_netArc->lemonArc(), cost);
}

double BcArc::cost()
{
  return _netArc->networkPtr()->arcCost(_netArc->lemonArc());
}

void BcArc::setElementaritySet(int elemSetId)
{
    _netArc->setElementaritySet(elemSetId);
}

void BcArc::setPackingSet(int packSetId)
{
    _netArc->setPackingSet(packSetId);
}

void BcArc::setCoveringSet(int covSetId)
{
    _netArc->setCoveringSet(covSetId);
}

void BcArc::addToElementaritySet(int elemSetId)
{
    _netArc->addToElementaritySet(elemSetId);
}

void BcArc::addToEnumerationOnlyElementaritySet(int elemSetId)
{
    _netArc->addToEnumerationOnlyElementaritySet(elemSetId);
}

void BcArc::addToPackingSet(int packSetId)
{
    _netArc->addToPackingSet(packSetId);
}

void BcArc::addToCoveringSet(int covSetId)
{
    _netArc->addToCoveringSet(covSetId);
}

void BcArc::setName(std::string name_)
{
    _netArc->name(std::move(name_));
}

void BcArc::addToMemoryOfElemSet(int elemSetId)
{
  _netArc->addToMemoryOfElemSet(elemSetId);
}

void BcArc::addBinaryResourceConsumption(int binaryResId, int consumption)
{
  _netArc->addBinaryResourceConsumption(binaryResId, consumption);
}

///BcNework methods

lemon::ListDigraph::Digraph & BcNetwork::getLemonDigraph()
{
  return _networkPtr->_digraph;
}

const lemon::ListDigraph::Digraph & BcNetwork::getLemonDigraph() const
{
  return _networkPtr->_digraph;
}

BcArc BcNetwork::getBcArc(lemon::ListDigraph::Arc arc)
{
  return BcArc(_networkPtr->_netArcPtsMap[arc]);
}

BcVertex BcNetwork::getBcVertex( lemon::ListDigraph::Node node)
{
  return BcVertex(_networkPtr->_netNodePtsMap[node]);
}

BcNetwork::BcNetwork(NetworkFlow * networkFlowPtr) :
    _networkPtr(networkFlowPtr)
{
}

BcNetwork::BcNetwork(BcFormulation & bcForm, int numElemSets, int numPackSets, int numCovSets) :
    _networkPtr(new NetworkFlow(numElemSets, numPackSets, numCovSets))
{
    bcForm.probConfPtr()->setNetworkFlowPtr(_networkPtr);
}

/// for backward compatibility
BcNetwork::BcNetwork(BcFormulation & bcForm, NetAttrMask optionalAttrMask, int numElemSets, int numPackSets,
                     int numCovSets) :
  _networkPtr(new NetworkFlow(numElemSets, numPackSets, numCovSets))
{
  bcForm.probConfPtr()->setNetworkFlowPtr(_networkPtr);
}

BcVertex BcNetwork::createVertex(bool isFake)
{
  BcVertex vertex(_networkPtr, isFake);
  return vertex;
}

BcArc BcNetwork::createArc(int tail, int head, double originalCost, bool isFake)
{
  BcArc arc(_networkPtr, getVertex(tail), getVertex(head), originalCost, isFake);
  return arc;
}

BcVertex BcNetwork::getVertex(int id) const
{
  lemon::ListDigraph::Node lemonVertex = lemon::ListDigraph::nodeFromId(id);
  return BcVertex(_networkPtr->netVertexPtr(lemonVertex));
}

BcArc BcNetwork::getArc(int id) const
{
  lemon::ListDigraph::Arc lemonArc = lemon::ListDigraph::arcFromId(id);
  return BcArc(_networkPtr->netArcPtr(lemonArc));
}

BcArc BcNetwork::getFirstArcBetweenVertices(int tailId, int headId) const
{
  lemon::ListDigraph::Node tailVertex = lemon::ListDigraph::nodeFromId(tailId);
  lemon::ListDigraph::Node headVertex = lemon::ListDigraph::nodeFromId(headId);
  lemon::ListDigraph::OutArcIt outArcIt(_networkPtr->digraph(), tailVertex);
  while ( (outArcIt != lemon::INVALID) && (_networkPtr->digraph().target(outArcIt) != headVertex) )
    {
      ++outArcIt;
    }
  return BcArc(_networkPtr->netArcPtr(outArcIt));
}

int BcNetwork::nVertices() const
{
  return lemon::countNodes(_networkPtr->digraph());
} //linear time complexity

int BcNetwork::nArcs() const
{
  return lemon::countArcs(_networkPtr->digraph());
} //linear time complexity

void BcNetwork::setPathSource(const BcVertex & vertex)
{
  _networkPtr->setPathSource(((NetworkVertex*) vertex)->lemonVertex());
}

void BcNetwork::setPathSink(const BcVertex & vertex)
{
  _networkPtr->setPathSink(((NetworkVertex*) vertex)->lemonVertex());
}

void BcNetwork::addToPackingSetCutNeighbourhood(int packSetId, int packSetIdToAdd)
{
  int numPackSets = (int)_networkPtr->packingSetPts().size();
  if ((packSetId >= numPackSets) || (packSetIdToAdd > numPackSets))
    {
      std::cerr << "BaPCod error : packing set id is too high in addToPackingSetCutNeighbourhood" << std::endl;
      exit(1);
    }
  NetworkSet * packSetPtrToAdd = _networkPtr->packingSetPts()[packSetIdToAdd];
  _networkPtr->packingSetPts()[packSetId]->addToCutSeparationNeighbourhood(packSetPtrToAdd);
}

void BcNetwork::setElemSetsDistanceMatrix(const std::vector<std::vector<double> > & matrix)
{
   int numSets = (int)_networkPtr->elementaritySetPts().size();
   if (numSets == 0)
   {
       std::cerr << "BaPCod error : distance matrix cannot be set if there is no elementary sets" << std::endl;
       exit(1);
   }
   if ((int)matrix.size() != numSets)
   {
       std::cerr << "BaPCod error : distance matrix size does not equal to the number of elementary sets" << std::endl;
       exit(1);
   }
   for (auto & row : matrix)
       if ((int) row.size() != numSets)
       {
           std::cerr << "BaPCod error : distance matrix size does not equal to the number of elementary sets"
                     << std::endl;
           exit(1);
       }
   _networkPtr->setElemSetsDistanceMatrix(matrix);
}

void BcNetwork::setBinaryResourceNonDisposable(const int binaryResId)
{
    _networkPtr->setSpecialResourceNonDisposable(binaryResId);
}

void BcNetwork::setSpecialResourceNonDisposable(const int binaryResId)
{
    _networkPtr->setSpecialResourceNonDisposable(binaryResId);
}

void BcNetwork::addPermanentRyanAndFosterConstraint(int firstPackSetId, int secondPackSetId, bool together)
{
    _networkPtr->addPermanentRyanAndFosterConstraint(firstPackSetId, secondPackSetId, together);
}

BcNetwork::~BcNetwork() = default;
