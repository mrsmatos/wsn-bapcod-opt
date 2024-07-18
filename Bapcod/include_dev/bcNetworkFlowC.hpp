/**
 *
 * This file bcNetworkFlowC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCNETWORKFLOWC_HPP
#define	BCNETWORKFLOWC_HPP


#include "bcVarConstrC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcBapcodInit.hpp"
#include <lemon/list_graph.h>

using namespace std;

class NetworkFlow;
class NetworkVertex;
class NetworkArc;

class BcNetwork;
class BcVertex;
class BcArc;

class Constraint;
class InstanciatedVar;

template <typename Type , typename CompleteOrder>
class MembershipMap : public std::map < Type *, Double, CompleteOrder >
{
};

class NetworkSet
{
  int _id;

  std::vector<NetworkArc *> _arcMemory;
  std::vector<NetworkVertex *> _vertexMemory;
  std::vector<NetworkSet *> _cutSeparationNeighbourhood;

public:

  NetworkSet(NetworkFlow * networkFlowPtr, int id) :
    _id(id)
  {
  }

  inline const int & id() const
  {
    return _id;
  }

  void addToMemory(NetworkArc * arcPtr)
  {
    _arcMemory.push_back(arcPtr);
  }

  void addToMemory(NetworkVertex * vertexPtr)
  {
    _vertexMemory.push_back(vertexPtr);
  }

  void addToCutSeparationNeighbourhood(NetworkSet * elemSetPtr)
  {
    _cutSeparationNeighbourhood.push_back(elemSetPtr);
  }

  inline const std::vector<NetworkArc *> & arcMemory() const
  {
    return _arcMemory;
  }

  inline const std::vector<NetworkVertex *> & vertexMemory() const
  {
    return _vertexMemory;
  }

  inline const std::vector<NetworkSet *> & cutSeparationNeighbourhood() const
  {
    return _cutSeparationNeighbourhood;
  }
};

class ScalableResource
{
    int _id;

    InstanciatedVar * _varPtr;

    lemon::ListDigraph::NodeMap<double> _vertexConsumptionLBMap;
    lemon::ListDigraph::NodeMap<double> _vertexConsumptionUBMap;
    lemon::ListDigraph::NodeMap<std::vector<int>> _vertexIncompatibleValuesForSelection;
    lemon::ListDigraph::ArcMap<double> _arcConsumptionMap;
    lemon::ListDigraph::ArcMap<double> _arcConsumptionLBMap;
    lemon::ListDigraph::ArcMap<double> _arcConsumptionUBMap;
    lemon::ListDigraph::ArcMap<double> _arcDepResWaitingTimeCoeff;
    lemon::ListDigraph::ArcMap<std::vector<std::pair<double, double>>> _arcConsumptionPWLfunctionMap;
    lemon::ListDigraph::ArcMap<std::vector<std::pair<double, double>>> _arcConsumptionPWLlowerBoundFunctionMap;
    lemon::ListDigraph::ArcMap<std::vector<std::pair<double, double>>> _arcConsumptionPWLupperBoundFunctionMap;
    lemon::ListDigraph::ArcMap<std::vector<int>> _arcIncompatibleValuesForSelection;

    double _step;
    bool _mainResource;
    bool _disposable;
    bool _baseResource;
    int _dependentWithBaseResId; /// >= 0 if this is a dependent resource
    int _nbValuesForSelection; /// if > 0, then this is a selection resource
    std::vector<std::pair<double, double> > _initAccumResConsFunction; /// only for dependent resources
    std::vector<std::pair<int, InstanciatedVar *>> _associatedVarsForSelection; /// only for selection resources

public:

    ScalableResource(NetworkFlow *networkFlowPtr, int id, int upperBoundOnUse = 0);

    virtual ~ScalableResource() = default;

    inline int id() const
    {
      return _id;
    }

    inline double step() const
    {
      return _step;
    }

    void step(double stepValue)
    {
      _step = stepValue;
    }

    void setAsMainResource()
    {
      _mainResource = true;
      _disposable = true;
      _baseResource = false;
      _dependentWithBaseResId = -1;
      _nbValuesForSelection = 0;
    }

    void disposableResource(bool disposable)
    {
      _disposable = disposable;
    }

    void setAsBaseResource()
    {
        _baseResource = true;
        _mainResource = false;
        _dependentWithBaseResId = -1;
        _nbValuesForSelection = 0;
    }

    void setAsDependentResource(int baseResourceId)
    {
        _dependentWithBaseResId = baseResourceId;
        _baseResource = false;
        _mainResource = false;
        _nbValuesForSelection = 0;
    }

    void setAsSelectionResource(int nbValues)
    {
        _nbValuesForSelection = nbValues;
        _baseResource = false;
        _mainResource = false;
        _dependentWithBaseResId = -1;
    }

    void setInitAccumResConsFunction(std::vector<std::pair<double, double> > & function)
    {
        _initAccumResConsFunction = function;
    }

    inline bool main() const
    {
        return _mainResource;
    }

    inline bool disposable() const
    {
        return _disposable;
    }

    inline bool base() const
    {
        return _baseResource;
    }

    inline bool dependent() const
    {
        return _dependentWithBaseResId >= 0;
    }

    inline int baseResourceId() const
    {
        return _dependentWithBaseResId;
    }

    inline bool selection() const
    {
        return _nbValuesForSelection > 0;
    }

    inline int nbValuesForSelection() const
    {
        return _nbValuesForSelection;
    }

    inline const std::vector<std::pair<double, double> > & initAccumResConsFunction() const
    {
        return _initAccumResConsFunction;
    }

    inline void associateVar(InstanciatedVar *varPtr)
    {
        _varPtr = varPtr;
    }

    inline void associateVarForSelectionValue(int value, InstanciatedVar * varPtr)
    {
        _associatedVarsForSelection.emplace_back(value, varPtr);
    }

    InstanciatedVar * associatedVar()
    {
        return _varPtr;
    }

    const std::vector<std::pair<int, InstanciatedVar *>> & associatedVarsForSelection() const
    {
        return _associatedVarsForSelection;
    }

    inline double arcConsumption(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionMap[lemonArc];
    }

    void arcConsumption(const lemon::ListDigraph::Arc & lemonArc, const double & consumption)
    {
        _arcConsumptionMap[lemonArc] = consumption;
    }

    inline double arcConsumptionLB(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionLBMap[lemonArc];
    }

    void arcConsumptionLB(const lemon::ListDigraph::Arc & lemonArc, const double & consumption)
    {
        _arcConsumptionLBMap[lemonArc] = consumption;
    }

    inline double arcConsumptionUB(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionUBMap[lemonArc];
    }

    double arcDepResWaitingTimeCoeff(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcDepResWaitingTimeCoeff[lemonArc];
    }

    const std::vector<std::pair<double, double>> &
        arcConsumptionPWLfunction(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionPWLfunctionMap[lemonArc];
    }

    const std::vector<std::pair<double, double>> &
        arcConsumptionPWLlowerBoundFunction(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionPWLlowerBoundFunctionMap[lemonArc];
    }

    const std::vector<std::pair<double, double>>
        & arcConsumptionPWLupperBoundFunction(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcConsumptionPWLupperBoundFunctionMap[lemonArc];
    }

    const std::vector<int> & incompatibleValuesForSelection(const lemon::ListDigraph::Node & lemonVertex) const
    {
        return _vertexIncompatibleValuesForSelection[lemonVertex];
    }

    const std::vector<int> & incompatibleValuesForSelection(const lemon::ListDigraph::Arc & lemonArc) const
    {
        return _arcIncompatibleValuesForSelection[lemonArc];
    }

    void arcConsumptionUB(const lemon::ListDigraph::Arc & lemonArc, const double & consumption)
    {
        _arcConsumptionUBMap[lemonArc] = consumption;
    }

    void arcDepResWaitingTimeCoeff(const lemon::ListDigraph::Arc & lemonArc, double coeff)
    {
        _arcDepResWaitingTimeCoeff[lemonArc] = coeff;
    }

    void arcConsumptionPWLfunction(const lemon::ListDigraph::Arc & lemonArc,
                                   const std::vector<std::pair<double, double>> & consFunction)
    {
        _arcConsumptionPWLfunctionMap[lemonArc] = consFunction;
    }

    void arcConsumptionPWLlowerBoundFunction(const lemon::ListDigraph::Arc & lemonArc,
                                             const std::vector<std::pair<double, double>> & consFunction)
    {
        _arcConsumptionPWLlowerBoundFunctionMap[lemonArc] = consFunction;
    }

    void arcConsumptionPWLupperBoundFunction(const lemon::ListDigraph::Arc & lemonArc,
                                             const std::vector<std::pair<double, double>> & consFunction)
    {
        _arcConsumptionPWLupperBoundFunctionMap[lemonArc] = consFunction;
    }

    inline double vertexConsumptionLB(const lemon::ListDigraph::Node & lemonVertex) const
    {
        return _vertexConsumptionLBMap[lemonVertex];
    }

    void vertexConsumptionLB(const lemon::ListDigraph::Node & lemonVertex, const double & consumptionLB)
    {
        _vertexConsumptionLBMap[lemonVertex] = consumptionLB;
    }

    inline double vertexConsumptionUB(const lemon::ListDigraph::Node & lemonVertex) const
    {
        return _vertexConsumptionUBMap[lemonVertex];
    }

    void vertexConsumptionUB(const lemon::ListDigraph::Node & lemonVertex, const double & consumptionUB)
    {
        _vertexConsumptionUBMap[lemonVertex] = consumptionUB;
    }

    inline void addIncompatibleValue(const lemon::ListDigraph::Node & lemonVertex, int value)
    {
        _vertexIncompatibleValuesForSelection[lemonVertex].push_back(value);
    }

    inline void addIncompatibleValue(const lemon::ListDigraph::Arc & lemonArc, int value)
    {
        _arcIncompatibleValuesForSelection[lemonArc].push_back(value);
    }

    virtual bool operator<(const ScalableResource & that)
    {
        return _id < that._id;
    }
};

class NetworkVertex
{
  friend class BcVertex;
  friend class BcArc;

private:
  NetworkFlow * _networkPtr;
  const lemon::ListDigraph & _digraph;
  lemon::ListDigraph::Node _lemonVertex;
  std::vector<NetworkSet *> _elementaritySetPts;
  std::vector<NetworkSet *> _enumElementaritySetPts;
  std::vector<NetworkSet *> _packingSetPts;
  std::vector<NetworkSet *> _coveringSetPts;
  std::vector<NetworkSet *> _inMemoryOfElemSets;
  std::map<int, std::pair<int,int>> _specResConsumptionBounds;
  std::string _name;
  bool _isFake; /// sometimes we need just to increase current vertexId without really creating the vertex

public:

  NetworkVertex(NetworkFlow * networkPtr,
                const lemon::ListDigraph & digraph,
                const lemon::ListDigraph::Node & lemonVertex, bool isFake) :
    _networkPtr(networkPtr), _digraph(digraph), _lemonVertex(lemonVertex), _elementaritySetPts(),
    _enumElementaritySetPts(), _packingSetPts(), _coveringSetPts(), _inMemoryOfElemSets(), _specResConsumptionBounds(),
    _isFake(isFake)
  {
  }

  virtual ~NetworkVertex() = default;
  
  inline const lemon::ListDigraph::Node & lemonVertex() const
  {
    return _lemonVertex;
  }
  
  inline void lemonVertex(lemon::ListDigraph::Node lemonVertex)
  {
    _lemonVertex = lemonVertex;
  }

  inline NetworkFlow * networkPtr()
  {
    return _networkPtr;
  }
  
  inline int id() const
  {
    return _digraph.id(_lemonVertex);
  }

  inline const lemon::ListDigraph & digraph()
  {
    return _digraph;
  }
  
  void setElementaritySet(int elemSetId);

  void addToElementaritySet(int elemSetId);

  void addToEnumerationOnlyElementaritySet(int elemSetId);

  inline const std::vector<NetworkSet *> & elementaritySetPts() const
  {
     return _elementaritySetPts;
  }

  inline const std::vector<NetworkSet *> & enumerationOnlyElementaritySetPts() const
  {
      return _enumElementaritySetPts;
  }

  void setPackingSet(int packSetId);

  void addToPackingSet(int packSetId);

  inline const std::vector<NetworkSet *> & packingSetPts() const
  {
    return _packingSetPts;
  }

  void setCoveringSet(int covSetId);

  void addToCoveringSet(int covSetId);

  inline const std::vector<NetworkSet *> & coveringSetPts() const
  {
    return _coveringSetPts;
  }
  
  void addToMemoryOfElemSet(int elemSetId);
  
  void setSpecialResourceConsumptionLB(int binaryResId, int lowerBound);
  void setSpecialResourceConsumptionUB(int binaryResId, int upperBound);

  void name(std::string name_)
  {
      _name = std::move(name_);
  }

  inline std::string & name()
  {
      return _name;
  }

  inline const std::vector<NetworkSet *> & inMemoryOfElemSets() const
  {
    return _inMemoryOfElemSets;
  }

  inline const std::map<int, std::pair<int, int> > & specResConsumptionBounds() const
  {
    return _specResConsumptionBounds;
  };

  inline bool isFake() const
  {
      return _isFake;
  }
};

class NetworkArc
{
  friend class BcVertex;
  friend class BcArc;

private:
  NetworkFlow * _networkPtr;
  const lemon::ListDigraph & _digraph;
  lemon::ListDigraph::Arc _lemonArc;
  std::vector<std::map<InstanciatedVar *, double>> _varToCoeffMaps;
  std::vector<NetworkSet *> _elementaritySetPts;
  std::vector<NetworkSet *> _enumElementaritySetPts;
  std::vector<NetworkSet *> _packingSetPts;
  std::vector<NetworkSet *> _coveringSetPts;
  std::vector<NetworkSet *> _inMemoryOfElemSets;
  std::map<int, int> _specResConsumption;
  std::string _name;
  bool _isFake; /// sometimes we need to add fake arc just to increase current arc id
  int _currentMapping;

public:

  NetworkArc(NetworkFlow * networkPtr,
             const lemon::ListDigraph & digraph,
             const lemon::ListDigraph::Arc & lemonArc, bool isFake) :
          _networkPtr(networkPtr), _digraph(digraph), _lemonArc(lemonArc), _varToCoeffMaps(1), _currentMapping(0),
          _isFake(isFake)
  {
  }

  virtual ~NetworkArc() = default;
  
  inline const lemon::ListDigraph::Arc & lemonArc() const
  {
    return _lemonArc;
  }
  
  inline void lemonArc(lemon::ListDigraph::Arc lemonArc)
  {
    _lemonArc = lemonArc;
  }

  inline NetworkFlow * networkPtr() const
  {
    return _networkPtr;
  }

  NetworkVertex * tailVertexPtr() const;
  NetworkVertex * headVertexPtr() const;

  inline const std::vector<std::map<InstanciatedVar *, double>> & varToCoeffMaps() const
  {
    return _varToCoeffMaps;
  }
    
  double cost() const;
  
  inline const int id() const
  {
    return _digraph.id(_lemonArc);
  }

  void addAlternativeMapping()
  {
      _currentMapping += 1;
      _varToCoeffMaps.emplace_back();
  }

    /// attention: not cumulative!
  void addVarAssociation(InstanciatedVar * instVarPtr, const double coeff = 1)
  {
    _varToCoeffMaps[_currentMapping][instVarPtr] = coeff;
    if (_currentMapping == 0)
        /// TO DO : we need to get rid of InstanciatedVar::_arcIdToCoeff,
        /// as it is used only in branchings which are not used anymore
        /// (pack.set. based assignment branching and paths-per-network branching)
        /// so that this call is not needed anymore
        instVarPtr->setArcMembership(id(), coeff);
  }

  /// backward compatibility
  void var(InstanciatedVar * instVarPtr)
  {
    addVarAssociation(instVarPtr, 1);
  }

    void setElementaritySet(int elemSetId);

    void addToElementaritySet(int elemSetId);

    void addToEnumerationOnlyElementaritySet(int elemSetId);

    inline const std::vector<NetworkSet *> & elementaritySetPts() const
    {
        return _elementaritySetPts;
    }

    inline const std::vector<NetworkSet *> & enumerationOnlyElementaritySetPts() const
    {
        return _enumElementaritySetPts;
    }

    void setPackingSet(int packSetId);

    void addToPackingSet(int packSetId);

    inline const std::vector<NetworkSet *> & packingSetPts() const
    {
        return _packingSetPts;
    }

    void setCoveringSet(int covSetId);

    void addToCoveringSet(int covSetId);

    inline const std::vector<NetworkSet *> & coveringSetPts() const
    {
        return _coveringSetPts;
    }

    void name(std::string name_);

  inline const std::string & name() const
  {
      return _name;
  }

  inline const lemon::ListDigraph & digraph()
  {
    return _digraph;
  }
  
  void addToMemoryOfElemSet(int elemSetId);
    
  inline const std::vector<NetworkSet *> & inMemoryOfElemSets() const
  {
    return _inMemoryOfElemSets;
  }

  void addBinaryResourceConsumption(int binaryResId, int consumption);

  inline const std::map<int, int> & specResConsumption() const
  {
    return _specResConsumption;
  };

  inline bool isFake() const
  {
    return _isFake;
  };
};

struct BcArcInfo;

#ifdef BCP_RCSP_IS_FOUND
namespace bcp_rcsp {
    struct GraphData;
}
#endif /* BCP_RCSP_IS_FOUND */

class NetworkFlow
{
  friend class BcArc;
  friend class BcVertex;
  friend class BcNetwork;
  friend class NetworkArc;
  friend class NetworkVertex;  
  
private:
  std::vector<BcArcInfo *> _arcIdToArcInfo;

  std::list<lemon::ListDigraph::Node> _sourceList;
  std::list<lemon::ListDigraph::Node> _sinkList;

  lemon::ListDigraph _digraph;

  //THE FOLLOWING MAPS HAVE O(1) ACCESS
  lemon::ListDigraph::NodeMap<NetworkVertex *> _netNodePtsMap;
  lemon::ListDigraph::ArcMap<NetworkArc *> _netArcPtsMap;
  lemon::ListDigraph::ArcMap<double> _arcCostMap;

  /// if one wants to record side constraint    
  std::list<ScalableResource *> _sideResourceConstrList;
  std::vector<NetworkSet *> _elementaritySetPts;
  std::vector<NetworkSet *> _packingSetPts;
  std::vector<NetworkSet *> _coveringSetPts;

  std::vector<std::vector<double> > _elemSetsDistanceMatrix;

  std::vector<std::tuple<int, int, bool> > _permRyanFosterConstrs;

  std::set<int> _nonDisposableSpecialResources;

public:

  NetworkFlow(const int numElemSets = 0, const int numPackSets = 0, const int numCovSets = 0) :
    _sourceList(), _sinkList(), _digraph(), _netNodePtsMap(_digraph), _netArcPtsMap(_digraph),
    _arcCostMap(_digraph), _sideResourceConstrList(), _elementaritySetPts(), _packingSetPts(), _coveringSetPts(),
    _elemSetsDistanceMatrix(), _permRyanFosterConstrs(), _nonDisposableSpecialResources()
  {
    for (int elemSetId = 0; elemSetId < numElemSets; ++elemSetId)
      _elementaritySetPts.push_back(new NetworkSet(this, elemSetId));
    for (int packSetId = 0; packSetId < numPackSets; ++packSetId)
      _packingSetPts.push_back(new NetworkSet(this, packSetId));
    for (int covSetId = 0; covSetId < numCovSets; ++covSetId)
      _coveringSetPts.push_back(new NetworkSet(this, covSetId));
  }
  
  virtual ~NetworkFlow();
  
  void generateArcInfo();

  void setSpecialResourceNonDisposable(const int binaryResId)
  {
    _nonDisposableSpecialResources.insert(binaryResId);
  }
  
  inline const BcArcInfo * getArcInfoPtr(const int arcId) const
  {
    return ((arcId >= 0) && (arcId < _arcIdToArcInfo.size()) ? _arcIdToArcInfo[arcId] : NULL);
  }

  inline const lemon::ListDigraph & digraph() const
  {
    return _digraph;
  }

  NetworkVertex * createVertex(bool isFake = false)
  {
    lemon::ListDigraph::Node lemonVertex = _digraph.addNode();
    NetworkVertex* netVertexPtr = new NetworkVertex(this, _digraph, lemonVertex, isFake);
    _netNodePtsMap[lemonVertex] = netVertexPtr;

    for (std::list<ScalableResource *>::iterator resPtrIt = _sideResourceConstrList.begin();
         resPtrIt != _sideResourceConstrList.end(); ++resPtrIt)
    {
      (*resPtrIt)->vertexConsumptionLB(lemonVertex, -BapcodInfinity);
      (*resPtrIt)->vertexConsumptionUB(lemonVertex, BapcodInfinity);
    }

    return netVertexPtr;
  }

    NetworkArc * createAndInsertArc(lemon::ListDigraph::Node tail, lemon::ListDigraph::Node head,
                                    double originalCost, bool isFake)
    {
        lemon::ListDigraph::Arc lemonArc = _digraph.addArc(tail, head);
        NetworkArc * netArcPtr = new NetworkArc(this, _digraph, lemonArc, isFake);

        _netArcPtsMap[lemonArc] = netArcPtr;

        for (std::list<ScalableResource *>::iterator resPtrIt = _sideResourceConstrList.begin();
             resPtrIt != _sideResourceConstrList.end(); ++resPtrIt)
        {
            (*resPtrIt)->arcConsumptionLB(lemonArc, (*resPtrIt)->vertexConsumptionLB(head));
            (*resPtrIt)->arcConsumptionUB(lemonArc, (*resPtrIt)->vertexConsumptionUB(head));
            (*resPtrIt)->arcConsumption(lemonArc, 0.0);
        }

        _arcCostMap[lemonArc] = originalCost;

        return netArcPtr;
    }

    inline NetworkArc * netArcPtr(const int & id) const
  {
    const lemon::ListDigraph::Arc & lemonArc = _digraph.arcFromId(id);
    if (lemonArc == lemon::INVALID)
      return NULL;

    return _netArcPtsMap[lemonArc];
  }
  
  inline NetworkVertex * netVertexPtr(const int & id) const
  {
    const lemon::ListDigraph::Node & lemonNode = _digraph.nodeFromId(id);
    if (lemonNode == lemon::INVALID)
      return NULL;

    return _netNodePtsMap[lemonNode];
  }
  
  inline NetworkArc * netArcPtr(const lemon::ListDigraph::Arc & lemonArc) const
  {
    return _netArcPtsMap[lemonArc];
  }
  
  inline NetworkVertex * netVertexPtr(const lemon::ListDigraph::Node & lemonVertex) const
  {
    return _netNodePtsMap[lemonVertex];
  }

  inline double arcCost(const lemon::ListDigraph::Arc & lemonArc) const
  {
    return arcCostMap()[lemonArc];
  }

  inline void arcCost(const lemon::ListDigraph::Arc & lemonArc, double cost)
  {
      _arcCostMap[lemonArc] = cost;
  }

  inline const lemon::ListDigraph & digraph()
  {
    return _digraph;
  }

  const lemon::ListDigraph::ArcMap<double> & arcCostMap() const
  {
      return _arcCostMap;
  }

  inline const std::list<lemon::ListDigraph::Node> & sourceList() const
  {
    return _sourceList;           
  }
  
  inline const std::list<lemon::ListDigraph::Node> & sinkList() const
  {
    return _sinkList;           
  }
  
  void addSideResource(ScalableResource * resourcePtr)
  {
    _sideResourceConstrList.push_back(resourcePtr);
  }
  
  inline const std::list<ScalableResource *> & sideResourceConstrList() const
  {
    return _sideResourceConstrList;
  }
  
  inline const std::vector<NetworkSet *> & elementaritySetPts() const
  {
    return _elementaritySetPts;
  }

  inline std::vector<NetworkSet *> & elementaritySetPts()
  {
    return _elementaritySetPts;
  }

  inline const std::vector<NetworkSet *> & packingSetPts() const
  {
    return _packingSetPts;
  }

  inline std::vector<NetworkSet *> & packingSetPts()
  {
    return _packingSetPts;
  }

  inline const std::vector<NetworkSet *> & coveringSetPts() const
  {
    return _coveringSetPts;
  }

  inline std::vector<NetworkSet *> & coveringSetPts()
  {
    return _coveringSetPts;
  }

  void setPathSource(const lemon::ListDigraph::Node & lemonVertex)
  {
      _sourceList.clear();
      _sourceList.push_back(lemonVertex);
  }

  void setPathSink(const lemon::ListDigraph::Node & lemonVertex)
  {
      _sinkList.clear();
      _sinkList.push_back(lemonVertex);
  }

  void setElemSetsDistanceMatrix(const std::vector<std::vector<double> > & matrix)
  {
    _elemSetsDistanceMatrix = matrix;
  }

  const std::set<int> & nonDisposableSpecialResourceIds() const
  {
    return _nonDisposableSpecialResources;
  }

  const std::vector<std::vector<double> > & elemSetsDistanceMatrix() const
  {
    return _elemSetsDistanceMatrix;
  }

  void addPermanentRyanAndFosterConstraint(int firstPackSetId, int secondPackSetId, bool together)
  {
      _permRyanFosterConstrs.push_back(std::make_tuple(firstPackSetId, secondPackSetId, together));
  }

  const std::vector<std::tuple<int, int, bool> > & permanentRyanFosterConstraints() const
  {
      return _permRyanFosterConstrs;
  }
};





#endif	/* BCNETWORKFLOWC_HPP */

