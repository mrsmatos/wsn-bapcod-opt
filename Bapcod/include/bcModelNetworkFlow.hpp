/**
 *
 * This file bcModelNetworkFlow.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELNETWORKFLOW_HPP
#define	BCMODELNETWORKFLOW_HPP

#include "bcModelVarC.hpp"
#include "bcModelConstrC.hpp"
#include <lemon/list_graph.h>

class BcNetwork;
class BcVertex;
class BcArc;

class NetworkSet;
class ScalableResource;
class NetworkArc;

class NetworkFlow;
class NetworkVertex;

class BcNetworkResource
{
  ScalableResource * _resourcePtr;
 
public:
  
  BcNetworkResource() : _resourcePtr(nullptr) {}
    
  BcNetworkResource(const BcNetwork & bcNetwork, int id, int upperBoundOnUse = 0);
  
  BcNetworkResource( const BcNetworkResource & that ) = default;
        
  int id() const;

  void setAsMainResource();

  void setAsNonDisposableResource(); /// by default disposable

  void setAsMainResourceWithStep(const double & stepValue);

  void setAsSelectionResourceWithNbValues(int nbValues);

  void setArcConsumption(const BcArc & bcArc, const double & consumption);

  void setArcConsumptionLB(const BcArc & bcArc, const double & consumptionLB);

  void setArcConsumptionUB(const BcArc & bcArc, const double & consumptionUB);

  void setVertexConsumptionLB(const BcVertex & bcVertex, const double & consumptionLB);

  void setVertexConsumptionUB(const BcVertex & bcVertex, const double & consumptionUB);

  void setVertexIncompatibleValueForSelection(const BcVertex & bcVertex, int value);

  void setAssociatedVar(const BcVar & bcVar);

  void setAssociatedVarForSelectionValue(int value, const BcVar & bcVar);

  void setAsBaseResource();

  void setAsBaseResourceWithStep(double stepValue);

  void setAsDependentResource(const BcNetworkResource & baseResource);

  void setInitResConsPWLfunction(std::vector<std::pair<double, double> > & function);

  void setInitResConsLinearFunction(double initValue, double slope = 0);

  void setWaitingTimeCoeff(const BcArc & bcArc, double coeff);

  void setArcConsumptionPWLfunction(const BcArc & bcArc, const std::vector<std::pair<double, double> > & function);

  void setArcConsumptionLinearFunction(const BcArc & bcArc, double initValue, double slope = 0);

  void setArcConsumptionLowerBoundPWLfunction(const BcArc & bcArc, std::vector<std::pair<double, double> > & function);

  void setArcConsumptionLowerBoundLinearFunction(const BcArc & bcArc, double initValue, double slope = 0);

  void setArcConsumptionUpperBoundPWLfunction(const BcArc & bcArc, std::vector<std::pair<double, double> > & function);

  void setArcConsumptionUpperBoundLinearFunction(const BcArc & bcArc, double initValue, double slope = 0);

  void setArcIncompatibleValueForSelection(const BcArc & bcArc, int value);

  explicit operator ScalableResource * () const
  {
    return _resourcePtr;
  }

  bool isDefined() const
  {
    return (_resourcePtr != nullptr);
  }
};


class BcVertex
{
  NetworkVertex * _netVertex;
            
public:

  lemon::ListDigraph::Node getLemonNode();
  
  explicit BcVertex(NetworkFlow * networkPtr, bool isFake = false); /// creates a new vertex

  explicit BcVertex(NetworkVertex * netVertex = nullptr): _netVertex(netVertex){}

  int ref() const;
  
  void setElementaritySet(int elemSetId);
  void setPackingSet(int packSetId);
  void setCoveringSet(int covSetId);

  void addToElementaritySet(int elemSetId);
  void addToEnumerationOnlyElementaritySet(int elemSetId);
  void addToPackingSet(int packSetId);
  void addToCoveringSet(int covSetId);

  void addToMemoryOfElemSet(int elemSetId);

  void setBinaryResourceConsumptionLB(int binaryResId, int lowerBound);
  void setBinaryResourceConsumptionUB(int binaryResId, int upperBound);
  void setSpecialResourceConsumptionLB(int binaryResId, int lowerBound); /// for backwards compatibility
  void setSpecialResourceConsumptionUB(int binaryResId, int upperBound); /// for backwards compatibility

  void setName(std::string name_);

  explicit operator NetworkVertex * () const
  { 
    return _netVertex;
  }
  
  /// return true if the BcVertex is defined
  bool defined();  
};

/// TO DO : BcArcInfo does not yet support the case when a vertex or an arc belongs to several elem., pack. or cov. sets
struct BcArcInfo
{
  /// constant fields (initialized in the constructor)
  int arcId;   /// lemon arc id
  int arcElemSetId;
  int arcPackSetId;
  int arcCovSetId;
  int tailVertId;
  int tailElemSetId;
  int tailPackSetId;
  int tailCovSetId;
  int headVertId;
  int headElemSetId;
  int headPackSetId;
  int headCovSetId;
  double * resConsumption;
  double * resConsumptionUB;
  std::string name;

  BcArcInfo(const NetworkArc * arcPtr, double * resConsumption_, double * resConsumptionUB_);
  ~BcArcInfo();
};

class BcArc
{
  NetworkArc * _netArc;
 
public:

  lemon::ListDigraph::Arc getLemonArc();

  BcArc(NetworkFlow * networkPtr, const BcVertex & tail, const BcVertex & head, double originalCost = 0,
        bool isFake = false);

  explicit BcArc(NetworkArc * netArc = nullptr): _netArc(netArc){}
  
  int ref() const;
     
  BcVertex tail() const;
  
  BcVertex head() const;
 
  /// backward compatibility
  void arcVar(const BcVar & newvar);
  
  void addVarAssociation(const BcVar & newvar, double coeff);

  void addAlternativeMapping();

  explicit operator NetworkArc * () const
  { 
    return  _netArc;
  }
  
  /// return true if the BcArc is defined
  bool defined();

  void cost(double val);
  double cost();

  void setElementaritySet(int elemSetId);
  void addToEnumerationOnlyElementaritySet(int elemSetId);
  void setPackingSet(int packSetId);
  void setCoveringSet(int covSetId);

  void addToElementaritySet(int elemSetId);
  void addToPackingSet(int packSetId);
  void addToCoveringSet(int covSetId);

  void setName(std::string name);

  void addToMemoryOfElemSet(int elemSetId);

  void addBinaryResourceConsumption(int binaryResId, int consumption);

};

/// For backward compatibility
#ifdef _WINDOWS
  enum NetAttrMask : int64_t
#else
  enum NetAttrMask
#endif
  {
     PathAttrMask                            = 0x01,
  };

class BcNetwork
{  
  NetworkFlow * _networkPtr;

public:
  BcNetwork() = delete;

  lemon::ListDigraph::Digraph & getLemonDigraph();
  const lemon::ListDigraph::Digraph & getLemonDigraph() const;
  BcArc getBcArc( lemon::ListDigraph::Arc );
  BcVertex getBcVertex( lemon::ListDigraph::Node);
  
  explicit BcNetwork(NetworkFlow * networkFlowPtr);

  /// for backward compatibility
  BcNetwork(BcFormulation & bcForm, NetAttrMask optionalAttrMask, int numElemSets = 0, int numPackSets = 0,
            int numCovSets = 0);

  explicit BcNetwork(BcFormulation & bcForm, int numElemSets = 0, int numPackSets = 0, int numCovSets = 0);

  BcVertex createVertex(bool isFake = false);

  BcArc createArc(int tail, int head, double originalCost = 0.0, bool isFake = false);

  BcVertex getVertex(int id) const;
  
  BcArc getArc(int id) const;
  BcArc getFirstArcBetweenVertices(int tailId, int headId) const;

  int nVertices() const; //linear time complexity
  int nArcs() const; //linear time complexity
  
  void setPathSource(const BcVertex & vertex);

  void setPathSink(const BcVertex & vertex);

  void addToPackingSetCutNeighbourhood(int packSetId, int packSetIdToAdd);

  void setElemSetsDistanceMatrix(const std::vector<std::vector<double> > & matrix);

  void setBinaryResourceNonDisposable(int binaryResId);
  void setSpecialResourceNonDisposable(int binaryResId); /// for backwards compatibility

  void addPermanentRyanAndFosterConstraint(int firstPackSetId, int secondPackSetId, bool together);

  explicit operator NetworkFlow * () const
  { 
    return _networkPtr;
  }    
  
  virtual ~BcNetwork();

};  
        


#endif	/* BCMODELNETWORKFLOW_HPP */
