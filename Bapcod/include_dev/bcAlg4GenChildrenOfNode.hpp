/**
 *
 * This file bcAlg4GenChildrenOfNode.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCGENCHILDNODESALG_HPP_
#define BCGENCHILDNODESALG_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcProbConfigC.hpp"


struct GenChildNodesInfo
{
  int numberOfNodes;

  GenChildNodesInfo() :
      numberOfNodes(0)
  {
  }

  virtual ~GenChildNodesInfo()
  {
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const
  {
    os << "GenChildNodesInfo with number of Nodes = " << numberOfNodes << std::endl;
    return os;
  }

};

class Alg4GenChildrenOfNode
{

protected:
  MasterCommons4GenChildNodesAlgorithm & _masterCommons;

  Node * _currentNodePtr;

  BapcodInit & bapcodInit() const
  {
    return _currentNodePtr->bapcodInit();
  }

  ControlParameters & param() const
  {
    return bapcodInit().param();
  }

  ProgStatus & progStatus() const
  {
    return bapcodInit().progStatus();
  }

public:

  Alg4GenChildrenOfNode(MasterCommons4GenChildNodesAlgorithm & masterCommons) :
      _masterCommons(masterCommons), _currentNodePtr(NULL)
  {
  }

  virtual ~Alg4GenChildrenOfNode()
  {
  }

  virtual bool setupAlgo(Node * nodePtr)
  {
    _currentNodePtr = nodePtr;
    return false;
  }

  virtual void run(int & globalTreatOrder)
  {
  }

  virtual void setDownAlgo()
  {
    _currentNodePtr = NULL;
  }

};

#endif /* BCGENCHILDNODESALG_HPP_ */
