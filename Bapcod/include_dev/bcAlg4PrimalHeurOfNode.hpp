/**
 *
 * This file bcAlg4PrimalHeurOfNode.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPrimalHeuristic_HPP_
#define BCPrimalHeuristic_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcProbConfigC.hpp"

class Alg4PrimalHeuristicOfNode
{
  virtual void runBody(int & globalTreatOrder){}

protected:
  MasterCommons4PrimalHeuristic & _masterCommons;

  std::vector<Node *> _generatedNodes;
  Node * _currentBaPNodePtr;
  Problem * _problemPtr;
  int _priority;

  BapcodInit & bapcodInit() const;

  ControlParameters & param() const
  {
    return bapcodInit().param();
  }

  ProgStatus& progStatus() const
  {
    return bapcodInit().progStatus();
  }

public:

  Alg4PrimalHeuristicOfNode(Problem * probPtr,
      MasterCommons4PrimalHeuristic & masterCommons) :
      _masterCommons(masterCommons), _generatedNodes(), _currentBaPNodePtr(NULL),
          _problemPtr(probPtr), _priority(1)
  {
  }

  virtual ~Alg4PrimalHeuristicOfNode();

  void setOptionPriority(int priority)
  {
    _priority = priority;
  }
  
  const int & priority() const 
  {
    return _priority;
  }

  void run(Node * nodePtr, int & globalTreatOrder);

};

struct PrimalHeurSort
{
  bool operator()(Alg4PrimalHeuristicOfNode * a, Alg4PrimalHeuristicOfNode * b) const
  {
    return a->priority() <= b->priority();
  }
};



#endif /* BCPrimalHeuristic_HPP_ */
