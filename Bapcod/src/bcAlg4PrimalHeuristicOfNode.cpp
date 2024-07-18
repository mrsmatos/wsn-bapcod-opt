/**
 *
 * This file bcAlg4PrimalHeuristicOfNode.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcModelC.hpp"

BapcodInit & Alg4PrimalHeuristicOfNode::bapcodInit() const
{
  return _currentBaPNodePtr->bapcodInit();
}

Alg4PrimalHeuristicOfNode::~Alg4PrimalHeuristicOfNode()
{
  for (std::vector<Node *>::iterator it = _generatedNodes.begin();
      it != _generatedNodes.end(); ++it)
    delete *it;
}

void Alg4PrimalHeuristicOfNode::run(Node * nodePtr, int & globalTreatOrder)
{
  _currentBaPNodePtr = nodePtr;

  runBody(globalTreatOrder);

  _currentBaPNodePtr = NULL;
}