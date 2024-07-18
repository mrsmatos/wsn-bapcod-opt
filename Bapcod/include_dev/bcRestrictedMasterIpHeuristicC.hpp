/**
 *
 * This file bcRestrictedMasterIpHeuristicC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCRESTRICTEDMASTERIPHEURISTIC_HPP_
#define BCRESTRICTEDMASTERIPHEURISTIC_HPP_

#include "bcUsefulHeadFil.hpp"
#include "bcNodeC.hpp"
#include "bcAlg4PreprocessingOfNode.hpp"
#include "bcAlg4GenChildrenOfNode.hpp"
#include "bcAlg4ProblemSetup.hpp"
#include "bcAlg4EvalByMip.hpp"

class RestrictedMasterIpHeuristic : public Alg4PrimalHeuristicOfNode
{
    bool _activateAllColumns;
    std::list<BranchingConstrBaseType *> _localNodeBrConstrList;

  Node * createRootNode();
  bool prepareNodeForTreatment(Node * nodePtr, const int globalNodesTreatOrder);
  virtual void runBody(int & globalTreatOrder);

public:

    void setOptionActivateAllColumns(const int value)
    {
        _activateAllColumns = value;
    }

  RestrictedMasterIpHeuristic(Problem * probPtr, MasterCommons4PrimalHeuristic & masterCommons);

  virtual ~RestrictedMasterIpHeuristic();
};

#endif /* BCRESTRICTEDMASTERIPHEURISTIC_HPP_ */
