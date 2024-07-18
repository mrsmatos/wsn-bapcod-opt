/**
 *
 * This file dcspCutSeparationRoutines.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dcspCutSeparationRoutines.hpp"
#include "dcspColGenSubproblem.hpp"
#include <vector>

using std::vector;



int CardinalityCutSeparationRoutine::operator() (BcFormulation formPtr,
						                         BcSolution & primalSol,
						                         double & maxViolation,
						                         std::list< BcConstr > & cutList)
{

  bool doPrint(false);

  if (doPrint) 
    std::cout <<"CardinalityCutSeparationRoutine ENTERING formPtr = "  << formPtr << std::endl;

  Double totalFractionalCardinality(0);

  std::set< BcVar > yvarList;
  primalSol.extractVar("Y", yvarList);

  for (std::set< BcVar >::const_iterator varIt = yvarList.begin(); varIt != yvarList.end(); ++varIt)
    {
      double varVal = (*varIt).solVal();
      totalFractionalCardinality += varVal;
      if (doPrint) 
	    std::cout << "CardinalityCutSeparationRoutine Y current value = "  << varVal << std::endl;
    }

  if (totalFractionalCardinality.fractional(0.001)) 
    {
      double lhs = Dceil(totalFractionalCardinality);

      BcCutConstrArray CardinalityConstr(formPtr, "CARD");

      if (CardinalityConstr.isDefinedAt(MultiIndex(0)))
	    {
	      std::cout <<"CardinalityCutSeparationRoutine: cut has already been added " << std::endl;
	      return 0; /// cut has already been added.
	    }
      BcConstr instantiatedCutConstr = CardinalityConstr(0);
      instantiatedCutConstr >= lhs;

      for (std::list< BcFormulation >::const_iterator formIt = formPtr.colGenSubProblemList().begin();
	       formIt != formPtr.colGenSubProblemList().end(); ++formIt)
	    {
	      BcVarArray Y(*formIt,"Y");
	      instantiatedCutConstr += Y[0];
	    }

      maxViolation = lhs - totalFractionalCardinality;

      if (doPrint) 
        std::cout <<"CardinalityCutSeparationRoutine violation = "  << maxViolation << std::endl;

      cutList.push_back(instantiatedCutConstr);

      return 1;

    }

  return 0;


}
