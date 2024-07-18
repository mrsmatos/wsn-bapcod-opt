/**
 *
 * This file dcspCutSeparationRoutines.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DcspCutSeparationRoutine_H_
#define DcspCutSeparationRoutine_H_

#include "dcspDataC.hpp"
#include "bcModelCutConstrC.hpp"

class CardinalityCutSeparationRoutine: public BcCutSeparationFunctor
{
  DcspData * _dataPtr;
public:
  CardinalityCutSeparationRoutine(DcspData * dataPtr) : _dataPtr(dataPtr) {}
  virtual ~CardinalityCutSeparationRoutine(){}
  virtual int operator() (BcFormulation formPtr,
			  BcSolution & primalSol,
			  double & maxViolation,
			  std::list< BcConstr > & cutList);
};


#endif /* DcspCutSeparationRoutine_H_ */
