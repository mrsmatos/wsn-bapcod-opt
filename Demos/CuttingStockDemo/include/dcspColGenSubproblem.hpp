/**
 *
 * This file dcspColGenSubproblem.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DCSPINTKNAPSACKSUBPROBLEM_H_
#define DCSPINTKNAPSACKSUBPROBLEM_H_

#include "bcKnapsackSolverC.hpp"
#include "dcspDataC.hpp"
#include "bcModelFormulationC.hpp"

bool GreedyKnap(const int & n, 
				const int & ni,
				const std::vector<double> & p,
				const std::vector<double> & w,
				std::vector<int> & ci,
                std::vector<KnpItem> & vit, 
                const double & minw,
                const double & cap,
                double &inc);

class IntKnapsackSubProb: public BcSolverOracleFunctor
{
  DcspData * _dataPtr;
  bool _firstCallInitializationDone;
  std::vector<bool> _flagOrderInUseV;
  std::vector<int> _boundV;
  std::vector<int> _maxciV;
  std::vector<double> _weightV;
  std::vector<double> _profitV;
  std::vector<KnpItem> _vitems;
  Double _minweight;

public:
  IntKnapsackSubProb(DcspData* dataPtr);
  
  virtual ~IntKnapsackSubProb(){}
  
  virtual bool operator()(BcFormulation spPtr,
                          int colGenPhase,
			              double & objVal,
			              double & dualBound,
			              BcSolution & primalSol);
};




#endif /* DCSPINTKNAPSACKSUBPROBLEM_H_ */
