/**
 *
 * This file bcModelBranchingConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGBranchingCONSTRC_H_
#define BCMODELINGBranchingCONSTRC_H_

#include "bcModelConstrC.hpp"
#include "bcModelFormulationC.hpp"
#include "bcModelParameterC.hpp"

class GenericConstr;
class InstanciatedConstr;
struct BcConstrIndex;
class BcFormulation;
class BcSolution;
class BcConstr;
class GenericBranchingConstr;
class GenAggrSubProbVarBranchingConstr;
class GenPathsPerNetworkBranchingConstr;
class GenPackSetAssignBranchingConstr;
class RyanAndFosterInstSubProbBranchConstr;


class BcRyanAndFosterBranchConstr
{
protected:
  RyanAndFosterInstSubProbBranchConstr * _iconstrPtr;
  BcVar _bcVarI;
  BcVar _bcVarJ;
public:
  BcRyanAndFosterBranchConstr(RyanAndFosterInstSubProbBranchConstr * iconstrPtr);
  const BcVar & bcVarI() const;
  const BcVar & bcVarJ() const;
  bool together() const;
};

class BcDisjunctiveBranchingConstrSeparationFunctor
{
public:
  BcDisjunctiveBranchingConstrSeparationFunctor(){}
  virtual ~BcDisjunctiveBranchingConstrSeparationFunctor(){}

  /// second element in the pair of returnBrConstrList is the string needed for
  /// 1) storing branching history for the same branching decision
  /// (thus, string for the same branching decision should be the same)
  /// 2) screen output
  virtual bool operator() (BcFormulation formPtr,
			               BcSolution & primalSol,
                           std::list<std::pair<double, BcSolution>> & columnsInSol,
			               const int & candListMaxSize,
			               std::list<std::pair<BcConstr, std::string> > & returnBrConstrList) = 0;
};

class BcVarBranchingPriorityFunctor
{
public:
  BcVarBranchingPriorityFunctor() = default;
  virtual ~BcVarBranchingPriorityFunctor() = default;
  virtual double operator() (const MultiIndex & index, const double & intPart, const double & fracPart) = 0;
};

class BcAggrSubProbVarBranching
{
  GenAggrSubProbVarBranchingConstr * _genAggrSubProbVarBranchingConstrPtr;

public:
  BcAggrSubProbVarBranching(const BcFormulation & formulation,
		                    const std::string & genVarName,
                            const double & highestPriorityFraction = 0.5,
		                    const double & priorityLevel = 1.0,
                            const int & numIgnoredIndices = 0,
		                    const bool & toBeUsedInPreprocessing = true);

  virtual ~BcAggrSubProbVarBranching();

  const BcAggrSubProbVarBranching & attach(BcVarBranchingPriorityFunctor * priorityRoutinePtr);
};

class BcBranchingConstrArray: public BcConstrArray
{
  GenericBranchingConstr * _genericBranchingConstrPtr;

public:
  BcBranchingConstrArray(const BcFormulation & formulation,
		                 const std::string & name,
		                 const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                         const double & priorityLevel = 1.0,
		                 const bool & toBeUsedInPreprocessing = true);

  ~BcBranchingConstrArray() override;

  GenericBranchingConstr * genericBranchingConstrPtr() const;

  const BcBranchingConstrArray & attach(BcDisjunctiveBranchingConstrSeparationFunctor * separationRoutinePtr);

};

std::ostream & operator<<(std::ostream& os, const BcBranchingConstrArray & that);


#endif /* BCMODELINGBranchingCONSTRC_H_ */
