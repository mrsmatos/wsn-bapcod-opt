/**
 *
 * This file bcModelSolutionC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef ModelSolutionC_h
#define ModelSolutionC_h

#include "bcUsefulHeadFil.hpp"
#include "bcModelVarC.hpp"
#include "bcModelConstrC.hpp"

class Solution;
class DualSolution;
class BcVar;
class BcModel;
class Problem;
class BcFormulation;

namespace bcp_rcsp {
    struct Solution;
}


class BcSolution
{
  friend class BcModel;
  friend class BcFormulation;
  friend class Problem;
public:

  Solution * _solutionPtr;

public:
  BcSolution(const BcFormulation & formulation);
  BcSolution(Solution * SolutionPointer);
  BcSolution clone(); ///< Creates a clone of this solution. Added by Boris
  ~BcSolution();
  /// includes var in the solution at its current value
  BcSolution & includeVarVal(const BcVar & var);
  /// includes var in the solution at its current value, if var is already in the solution then its value
  /// is updated by the current value
  BcSolution & updateVarVal(const BcVar & var);
  BcSolution & updateVarVal(const BcVar & var, const double & val);
  /// appends solution in the chain defined by next
  BcSolution & appendSol(BcSolution & sol);
  BcSolution & includeVarVal(BcVarIndex & varInd) {return includeVarVal((BcVar) varInd);}
  BcSolution & updateVarVal(BcVarIndex & varInd) {return updateVarVal((BcVar) varInd);}

#ifdef BCP_RCSP_IS_FOUND
  void setRCSPsolution(bcp_rcsp::Solution * rcspSolPtr);
  const bcp_rcsp::Solution * getRCSPsolution();
#endif

  void initializeOrderedSolution(const std::vector<double> & resConsumption);
  void addToOrderedSolution(const int & id, const std::vector<double> & resConsumption,
                            const bool & doNotAddSameConsecutiveIds = false);
  const std::vector<int> & orderedIds() const;
#ifdef BCP_RCSP_IS_FOUND
  const std::vector<double> & waitingTimes() const;
#endif
  const std::vector<std::vector<double> > & resConsumption() const;
  void setEnumeratedFlag(const bool flag);
  bool defined() const;

  BcFormulation formulation() const;
  void formulation(const BcFormulation & formulation);
  /// returns how many times this subproblem solution is repeated in the master solution
  int getMultiplicity() const;
  /// returns the cost the solution
  double cost() const;
  double resetCost();
  /// returns the next solution in the chain of subproblem solutions that define a master solution
  BcSolution next() const;
  /// returns the list of non zero variables in the solution whose value can be obtained using solVal()
  void getVar(std::set< BcVar > & varSet);
  /// returns the list of non zero variables in the solution chain whose value can be obtained using solVal()
  void extractVar(std::set< BcVar > & varSet);
  /// returns the list of non-zero variables in an aggregated solution obtained by projecting a master solution
  /// in the original variable place; the list is filtered to retain only those associated with a given generic name
  void extractVar(const std::string & genericName, std::set< BcVar > & varSet);
  std::set< BcVar > extractVar(const std::string & genericName = ""); /// for backwards compatibility
  /// idem but the list is filtered to retain only those associated with a given generic name and formulation
  void extractVar(const std::string & genericName, const BcFormulation & formulation, std::set< BcVar > & varSet);
  /// idem but the list is filtered to retain only those associated with a given generic name and first index
  void extractVar(const std::string & genericName, int firstIndex, std::set< BcVar > & varSet);
  /// returns the value of a variable in the solution
  double getVarVal(BcVar vptr);
  void clear();
  void deleteSolution();
  void deleteSolutionsChain();
  operator Solution *() const;
  std::ostream & print(std::ostream& os) const;
  void printOrderedSolution(std::ostream& os) const;

  const Solution * solutionPtr() const;
  Solution * solutionPtr();

  BcSolution & operator= (const BcVar & var) {return includeVarVal(var);}
  BcSolution & operator+= (const BcVar & var) {return updateVarVal(var);}
  BcSolution & operator= (BcVarIndex & varInd) {return includeVarVal(varInd);}
  BcSolution & operator+= (BcVarIndex & varInd) {return updateVarVal(varInd);}
  BcSolution & operator+= (BcSolution & sol) {return appendSol(sol);}

  friend std::ostream& operator<<(std::ostream& os, const BcSolution & that) {return that.print(os);}
};

class BcDualSolution
{
    friend class Problem;
public:

    DualSolution * _dualSolutionPtr;
public:
    BcDualSolution(DualSolution * SolutionPointer);
    BcDualSolution(const BcFormulation & formulation);
    operator DualSolution *() const;
    BcDualSolution & includeConstr(const BcConstr & constr);
    BcDualSolution & updateConstrVal(const BcConstr & constr);
    BcDualSolution & appendSol(BcDualSolution & sol); /// appends solution in the chain defined by next
    DualSolution * dualSolPtr() const{return _dualSolutionPtr;} /// added by Ruslan, could not use operator DualSolution *() const;
    bool isDefined() const;

    BcDualSolution & includeConstr(const BcConstrIndex & constrInd)
    {
        return includeConstr((const BcConstr&) constrInd);
    }
    BcDualSolution & updateConstrVal(const BcConstrIndex & constrInd)
    {
        return updateConstrVal((const BcConstr&) constrInd);
    }
    BcDualSolution & operator= (const BcConstr & constr) {return includeConstr(constr);}
    BcDualSolution & operator+= (const BcConstr & constr) {return updateConstrVal(constr);}
    BcDualSolution & operator= (const BcConstrIndex & constrInd) {return includeConstr(constrInd);}
    BcDualSolution & operator+= (const BcConstrIndex & constrInd) {return updateConstrVal(constrInd);}
    BcDualSolution & operator+= (BcDualSolution & sol) {return appendSol(sol);}

    std::ostream & print(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& os, const BcDualSolution & that) {return that.print(os);}
};


#endif // ModelSolutionC_h
