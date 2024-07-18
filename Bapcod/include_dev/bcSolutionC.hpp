/**
 *
 * This file bcSolutionC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef SolutionClass_h
#define SolutionClass_h

#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcVarConstrC.hpp"

class DualSolution;
class ColGenSpConf;
class BcVar;

#ifdef BCP_RCSP_IS_FOUND
namespace bcp_rcsp {
    struct Solution;
}
#endif /* BCP_RCSP_IS_FOUND */

class Solution
{
  friend class BcFormulation;
 protected:
  ProbConfig * _probConfPtr;
  int _ref;
  Double _cost;

  /**
   * the multiplicity is number of copy of this solution (or is associated master column)
   * that is used in defing a global solution to the master
   */
  int _multiplicity;

  Solution * _previousSolPtr;
  Solution * _nextSolPtr;

  VarPtr2DoubleMap _solVarValMap;
  /// this vector here is for keeping an additional information about the pricing problem solution
  /// in the case of paths, we keep the ids of arcs of the path
  std::vector<int> _orderedIds;
  std::vector<std::vector<double> > _resConsumption;

#ifdef BCP_RCSP_IS_FOUND
  const bcp_rcsp::Solution * _rcspSolPtr;
#endif /* BCP_RCSP_IS_FOUND */

  bool _enumeratedFlag;

 public:
  Solution(ProbConfig * probConfPtr = NULL, Solution * previousSolPtr = NULL);
  Solution(const Solution & that);
  Solution(ProbConfig * probConfPtr, const VarPtr2DoubleMap  & solVarValMap);
  virtual ~Solution();
  ProbConfig * probConfPtr() const {return _probConfPtr;}
  int ref() const {return _ref;}
  void probConfigPtr(ProbConfig * probConfPtr)  {_probConfPtr =  probConfPtr;}

  void deleteSolutionsChain();
  int solutionsChainSize();

  Solution * previousSolPtr() const {return _previousSolPtr;}
  Solution * nextSolPtr() const {return _nextSolPtr;}

  /// Virtual copy constructor
  virtual Solution * clone() const;
  virtual void includeVarSet(const VarPtrSet & varPtrSet);
  void getVar(std::set< BcVar > & varSet);
  void extractVar(std::set< BcVar > & varSet);
  void extractVarWithGenericName(const std::string & name, std::set< BcVar > & varSet);
  void extractVarWithGenericName(const std::string & name, ProbConfig * probConfPtr, std::set< BcVar > & varSet);
  void extractVarWithGenericName(const std::string & name, int firstIndex, std::set< BcVar > & varSet);
  const Double & solVal(Variable * vptr);
  void includeVars(const VarPtr2DoubleMap & varValMap, const bool & cumulativeVal);
  void includeVar(Variable * varPtr, const Double & val, const bool & cumulativeVal);
  virtual const VarPtr2DoubleMap & solVarValMap() const;
  void clear();
  int size() {return _solVarValMap.size();}
  bool empty() {return _solVarValMap.empty();}
  void addToOrderedIds(const int & vertexId)
  {
    _orderedIds.push_back(vertexId);
  }
  void addToResConsumption(const std::vector<double> & resConsumption)
  {
    _resConsumption.push_back(resConsumption);
  }
  inline const std::vector<int> & orderedIds() const
  {
    return _orderedIds;
  }
  inline const std::vector<std::vector<double> > & resConsumption() const
  {
    return _resConsumption;
  }

#ifdef BCP_RCSP_IS_FOUND
  void setRcspSolPtr(const bcp_rcsp::Solution * solPtr);
  inline const bcp_rcsp::Solution * rcspSolPtr() const {return _rcspSolPtr;}
  const bcp_rcsp::Solution * copyRcspSolPtr() const;
#endif /*  BCP_RCSP_IS_FOUND */

  inline bool enumeratedFlag() const {return _enumeratedFlag;}
  void enumeratedFlag(bool flag);
  virtual const Double & cost() const {return _cost;}
  virtual const Double & resetCost();
  virtual const int & multiplicity()  const {return _multiplicity;}
  virtual void multiplicity(const int & mult) {_multiplicity = mult;}
  virtual void cost(const Double & c)  {_cost = c;}
  virtual void previousSolPtr(Solution * prevSolPtr);
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void shortPrint(std::ostream& os = std::cout) const;
  virtual void printOrderedSolution(std::ostream& os = std::cout) const;
  virtual std::ostream& printVar(std::ostream& os = std::cout) const;
  virtual bool operator<(const Solution  & that) const;
  virtual void printDetailedSolution(std::ostream& os = std::cout) const;

};

inline std::ostream& operator<<(std::ostream& os, const Solution & m)
{return m.print(os);}

class DualSolution
{
 protected:
  ProbConfig * _probConfPtr;
  int _ref;
  Double _rhs;
  ConstrPtr2DoubleMap _dualSolConstrValMap;
  DualSolution * _previousSolPtr;
  DualSolution * _nextSolPtr;

 public:
  DualSolution(ProbConfig * probConfPtr = NULL);
  DualSolution(const DualSolution & that);
  DualSolution(const ConstrPtr2DoubleMap & dualSolConstrValMap);
  virtual ~DualSolution() {}
  ProbConfig * probConfPtr() const {return _probConfPtr;}
  int ref() const {return _ref;}
  void probConfigPtr(ProbConfig * probConfigPtr) {_probConfPtr = probConfigPtr;}
  const ConstrPtr2DoubleMap & dualSolConstrValMap() 
  {
    return _dualSolConstrValMap;
  }
  virtual DualSolution * clone() const;
  const Double & computeTrueRhs();
  virtual const Double & rhs() const {return _rhs;}
  virtual void rhs(const Double & c)  {_rhs = c;}
  virtual void includeConstr(const ConstrPtrSet & constrPtrSet);
  virtual void includeConstr(Constraint * constrPtr, const Double & val, const bool & cumulativeVal);
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual std::ostream& printConstr(std::ostream& os = std::cout) const;
  DualSolution * previousSolPtr() const {return _previousSolPtr;}
  DualSolution * nextSolPtr() const {return _nextSolPtr;}
  virtual void previousSolPtr(DualSolution * prevSolPtr);
  const int size() const {return _dualSolConstrValMap.size();} //(hsen) to calculate dualSol density
};

#endif // SolutionClass_h
