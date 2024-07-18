/**
 *
 * This file bcModelVarC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGVARC_H_
#define BCMODELINGVARC_H_

#include <map>
#include <vector>

#include "bcMultiIndexC.hpp"
#include "bcModelParameterC.hpp"

#include <stdarg.h>

class GenericVar;
class ProbConfig;
class InstanciatedVar;
class BcModel;
class BcConstr;
class BcObjective;
class BcSolution;
class BcVarArray;
class BcConstrArray;
class BcFormulation;
struct BcVarIndex;
class BcRowExpression;

class BcVar
{
  friend class BcConstr;
  friend class BcObjective;
  friend class BcSolution;
  friend class BcVarArray;
  friend class BcConstrArray;  
  friend class BcArc;  
  friend class NetworkArc;
  friend class BcNetworkFlowArcVarArray;
protected:
  InstanciatedVar * _ivarPtr;
public:
  BcVar(InstanciatedVar * const ivarPtr);
  BcVar(BcVarIndex & varInd);
  BcVar();
  ~BcVar();
  bool isDefined() const;
  const std::string & name() const;
  const std::string & genericName() const;
  const BcFormulation formulation() const;
  const MultiIndex & id() const;
  double curLb() const;
  double curUb() const;
  /// modified cost in the current problem (could represent a reduced cost)
  virtual double curCost() const;
  /// originally defined cost for this variable.
  double originalCost() const;

  /// true if the variable is active in the current formulation
  bool inCurForm() const;

  /// to assign variable value
  void curVal(double value);
  double curVal() const;
  inline const BcVar & operator=(double value)
  {
      curVal(value);
      return *this;
  }
  /// variable value in the currently recorded solution
  double solVal() const;

  /// lb local to the subproblem formulation
  BcVar & localLb(double lb);
  /// ub local to the subproblem formulation
  BcVar & localUb(double ub);
  /// ub local to the subproblem formulation
  inline BcVar & operator<=(double ub)
  {
    return localUb(ub);
  }
  /// lb local to the subproblem formulation
  inline BcVar & operator>=(double lb)
  {
    return localLb(lb);
  }

  /// global ub valid in the master on the sum of such subproblem variables
  BcVar & globalLb(double lb);
  /// global ub valid in the master on the sum of such subproblem variables
  BcVar & globalUb(double ub);

//  inline const BcVar & operator==(double value)
//  {
//    localLb(value);
//    localUb(value);
//    type('F');
//    return *this;
//  }

  /// type 'I' = integer (= default)
  /// type 'B' = binary,  
  /// type 'C' = continuous,
  void type(const char & flag);

  /// sense 'P' = positive (= default)
  /// sense  'N' = negative
  /// sense 'F' = free
  void sense(const char & flag);

  void setImplicit();

  /// A higher priority means that variable is selected first for branching
  BcVar & branchingPriority(double priority);

  /// 'U' or 'D' branch treated in priority
  BcVar & branchingDirective(char directive);

  /// to create a pair < Variable,  coefficient >
  BcVarCoef addCoef(double coef);
  BcVarCoef operator*(double coef);

  /// to add a  pair < Variable,  coefficient = 1 > in an expression
  BcRowExpression addCoef(const BcVar & var);
  BcRowExpression operator+(const BcVar & var);
 
  /// to add a  pair < Variable,  coefficient > in an expression
  BcRowExpression addCoef(const BcVarCoef & varCoef);
  BcRowExpression operator+(const BcVarCoef & varCoef);


  /// to add a  pair < Variable,  coefficient = 1 > in an expression
  BcRowExpression addCoef(BcVarIndex & varIndex);
  BcRowExpression operator+(BcVarIndex & varIndex);

  /// to cast a BcVar into InstanciatedVar * for advanced users
  operator InstanciatedVar *() const
  {
    return _ivarPtr;
  }

  /// to sort BcVar
  bool operator<(const BcVar & that) const;
};


class BcVarArray
{
  GenericVar * _genericVarPtr;
  /// cache memory for operator()
  BcVar _curInstVarPtr;
  
public:

  //If the max indices are provided direct access will be used for instanciedVars instead of map.
  BcVarArray(const BcFormulation & formulation, const std::string & name, int firstIndexMax = -1,
             int secondIndexMax = -1, int thirdIndexMax = -1);
  
  BcVarArray(GenericVar * genericVarPtr = NULL);

  virtual ~BcVarArray();

  const std::string & genericName() const;
  GenericVar * genericVarPtr() const
  {
    return _genericVarPtr;
  }
  /// returns the formulation to which the BcVarArray belongs
  BcFormulation formulation() const;
  /// returns the model to which the BcVarArray belongs
  const BcModel& model() const;
  /// returns the model to which the BcVarArray belongs
  BcModel & model();
  
  // the last accessed BcVar.
  BcVar & curInstVarPtr()
  {
    return _curInstVarPtr;
  }

  BcVar & getElement(const MultiIndex & indexArray);
  /// call getElement, but creates it if does not exists
  BcVar & createElement(const MultiIndex & indexArray);
  /// creations  calls addVarFunctor if defined
  
  inline BcVar & createElement(int firstIndex = 0, int secondIndex = -1, int thirdIndex = -1,
                               int fourthIndex = -1, int fifthIndex = -1, int sixthIndex = -1,
                               int seventhIndex = -1)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }
  
  inline BcVar & operator()(int firstIndex)
  {
    return createElement(MultiIndex(firstIndex));
  }
  
  inline BcVar & operator()(int firstIndex, int secondIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex));
  }
  
  inline BcVar & operator()(int firstIndex, int secondIndex, int thirdIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex));
  }

  inline BcVar & operator()(int firstIndex, int secondIndex,
                            int thirdIndex, int fourthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex));
  }
  
  inline BcVar & operator()(int firstIndex, int secondIndex,
                            int thirdIndex, int fourthIndex,
                            int fifthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex));
  }
  
  inline BcVar & operator()(int firstIndex, int secondIndex,
                            int thirdIndex, int fourthIndex,
                            int fifthIndex, int sixthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex));
  }
  
  inline BcVar & operator()(int firstIndex, int secondIndex,
                            int thirdIndex, int fourthIndex,
                            int fifthIndex, int sixthIndex,
                            int seventhIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }

  inline BcVar & operator()(const MultiIndex & indexArray)
  {
    return createElement(indexArray);
  }

  /// call getElement in the end of the sequence of operator[]
  BcVarIndex operator[](const int & index);
  
  /// checks whether the array has an element defined at index.
  bool isDefinedAt(const MultiIndex & indexArray);  

  /// set a branching rule common to all elements of this BcVarArray
  const BcVarArray & ruleForBranchingVarSelection(const SelectionStrategy & branchingPriorityRule);

  /// Default priority level for variables of this class for componentBoundBranching: 
  /// a higher priority means that variable is selected first for branching
  /// a negative (or null) priority means that variable is not selected for branching in subproblem
  const BcVarArray & priorityForMasterBranching(double priorityValue);

  const BcVarArray & priorityForSubproblemBranching(double priorityValue);

  const BcVarArray & priorityForRyanFosterBranching(double priorityValue);

  /// Default branching direction for variables of this class: 'U' or 'D' branch treated in priority
  const BcVarArray & branchingDirective(const char & defaultDirective);

  /**
   * Set a new index reference letter.
   * @param multiIndexNames index reference lette
   * @return the current BcVarArray with the new index.
   */
  const BcVarArray & defineIndexNames(const MultiIndexNames & multiIndexNames);
  const BcVarArray & defineIndexName(const MultiIndexNames & multiIndexNames) ///backward compatibility
  {
      return defineIndexNames(multiIndexNames);
  }

  /// global ub valid in the master on the sum of such subproblem variables
  const BcVarArray & globalUb(double vdefaultGlobUb);

  /// global ub valid in the master on the sum of such subproblem variables
  const BcVarArray & globalLb(double vdefaultGlobLb);

  /// ub local to the subproblem formulation
  const BcVarArray & localUb(double defaultLocUb);
  /// lb local to the subproblem formulation
  const BcVarArray & localLb(double defaultLocLb);
  inline const BcVarArray & operator<=(double vdefaultUb)
  {
    return localUb(vdefaultUb);
  }
  inline const BcVarArray & operator>=(double vdefaultLb)
  {
    return localLb(vdefaultLb);
  }

  /// assigned a default value to all elements of this BcVarArray
  const BcVarArray & defaultVal(double defaultValue);
//  inline const BcVarArray & operator==(double defaultFixVal)
//  {
//    localLb(defaultFixVal);
//    localUb(defaultFixVal);
//    return type('F');
//  }
  const BcVarArray & type(const char & defaultType);
  const BcVarArray & sense(const char & defaultSense);
  const BcVarArray & setImplicit();

};

std::ostream& operator<<(std::ostream& os, const BcVarArray & that);

struct BcVarIndex: public std::pair<BcVarArray, MultiIndex>
{
  BcVarIndex(const BcVarArray & modvar, const MultiIndex & indexArray);
  BcVarIndex(const BcVarArray & modvar);

  /// return true if the BcConstrIndex is
  bool isDefined();
  /// call getElement in the end of the sequence of operator[]
  BcVarIndex & operator[](const int & index);
  const BcVar & operator=(double value);
  BcVar & operator<=(double ub);
  BcVar & operator>=(double lb);

  /// type 'I' = integer (= default)
  /// type 'B' = binary,  
  /// type 'C' = continuous,
  /// sense 'P' = positive (= default)
  /// sense  'N' = negative
  /// sense 'F' = free 
  void type(const char & flag);

  BcVarCoef addCoef(double coef);
  BcVarCoef operator*(double coef);

  BcRowExpression addCoef(const BcVar & var);
  BcRowExpression operator+(const BcVar & var);

  BcRowExpression addCoef(const BcVarCoef & varCoef);
  BcRowExpression operator+(const BcVarCoef & varCoef);

  BcRowExpression addCoef(BcVarIndex & varIndex);
  BcRowExpression operator+(BcVarIndex & varIndex);

  const MultiIndex & id();
  double curLb();
  double curUb();
  double curCost();
  /// variable value in the current state of the algorithm
  double curVal();
  /// variable value in the currently recorded solution
  double solVal();
  bool inCurForm();

  //friend std::ostream& operator<<(std::ostream& os, BcVarIndex & that);

  BcVar & branchingPriority(double priority);

  /// 'U' or 'D' branch treated in priority
  BcVar & branchingDirective(char directive);

  operator BcVar();
};



struct BcVarCoef: public std::pair<BcVar, double>
{
  BcVarCoef(const BcVar & var, double coef);
  BcVarCoef(BcVarIndex & modVarInd, double coef);
  BcVarCoef(const BcVar & var);
  BcVarCoef(BcVarIndex & modVarInd);

  /// to add a  pair < Variable,  coefficient = 1 > in an expression
  BcRowExpression addCoef(const BcVar & var);
  BcRowExpression operator+(const BcVar & var);

  /// to add a  pair < Variable,  coefficient> in an expression
  BcRowExpression addCoef(const BcVarCoef & varCoef);
  BcRowExpression operator+(const BcVarCoef & varCoef);

  /// to add a  pair < Variable,  coefficient = 1 > in an expression
  BcRowExpression addCoef(BcVarIndex & varIndex);
  BcRowExpression operator+(BcVarIndex & varIndex);

};

class BcRowExpression: public std::list<BcVarCoef>
{
  friend class BcConstr;
  friend class BcObjective;
  double _mult;
public:
  BcRowExpression();
  const BcRowExpression & addCoef(const BcVarCoef & coef);
  const BcRowExpression & removeCoef(BcVarCoef & coef);

  inline const BcRowExpression & operator+=(const BcVarCoef & coef)
  {
    return addCoef(coef);
  }
  inline const BcRowExpression & operator-=(BcVarCoef & coef)
  {
    return removeCoef(coef);
  }
  const BcRowExpression & operator=(const BcVarCoef & coef)
  {
    return addCoef(coef);
  }
  inline const BcRowExpression & operator+(const BcVarCoef & coef)
  {
    return addCoef(coef);
  }
  inline const BcRowExpression & operator-(BcVarCoef & coef)
  {
    return removeCoef(coef);
  }

  const BcRowExpression & addCoef(const BcVar & coef);
  const BcRowExpression & removeCoef(const BcVar & coef);
  inline const BcRowExpression & operator=(const BcVar & var)
  {
    return addCoef(var);
  }
  const BcRowExpression & operator+(const BcVar & var)
  {
    return addCoef(var);
  }
  const BcRowExpression & operator-(const BcVar & var)
  {
    return removeCoef(var);
  }
  const BcRowExpression & operator+=(const BcVar & var)
  {
    return addCoef(var);
  }
  const BcRowExpression & operator-=(const BcVar & var)
  {
    return removeCoef(var);
  }

  const BcRowExpression & addCoef(BcVarIndex & coef);
  const BcRowExpression & removeCoef(BcVarIndex & coef);

  const BcRowExpression & operator=(BcVarIndex & modVarInd)
  {
    return addCoef(modVarInd);
  }
  inline const BcRowExpression & operator+(BcVarIndex & modVarInd)
  {
    return addCoef(modVarInd);
  }
  inline const BcRowExpression & operator-(BcVarIndex & modVarInd)
  {
    return removeCoef(modVarInd);
  }
  inline const BcRowExpression & operator+=(BcVarIndex & modVarInd)
  {
    return addCoef(modVarInd);
  }
  inline const BcRowExpression & operator-=(BcVarIndex & modVarInd)
  {
    return removeCoef(modVarInd);
  }

  BcRowExpression operator*(const double mult) const;
  friend BcRowExpression operator*(const double mult, const BcRowExpression & exp);
};

inline BcRowExpression operator+(const BcVar & var1, const BcVar & var2)
{
  return var1 + var2;
}

inline BcRowExpression operator+(const BcVarCoef & var1, const BcVarCoef & var2)
{
  return var1 + var2;
}
inline BcVarCoef BcVar::operator*(double coef)
{
  return addCoef(coef);
}

inline BcRowExpression BcVar::operator+(const BcVar & var)
{
  return addCoef(var);
}
inline BcRowExpression BcVar::operator+(const BcVarCoef & varCoef)
{
  return addCoef(varCoef);
}

inline BcRowExpression BcVar::operator+(BcVarIndex & varIndex)
{
  return addCoef(varIndex);
}

inline BcRowExpression BcVarCoef::operator+(const BcVar & var)
{
  return addCoef(var);
}

inline BcRowExpression BcVarCoef::operator+(const BcVarCoef & varCoef)
{
  return addCoef(varCoef);
}


inline BcRowExpression BcVarCoef::operator+(BcVarIndex & varIndex)
{
  return addCoef(varIndex);
}

inline BcVarCoef operator*(double coef, const BcVar & var)
{
  return BcVarCoef(var, coef);
}

inline BcVarCoef operator*(double coef, BcVarIndex & modVarInd)
{
  return BcVarCoef(modVarInd, coef);
}


#endif /* BCMODELINGVARC_H_ */
