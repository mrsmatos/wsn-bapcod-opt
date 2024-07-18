/**
 *
 * This file bcModelConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGCONSTRC_H_
#define BCMODELINGCONSTRC_H_

#include "bcModelVarC.hpp"

class GenericConstr;
class InstanciatedConstr;
struct BcConstrIndex;
class BcFormulation;
class BcConstr;

#ifdef _MSC_VER
#define EXPORTED_CONSTR  __declspec( dllexport )
#else
#define EXPORTED_CONSTR __attribute__ ((visibility ("default")))
#endif

class EXPORTED_CONSTR BcConstr
{
  friend class BcVar;
  friend class BcVarArray;
  friend class BcConstrArray;
  friend class BcCutConstrArray;
  friend class BcBranchingConstrArray;
  friend class BcDualSolution;
protected:
  InstanciatedConstr * _iconstrPtr;
public:
  BcConstr() = delete;
  BcConstr(InstanciatedConstr * iconstrPtr);
  /// return true if the Constr is not NULL
  bool isDefined() const;
  /// name of the class of constraint
  const std::string & genericName() const;
  /// print the constraint to the standart output in a compact form
  void nicePrint() const;
  /// name of the specific constraint within the class
  const std::string & name() const;
  /// multiple index identification
  const MultiIndex & id() const;
  /// current RHS value
  double curRhs() const;
  /// current dual solution value
  double curDualVal() const;
  const BcFormulation formulation() const;
  /// whether or not the constraint is currently in the active formulation
  bool inCurProb() const;
  const BcConstr & add(const BcVarCoef & coef);
  const BcConstr & remove(const BcVarCoef & coef);
  const BcConstr & add(const BcRowExpression & exp);
  const BcConstr & remove(const BcRowExpression & exp);
  const BcConstr & operator+=(const BcRowExpression & exp)
  { return add(exp);}
  const BcConstr & operator=(const BcRowExpression & exp)
  { return add(exp);}
  const BcConstr & operator+(const BcRowExpression & exp)
  { return add(exp);}
  const BcConstr & operator-=(const BcRowExpression & exp)
  { return remove(exp);}
  const BcConstr & operator-(const BcRowExpression & exp)
  { return remove(exp);}
  const BcConstr & operator+(const BcVarCoef & coef)
  { return add(coef);}
  const BcConstr & operator-(const BcVarCoef & coef)
  { return remove(coef);}
  const BcConstr & operator+=(const BcVarCoef & coef)
  { return add(coef);}
  const BcConstr & operator-=(const BcVarCoef & coef)
  { return remove(coef);}
  const BcConstr & operator+(const BcVar & var)
  { return add(BcVarCoef(var,1));}
  const BcConstr & operator-(const BcVar & var)
  { return add(BcVarCoef(var,-1));}
  const BcConstr & operator+=(const BcVar & var)
  { return add(BcVarCoef(var,1));}
  const BcConstr & operator-=(const BcVar & var)
  { return add(BcVarCoef(var,-1));}
  const BcConstr & operator+(BcVarIndex & modVarInd);
  const BcConstr & operator-(BcVarIndex & modVarInd);
  const BcConstr & operator+=(BcVarIndex & modVarInd);
  const BcConstr & operator-=(BcVarIndex & modVarInd);

  /**
   * @param senseV 'L' (Less or equal), 'G' (Greater or equal) or 'E' (equal)
   */
  void sense(const char & sense);
  void rhs(double rhs);

  const BcConstr & operator<=(double rhs)
  {
    this->sense('L');
    this->rhs(rhs);
    return *this;
  }
  const BcConstr & operator<=(Double rhs)
  {
    return operator<=(rhs.val());
  }

  const BcConstr & operator<=(int rhs)
  {
    return operator<=((double) rhs);
  }

  const BcConstr & operator>=(double rhs)
  {
    this->sense('G');
    this->rhs(rhs);
    return *this;
  }

  const BcConstr & operator>=(Double rhs)
  {
    return operator>=(rhs.val());
  }

  const BcConstr & operator>=(int rhs)
  {
    return operator>=((double) rhs);
  }

  const BcConstr & operator==(double rhs)
  {
    this->sense('E');
    this->rhs(rhs);
    return *this;
  }

  const BcConstr & operator==(Double rhs)
  {
     return operator==(rhs.val());
  }

  const BcConstr & operator==(int rhs)
  {
    return operator==((double) rhs);
  }

  /// type = 'I' for constraints defining a implicit constraint used for preprocessing but not place in the formulation
  /// type = 'E' for constraints defining a explicit constraint inserted in the formulation
  /// type = 'C' for core (required for the IP formulation,
  /// type = 'F' for facultative (only helpfull for the LP formulation),
  /// type = 'S' for constraints defining a subsystem in column generation for extended formulation approach
  /// type = 'M' for constraints defining a pure master constraint
  /// type = 'X' for constraints defining a subproblem convexity constraint in the master
  void type(const char & type);

  void dualVal(double dualVal);

  operator InstanciatedConstr *() const;
  bool operator<(const BcConstr & that) const;
};

std::ostream& operator<<(std::ostream& os, const BcConstr & that);

class BcAddConstrFunctor
{
public:
    BcAddConstrFunctor()
    {
    }
    virtual ~BcAddConstrFunctor()
    {
    }

    virtual void operator()(const MultiIndex & indexArray);
};

class BcConstrArray
{
protected:
  GenericConstr * _genericConstrPtr;
  BcConstr _curInstConstrPtr; /// cache memory for operator()

  BcConstrArray();
public:
  explicit BcConstrArray(const BcFormulation & formulation, const std::string & name = "");

  virtual ~BcConstrArray();
  const std::string & genericName() const;
  const BcFormulation formulation() const;
  const BcModel& model() const;
  BcModel& model();
  
  /// checks whether the array has an element defined at index.
  bool isDefinedAt(const MultiIndex & indexArray);  

  virtual BcConstr & getElement(const MultiIndex & indexArray);
  virtual BcConstr & createElement(const MultiIndex & indexArray);
  inline  BcConstr & createElement(int firstIndex = 0, int secondIndex = -1,
                                   int thirdIndex = -1, int fourthIndex = -1,
                                   int fifthIndex = -1, int sixthIndex = -1,
                                   int seventhIndex = -1)
  {
      return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                      fifthIndex, sixthIndex, seventhIndex));
  }

   inline BcConstr & operator()(int firstIndex)
  {
    return createElement(MultiIndex(firstIndex));
  }
  
  inline BcConstr & operator()(int firstIndex, int secondIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex));
  }
  
  inline BcConstr & operator()(int firstIndex, int secondIndex, int thirdIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex));
  }

  inline BcConstr & operator()(int firstIndex, int secondIndex,
                               int thirdIndex, int fourthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex));
  }
  
  inline BcConstr & operator()(int firstIndex, int secondIndex,
                               int thirdIndex, int fourthIndex,
                               int fifthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex, fifthIndex));
  }
  
  inline BcConstr & operator()(int firstIndex, int secondIndex,
                               int thirdIndex, int fourthIndex,
                               int fifthIndex, int sixthIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex));
  }
  
  inline BcConstr & operator()(int firstIndex, int secondIndex,
                            int thirdIndex, int fourthIndex,
                            int fifthIndex, int sixthIndex,
                            int seventhIndex)
  {
    return createElement(MultiIndex(firstIndex, secondIndex, thirdIndex, fourthIndex,
                                    fifthIndex, sixthIndex, seventhIndex));
  }

  BcConstr & operator()(const MultiIndex & indexArray)
  {
    return createElement(indexArray);
  }

  BcConstrIndex operator[](const int & index); // call getConstr in the end

  /// to set index reference letter
  const BcConstrArray & defineIndexNames(const MultiIndexNames & multiIndexNames);
  const BcConstrArray & defineIndexName(const MultiIndexNames & multiIndexNames) ///backward compatibility
  {
      return defineIndexNames(multiIndexNames);
  }

  void sense(char sense);
  void rhs(double rhs);

  const BcConstrArray & operator<=(double defaultRhs)
  {
    this->sense('L');
    this->rhs(defaultRhs);
    return *this;
  }

  const BcConstrArray & operator<=(int defaultRhs)
  {
    this->sense('L');
    this->rhs(defaultRhs);
    return *this;
  }

  const BcConstrArray & operator>=(double defaultRhs)
  {
    this->sense('G');
    this->rhs(defaultRhs);
    return *this;
  }

  const BcConstrArray & operator>=(int defaultRhs)
  {
    this->sense('G');
    this->rhs(defaultRhs);
    return *this;
  }

  const BcConstrArray & operator==(double defaultRhs)
  {
    this->sense('E');
    this->rhs(defaultRhs);
    return *this;
  }

  const BcConstrArray & operator==(int defaultRhs)
  {
    this->sense('E');
    this->rhs(defaultRhs);
    return *this;
  }

  /// by default a constraint is propagated in proprocessing; the user can set this off
  void toBeUsedInPreprocessing(bool flag);
  void considerAsEqualityInPreprocessing(bool flag);

  void type(char type);
  void flag(char flag);

  void dualVal(double defaultDualVal);
  std::ostream & print(std::ostream& os = std::cout) const;

  const BcConstrArray & attach(BcAddConstrFunctor * addConstrRoutinePtr);
};

std::ostream& operator<<(std::ostream& os, const BcConstrArray & that);

struct BcConstrIndex: public std::pair<BcConstrArray, MultiIndex> // BcConstrIndex
{
  BcConstrIndex(const BcConstrArray & modConstr, const MultiIndex & indexArray);
  //BcConstrIndex(const BcConstrArray & modConstr);
  BcConstrIndex & operator[](const int & index);
  const MultiIndex & id();

  bool isDefined();
  const BcConstr & operator=(const BcRowExpression & exp);
  const BcConstr & operator+(const BcRowExpression & exp);
  const BcConstr & operator+=(const BcRowExpression & exp);
  const BcConstr & operator-(const BcRowExpression & exp);
  const BcConstr & operator-=(const BcRowExpression & exp);
  const BcConstr & operator+(const BcVarCoef & coef);
  const BcConstr & operator-(const BcVarCoef & coef);
  const BcConstr & operator+=(const BcVarCoef & coef);
  const BcConstr & operator-=(const BcVarCoef & coef);
  const BcConstr & operator+(const BcVar & var);
  const BcConstr & operator-(const BcVar & var);
  const BcConstr & operator+=(const BcVar & var);
  const BcConstr & operator-=(const BcVar & var);
  const BcConstr & operator+(BcVarIndex & modVarInd);
  const BcConstr & operator-(BcVarIndex & modVarInd);
  const BcConstr & operator+=(BcVarIndex & modVarInd);
  const BcConstr & operator-=(BcVarIndex & modVarInd);

  void sense(const char & sense);
  void rhs(double rhs);

  const BcConstrIndex & operator<=(double rhs)
  {
    this->sense('L');
    this->rhs(rhs);
    return *this;
  }

  const BcConstrIndex & operator<=(int rhs)
  {
    this->sense('L');
    this->rhs(rhs);
    return *this;
  }

  const BcConstrIndex & operator>=(double rhs)
  {
    this->sense('G');
    this->rhs(rhs);
    return *this;
  }

  const BcConstrIndex & operator>=(int rhs)
  {
    this->sense('G');
    this->rhs(rhs);
    return *this;
  }

  const BcConstrIndex & operator==(double rhs)
  {
    this->sense('E');
    this->rhs(rhs);
    return *this;
  }

  const BcConstrIndex & operator==(int rhs)
  {
    this->sense('E');
    this->rhs(rhs);
    return *this;
  }

  /// type = 'I' for constraints defining a implicit constraint used for preprocessing but not place in the formulation
  /// type = 'E' for constraints defining a explicit constraint inserted in the formulation
  /// type = 'C' for core (required for the IP formulation,
  /// type = 'F' for facultative (only helpful for the LP formulation),
  /// type = 'S' for constraints defining a subsystem in column generation for extended formulation approach
  /// type = 'X' for constraints defining a subproblem convexity constraint in the master
  void type(const char & type);

  void dualVal(double dualVal);

  bool inCurProb();

  explicit operator BcConstr();
};

#endif /* BCMODELINGCONSTRC_H_ */
