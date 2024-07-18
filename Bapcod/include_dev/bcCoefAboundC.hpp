/**
 *
 * This file bcCoefAboundC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef coefAboundC_h
#define coefAboundC_h
#include "bcUsefulHeadFil.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"

class Variable;
class InstanciatedVar;
class ColGenSpConf;
class Solution;
class CompBoundSetGenBranchConstr;
class BranchingConstrGenerator;

/**
 * @brief Describe ProbCoef  here
 *
 */
struct ProbCoef
{
  int rowRef;
  int colRef;
  Double coef;

  /**
   *
   *
   *
   * @return
   */
  ProbCoef();

  /**
   *
   *
   * @param pc
   *
   * @return
   */
  ProbCoef(const ProbCoef & pc);

  /**
   *
   *
   * @param rowNb
   * @param colNb
   * @param coefVal
   *
   * @return
   */
  ProbCoef(const int & rowNb, const int & colNb, const Double & coefVal);

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator<(const ProbCoef & b) const;

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator==(const ProbCoef & b) const;

  /**
   *
   *
   * @param os
   *
   * @return
   */
  std::ostream & print(std::ostream& os = std::cout) const;

};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ProbCoef & that)
{
  return that.print(os);
}


/**
 * @brief Describe ProbCoefColSmallerThan here
 *
 */
struct ProbCoefColSmallerThan
{
  /**
   *
   *
   *
   * @return
   */
  bool operator()(const ProbCoef  & a, const ProbCoef & b) const
  {
    return (a < b);
  }
};

/**
 * @brief Describe ProbCoefRowSmallerThan here
 *
 */
struct ProbCoefRowSmallerThan
{

  /**
   *
   *
   *
   * @return
   */
  bool operator()(const ProbCoef  & a, const ProbCoef & b) const
  {
    if (a.rowRef < b.rowRef) return true;
    if (a.rowRef > b.rowRef) return false;
    if (a.colRef < b.colRef) return true;
    return false;
  }
};

/**
 * @brief Describe ProbSetCoef
 *
 */
struct ProbSetCoef
{
  int setRef;
  char setType;
  int colRef;
  Double coef;
  Double setPriorityIndex;

  /**
   *
   *
   *
   * @return
   */
  ProbSetCoef();

  /**
   *
   *
   * @param setNb
   * @param type
   * @param colNb
   * @param coefVal
   * @param pri
   *
   * @return
   */
  ProbSetCoef(const int & setNb, const char & type, const int & colNb, const Double & coefVal, const Double & pri);

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator<(const ProbSetCoef & b) const;

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator==(const ProbSetCoef & b) const;

  /**
   *
   *
   * @param os
   *
   * @return
   */
  std::ostream & print(std::ostream& os = std::cout) const;
};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ProbSetCoef & that)
{
  return that.print(os);
}


/**
 * @brief Describe ProbBound here
 *
 */
struct ProbBound
{
  int ref;
  char sense;
  Double bound;
public:
  ProbBound();
  /**
   *
   *
   * @param Lref
   * @param Lsense
   * @param Lbound
   *
   * @return
   */
  ProbBound(const int & Lref, const char & Lsense, const Double & Lbound);

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator<(const ProbBound & b) const;

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator==(const ProbBound & b) const;

  /**
   *
   *
   * @param os
   *
   * @return
   */
  std::ostream & print(std::ostream& os = std::cout) const;
};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ProbBound & that)
{
  return that.print(os);
}


/**
 * @brief Describe ProbType here
 *
 */
struct ProbType
{
  int ref;
  char type;
public:
  ProbType();

  /**
   *
   *
   * @param Lref
   * @param Ltype
   *
   * @return
   */
  ProbType(const int & Lref, const char & Ltype);
  bool operator<(const ProbType & b) const;
  bool operator==(const ProbType & b) const;
  std::ostream & print(std::ostream& os = std::cout) const;
};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ProbType & that)
{
  return that.print(os);
}


/**
 * @brief Describe ProbIntC here
 *
 */
struct ProbIntC
{
  int ref;
  char type;
  Double val;
public:
  ProbIntC();

  /**
   *
   *
   * @param Lref
   * @param Ltype
   * @param Lval
   *
   * @return
   */
  ProbIntC(const int & Lref, const char & Ltype, const Double & Lval);

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator<(const ProbIntC & b) const;

  /**
   *
   *
   * @param b
   *
   * @return
   */
  bool operator==(const ProbIntC & b) const;

  /**
   *
   *
   * @param os
   *
   * @return
   */
  std::ostream & print(std::ostream& os = std::cout) const;
};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ProbIntC & that)
{
  return that.print(os);
}

/**
 * Describe ComponentBound here
 *
 */
class ComponentBound
{
 public:
  ///variable concerned by the component bound
  InstanciatedVar * _varPtr;
  ///value of the bound
  Double _val;
  /// 'G' for imposing a  ">= " bound, i.e.,  a lower bound; 'L' otherwise
  char _sign;
  /// number of col satisfying the component bound
  Double _cardinality;
  /// number of col not satisfying the component bound
  Double _complementCard;
  /// holds current branching priority level
  //Double _priorityFactor;

 public:
 ComponentBound(): _varPtr(NULL), _sign(' ') {}

  /**
   *
   *
   * @param varPtr
   * @param val
   * @param sign
   * @param card
   * @param compCard
   */
  ComponentBound(InstanciatedVar * varPtr, const Double & val = 1.0, const char & sign = 'G', const Double & card = 0, const Double & compCard = 0);

  /**
   *
   *
   * @param varPtr
   * @param val
   * @param sign
   * @param card
   * @param compCard
   */
  ComponentBound(Variable * varPtr, const Double & val = 1.0, const char & sign = 'G', const Double & card = 0, const Double & compCard = 0);

  virtual ~ComponentBound(){}

  /**
   *
   *
   */
  void complement();

  /**
   *
   *
   *
   * @return
   */
  InstanciatedVar * varPtr()  const {return  _varPtr;}

  /**
   *
   *
   *
   * @return
   */
  const Double & val() const {return _val;}

  /**
   *
   *
   *
   * @return
   */
  Double roundedVal() const {return (_sign == 'G'? _val:_val+1);}

  /**
   *
   *
   * @param varPtr
   * @param fracVal
   *
   * @return
   */
  bool satisfiedBy(Variable * varPtr, const Double & fracVal) const;

  /**
   *
   *
   *
   * @return
   */
  char sign() const {return _sign;}

  /**
   *
   *
   * @param val
   *
   * @return
   */
  bool satisfiedBy(const Double & val) const;

  /**
   *
   *
   * @param c
   */
  void cardinality(const Double & c) { _cardinality = c;}

  /**
   *
   *
   *
   * @return
   */
  const Double & cardinality() const {return _cardinality;}

  /**
   *
   *
   * @param c
   */
  void complementCard(const Double & c) { _complementCard = c;}

  /**
   *
   *
   *
   * @return
   */
  const Double & complementCard() const {return _complementCard;}

  /**
   *
   *
   * @param os
   *
   * @return
   */
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual bool operator==(const ComponentBound & that) const;
  virtual bool operator!=(const ComponentBound & that) const {return !this->operator==(that);}
};

/**
 *
 *
 * @param os
 * @param that
 *
 * @return
 */
inline std::ostream& operator<<(std::ostream& os, const ComponentBound & that)
{
  return that.print(os);
}




class ComponentSequence: public std::vector <ComponentBound>
{
  ComponentSequence() = delete;

public:
  ColGenSpConf * _cgSpConfPtr;
  bool _allBranchesGenerated;
  Double _fracWeight;
  // _fracWeightRU = Rounded up of fract weight
  Double _fracWeightRU;
  /// _fracWeightRD; = Rounded down of fract weight
  Double _fracWeightRD;
  /// Stands for the rounding sense of the last component bound: 'L' or 'G'
  /**
   * Default is 'G' for empty class sequence
   * Default is 'G' for all non-empty class sequence
   */
  char _activeSense;
  /// It is the sense of the branching constraint inequality only when comp. seq. is empty
  ComponentSequence *_directPred;

  int _additionalNbOfCpBd;
  Double _priorityFactor;
  //Double _colClassFractionality;
  //Double _totalCardinalityOfFractCol;

 public:
   ComponentSequence(ColGenSpConf * cgSpConfPtr);
   ComponentSequence(const ComponentSequence & cs);
   void reset();
  // ComponentSequence(const ComponentSequence::iterator & begin, const ComponentSequence::iterator & end);

  /**
   *
   *
   *
   * @return
   */
  virtual ~ComponentSequence(){}

  /**
   *
   *
   *
   * @return
   */
  const Double & classCardinality() const;

  /**
   *
   *
   *
   * @return
   */
  const Double & fracWeight() const;
  const Double & fracWeightRU() const {return _fracWeightRU;}
  const Double & fracWeightRD() const {return _fracWeightRD;}
  const int & additionalNbOfCpBd() const {return _additionalNbOfCpBd;}
 
  /**
   *
   *
   *
   * @return
   */
  const ComponentBound & lastCP() const
  {
    //  return *((((const ComponentBound)(*(this->rbegin()))).ptr));
    return *this->rbegin();
  }

  bool satisfiedBy(Solution * solPtr) const;

  /**
   *
   *
   *
   * @return
   */
  InstanciatedVar * lastCompLbVarPtr() const;

  void allCompLbVarPts(std::vector<InstanciatedVar *> & varPts) const;

  /**
   *
   *
   *
   * @return
   */
  ComponentSequence * directPred() const;

  /**
   *
   *
   * @param fw
   */
  void fracWeight(const Double & fw);

  /**
   *
   *
   */
  void roundFracWeight();

  void complement();

  void additionalNbOfCpBd(const int & nb)
  { _additionalNbOfCpBd = nb;}

  void priorityFact(const Double & priority)
  { _priorityFactor = priority;}

  char activeSense() const;
  ColGenSpConf * cgSpConfPtr() const {return _cgSpConfPtr;}

  const bool & allBranchesGenerated() const {return _allBranchesGenerated;}

  void allBranchesGenerated(const bool & flag) {_allBranchesGenerated = flag;}

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  enum InclusionStatus {different = -1, identical = 0, superclass = 1, subclass = 2};
  friend InclusionStatus compareCbS(const ComponentSequence & listA, const ComponentSequence & listB);
};

inline std::ostream& operator<<(std::ostream& os, const ComponentSequence & that)
{
    return that.print(os);
}

ComponentSequence::InclusionStatus compareCbS(const ComponentSequence & listA, const ComponentSequence & listB);


#endif	// coefAboundC_h
