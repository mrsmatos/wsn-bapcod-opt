/**
 *
 * This file bcOvfVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef OvfVarConstrClasses_h
#define OvfVarConstrClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcDoubleC.hpp"
#include "bcGenVarConstrC.hpp"
#include "bcInstanciatedVarConstrC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"

class Node;
class ProbConfig;
class SubProbVariable;
class MasterConf;
class ColGenSpConf;

class OvfVar: public Variable
{
 protected:
  ProbConfig * _originatingPconfPtr;
  Variable * _originatingVarPtr;
  int _cnt;
 public:
  OvfVar(ProbConfig * originatingPconfPtr,
	 Variable * originatingVarPtr,
	 int cnt = 0);

  OvfVar(ProbConfig * originatingPconfPtr, int cnt = 0);

  virtual ~OvfVar()
    {
    }

  ProbConfig * originatingPconfPtr() const
  {
    return _originatingPconfPtr;
  }

  Variable * originatingVarPtr() const
  {
    return _originatingVarPtr;
  }
  void  originatingVarPtr(Variable * ofvPtr)
  {
    _originatingVarPtr = ofvPtr;
  }

  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);
  virtual void setMembership();
  int cnt() const
  {
    return _cnt;
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class SpSetupOvfVar: public OvfVar
{
 public:
  SpSetupOvfVar(ProbConfig * originatingPconfPtr, int cnt);
  virtual ~SpSetupOvfVar() {}
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class OvfConstr: public Constraint
{
 protected:
  ProbConfig * _originatingPconfPtr;
  Constraint * _originatingConstrPtr;
  int _cnt;
 public:
  OvfConstr(ProbConfig * originatingPconfPtr,
	    Constraint * originatingConstrPtr = NULL,
	    int cnt = 0);

  virtual ~OvfConstr()
    {
    }

  ProbConfig * originatingPconfPtr() const
  {
    return _originatingPconfPtr;
  }

  Constraint * originatingConstrPtr() const
  {
    return _originatingConstrPtr;
  }

  void  originatingConstrPtr(Constraint * ofcPtr)
  {
    _originatingConstrPtr = ofcPtr;
  }

  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual bool membCount(ConstVarConstrConstPtr vcPtr);
  virtual const Double & membCoef(ConstVarConstrConstPtr vcPtr);
  virtual void setMembership();
  virtual const Double curRhs() const;
  int cnt() const
  {
    return _cnt;
  }

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

/**
 * Common base for Constraints that are not
 * derived from master or subproblem constraints
 */
class SpecificOvfConstr: public OvfConstr
{
 public:
 SpecificOvfConstr(ProbConfig * originatingPconfPtr,  int cnt):
    OvfConstr(originatingPconfPtr, NULL, cnt)
    {
    }

 /**
  * Test if the current class is type of vcIdentifier.
  * @param vcIdentifier
  * @return true if current class is type of vcIdentifier, false otherwise.
  */
 virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class SpVarLbOvfConstr: public SpecificOvfConstr
{
 protected:
  Variable * _originatingVarPtr;
 public:
  SpVarLbOvfConstr(ProbConfig * originatingPconfPtr, Variable * _originatingVarPtr, int cnt);

  virtual ~SpVarLbOvfConstr() {}

  Variable * originatingVarPtr() const
  {
    return _originatingVarPtr;
  }

  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class SpVarUbOvfConstr: public SpecificOvfConstr
{
 protected:
  Variable * _originatingVarPtr;
 public:
  SpVarUbOvfConstr(ProbConfig * originatingPconfPtr,
		   Variable * _originatingVarPtr,
		   int cnt);

  virtual ~SpVarUbOvfConstr()
    {
    }

  Variable * originatingVarPtr() const
  {
    return _originatingVarPtr;
  }

  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

class SpLbOvfConstr: public SpecificOvfConstr
{
 public:
  SpLbOvfConstr(ProbConfig * originatingPconfPtr);
  virtual ~SpLbOvfConstr() {}
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr  vcPtr);
  virtual void setMembership();
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

#endif // OvfVarConstrClasses_h
