/**
 *
 * This file bcInstanciatedVarConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef InstanciatedVarConstrClasses_h
#define InstanciatedVarConstrClasses_h

#include "bcUsefulHeadFil.hpp"
#include "bcBoundLevC.hpp"
#include "bcCoefAboundC.hpp"
#include "bcDoubleC.hpp"
#include "bcPrintC.hpp"
#include "bcSolutionC.hpp"
#include "bcVarConstrC.hpp"

class GenericVarConstr;
class GenericVar;
class GenericConstr;
class GenericBranchingConstr;
class GenVarGenBranchConstr;
class NonLinearGenericConstr;
class BranchingConstrBaseType;
class InstanciatedConstr;
class MastColumn;
class Node;
class ProbConfig;
class OvfConf;
class MasterConf;
class Model;
class CompSetInstMastBranchConstr;


/********************************************************************/

/**
 * @brief Common base for Entities that are instanciated from a Generic class (Variables and Constraints)
 *
 */
class InstanciatedVarConstr 
{
 protected:
  IndexCell _id; /// data pointer id's
  GenericVarConstr * _genericVarConstrPtr;
  ProbConfig * _probConfigPtr;
 public:
  InstanciatedVarConstr(const IndexCell & id, GenericVarConstr* genericVarConstrPtr, ProbConfig* probConfigPtr);
  InstanciatedVarConstr(InstanciatedVarConstr * imcPtr, GenericVarConstr * genericVarConstrPtr);
  InstanciatedVarConstr(const InstanciatedVarConstr & that);
  virtual ~InstanciatedVarConstr(){}
  virtual const IndexCell & id() const;
  virtual const MultiIndex & multiIndex() const {return id();}
  virtual const Double & costrhs() const = 0;
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr) = 0;
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr) = 0;
  virtual std::ostream & print(std::ostream& os = std::cout) const = 0;
  virtual bool consecutive2varWhenBrOnCBS(InstanciatedVar * varPtr) {return false;}

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};

struct InstVarConstrSortByIndexCell
{
  bool operator()(InstanciatedVarConstr * a, InstanciatedVarConstr * b) const
  {
    return (a->id() < b->id());
  }
};


/********************************************************************/


/**
 * @brief This class model Variables that are Instanciated from a Generic Variable 
 *
 * 
 */
class InstanciatedVar: public Variable, public InstanciatedVarConstr 
{
  GenericVar * _genVarPtr;
  std::map< int, double> _arcIdToCoeff;
 public:
  /**
   *
   * @param id Identifier of the Variable: Made of pointers to the date structure associated to the variable indices  
   * @param genericVarConstrPtr Pointer to the generic mathematical variable from which this model derives  
   * @param probConfigPtr Pointer to the Problem configuratin to which this variable belongs  
   * @param name name of the variable model  
   * @param costrhs  Cost of the variable model  
   * @param sense  'P' = positive ,  'N' = negative , 'F' = free ;  
   * @param type 'C' = continuous, 'B' = binary, 'I' = integer;  
   * @param kind 'E' = explicitly in the formulation, 'I' = implicitly in the formulation, 'R' = explicit but delayed/retarded inclusion in the formulation   
   * @param upperBound local upper bound  
   * @param lowerBound local lower bound  
   * @param flag 's' = static belonging to the problem and erased when the problem is erased;  'd' = dynamic VarConstr not belonging to the problem  
   * @param directive 'U' = priority do rounding the variable up when branching; 'D' = down  
   * @param priority level of the variable from 1 to 500 the higher the most priority is given to this varaible when branching  
   * @param val The intital value of the variable  
   * @param globalUb global upper bound can differ from local upper bound for subproblem var. only  
   * @param globalLb global lower bound can differ from local lower bound for subproblem var. only
   * @param presetMembership == true if the membership of this variable in any constraint is ENTIRELY modeled by functions "buildmembership" in proper constraints
   */
  InstanciatedVar(const IndexCell & id,                    
                  GenericVar * genericVarConstrPtr,  
                  ProbConfig * probConfigPtr,              
                  const std::string & name = "",                
                  const Double & costrhs = 0,              
                  const char & sense = 'P',                
                  const char & type = 'I',                 
                  const char & kind = 'E',                 
                  const Double & upperBound = BapcodInfinity,
                  const Double & lowerBound = 0,           
                  const char & flag = 's',                 
                  const char & directive = 'U',            
                  const Double & priority = 0.0,           
                  const Double & val = 0,                  
                  const Double & globalUb = BapcodInfinity,
                  const Double & globalLb = 0,            
		          const bool & presetMembership = true);
  InstanciatedVar(const InstanciatedVar & ivar);
  virtual ~InstanciatedVar();
  virtual ProbConfig * probConfPtr() const {return _probConfigPtr;}
  virtual GenericVar * genVarPtr() const {return _genVarPtr;}

  virtual const Double & costrhs() const;

  virtual void costrhs(const Double & newCostRhs);
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void setMembership();
  virtual GenericVarConstr * genVarConstrPtr() const;
  virtual void enumerativeSetMembership(); // set membership by enumerating problem constraints
  virtual bool consecutive2varWhenBrOnCBS(InstanciatedVar * varPtr);
  virtual const std::string & genericName();
  virtual const Double fracPartRelativeToTarget() const;

  void createStabInfo(const BcObjStatus::MinMaxIntFloat & minmax);

  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
  void setArcMembership(const int arcId, const double value);
  inline virtual const std::map<int, double> & arcIdToCoeff() const { return _arcIdToCoeff; }
};
inline std::ostream& operator<<(std::ostream& os, const InstanciatedVar & that)
{return that.print(os);}

/********************************************************************/

/**
 * @brief This class models Constraints that are Instanciated from a Generic Constraint 
 * 
 */

class InstanciatedConstr: public Constraint, public InstanciatedVarConstr /// @brief Constraints that are Instanciated from a Generic Constraint  
{
  GenericConstr * _genConstrPtr;

 public:
  /**
   * @param id Identifier of the Variable made of pointers to the date structure associated to the variable indices
   * @param genericVarConstrPtr  Pointer to the generic mathematical variable from which this instanciated constraint derives
   * @param probConfigPtr  Pointer to the Problem configuratin to which this variable belongs
   * @param name  name of the constraint
   * @param rhs   Righ-Hand-Size of the constraint
   * @param sense  'G' for >= , 'L' for >= , 'E' for == constraint
   * @param type   'C' for core (required for the IP formulation, 'F' for facultative (only helpfull for the LP formulation)
   * @param kind  'E' = explicit, 'I' = implicit, 'R' = delayed constraints added only after solving the LP root node without them
   * @param flag  's' = static belonging to the problem and erased when the problem is erased;  'd' = dynamic VarConstr not belonging to the problem
   * @param val  The intital value of the dual variable associated to this constraint
   * @param upperBound  local upper bound on the  dual variable associated to this constraint
   * @param lowerBound  local lower bound on the  dual variable associated to this constraint
   * @param directive ('U' or 'D') is used for disjunctive branching on constraints only, it says whether the first branch to be explore is the rounding up of the rhs 'U' or its rounding down 'D'
   * @param  priority  is used for disjunctive branching on constraints only, the higher value the more chance this branching on constraints is chose in priority for branching
   * @paramy toBeUsedInPreprocessing == true if constraint can be used for preprocesing purposes by bapcod, == false otherwise
   * @paramy presetMembership == true if the membership of this variable in any constraint is ENTIRELY modeled by functions "buildmembership" in proper variables, == false otherwise
   */
  //TODO remove  presetMembership from parameter list
  InstanciatedConstr(const IndexCell& id, 
                     GenericConstr* genericVarConstrPtr, 
                     ProbConfig* probConfigPtr, 
                     const std::string & name = "", 
                     const Double & rhs = 0, 
                     const char & sense= 'L', 
                     const char & type = 'C', 
                     const char & kind = 'E', 
                     const char & flag = 's', 
                     const Double & val = 0, 
                     const Double & upperBound = BapcodInfinity,
                     const Double & lowerBound = - BapcodInfinity,
		             const char & directive = 'U',
		             const Double & priority = 1.0,
		             const bool & presetMembership = true,
		             const bool & toBeUsedInPreprocessing = true,
		             const bool & considerAsEqualityInPreprocessing = false);

  InstanciatedConstr(const InstanciatedConstr & that);

  InstanciatedConstr(InstanciatedConstr * imcPtr,  
		             GenericConstr * genBrConstrPtr,
		             const std::string & name,
		             const Double & rhs,
		             const char & sense,
		             const char & type,
		             const char & kind,
		             const char & flag);

  virtual ~InstanciatedConstr();

  virtual GenericConstr * genConstrPtr() const {return _genConstrPtr;}

  virtual const Double & costrhs() const;
  virtual void costrhs(const Double & newCostRhs)
  {
    Constraint::costrhs(newCostRhs);
     return;
   }
  virtual ProbConfig * probConfPtr() const {return _probConfigPtr;}
  virtual GenericVarConstr * genVarConstrPtr() const;
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual void setMembership();
  virtual void enumerativeSetMembership(); // set membership by enumerating problem variables
  virtual const std::string & genericName();
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void nicePrint(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};
inline std::ostream& operator<<(std::ostream& os, const InstanciatedConstr & that)
{return that.print(os);}


/********************************************************************/

class NonLinearInstConstr: public InstanciatedConstr, public Base4NonLinearConstraint
{
 public:
  NonLinearInstConstr(const IndexCell& id,  
                     GenericConstr* genConstrPtr, 
                     ProbConfig* probConfigPtr,
                     const std::string& name  = "", 
                     const Double& costrhs = 0, 
                     const char& sense= 'L', 
                     const char& type = ' ', 
                     const char& kind = 'E',
                     const char& flag = 's', 
                     const Double& val = 0,
                     const Double& upperBound = BapcodInfinity,
                     const Double& lowerBound = - BapcodInfinity,
		             const char & directive = 'U', // for branching on constraints
		             const Double & priority = 1.0);  // for branching on constraints
  virtual ~NonLinearInstConstr(){}
  virtual bool computeCount(ConstVarConstrConstPtr vcPtr);
  virtual const LpCoef computeCoef(ConstVarConstrConstPtr vcPtr);
  virtual std::ostream & print(std::ostream& os = std::cout) const;

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};
inline std::ostream& operator<<(std::ostream& os, const NonLinearInstConstr & that)
{return that.print(os);}



/********************************************************************/


class BranchingConstrBaseType
{
 protected:
  int _depthWhenGenerated;
  std::set<ProbConfig *> _probConfigPtrSet;
 public:
  BranchingConstrBaseType(ProbConfig * probConfPtr) {_probConfigPtrSet.insert(probConfPtr);}
  BranchingConstrBaseType(const std::set<ProbConfig *> & probConfigPtrSet):
    _probConfigPtrSet(probConfigPtrSet) {}
  virtual ~BranchingConstrBaseType(){}
  virtual const std::set<ProbConfig *> & probConfigPtrSet() const {return _probConfigPtrSet;}
  virtual std::ostream & print(std::ostream& os = std::cout) const;
  virtual void shortPrint(std::ostream& os = std::cout) const {}
  virtual std::vector<std::string> forDotPrint() const {std::vector<std::string> a;return a;} // All elements of this string must have the same size.
  void append2name(const std::string & ad);
  void depthWhenGenerated(const int depth)
  {
    _depthWhenGenerated = depth;
  }

  /**
   * Test if the current class is type of vcIdentifier.
   * @param vcIdentifier
   * @return true if current class is type of vcIdentifier, false otherwise.
   */
  virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
};
inline std::ostream& operator<<(std::ostream& os, const BranchingConstrBaseType & that)
{return that.print(os);}

#endif // InstanciatedVarConstrClasses_h
