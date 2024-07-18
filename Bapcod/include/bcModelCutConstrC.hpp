/**
 *
 * This file bcModelCutConstrC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCMODELINGCUTCONSTRC_H_
#define BCMODELINGCUTCONSTRC_H_

#include "bcModelConstrC.hpp"
#include "bcModelFormulationC.hpp"
#include "bcModelParameterC.hpp"

class GenericConstr;
class InstanciatedConstr;
struct BcConstrIndex;
class BcFormulation;
class BcSolution;
class BcConstr;
class BcNetwork;
class GenericCutConstr;
class LimMemRankOneCut;
class ExtendedArcCut;

class BcCutSeparationFunctor
{
private:
  char type;
public:
  BcCutSeparationFunctor():type('F'){}
  BcCutSeparationFunctor(const char & type):type(type){}
  virtual ~BcCutSeparationFunctor(){}
  
  char getType() { return type; }

  #define IN_ // Indicates what may be used for input
  #define OUT_ // Indicates what may be used for output
  virtual int operator() (IN_ BcFormulation formPtr,
                          IN_ BcSolution & primalSol,
			              OUT_ double & maxViolation,
			              OUT_ std::list< BcConstr > & cutList);
};

class CliqueCut;
class LimMemRankOneCut;
class ResConsKnapsackCut;
class ExtendedArcCut;
class CustomNonLinearCut;

class BcCutConstrArray: public BcConstrArray
{
protected:
  GenericCutConstr * _genericCutConstrPtr;

  BcCutConstrArray();
public:
  /**
   *
   * @param formulation
   * @param name
   * @param type 'C' for core (required for the IP formulation), 'F' for facultative (only helpfull for the LP formulation),
   * @param priorityRule
   * @param priorityLevel
   * @param toBeUsedInPreprocessing
   */
  BcCutConstrArray(const BcFormulation & formulation,
                   const std::string & name,
                   const char & type = 'F',
                   const double & rootPriorityLevel = 1.0,
                   const double & nonRootPriorityLevel = 1.0,
                   const bool & toBeUsedInPreprocessing = true);

  ~BcCutConstrArray() override;

  const std::list< BcConstr > elementList() const;

  void setRootPriorityLevel(const double & rootPriorityLevel_);

  void attach(BcCutSeparationFunctor * separationRoutinePtr);

};

struct BcCustomNonLinearCutInfo
{
  BcCustomNonLinearCutInfo()
  {
  }
  virtual ~BcCustomNonLinearCutInfo()
  {
  }
};

class CustomNonLinearCut;

class BcCustomNonLinearCut: public BcConstr
{
friend class BcCustomNonLinearCutArrayFunctor;
protected:
  CustomNonLinearCut * _custNonLinCutPtr;
public:
  BcCustomNonLinearCut(CustomNonLinearCut * cutPtr);
  const BcCustomNonLinearCutInfo * cutInfoPtr() const;
};

class BcCustomNonLinearCutArrayFunctor: public BcCutConstrArray
{
  int _numberOfCuts;
public:
  BcCustomNonLinearCutArrayFunctor(const BcFormulation & formulation,
                                   const std::string & name,
                                   const char & type = 'F',
                                   const SelectionStrategy & priorityRule = SelectionStrategy::MostViolated,
                                   const double & priorityLevel = 1.0);
  BcCustomNonLinearCut createNewCut(BcCustomNonLinearCutInfo * cutInfoPtr);
  virtual ~BcCustomNonLinearCutArrayFunctor();
  virtual int cutSeparationRoutine(BcFormulation formPtr,
				     	           BcSolution & projectedSol,
                                   std::list<std::pair<double, BcSolution> > & columnsInSol,
					               const double & violationTolerance,
					               std::list<BcCustomNonLinearCut> & cutList);
  virtual double getCoefficient(const BcCustomNonLinearCut & cut,
                                const BcSolution & spSol);
};

class SoftConflictsCut;

class BcSoftConflictsCut: public BcConstr
{
    friend class BcSoftConflictsCutArrayFunctor;
protected:
    SoftConflictsCut * _cutPtr;
public:
    BcSoftConflictsCut(SoftConflictsCut * cutPtr);
    double getDualVal() const;
    int type() const;
    const std::vector<std::pair<BcVar, BcVar> > conflicts() const;
};

class BcSoftConflictsCutArrayFunctor: public BcCutConstrArray
{
    int _numberOfCuts;
public:
    BcSoftConflictsCutArrayFunctor(const BcFormulation & formulation,
                                   const std::string & name,
                                   const char & type = 'C',
                                   const SelectionStrategy & priorityRule = SelectionStrategy::MostViolated,
                                   const double & priorityLevel = 1.0);

    /// type = 0 : coefficient of the column in the cut is equal to the number of conflicts in the column
    ///            (in the subproblem, we need an indication variable per conflict)
    /// type = 1 : coefficient of the column in the cut is 1 if there are conflicts in the column,
    ///            0 otherwise (in the subproblem, we need an indicator variable per cut)
    BcConstr createNewCut(const double & rhs, std::vector<std::pair<BcVar, BcVar> > & conflicts, const int cutType = 0);

    virtual ~BcSoftConflictsCutArrayFunctor();
    virtual int cutSeparationRoutine(BcFormulation formPtr,
                                     std::list<std::pair<double, BcSolution> > & columnsInFixedSol,
                                     std::list<std::pair<double, BcSolution> > & columnsInSol,
                                     std::list<BcConstr> & cutList);
    virtual int cutSeparationBasedOnFixedSol(BcFormulation formPtr,
                                             std::list<std::pair<double, BcSolution> > & columnsInOldFixedSol,
                                             std::list<std::pair<double, BcSolution> > & columnsInNewFixedSol,
                                             std::list<BcConstr> & cutList);
};

#endif /* BCMODELINGCUTCONSTRC_H_ */
