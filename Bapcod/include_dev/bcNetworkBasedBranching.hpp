/**
 *
 * This file bcNetworkBasedBranching.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

//
//  bcNetworkBasedBranching.hpp
//  Project
//
//  Created by Ruslan Sadykov on 6/10/2017.
//
//

#ifndef Project_bcNetworkBasedBranching_hpp
#define Project_bcNetworkBasedBranching_hpp

#include "bcGenBranchingConstrC.hpp"
#include "bcSpVarConstrC.hpp"

class GenPackSetAssignBranchingConstr: public GenericBranchingConstr
{
  int _numGeneratedBrConstrs;
  int _numElemSets;
  /// variables corresponding to arcs incident to vertices belonging to an elem. set
  std::map<const ProbConfig *, std::vector<std::set<InstanciatedVar *> > > _incidBranchingVars;
  /// variables corresponding to arcs belonging to an elem. set
  std::map<const ProbConfig *, std::vector<std::set<InstanciatedVar *> > > _belongBranchingVars;

  void augmentLhs(const ProbConfig * probConfPtr, std::map<InstanciatedVar *, double> & varValueMap,
                  std::vector<double> & lhsVector);

public:
  GenPackSetAssignBranchingConstr(Model * modelPtr,
			                      ProbConfig * probConfPtr,
			                      const std::string & name,
                                  const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                                  const Double & nonRootPriorityLevel = 1.0,
                                  const Double & rootPriorityLevel = 1.0,
			                      const bool & toBeUsedInPreprocessing = true);

  virtual ~GenPackSetAssignBranchingConstr();
  
  virtual bool prepareSeparation();
  
  virtual void branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                 const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);
};

class GenPathsPerNetworkBranchingConstr: public GenericBranchingConstr
{
  int _numGeneratedBrConstrs;
  std::map<const ProbConfig *, std::set<InstanciatedVar *> > _branchingVars;

public:
  GenPathsPerNetworkBranchingConstr(Model * modelPtr,
			                        ProbConfig * probConfPtr,
			                        const std::string & name,
                                    const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
                                    const Double & nonRootPriorityLevel = 1.0,
                                    const Double & rootPriorityLevel = 1.0,
			                        const bool & toBeUsedInPreprocessing = true);

  virtual ~GenPathsPerNetworkBranchingConstr();
  
  virtual bool prepareSeparation();
  
  virtual void branchingSeparationFindCandidates(const MasterVarSolution & curListOfMastAndSubprobVar,
                                                 const MasterColSolution & curListOfMasterCol,
                                                 const int & maxNumOfCandidates,
					                             BranchGeneratorsSet & generatedBrConstrGeneratorSet);
};

#ifdef BCP_RCSP_IS_FOUND

/**
 *   Elementary set resource consumption branching
 */

namespace bcp_rcsp
{
    struct AccumResConsBranchGenerator;
    class AccumResConsBranchingInterface;
    struct AccumResConsBranchConstraint;
    struct FractionalMasterSolution;
}

class PackSetResConsGenBranchConstr;

class PackSetResConsBranchConstrGenerator: public BranchingConstrGenerator
{
	/// Specific instanciated first var on which we branch
    const bcp_rcsp::AccumResConsBranchGenerator * _rcspGeneratorPtr;
    long int _generatorId;
	PackSetResConsGenBranchConstr * _ESRCgenBrGenPtr;

	void instanciateBrConstrs(const int & parentNodeNb, const int & parentNodeTreatId,
							  const int & childNb, const bool & greateOrEqual,
							  std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);


public:
	PackSetResConsBranchConstrGenerator(PackSetResConsGenBranchConstr * gbcPtr,
                                        const bcp_rcsp::AccumResConsBranchGenerator * generator,
										long int generatorId, const Double & candidateLhs,
										const char & priorityDir = 'U');
	PackSetResConsBranchConstrGenerator(const PackSetResConsBranchConstrGenerator & that);

	virtual ~PackSetResConsBranchConstrGenerator();

	virtual BranchingConstrGenerator * clone() const
	{
		return new PackSetResConsBranchConstrGenerator(*this);
	}

	virtual bool nextNodeBrConstr(Node * parentNodePtr,
								  std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
								  const ConstrPtrSet & existingMasterBranchingConstr);

	virtual void computeLhs(const SolutionVarInfoPtrList & curSol);

	virtual std::ostream & print(std::ostream& os = std::cout) const;
	virtual void nicePrint(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const PackSetResConsBranchConstrGenerator & that)
{
	return that.print(os);
}

class PackSetResConsInstMastBranchConstr;

class PackSetResConsGenBranchConstr: public GenericBranchingConstr, public Base4NonLinearGenericConstr
{
    long int _numGenerators;
    bcp_rcsp::AccumResConsBranchingInterface * _interfacePtr;

public:
	PackSetResConsGenBranchConstr(Model * modelPtr, ProbConfig * probConfPtr, const std::string & name,
								  const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
								  const Double & priorityLevel = 1.0);

	virtual ~PackSetResConsGenBranchConstr();

	virtual bool prepareSeparation();

	virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
												   const int & maxNumOfCandidates,
												   BranchGeneratorsSet & generatedBrConstrGeneratorSet);

    virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;

    virtual void buildMembership(InstanciatedConstr * iconstrPtr);

	const LpCoef getMastColumnCoeff(PackSetResConsInstMastBranchConstr * constrPtr, MastColumn * colPtr) const;

    double violation(const bcp_rcsp::FractionalMasterSolution & rcspFracSolution,
                     const bcp_rcsp::AccumResConsBranchGenerator & generator);

    virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const PackSetResConsGenBranchConstr * that)
{
	return that->print(os);
}

class PackSetResConsInstMastBranchConstr:  public InstMasterBranchingConstr, public Base4NonLinearConstraint
{
	friend class BcPackSetResConsBranchConstr;
	friend class PackSetResConsGenBranchConstr;

	int _treatOrderId; /// treat order of the node where this constraint has been generated
    bcp_rcsp::AccumResConsBranchConstraint * _rcspConstrPtr;

	PackSetResConsGenBranchConstr * _PackSetResConstrGenBranchConstrPtr;
public:

	PackSetResConsInstMastBranchConstr(const IndexCell & id,
			                           PackSetResConsGenBranchConstr * PackSetResConstrGenBranchConstrPtr,
			                           ProbConfig * probConfigPtr,
									   const std::string & name,
									   int treatOrderId,
									   bcp_rcsp::AccumResConsBranchGenerator generator,
							           bool greaterOrEqual);
	virtual ~PackSetResConsInstMastBranchConstr() {}

	virtual void setMembership();
	virtual std::ostream & print(std::ostream& os = std::cout) const;
	virtual void shortPrint(std::ostream& os) const;

    const bcp_rcsp::AccumResConsBranchConstraint * getRcspConstrPtr() const {return _rcspConstrPtr;}

    virtual bool isTypeOf(const VcId::VcIdentifier& vcIdentifier) const;
	virtual std::vector<std::string> forDotPrint() const;
};

/**
 *   Packing sets resource consumption branching
 */

namespace bcp_rcsp
{
    struct RyanFosterBranchGenerator;
    class RyanFosterBranchingInterface;
    struct RyanFosterBranchConstraint;
}

class PackSetRyanFosterGenBranchConstr;

class PackSetRyanFosterBranchConstrGenerator: public BranchingConstrGenerator
{
	const bcp_rcsp::RyanFosterBranchGenerator * _rcspGeneratorPtr;
	long int _generatorId;
	PackSetRyanFosterGenBranchConstr * _PSRFgenBrGenPtr;

	void instanciateBrConstrs(const int & parentNodeNb, const int & parentNodeTreatId,
							  const int & childNb, const bool & together,
							  std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList);


public:
	PackSetRyanFosterBranchConstrGenerator(PackSetRyanFosterGenBranchConstr * gbcPtr,
                                           const bcp_rcsp::RyanFosterBranchGenerator * rcspGeneratorPtr,
										   long int generatorId, const Double & candidateLhs,
										   const char & priorityDir = 'U');

	PackSetRyanFosterBranchConstrGenerator(const PackSetRyanFosterBranchConstrGenerator & that);

	virtual ~PackSetRyanFosterBranchConstrGenerator();

	virtual BranchingConstrGenerator * clone() const
	{
		return new PackSetRyanFosterBranchConstrGenerator(*this);
	}

	virtual bool nextNodeBrConstr(Node * parentNodePtr,
								  std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
								  const ConstrPtrSet & existingMasterBranchingConstr);

	virtual void computeLhs(const SolutionVarInfoPtrList & curSol);

	virtual std::ostream & print(std::ostream& os = std::cout) const;
	virtual void nicePrint(std::ostream& os = std::cout) const;
};

class PackSetRyanFosterInstMastBranchConstr;

class PackSetRyanFosterGenBranchConstr: public GenericBranchingConstr, public Base4NonLinearGenericConstr
{
	long int _numGenerators;
    bool _usePackingSets;
    bcp_rcsp::RyanFosterBranchingInterface * _interfacePtr;

public:
	PackSetRyanFosterGenBranchConstr(Model * modelPtr, ProbConfig * probConfPtr, const std::string & name,
								     const SelectionStrategy & priorityRule = SelectionStrategy::MostFractional,
								     const Double & priorityLevel = 1.0, const bool & usePackingSets = true);

	virtual ~PackSetRyanFosterGenBranchConstr();

	virtual bool prepareSeparation();

	virtual void branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
												   const int & maxNumOfCandidates,
												   BranchGeneratorsSet & generatedBrConstrGeneratorSet);

	virtual const LpCoef genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const;

	virtual void buildMembership(InstanciatedConstr * iconstrPtr);

	const LpCoef getMastColumnCoeff(PackSetRyanFosterInstMastBranchConstr * constrPtr, MastColumn * colPtr) const;

    double violation(const bcp_rcsp::FractionalMasterSolution & rcspFracSolution,
                     const bcp_rcsp::RyanFosterBranchGenerator & generator);

    virtual std::ostream & print(std::ostream& os = std::cout) const;
};

inline std::ostream& operator<<(std::ostream& os, const PackSetRyanFosterGenBranchConstr * that)
{
	return that->print(os);
}

class PackSetRyanFosterInstMastBranchConstr:  public InstMasterBranchingConstr, public Base4NonLinearConstraint
{
	friend class BcPackSetRyanFosterBranchConstr;
	friend class PackSetRyanFosterGenBranchConstr;

	int _treatOrderId; /// treat order of the node where this constraint has been generated
	bcp_rcsp::RyanFosterBranchConstraint * _rcspConstrPtr;

	PackSetRyanFosterGenBranchConstr * _packSetRyanFosterGenBranchConstrPtr;
public:

	PackSetRyanFosterInstMastBranchConstr(const IndexCell & id,
										  PackSetRyanFosterGenBranchConstr * packSetRyanFosterGenBranchConstrPtr,
									      ProbConfig * probConfigPtr,
									      const std::string & name,
                                          int treatOrderId,
                                          bcp_rcsp::RyanFosterBranchGenerator generator,
 									      bool together);

	virtual ~PackSetRyanFosterInstMastBranchConstr();

    const bcp_rcsp::RyanFosterBranchConstraint * getRcspConstrPtr() const {return _rcspConstrPtr;}

	virtual void setMembership();
	virtual std::ostream & print(std::ostream & os = std::cout) const;
	virtual void shortPrint(std::ostream & os) const;

	virtual bool isTypeOf(const VcId::VcIdentifier & vcIdentifier) const;
	virtual std::vector<std::string> forDotPrint() const;
};

#endif /* BCP_RCSP_IS_FOUND */

#endif /* Project_bcNetworkBasedBranching_hpp */
