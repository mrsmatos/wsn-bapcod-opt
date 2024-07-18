/**
 *
 * This file bcPackSetResConsBranching.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcPackSetResConsBranching.cpp
//  Project
//
//  Created by Ruslan Sadykov on 24/09/2018.
//
//

#include "bcNetworkBasedBranching.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"

#ifdef BCP_RCSP_IS_FOUND

#include "rcsp_interface.hpp"

#define PackSetResConsBranchingPrintLevel 3

//////////////////////////////////////////////////////////////////////
/********************************************************************/
/****** Methods of class PackSetResConsBranchConstrGenerator ********/
/********************************************************************/
//////////////////////////////////////////////////////////////////////

PackSetResConsBranchConstrGenerator
::PackSetResConsBranchConstrGenerator(PackSetResConsGenBranchConstr * gbcPtr,
                                      const bcp_rcsp::AccumResConsBranchGenerator * generatorPtr, long int generatorId,
                                      const Double & candidateLhs, const char & priorityDir) :
        BranchingConstrGenerator(gbcPtr, priorityDir, candidateLhs), _rcspGeneratorPtr(generatorPtr),
        _generatorId(generatorId), _ESRCgenBrGenPtr(gbcPtr)
{
  if (_rcspGeneratorPtr == nullptr)
      return;
  std::stringstream ss;
  ss << "ES_" << generatorPtr->packSetId << "_rID_" << generatorPtr->resId << "_" << generatorPtr->threshold;
  _description = ss.str();
}



PackSetResConsBranchConstrGenerator
::PackSetResConsBranchConstrGenerator(const PackSetResConsBranchConstrGenerator & that):
        BranchingConstrGenerator(that), _rcspGeneratorPtr(that._rcspGeneratorPtr), _generatorId(that._generatorId),
        _ESRCgenBrGenPtr(that._ESRCgenBrGenPtr)
{
    if (that._rcspGeneratorPtr != nullptr)
        _rcspGeneratorPtr = new bcp_rcsp::AccumResConsBranchGenerator(*that._rcspGeneratorPtr);
}

PackSetResConsBranchConstrGenerator::~PackSetResConsBranchConstrGenerator()
{
    delete _rcspGeneratorPtr;
}


void PackSetResConsBranchConstrGenerator
     ::instanciateBrConstrs(const int & parentNodeNb, const int & parentNodeTreatId, const int & childNb,
                            const bool & greaterOrEqual,
                            std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
  std::string name("ESRCBC");
  name = name + "_" + _rcspGeneratorPtr->packSetId + "_" + _rcspGeneratorPtr->resId;

  if (printL(5))
    std::cout << "PackSetResConsBranchConstrGenerator::instanciateBrConstr() " << name << std::endl;


  BranchingConstrBaseType *constrPtr = new PackSetResConsInstMastBranchConstr(MultiIndex(_generatorId, childNb),
                                                                              _ESRCgenBrGenPtr,
                                                                              _ESRCgenBrGenPtr->probConfPtr(),
                                                                              name + "g" + _generatorId + "c" + childNb,
                                                                              parentNodeTreatId, *_rcspGeneratorPtr,
                                                                              greaterOrEqual);
  if (printL(5))
      constrPtr->print(cout);

  nextBranchingConstrPtrList.push_back(constrPtr);

  return;
}

bool PackSetResConsBranchConstrGenerator
     ::nextNodeBrConstr(Node * parentNodePtr, std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
                        const ConstrPtrSet & existingMasterBranchingConstr)
{
  /// Node branching constraint defined last is treated first (LIFO)
  nextNodeBranchingConstrPtrList.clear();
  bool success(true);

  int ancestorNodeRef(-1);
  int treatOrderId(-1);
  if (parentNodePtr != NULL)
  {
    ancestorNodeRef = parentNodePtr->ref();
    treatOrderId = parentNodePtr->treatOrder();
  }

  switch (_direction)
  {
    case 'U':
    {
      if (_childNbCounter == 0)
        instanciateBrConstrs(ancestorNodeRef, treatOrderId, ++_childNbCounter, true, nextNodeBranchingConstrPtrList);
      else if (_childNbCounter == 1)
        instanciateBrConstrs(ancestorNodeRef, treatOrderId, ++_childNbCounter, false, nextNodeBranchingConstrPtrList);
      else
        success = false;
      break;
    }
    default:
    {
      if (_childNbCounter == 0)
        instanciateBrConstrs(ancestorNodeRef, treatOrderId, ++_childNbCounter, false, nextNodeBranchingConstrPtrList);
      else if (_childNbCounter == 1)
        instanciateBrConstrs(ancestorNodeRef, treatOrderId, ++_childNbCounter, true, nextNodeBranchingConstrPtrList);
      else
        success = false;
      break;
    }
  }
  return success;

}

void PackSetResConsBranchConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curSol)
{
  bcp_rcsp::FractionalMasterSolution rcspFracSolution;
  rcspFracSolution.solPts.reserve(curSol.size());
  rcspFracSolution.values.reserve(curSol.size());
  for (SolutionVarInfoPtrList::const_iterator infoIt = curSol.begin(); infoIt != curSol.end(); infoIt++)
  {
    if (!(*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
      continue;

    MastColumn * colPtr = static_cast<MastColumn *>((*infoIt)->varPtr);

    auto pathPtr = colPtr->spSol()->rcspSolPtr();
    if (pathPtr != nullptr)
    {
      rcspFracSolution.solPts.push_back(pathPtr);
      rcspFracSolution.values.push_back((*infoIt)->value);
    }
  }

  if (_rcspGeneratorPtr == nullptr)
  {
      _candidateLhs = 0.0;

  }
  else {
      _candidateLhs = _ESRCgenBrGenPtr->violation(rcspFracSolution, *_rcspGeneratorPtr);

      if (printL(5))
          std::cout << "Accum. res. cons. branching generator (packSetId = " << _rcspGeneratorPtr->packSetId
                    << ",  resId = " << _rcspGeneratorPtr->resId << ", threshold = " << _rcspGeneratorPtr->threshold
                    << "); violation = " << _candidateLhs << std::endl;
  }

  return;

}

std::ostream & PackSetResConsBranchConstrGenerator::print(std::ostream & os) const
{
  BranchingConstrGenerator::print(os);
  os << "PackSetResConsBranchConstrGenerator" << std::endl;
  if (_rcspGeneratorPtr != nullptr)
  {
      os << "   packSetId = " << _rcspGeneratorPtr->packSetId << std::endl;
      os << "   resId = " << _rcspGeneratorPtr->resId << std::endl;
      os << "   accum. res. cons. threshold = " << _rcspGeneratorPtr->threshold << std::endl;
  }
  os << "   candidateLhs = " <<  _candidateLhs << std::endl;

  return(os);
}

void PackSetResConsBranchConstrGenerator::nicePrint(std::ostream& os) const
{
  if (_rcspGeneratorPtr != nullptr)
    os << "PackSetId " << _rcspGeneratorPtr->packSetId << " ResId " << _rcspGeneratorPtr->resId
       << " Thr. " << _rcspGeneratorPtr->threshold << " (lhs=" << _candidateLhs << ")";
}


//////////////////////////////////////////////////////////////////////
/********************************************************************/
/********** Methods of class PackSetResConsGenBranchConstr **********/
/********************************************************************/
//////////////////////////////////////////////////////////////////////


PackSetResConsGenBranchConstr::PackSetResConsGenBranchConstr(Model * modelPtr,
                                                             ProbConfig * probConfPtr,
                                                             const std::string & name,
                                                             const SelectionStrategy & priorityRule,
                                                             const Double & priorityLevel):
        GenericBranchingConstr(modelPtr, probConfPtr, name, priorityRule, priorityLevel, priorityLevel, false),
        Base4NonLinearGenericConstr(NULL), _numGenerators(0), _interfacePtr(NULL)
{
}

PackSetResConsGenBranchConstr::~PackSetResConsGenBranchConstr()
{
}

bool PackSetResConsGenBranchConstr::prepareSeparation()
{
    std::vector<const bcp_rcsp::GraphData *> graphs;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() != NULL)
            graphs.push_back((*cgSpConfPtrIt)->rcspGraphPtr());
    }
    int resId = 0; /// TO DO : take resId as a parameter
    _interfacePtr = bcp_rcsp::createAndPrepareAccumResConsBranching(graphs, resId);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "BaPCod error : could not prepare. accum. res. cons. branching " << std::endl;
    }

    return true;
}

/// for the moment, the separtion is implemented only for the resource whose algId is equal to zero,
/// i.e. for the most critical resource
void PackSetResConsGenBranchConstr
     ::branchingSeparationFindCandidates(const MasterColSolution & curListOfMasterCol,
                                         const int & maxNumOfCandidates,
                                         BranchGeneratorsSet & generatedBrConstrGeneratorSet)
{
    if (_interfacePtr == NULL)
        return;

    bcp_rcsp::FractionalMasterSolution rcspFracSolution;
    rcspFracSolution.solPts.reserve(curListOfMasterCol.size());
    rcspFracSolution.values.reserve(curListOfMasterCol.size());
    for (MasterColSolution::const_iterator colIt = curListOfMasterCol.begin(); colIt != curListOfMasterCol.end(); ++colIt)
    {
        auto rcspSolPtr = colIt->first->spSol()->rcspSolPtr();
        if (rcspSolPtr != nullptr)
        {
            rcspFracSolution.solPts.push_back(rcspSolPtr);
            rcspFracSolution.values.push_back(colIt->second._value);
        }
    }

    std::vector<const bcp_rcsp::AccumResConsBranchGenerator *> generatorPts;
    if (!_interfacePtr->separate(rcspFracSolution, maxNumOfCandidates, generatorPts))
        return;

    /// we create the candidates
    for (std::vector<const bcp_rcsp::AccumResConsBranchGenerator *>::iterator genIt = generatorPts.begin();
         genIt != generatorPts.end(); ++genIt)
    {
        Double generatorLhs(_interfacePtr->violation(rcspFracSolution, **genIt));
        generatedBrConstrGeneratorSet.insert(new PackSetResConsBranchConstrGenerator(this, *genIt, _numGenerators++,
                                                                                     generatorLhs));
    }
}

const LpCoef PackSetResConsGenBranchConstr::genericMastColumnCoef(InstanciatedConstr * icPtr, MastColumn * colPtr) const
{
  if (!icPtr->isTypeOf(VcId::PackSetResConsInstMastBranchConstrMask))
    return LpCoef(false, 0.0);

  return getMastColumnCoeff(static_cast<PackSetResConsInstMastBranchConstr *>(icPtr), colPtr);
}

void PackSetResConsGenBranchConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

const LpCoef PackSetResConsGenBranchConstr
             ::getMastColumnCoeff(PackSetResConsInstMastBranchConstr * constrPtr, MastColumn * colPtr) const
{
    LpCoef lpCoeff(false, 0.0);

    if (_interfacePtr == NULL)
        return lpCoeff;

    /// if the column has been generated after the constraint, then the coefficient should be zero,
    /// as the constraint changes the resource consumption time window and this new window is imposed in the pricing
    /// (we need to verify it, as a path can be feasible for both branches, i.e. the below verification does not
    ///  exactly correspond to the pricing constraint, notably for the backward part of the graph)
    if (!colPtr->spSol()->enumeratedFlag() && (colPtr->treatOrderId() > constrPtr->_treatOrderId))
        return lpCoeff;

    if (!_interfacePtr->satisfies(colPtr->spSol()->rcspSolPtr(), constrPtr->_rcspConstrPtr))
        lpCoeff.second = 1;

    lpCoeff.second.Cfloor();
    if (!lpCoeff.second.isZero())
        lpCoeff.first = true;

    return lpCoeff;
}

double PackSetResConsGenBranchConstr::violation(const bcp_rcsp::FractionalMasterSolution & rcspFracSolution,
                                                const bcp_rcsp::AccumResConsBranchGenerator & generator)
{
    return (_interfacePtr != NULL) ? _interfacePtr->violation(rcspFracSolution, generator) : 0.0;
}

std::ostream & PackSetResConsGenBranchConstr::print(std::ostream & os) const
{
  return os << "PackSetResConsGenBranchConstr" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/************* Methods of class PackSetResConsInstMastBranchConstr ***********/
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////


PackSetResConsInstMastBranchConstr
::PackSetResConsInstMastBranchConstr(const IndexCell & id,
                                     PackSetResConsGenBranchConstr * PackSetResConstrGenBranchConstrPtr,
                                     ProbConfig * probConfigPtr,
                                     const std::string & name,
                                     int treatOrderId,
                                     bcp_rcsp::AccumResConsBranchGenerator generator,
                                     bool greaterOrEqual) :
        InstMasterBranchingConstr(id, PackSetResConstrGenBranchConstrPtr, probConfigPtr, name, 0, 'L', ' ', 'E'),
        Base4NonLinearConstraint(), _treatOrderId(treatOrderId), _rcspConstrPtr(nullptr),
        _PackSetResConstrGenBranchConstrPtr(PackSetResConstrGenBranchConstrPtr)
{
    _rcspConstrPtr = new bcp_rcsp::AccumResConsBranchConstraint(id.first(), generator, greaterOrEqual);
}

void PackSetResConsInstMastBranchConstr::setMembership()
{
  if(!buildMembershipHasBeenPerformed())
  {
    genVarConstrPtr()->buildMembership(this);
    buildMembershipHasBeenPerformed(true);
  }

  bool cumulativeCoef(false);

  VarIndexManager::const_iterator it;
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Active, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Active, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _PackSetResConstrGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _PackSetResConstrGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  /// branching constraint is created before activation of columns, so we need to calculate the coefficients
  /// for unsuitable columns too
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _PackSetResConstrGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }

  Constraint::setMembership();

  return;
}

std::ostream & PackSetResConsInstMastBranchConstr::print(std::ostream& os) const
{
  os << "PackSetResConsInstMastBranchConstr" << std::endl;
  os << "   packSetId = " <<  _rcspConstrPtr->generator.packSetId << std::endl;
  os << "       resId = " <<  _rcspConstrPtr->generator.resId << std::endl;
  os << "  constraint = " << (_rcspConstrPtr->greaterOrEqual ? ">=" : "<") <<  _rcspConstrPtr->generator.threshold
     << std::endl;

  InstMasterBranchingConstr::print(os);

  return(os);
}

void PackSetResConsInstMastBranchConstr::shortPrint(std::ostream& os) const
{
  os << "PackSet " << _rcspConstrPtr->generator.packSetId << " rc" << (_rcspConstrPtr->greaterOrEqual ? ">=" : "<")
     <<  _rcspConstrPtr->generator.threshold;
}

bool PackSetResConsInstMastBranchConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::PackSetResConsInstMastBranchConstrMask, vcIdentifier);
}

std::vector<std::string> PackSetResConsInstMastBranchConstr::forDotPrint() const
{
  std::stringstream sstream;
  shortPrint(sstream);
  std::string s = sstream.str();

  s = s.substr(0, s.size()); //to remove the tailing one space.

  vector<string> ret;
  ret.push_back(s);
  return ret;
}

#endif /* BCP_RCSP_IS_FOUND */
