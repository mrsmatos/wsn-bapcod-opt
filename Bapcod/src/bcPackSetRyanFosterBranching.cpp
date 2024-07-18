/**
 *
 * This file bcPackSetRyanFosterBranching.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

#ifdef BCP_RCSP_IS_FOUND

#include "bcNetworkBasedBranching.hpp"
#include "bcNetworkFlowC.hpp"
#include "bcColGenSpConfC.hpp"
#include "bcModelC.hpp"
#include "rcsp_interface.hpp"

//////////////////////////////////////////////////////////////////////
/********************************************************************/
/****** Methods of class PackSetRyanFosterBranchConstrGenerator *****/
/********************************************************************/
//////////////////////////////////////////////////////////////////////

PackSetRyanFosterBranchConstrGenerator
::PackSetRyanFosterBranchConstrGenerator(PackSetRyanFosterGenBranchConstr * gbcPtr,
                                         const bcp_rcsp::RyanFosterBranchGenerator * rcspGeneratorPtr,
                                         long int generatorId, const Double & candidateLhs,
                                         const char & priorityDir) :
        BranchingConstrGenerator(gbcPtr, priorityDir, candidateLhs), _rcspGeneratorPtr(rcspGeneratorPtr),
        _generatorId(generatorId), _PSRFgenBrGenPtr(gbcPtr)
{
  std::stringstream ss;
  if (_rcspGeneratorPtr != nullptr)
  ss << "PS_" << _rcspGeneratorPtr->firstSetId << "_" << _rcspGeneratorPtr->secondSetId;
  _description = ss.str();
}



PackSetRyanFosterBranchConstrGenerator
::PackSetRyanFosterBranchConstrGenerator(const PackSetRyanFosterBranchConstrGenerator & that):
        BranchingConstrGenerator(that), _rcspGeneratorPtr(nullptr), _generatorId(that._generatorId),
        _PSRFgenBrGenPtr(that._PSRFgenBrGenPtr)
{
    if (that._rcspGeneratorPtr != nullptr)
        _rcspGeneratorPtr = new bcp_rcsp::RyanFosterBranchGenerator(*that._rcspGeneratorPtr);
}

PackSetRyanFosterBranchConstrGenerator::~PackSetRyanFosterBranchConstrGenerator()
{
}

void PackSetRyanFosterBranchConstrGenerator
     ::instanciateBrConstrs(const int & parentNodeNb, const int & parentNodeTreatId, const int & childNb,
                            const bool & together, std::list<BranchingConstrBaseType *> & nextBranchingConstrPtrList)
{
  std::string name("PSRFBC");
  if (_rcspGeneratorPtr != nullptr)
    name = name + "_" + _rcspGeneratorPtr->firstSetId + "_" + _rcspGeneratorPtr->secondSetId;

  if (printL(5))
    std::cout << "PackSetRyanFosterBranchConstrGenerator::instanciateBrConstr() " << name << std::endl;


  BranchingConstrBaseType *constrPtr
          = new PackSetRyanFosterInstMastBranchConstr(MultiIndex(_generatorId, childNb),
                                                      _PSRFgenBrGenPtr, _PSRFgenBrGenPtr->probConfPtr(),
                                                      name + "g" + _generatorId + "c" + childNb,
                                                      parentNodeTreatId, *_rcspGeneratorPtr, together);
  if (printL(5))
      constrPtr->print(cout);

  nextBranchingConstrPtrList.push_back(constrPtr);

  return;
}

bool PackSetRyanFosterBranchConstrGenerator
     ::nextNodeBrConstr(Node * parentNodePtr, std::list<BranchingConstrBaseType *> & nextNodeBranchingConstrPtrList,
                        const ConstrPtrSet & existingMasterBranchingConstr)
{
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

void PackSetRyanFosterBranchConstrGenerator::computeLhs(const SolutionVarInfoPtrList & curSol)
{
  bcp_rcsp::FractionalMasterSolution rcspFracSolution;
  rcspFracSolution.solPts.reserve(curSol.size());
  rcspFracSolution.values.reserve(curSol.size());
  for (SolutionVarInfoPtrList::const_iterator infoIt = curSol.begin(); infoIt != curSol.end(); infoIt++)
  {
    if (!(*infoIt)->varPtr->isTypeOf(VcId::MastColumnMask))
      continue;

    MastColumn * colPtr = static_cast<MastColumn *>((*infoIt)->varPtr);

    auto rcspSolPtr = colPtr->spSol()->rcspSolPtr();
    if (rcspSolPtr != nullptr)
    {
        rcspFracSolution.solPts.push_back(rcspSolPtr);
        rcspFracSolution.values.push_back((*infoIt)->value);
    }
  }

  if (_rcspGeneratorPtr == nullptr)
  {
      _candidateLhs = 0.0;
  }
  else
  {
    _candidateLhs = _PSRFgenBrGenPtr->violation(rcspFracSolution, *_rcspGeneratorPtr);

    if (printL(5))
        std::cout << "PackSetRyanFosterBranchConstrGenerator (firstPackSetId = " << _rcspGeneratorPtr->firstSetId
                  << ", secondPackSetId = " << _rcspGeneratorPtr->secondSetId << "); _candidateLhs = "
                  << _candidateLhs << std::endl;

  }

  return;
}

std::ostream & PackSetRyanFosterBranchConstrGenerator::print(std::ostream & os) const
{
  BranchingConstrGenerator::print(os);
  os << "PackSetRyanFosterBranchConstrGenerator" << std::endl;
  if (_rcspGeneratorPtr != nullptr)
  {
      os << "   firstPackSetId = " << _rcspGeneratorPtr->firstSetId << std::endl;
      os << "   secondPackSetId = " << _rcspGeneratorPtr->secondSetId << std::endl;
  }
  os << "   candidateLhs = " <<  _candidateLhs << std::endl;

  return(os);
}

void PackSetRyanFosterBranchConstrGenerator::nicePrint(std::ostream& os) const
{
    if (_rcspGeneratorPtr != nullptr)
        os << "Ryan&Foster pack.set pair " << _rcspGeneratorPtr->firstSetId << " and "
           << _rcspGeneratorPtr->secondSetId << " (lhs=" << _candidateLhs << ")";
}


//////////////////////////////////////////////////////////////////////
/********************************************************************/
/******** Methods of class PackSetRyanFosterGenBranchConstr *********/
/********************************************************************/
//////////////////////////////////////////////////////////////////////


PackSetRyanFosterGenBranchConstr::PackSetRyanFosterGenBranchConstr(Model * modelPtr,
                                                                   ProbConfig * probConfPtr,
                                                                   const std::string & name,
                                                                   const SelectionStrategy & priorityRule,
                                                                   const Double & priorityLevel,
                                                                   const bool & usePackingSets):
        GenericBranchingConstr(modelPtr, probConfPtr, name, priorityRule, priorityLevel, priorityLevel, false),
        Base4NonLinearGenericConstr(NULL), _interfacePtr(NULL), _numGenerators(0),
        _usePackingSets(usePackingSets)
{
}

PackSetRyanFosterGenBranchConstr::~PackSetRyanFosterGenBranchConstr()
{
    delete _interfacePtr;
}

bool PackSetRyanFosterGenBranchConstr::prepareSeparation()
{
    std::vector<const bcp_rcsp::GraphData *> graphs;
    for (std::vector<ColGenSpConf *>::const_iterator cgSpConfPtrIt = probConfPtr()->colGenSubProbConfPts().begin();
         cgSpConfPtrIt != probConfPtr()->colGenSubProbConfPts().end(); ++cgSpConfPtrIt)
    {
        if ((*cgSpConfPtrIt)->rcspGraphPtr() != NULL)
            graphs.push_back((*cgSpConfPtrIt)->rcspGraphPtr());
    }
    _interfacePtr = bcp_rcsp::createAndPrepareRyanFosterBranching(graphs, _usePackingSets);
    if (_interfacePtr == nullptr)
    {
        std::cerr << "BaPCod error : could not prepare pack. set. based Ryan&Foster branching" << std::endl;
        return false;
    }
    return true;
}

void PackSetRyanFosterGenBranchConstr
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

    std::vector<const bcp_rcsp::RyanFosterBranchGenerator *> generatorPts;
    if (!_interfacePtr->separate(rcspFracSolution, maxNumOfCandidates, generatorPts))
        return;

    /// we create the candidates
    for (std::vector<const bcp_rcsp::RyanFosterBranchGenerator * >::iterator genIt = generatorPts.begin();
         genIt != generatorPts.end(); ++genIt)
    {
        Double generatorLhs(_interfacePtr->violation(rcspFracSolution, **genIt));
        generatedBrConstrGeneratorSet.insert(new PackSetRyanFosterBranchConstrGenerator(this, *genIt, _numGenerators++,
                                                                                        generatorLhs));
    }

}

const LpCoef PackSetRyanFosterGenBranchConstr::genericMastColumnCoef(InstanciatedConstr * icPtr,
                                                                     MastColumn * colPtr) const
{
  if (!icPtr->isTypeOf(VcId::PackSetRyanFostInstMastBranchConstrMask))
    return LpCoef(false, 0.0);

  return getMastColumnCoeff(static_cast<PackSetRyanFosterInstMastBranchConstr *>(icPtr), colPtr);
}

void PackSetRyanFosterGenBranchConstr::buildMembership(InstanciatedConstr * iconstrPtr)
{
  iconstrPtr->presetMembership(true);
  return;
}

const LpCoef PackSetRyanFosterGenBranchConstr
             ::getMastColumnCoeff(PackSetRyanFosterInstMastBranchConstr * constrPtr, MastColumn * colPtr) const
{
    LpCoef lpCoeff(false, 0.0);

    if (_interfacePtr == NULL)
        return lpCoeff;

    if (!_interfacePtr->satisfies(colPtr->spSol()->rcspSolPtr(), constrPtr->_rcspConstrPtr))
        lpCoeff.second = 1;

    lpCoeff.second.Cfloor();
    if (!lpCoeff.second.isZero())
        lpCoeff.first = true;

    return lpCoeff;
}

double PackSetRyanFosterGenBranchConstr::violation(const bcp_rcsp::FractionalMasterSolution & rcspFracSolution,
                                                   const bcp_rcsp::RyanFosterBranchGenerator & generator)
{
    return (_interfacePtr != NULL) ? _interfacePtr->violation(rcspFracSolution, generator) : 0.0;
}

std::ostream & PackSetRyanFosterGenBranchConstr::print(std::ostream & os) const
{
  return os << "PackSetRyanFosterGenBranchConstr" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/*********** Methods of class PackSetRyanFosterInstMastBranchConstr **********/
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////


PackSetRyanFosterInstMastBranchConstr
::PackSetRyanFosterInstMastBranchConstr(const IndexCell & id,
                                        PackSetRyanFosterGenBranchConstr * packSetRyanFosterGenBranchConstrPtr,
                                        ProbConfig * probConfigPtr,
                                        const std::string & name,
                                        int treatOrderId,
                                        bcp_rcsp::RyanFosterBranchGenerator generator,
                                        bool together) :
        InstMasterBranchingConstr(id, packSetRyanFosterGenBranchConstrPtr, probConfigPtr, name, 0, 'L', ' ', 'E'),
        Base4NonLinearConstraint(), _treatOrderId(treatOrderId), _rcspConstrPtr(nullptr),
        _packSetRyanFosterGenBranchConstrPtr(packSetRyanFosterGenBranchConstrPtr)
{
    _rcspConstrPtr = new bcp_rcsp::RyanFosterBranchConstraint(id.first(), generator, together);
}

PackSetRyanFosterInstMastBranchConstr::~PackSetRyanFosterInstMastBranchConstr()
{
    delete _rcspConstrPtr;
}

void PackSetRyanFosterInstMastBranchConstr::setMembership()
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
      LpCoef lpCoeff = _packSetRyanFosterGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Inactive, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Inactive, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _packSetRyanFosterGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }
  /// branching constraint is created before activation of columns, so we need to calculate the coefficients
  /// for unsuitable columns too
  for (it = problemPtr()->probVarSet().begin(VcIndexStatus::Unsuitable, 'd');
       it != problemPtr()->probVarSet().end(VcIndexStatus::Unsuitable, 'd'); ++it)
    if ((*it)->isTypeOf(VcId::MastColumnMask))
    {
      LpCoef lpCoeff = _packSetRyanFosterGenBranchConstrPtr->getMastColumnCoeff(this, static_cast<MastColumn *>(*it));
      if (lpCoeff.first)
        includeMember(*it, lpCoeff.second, cumulativeCoef);
    }

  Constraint::setMembership();

  return;
}

std::ostream & PackSetRyanFosterInstMastBranchConstr::print(std::ostream& os) const
{
  os << "PackSetRyanFosterInstMastBranchConstr" << std::endl;
  os << "  firstPackSetId = " << _rcspConstrPtr->generator.firstSetId << std::endl;
  os << "  secondPackSetId = " << _rcspConstrPtr->generator.secondSetId << std::endl;
  os << "  together = " << _rcspConstrPtr->together << std::endl;

  InstMasterBranchingConstr::print(os);

  return(os);
}

void PackSetRyanFosterInstMastBranchConstr::shortPrint(std::ostream& os) const
{
  os << "Pack.sets " << _rcspConstrPtr->generator.firstSetId << "," << _rcspConstrPtr->generator.secondSetId
     << ((_rcspConstrPtr->together) ? "" : " not") << " together";
}

bool PackSetRyanFosterInstMastBranchConstr::isTypeOf(const VcId::VcIdentifier& vcIdentifier) const
{
  return compareIdentifier(VcId::PackSetRyanFostInstMastBranchConstrMask, vcIdentifier);
}

std::vector<std::string> PackSetRyanFosterInstMastBranchConstr::forDotPrint() const
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
