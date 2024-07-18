/**
 *
 * This file bcNonPublicHomExtCapCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcHomExtCapCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 17/05/2016.
//
//

#ifdef BCP_RCSP_IS_FOUND

#include "bcUsefulHeadFil.hpp"
#include "bcNonPublicCuts.hpp"
#include "bcProblemC.hpp"
#include "bcMastColumnC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"

#define HomExtCapCutsPrintLevel 1

#include "rcsp_interface.hpp"

#ifdef RHECC_SEP_IS_FOUND

#include "cutInterface.h"
#include "CutGenerator2.h"
#include "CutGenerator.h"

HomogExtendedCapacityCut::HomogExtendedCapacityCut(GenericExtendedArcCutConstr * genConstrPtr,
                                                   ProbConfig * probConfigPtr,
                                                   const std::string & name,
                                                   const Double & rhs,
                                                   const int & numerator,
                                                   const int & denominator,
                                                   const int & resourceId,
                                                   const std::set<int> verticesSet,
                                                   const int & sumDemands):
  ExtendedArcCut(genConstrPtr, probConfigPtr, name, rhs, 'G'), _useNewSep(false),
  _numerator(numerator), _denominator(denominator), _resourceId(resourceId),
  _sumDemands(sumDemands), _verticesSet(verticesSet)
{
  int max = 0;
  for (std::set<int>::iterator it = verticesSet.begin(); it != verticesSet.end(); it++)
    {
      if (max < *it) max = *it;
    }
  _verticesBitmap.resize(max+1, false);
  for (std::set<int>::iterator it = verticesSet.begin(); it != verticesSet.end(); it++)
    {
      _verticesBitmap[*it] = true;
    }
}

HomogExtendedCapacityCut::~HomogExtendedCapacityCut()
{
}

void HomogExtendedCapacityCut::nicePrint(std::ostream& os)
{
  os << "HECC "<< name() << ": rows = (";
  std::set<int>::iterator setIt = _verticesSet.begin();
  if (setIt != _verticesSet.end())
    {
      std::cout << *setIt;
      for (++setIt; setIt != _verticesSet.end(); ++setIt)
        os << ", " << *setIt;
    }
  os << "), multiplier = " << _numerator << "/" << _denominator;
}

double HomogExtendedCapacityCut::getArcCoefficient(const int & tailVertId, const int & headVertId,
                                                   const double * tailResCons) const
{
  /// test whether head and tail belong to S
  bool tailIn = (tailVertId >= int(_verticesBitmap.size()))? false: _verticesBitmap[tailVertId];
  bool headIn = (headVertId >= int(_verticesBitmap.size()))? false: _verticesBitmap[headVertId];
  if (tailIn == headIn)
    return 0.0;

     if(_useNewSep)
     {
	  double fCoeff, alpha;

	  // check the integrality of resource constraints (required for HECCs)
	  double intTol = 1e-5;
	  if ((int(tailResCons[_resourceId] - intTol) == int(tailResCons[_resourceId] + intTol)) &&
		    (int(tailResCons[_resourceId] + intTol) > 0))
	  {
	       std::cerr << "ERROR: HomogExtendedCapacityCut::getArcCoefficient called with tailResCons[" << _resourceId
		             << "]=" << tailResCons[_resourceId] << std::endl;
	       exit(1);
	  }
	  int d = int(tailResCons[_resourceId] + intTol);

	  if (tailIn && !headIn)
	  {
	       // Set coefficients for variables leaving S
	       fCoeff = arturnewsep::CutGenerator2::eccfrac((long long)d * (long long)_numerator, _denominator);
	       if( fCoeff == 0 )
	       {
		    return 0.0;
	       }
	       else
	       {
		    alpha = 0.0;// eccfrac((long long)demand * (long long)num, den);
		    if( fCoeff >= alpha - 0.0001)
			 return 1.0 - fCoeff;
		    else
			 return fCoeff/alpha - fCoeff;
	       }
	  }
      else /// (!tailIn && headIn)
	  {
	       // Set coefficients for variables entering S
	       fCoeff = arturnewsep::CutGenerator2::eccfrac(-(long long)d * (long long)_numerator, _denominator);
	       if( fCoeff == 0 )
		    return 0.0;
	       else
	       {
		    alpha = 0.0; //eccfrac((long long)demand * (long long)num, den);
		    if( fCoeff >= alpha - 0.0001)
			 return 1.0 - fCoeff;
		    else
			 return fCoeff/alpha - fCoeff;
	       }
	  }
     }
     else
     {
	  double fCoeff, alpha;

	  // check the integrality of resource constraints (required for HECCs)
	  double intTol = 1e-5;
	  if ((int(tailResCons[_resourceId] - intTol) == int(tailResCons[_resourceId] + intTol)) &&
		    (int(tailResCons[_resourceId] + intTol) > 0))
	  {
	       std::cerr << "ERROR: HomogExtendedCapacityCut::getArcCoefficient called with tailResCons[" << _resourceId
		             << "]=" << tailResCons[_resourceId] << std::endl;
	       exit(1);
	  }
	  int d = int(tailResCons[_resourceId] + intTol);

	  if (tailIn && !headIn)
	  {
	       // Set coefficients for variables leaving S
	       fCoeff = CutGenerator::eccfrac(d * _numerator, _denominator);
	       if( fCoeff == 0 )
	       {
		    return (double(d * _numerator) / double(_denominator));
	       }
	       else
	       {
		    alpha = CutGenerator::eccfrac(_sumDemands * _numerator, _denominator);
		    if( fCoeff >= alpha - 0.0001)
			 return (double)CutGenerator::eccceil(d * _numerator, _denominator);
		    else
			 return ((double)CutGenerator::eccceil(d * _numerator, _denominator)
				   -1.0 + fCoeff/alpha); // lifting of gomory cuts
	       }
	  }
      else /// (!tailIn && headIn)
	  {
	       // Set coefficients for variables entering S
	       fCoeff = CutGenerator::eccfrac(-d * _numerator, _denominator);
	       if( fCoeff == 0 )
	       {
		    return ((double)(-d * _numerator) / double(_denominator));
	       }
	       else
	       {
		    alpha = CutGenerator::eccfrac(_sumDemands * _numerator, _denominator);
		    if( fCoeff >= alpha - 0.0001)
			 return (double)CutGenerator::eccceil(-d * _numerator, _denominator);
		    else
			 return ((double)CutGenerator::eccceil(-d * _numerator, _denominator)
				   -1.0 + fCoeff/alpha); // lifting of gomory cuts

	       }
	  }
     }
}

double HomogExtendedCapacityCut::getRouteCoefficient(const std::vector<int> & routeVertIds,
	  const std::vector<std::vector<double> > & routeResCons) const
{
     //if new separation is used, route coefficient equals the sum of arcs coeffs, no lifting for now
     if(_useNewSep)
     {
	  const std::vector<int> & route = routeVertIds;
	  const std::vector<std::vector<double> > & resCons = routeResCons;

	  double coeff = 0.0;

	  std::vector<int>::const_iterator nextVertIt, vertIt = route.begin();
	  std::vector<std::vector<double> >::const_iterator resConsIt = resCons.begin();
	  if (vertIt != route.end())
	  {
	       nextVertIt = vertIt;
	       ++nextVertIt;
	       while (nextVertIt != route.end())
	       {
		    coeff += getArcCoefficient(*vertIt, *nextVertIt, &(*resConsIt)[0]);
		    ++vertIt;
		    ++nextVertIt;
		    ++resConsIt;
	       }
	  }

	  return coeff;
     }
     else
     {
	  double fCoeff, alpha;

	  // check the integrality of resource constraints (required for HECCs)
	  double intTol = 1e-5;
	  std::vector<int>::const_iterator vIt;
	  std::vector<std::vector<double> >::const_iterator rcIt;
	  for (rcIt = routeResCons.begin(); rcIt != routeResCons.end(); rcIt++)
	       if ((int((*rcIt)[_resourceId] - intTol) == int((*rcIt)[_resourceId] + intTol)) &&
			 (int((*rcIt)[_resourceId] + intTol) > 0))
	       {
		    std::cerr << "ERROR: HomogExtendedCapacityCut::getRouteCoefficient called with routeResCons["
			          << (rcIt - routeResCons.begin()) << "][" << _resourceId
			          << "]=" << (*rcIt)[_resourceId] << std::endl;
		    exit(1);
	       }

	  // calculate the total demand served by the route inside the set S that defines the HECC
	  double dd = 0.0;
	  double prev = 0.0;
	  rcIt = routeResCons.begin();
	  for (vIt = routeVertIds.begin(); vIt != routeVertIds.end(); vIt++)
	  {
	       // test whether the current vertex belongs to S
	       if ((*vIt < int(_verticesBitmap.size())) && _verticesBitmap[*vIt])
		    dd += (*rcIt)[_resourceId] - prev;

	       // next vertex
	       prev = (*rcIt)[_resourceId];
	       rcIt++;
	  }
	  int d = int(dd + intTol);

	  // the coefficient of the route in the lifted cut is simply the rounded scaled served demand
	  fCoeff = CutGenerator::eccfrac(d * _numerator, _denominator);
	  if( fCoeff == 0 )
	  {
	       return (double(d * _numerator) / double(_denominator));
	  }
	  else
	  {
	       alpha = CutGenerator::eccfrac(_sumDemands * _numerator, _denominator);
	       if( fCoeff >= alpha - 0.0001)
		    return (double)CutGenerator::eccceil(d * _numerator, _denominator);
	       else
		    return ((double)CutGenerator::eccceil(d * _numerator, _denominator)
			      -1.0 + fCoeff/alpha); // lifting of gomory cuts is still valid
	  }
	  return 0.0; // never reach this point
     }
}

GenericHomExtCapCutConstr::GenericHomExtCapCutConstr(Model * modelPtr,
                                                     ProbConfig * probConfPtr,
                                                     const std::string & name,
                                                     const Double & nonRootPriorityLevel,
                                                     const Double & rootPriorityLevel,
                                                     const int & maxNumCutsPerRound,
                                                     const std::vector<int> & demands,
                                                     const std::vector<int> & capacities,
                                                     const int & resourceId):
  GenericExtendedArcCutConstr(modelPtr, probConfPtr, name, nonRootPriorityLevel, rootPriorityLevel),
  _useNewSep(false), _maxNumCutsPerRound(maxNumCutsPerRound), _resourceId(resourceId),
  _demands(demands), _capacities(capacities), cutGenObj(NULL), cutGenObjNew(NULL)
{
}

bool GenericHomExtCapCutConstr::prepareSeparation()
{
  if (_useNewSep)
	{
		// initialize the instance information
		arturnewsep::InstanceInfo instNew;
		instNew.capacity = 0;
		std::vector<int> tmpDemands(_demands.size()+1, 0);
		instNew.demand = &tmpDemands[0];
		for (int i = 1; i < int(_demands.size()+1); i++)
		{
			instNew.demand[i] = _demands[i-1];
		}
		instNew.nrootbranches = 0; // not used
		instNew.numNodes = _demands.size()+1;

		// calculate the maximum capacity of a subproblem
		for (int i = 0; i < int(_capacities.size()); i++)
		{
			int c = _capacities[i];
			if (c > instNew.capacity) instNew.capacity = c;
		}

		// create the cut generator object
		cutGenObjNew = new arturnewsep::CutGenerator2(&instNew);
	}
	else
	{
		InstanceInfo inst;
		inst.capacity = 0;
		std::vector<int> tmpDemands(_demands.size()+1, 0);
		inst.demand = &tmpDemands[0];
		for (int i = 1; i < int(_demands.size()+1); i++)
		{
			inst.demand[i] = _demands[i-1];
		}
		inst.nrootbranches = 0; // not used
		inst.numNodes = _demands.size()+1;

		// calculate the maximum capacity of a subproblem
		for (int i = 0; i < int(_capacities.size()); i++)
		{
			int c = _capacities[i];
			if (c > inst.capacity) inst.capacity = c;
		}

		// create the cut generator object
		cutGenObj = new CutGenerator(&inst);
	}

    return GenericExtendedArcCutConstr::prepareSeparation();
}

GenericHomExtCapCutConstr::~GenericHomExtCapCutConstr()
{
	if(_useNewSep)
		delete cutGenObjNew;
	else
		delete cutGenObj;
}

void GenericHomExtCapCutConstr
     ::cutSeparationRoutine(const VarPtrSet & curSol,
                            std::multiset < InstanciatedConstr * , CutSeparationPriorityComp > & generatedCutConstrSet)
{
	double inteps = 1e-5;
	/// this is the cut separation routine

	// initialize data structures to store the fractional solution information
	std::map<std::pair<std::pair<int,int>,int>,double> capArcVals;
	std::map<std::pair<std::pair<int,int>,int>,double>::iterator capArcIt;
	std::map<std::pair<int,int>,double> arcVals;
	std::map<std::pair<int,int>,double>::iterator arcIt;

	/// we loop over all the columns participating in the fractional solution (with positive value)
	/// and for each column we find a route through elementarity sets
	int i, j, d;
	for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
	{
		if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
			continue;
		MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);

		if ((colPtr->spSol() == nullptr) || (colPtr->spSol()->rcspPathPtr() == nullptr))
		    continue;

		int graphId = colPtr->spSol()->rcspPathPtr()->graphId;
		//const std::vector<int> & orderedArcIds = colPtr->spSol()->rcspPathPtr()->arcIds;
		const std::vector<std::vector<double> > & resCons = colPtr->spSol()->rcspPathPtr()->resConsumption;

        /// for the moment, uses the transformation of the sequence of arc ids to the sequence of vertex ids
        /// TO DO : functions getRouteCoefficient and getArcCoefficient should be rewritten to deal directly
        ///         with the sequence of arc ids to increase the performance

        std::vector<int> vertexIds;
        _graphPts[graphId]->obtainVertexIds(colPtr->spSol()->rcspPathPtr(), vertexIds);

        std::vector<int>::const_iterator vertIt;
		std::vector<std::vector<double> >::const_iterator resConsIt;
		for (vertIt = vertexIds.begin(), resConsIt = resCons.begin(); vertIt != vertexIds.end(); ++vertIt, ++resConsIt)
		{
			/// current vertex in the route is *vertIt (0 is depot, client first index is 1)
			/// with cumulative demand (*resConsIt)[0] (including current vertex)
			/// value of column is colPtr->val()
          
            /// skip additional auxiliary nodes
            if (*vertIt > _demands.size() + 1)
              continue;

			// skip the first visit to the depot
			if (*vertIt == 0)
			{
				i = 0;
				d = 0;
				continue;
			}

			// accumulate the value of the current capacitated arc (i,j,d) and arc (i,j)
			j = *vertIt;
			if (j == (_demands.size()+1)) j = 0;
			capArcIt = capArcVals.find(std::pair<std::pair<int,int>,int>(std::pair<int,int>(i,j),d));
			if (capArcIt != capArcVals.end()) capArcIt->second += colPtr->val();
			else capArcVals[std::pair<std::pair<int,int>,int>(std::pair<int,int>(i,j),d)] = colPtr->val();
			arcIt = arcVals.find(std::pair<int,int>(i,j));
			if (arcIt != arcVals.end()) arcIt->second += colPtr->val();
			else arcVals[std::pair<int,int>(i,j)] = colPtr->val();

			// the current head is the next tail, and the demand on it is the next arc demand
			i = j;
			d = int((*resConsIt)[0] + inteps);
		}
		if (printL(HomExtCapCutsPrintLevel))
		{
			std::cout << "Route from sp." << colPtr->cgSpConfPtr()->id().first()
	                    				 << " with val = " << colPtr->val()
										 << " and cost = " << colPtr->costrhs() << " : ";
			colPtr->spSol()->printOrderedSolution();
		}
	}

	// build the fractional solution information (inverting the arcs to follow the convention of the CutGenerator)
	if (_useNewSep)
	{
		std::vector<arturnewsep::CUTS_ArcCapVariable> tmpCapArc(capArcVals.size());

		int k;
		for (k = 0, capArcIt = capArcVals.begin(); k < int(tmpCapArc.size()); k++, capArcIt++)
		{
			tmpCapArc[k].j = capArcIt->first.first.first;
			tmpCapArc[k].i = capArcIt->first.first.second;
			tmpCapArc[k].d = capArcIt->first.second;
			tmpCapArc[k].value = capArcIt->second;
		}
		std::vector<arturnewsep::CUTS_ArcVariable> tmpArc(arcVals.size());
		for (k = 0, arcIt = arcVals.begin(); k < int(tmpArc.size()); k++, arcIt++)
		{
			tmpArc[k].j = arcIt->first.first;
			tmpArc[k].i = arcIt->first.second;
			tmpArc[k].value = arcIt->second;
		}

		arturnewsep::LpSolution sol;
		sol.arcs = &tmpArc[0];
		sol.arcsCap = &tmpCapArc[0];
		sol.numArcs = tmpArc.size();
		sol.numArcsCap = tmpCapArc.size();

		/// we need to generate up to _maxNumCutsPerRound cuts
		/// with minimum violation param().BapCodCutViolationTolerance()
		/// using demands _demands (client first index is 0)
		/// and vehicle capacities _capacities
		/// for a cut with vertices set vertSet (std::set<int>), sum of vertices demands sumDem (int)
		/// multiplier (numerator (int) / denominator(int)), and right-hand side rhs (double)
		/// one should writeta = auxdata

		///  std::string name("HECC");
		///  MultiIndex newCutId(_numGeneratedCuts++);
		///  newCutId.appendRef2name(name, multiIndexNames());
		///  generatedCutConstrSet.insert(new HomogExtendedCapacityCut(newCutId, this, probConfPtr(), name, rhs, numerator, denominator, _resourceId, vertSet, sumDem));

		// Allocate a cut list
		arturnewsep::CutList cuts;
		std::vector<arturnewsep::Cut> cutVector;
		cutVector.resize( _maxNumCutsPerRound );
		cuts.cuts = &(cutVector[0]);
		cuts.numCuts = 0;

		// Call the separation routine that uses heuristic
		arturnewsep::SeparateECCbyHeuristic( cutGenObjNew, &sol, &cuts,
				_maxNumCutsPerRound, arturnewsep::PROB_PATH, 2);
		// Add the cuts to the problem
		for (i = 0; i < cuts.numCuts; i++)
		{
			arturnewsep::ExtCCData* eccData = (arturnewsep::ExtCCData*)cuts.cuts[i].data;
			std::set<int> vertSet;
			for (j = 0; j < int(eccData->S.size()); j++)
				if (eccData->S[j]) {vertSet.insert(j); };
			generatedCutConstrSet.insert(new HomogExtendedCapacityCut(this, probConfPtr(), "HECC",
					cuts.cuts[i].rhs, eccData->num, eccData->den, _resourceId,
					vertSet, eccData->weightS));
			// Release the cut data
			cutVector[i].DestroyData( cutVector[i].data );
		}
		if (printL(0))
		  std::cout << "ECC added by new heuristic. Total cuts: " << cuts.numCuts << std::endl;
	}
	else
	{
		// build the fractional solution information (inverting the arcs to follow the convention of the CutGenerator2)
		std::vector<ArcCapVariable> tmpCapArc(capArcVals.size());
		int k;
		for (k = 0, capArcIt = capArcVals.begin(); k < int(tmpCapArc.size()); k++, capArcIt++)
		{
			tmpCapArc[k].j = capArcIt->first.first.first;
			tmpCapArc[k].i = capArcIt->first.first.second;
			tmpCapArc[k].d = capArcIt->first.second;
			tmpCapArc[k].value = capArcIt->second;
		}

		std::vector<ArcVariable> tmpArc(arcVals.size());
		for (k = 0, arcIt = arcVals.begin(); k < int(tmpArc.size()); k++, arcIt++)
		{
			tmpArc[k].j = arcIt->first.first;
			tmpArc[k].i = arcIt->first.second;
			tmpArc[k].value = arcIt->second;
		}

		LpSolution sol;
		sol.arcs = &tmpArc[0];
		sol.arcsCap = &tmpCapArc[0];
		sol.numArcs = tmpArc.size();
		sol.numArcsCap = tmpCapArc.size();

		/// we need to generate up to _maxNumCutsPerRound cuts
		/// with minimum violation param().BapCodCutViolationTolerance()
		/// using demands _demands (client first index is 0)
		/// and vehicle capacities _capacities
		/// for a cut with vertices set vertSet (std::set<int>), sum of vertices demands sumDem (int)
		/// multiplier (numerator (int) / denominator(int)), and right-hand side rhs (double)
		/// one should writeta = auxdata

		///  std::string name("HECC");
		///  MultiIndex newCutId(_numGeneratedCuts++);
		///  newCutId.appendRef2name(name, multiIndexNames());
		///  generatedCutConstrSet.insert(new HomogExtendedCapacityCut(newCutId, this, probConfPtr(), name, rhs, numerator, denominator, _resourceId, vertSet, sumDem));

		// Allocate a cut list
		CutList cuts;
		std::vector<Cut> cutVector;
		cutVector.resize( _maxNumCutsPerRound );
		cuts.cuts = &(cutVector[0]);
		cuts.numCuts = 0;

		// Call the separation routine that uses heuristic
		int depth = 6;
		cutGenObj->setLpSolution( &sol );
		cutGenObj->setCutBatch( _maxNumCutsPerRound );
		cutGenObj->setProbType( PROB_PATH );
		cutGenObj->setMinSetSize( depth+1 );
		cutGenObj->extCapCutGenByHeur();
		cutGenObj->getExtCapCuts( &cuts );
		cutGenObj->getECCsInOldSets( &cuts );

		// Add the cuts to the problem
		for (i = 0; i < cuts.numCuts; i++)
		{
			ExtCCData* eccData = (ExtCCData*)cuts.cuts[i].data;
			std::set<int> vertSet;
			for (j = 0; j < int(eccData->S.size()); j++) {
				if (eccData->S[j])
					vertSet.insert(j);
			}
			generatedCutConstrSet.insert(new HomogExtendedCapacityCut(this, probConfPtr(), "HECC",
					cuts.cuts[i].rhs, eccData->num, eccData->den, _resourceId,
					vertSet, eccData->weightS));
			// Release the cut data
			cutVector[i].DestroyData( cutVector[i].data );
		}
		if (printL(0))
			std::cout << "ECC added by old heuristic. Total cuts: " << cuts.numCuts << std::endl;


		// if the number of violated cuts is smaller than the maximum and the enumeration depth is at least 6
		// (since 5-SRC make sets smaller than 6 useless)
		if ((depth >= 6) && (cuts.numCuts < _maxNumCutsPerRound))
		{
			// Call the separation routine that uses enumeration
			cutGenObj->setLpSolution( &sol );
			cutGenObj->setCutBatch( _maxNumCutsPerRound - cuts.numCuts );
			cutGenObj->setProbType( PROB_PATH );
			cutGenObj->setMaxDepth( depth );
			cutGenObj->setMaxCalls( 0 );
			cutGenObj->setMinSetSize( 6 );
			cutGenObj->extCapCutGenByEnum();
			cutGenObj->getExtCapCuts( &cuts );
			// Add the cuts to the problem
			for (i = 0; i < cuts.numCuts; i++)
			{
				ExtCCData* eccData = (ExtCCData*)cuts.cuts[i].data;
				std::set<int> vertSet;
				for (j = 0; j < int(eccData->S.size()); j++)
					if (eccData->S[j]) vertSet.insert(j);
				generatedCutConstrSet.insert(new HomogExtendedCapacityCut(this, probConfPtr(), "HECC",
						cuts.cuts[i].rhs, eccData->num, eccData->den, _resourceId,
						vertSet, eccData->weightS));

				// Release the cut data
				cutVector[i].DestroyData( cutVector[i].data );
			}
			if (printL(0))
  			  std::cout << "ECC added by enumeration. Total cuts: " << cuts.numCuts << std::endl;
		}
	}
  // Artur: FOR DEBUGGING ONLY
  // recalculate every cut violation
  /*std::multiset<InstanciatedConstr*,CutSeparationPriorityComp>::iterator it;
  for (it = generatedCutConstrSet.begin(); it != generatedCutConstrSet.end(); it++)
    {
      HomogExtendedCapacityCut* h = dynamic_cast<HomogExtendedCapacityCut*> (*it);
      if (h == 0) continue;
      double lhs = 0.0;
      for (VarPtrSet::const_iterator varPtrIt = curSol.begin(); varPtrIt != curSol.end(); ++varPtrIt)
        {
          if (!(*varPtrIt)->isTypeOf(VcId::MastColumnMask))
            continue;
          MastColumn * colPtr = static_cast<MastColumn *>(*varPtrIt);
          lhs += getMastColumnCoeff(h, colPtr).second * colPtr->val();
        }
      h->nicePrint(std::cout);
      std::cout << "  violation = " << h->rhs() - lhs << std::endl;
      i = 0;
  }

  // check if the cuts invalidate the following manually entered solution
  const int numRtSol = 6;
  const int testSol[][10] =
  { {7, 25, 35, 16, 0, 0, 0, 0, 0, 0},
    {18, 31, 19, 9, 21, 26, 0, 0, 0, 0},
    {14, 6, 36, 29, 24, 0, 0, 0, 0, 0},
    {33, 2, 28, 23, 22, 12, 11, 10, 4, 0},
    {13, 30, 15, 32, 27, 0, 0, 0, 0, 0},
    {20, 8, 5, 3, 1, 34, 17, 0, 0, 0} };
  std::multiset<InstanciatedConstr*,CutSeparationPriorityComp>::iterator it;
  for (it = generatedCutConstrSet.begin(); it != generatedCutConstrSet.end(); it++)
    {
      HomogExtendedCapacityCut* h = dynamic_cast<HomogExtendedCapacityCut*> (*it);
      if (h == 0) continue;
      double lhs = 0.0;
      for (int k = 0; k < numRtSol; k++)
        {
          Solution s;
          s.addToOrderedIds(0);
          for (int i = 0; testSol[k][i] != 0; i++)
            s.addToOrderedIds(testSol[k][i]);
          s.addToOrderedIds(_demands.size()+1);
          int d = 0;
          s.addToResConsumption(std::vector<double>(1, d));
          for (int i = 0; testSol[k][i] != 0; i++)
            {
              d += _demands[testSol[k][i]-1];
              s.addToResConsumption(std::vector<double>(1, d));
            }
          s.addToResConsumption(std::vector<double>(1, d));
          MastColumn c(h->masterConfPtr(), sp, &s);  // create a dummy column
          lhs += getMastColumnCoeff(h, &c).second;
        }
      if (lhs + 1e-7 < h->rhs())
        {
          std::cout << "ERROR: RHECC cuts optimal solution!" << std::endl;
        }
    } */
}

#endif /* RHECC_SEP_IS_FOUND */

#endif //BCP_RCSP_IS_FOUND
