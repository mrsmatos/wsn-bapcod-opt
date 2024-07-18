/**
 *
 * This file bcNonPublicOverlElimCuts.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
//  bcOverlElimCuts.cpp
//  Project
//
//  Created by Ruslan Sadykov on 12/07/2016.
//
//

#ifdef BCP_RCSP_IS_FOUND
#include "bcUsefulHeadFil.hpp"
#include "bcNonPublicCuts.hpp"
#include "bcProblemC.hpp"
#include "bcMastColumnC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"

#define OverlEliminCutsPrintLevel 1

#ifdef RHECC_SEP_IS_FOUND

#include "CutGenerator2.h"
#include "CutGenerator.h"


OverloadEliminationCut::OverloadEliminationCut(GenericExtendedArcCutConstr * genConstrPtr,
                                               ProbConfig * probConfigPtr,
                                               const std::string & name,
                                               const Double & rhs,
                                               const int & resourceId,
                                               const int & resConsumption,
                                               const std::set<int> verticesSet,
                                               const int & sumDemands,
											   const int & capacity,
											   const int & nbMachines) :
  ExtendedArcCut(genConstrPtr, probConfigPtr, name, rhs, 'G'),
  _resourceId(resourceId), _resConsumption(resConsumption), _sumDemands(sumDemands),
  _capacity(capacity), _nbMachines(nbMachines), _verticesSet(verticesSet)
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

OverloadEliminationCut::~OverloadEliminationCut()
{
}

void OverloadEliminationCut::nicePrint(std::ostream& os)
{
	os << "OEC "<< name() << ": rows = (";
	std::set<int>::iterator setIt = _verticesSet.begin();
	if (setIt != _verticesSet.end())
	{
		std::cout << *setIt;
		for (++setIt; setIt != _verticesSet.end(); ++setIt)
			os << ", " << *setIt;
	}
	os << "), resourceConsumption = " << _resConsumption;
}

double OverloadEliminationCut::getArcCoefficient(const int & tailVertId, const int & headVertId,
		                                         const double * tailResCons) const
{
	/// test whether head and tail belong to S
	bool tailIn = (tailVertId >= int(_verticesBitmap.size()))? false: _verticesBitmap[tailVertId];
	bool headIn = (headVertId >= int(_verticesBitmap.size()))? false: _verticesBitmap[headVertId];
	if (tailIn == headIn)
      return 0.0;

	// check the integrality of resource constraints (required for HECCs)
	double intTol = 1e-5;
	if ((int(tailResCons[_resourceId] - intTol) == int(tailResCons[_resourceId] + intTol)) &&
			(int(tailResCons[_resourceId] + intTol) > 0))
	{
		std::cerr << "ERROR: OverloadEliminationCut::getArcCoefficient called with tailResCons[" << _resourceId
				  << "]=" << tailResCons[_resourceId] << std::endl;
		exit(1);
	}
	int d = int(tailResCons[_resourceId] + intTol);

	int t1 = _sumDemands - _resConsumption - (_nbMachines-2)*(_resConsumption-1);
	// Set coefficients for variables leaving S
	if (tailIn && !headIn &&  (d >= _resConsumption))
	{
		if(d <= t1)
			return 1.0;
		else
			return 2.0;
	}

	if (!tailIn && headIn)
	{
		int t2 = _capacity - _sumDemands + _nbMachines*(_resConsumption - 1) + 1;
		int tz = t1 > t2 ? t1 : t2;
		if(d >= tz && d < _capacity)
			return -1.0;
		else
			return 0.0;
	}

	return 0.0;
}

double OverloadEliminationCut::getRouteCoefficient(const std::vector<int> & routeVertIds,
                                                   const std::vector<std::vector<double> > & routeResCons) const
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

GenericOverlEliminCutConstr::GenericOverlEliminCutConstr(Model * modelPtr,
                                                         ProbConfig * probConfPtr,
                                                         const std::string & name,
                                                         const Double & nonRootPriorityLevel,
                                                         const Double & rootPriorityLevel,
                                                         const int & maxNumCutsPerRound,
                                                         const std::vector<int> demands,
                                                         const std::vector<int> capacities,
                                                         const int & resourceId) :
  GenericExtendedArcCutConstr(modelPtr, probConfPtr, name, nonRootPriorityLevel, rootPriorityLevel),
  _maxNumCutsPerRound(maxNumCutsPerRound), _resourceId(resourceId), _demands(demands), _capacities(capacities)
{
	arturnewsep::InstanceInfo inst;
	inst.capacity = 0;
	std::vector<int> tmpDemands(_demands.size()+1, 0);
	inst.demand = &tmpDemands[0];
	for (int i = 1; i < int(_demands.size()+1); i++)
	{
		inst.demand[i] = _demands[i-1];
	}
	inst.nrootbranches = capacities.size();
	inst.numNodes = _demands.size()+1;

	// calculate the maximum capacity of a subproblem
	for (int i = 0; i < int(_capacities.size()); i++)
	{
		int c = _capacities[i];
		if (c > inst.capacity) inst.capacity = c;
	}

	// create the cut generator object
	cutGenObj = new arturnewsep::CutGenerator2(&inst);
}
    
GenericOverlEliminCutConstr::~GenericOverlEliminCutConstr()
{
     delete cutGenObj;
}

void GenericOverlEliminCutConstr
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
		const std::vector<int> & orderedArcIds = colPtr->spSol()->orderedIds();
		const std::vector<std::vector<double> > & resCons = colPtr->spSol()->resConsumption();

        /// for the moment, uses the transformation of the sequence of arc ids to the sequence of vertex ids
        /// TO DO : functions getRouteCoefficient and getArcCoefficient should be rewritten to deal directly
        ///         with the sequence of arc ids to increase the performance
        const NetworkFlow * netFlowPtr = colPtr->cgSpConfPtr()->networkFlowPtr();
        std::vector<int> vertexIds;
        vertexIds.reserve(orderedArcIds.size() + 1);
        std::vector<int>::const_iterator arcIdIt = orderedArcIds.begin();
        vertexIds.push_back(netFlowPtr->netArcPtr(*arcIdIt)->tailVertexPtr()->id());
        while (arcIdIt != orderedArcIds.end())
          {
            vertexIds.push_back(netFlowPtr->netArcPtr(*arcIdIt)->headVertexPtr()->id());
            ++arcIdIt;
          }

		std::vector<int>::const_iterator vertIt;
		std::vector<std::vector<double> >::const_iterator resConsIt;
		for (vertIt = vertexIds.begin(), resConsIt = resCons.begin(); vertIt != vertexIds.end(); ++vertIt, ++resConsIt)
		{
			/// current vertex in the route is *vertIt (0 is depot, client first index is 1)
			/// with cumulative demand (*resConsIt)[0] (including current vertex)
			/// value of column is colPtr->val()

            // skip auxiliary vertices
            if (*vertIt > _demands.size()+1)
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
	}

	// build the fractional solution information (inverting the arcs to follow the convention of the CutGenerator2)
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

	// Allocate a cut list
	std::vector<arturnewsep::Cut> cutVector;

	arturnewsep::SeparateOECbyHeuristic( cutGenObj, &sol, cutVector,
			_maxNumCutsPerRound, arturnewsep::PROB_PATH, 2);

	// Add the cuts to the problem
	for (i = 0; i < cutVector.size(); i++)
	{
		arturnewsep::ExtCCData* eccData = (arturnewsep::ExtCCData*)cutVector[i].data;
		std::set<int> vertSet;
		for (j = 0; j < int(eccData->S.size()); j++)
		{
			if (eccData->S[j]){ vertSet.insert(j); }
		}
		generatedCutConstrSet.insert(new OverloadEliminationCut(this, probConfPtr(), "OEC",
				cutVector[i].rhs, _resourceId, eccData->t, vertSet, eccData->weightS, eccData->capacity, eccData->m));
		// Release the cut data
		cutVector[i].DestroyData( cutVector[i].data );
	}
	if (printL(0))
	std::cout << "OEC added by heuristic. Total cuts: " << cutVector.size() << std::endl;

	// Artur: FOR DEBUGGING ONLY
	  // recalculate every cut violation
	 /* std::multiset<InstanciatedConstr*,CutSeparationPriorityComp>::iterator it;
	  for (it = generatedCutConstrSet.begin(); it != generatedCutConstrSet.end(); it++)
	    {
	      OverloadEliminationCut* h = dynamic_cast<OverloadEliminationCut*> (*it);
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

	  // check if the cuts invalidate the following manually entered solution for the instance wt40-2m-71
	const int numRtSol = 2;
	const int testSol[][22] = {
	 {34, 24, 10,5,19,6,37,7,14,4,31,2,40,15,13,11,36,26,12,0,0,0},
			{35,30,29,17,3,23,21,28,16,38,32,9,8,22,39,1,25,33,27,20,18,0} };

	std::multiset<InstanciatedConstr*,CutSeparationPriorityComp>::iterator it;
	for (it = generatedCutConstrSet.begin(); it != generatedCutConstrSet.end(); it++)
	{
		OverloadEliminationCut* h = dynamic_cast<OverloadEliminationCut*> (*it);
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
			std::cout << "ERROR: OEC cuts optimal solution!" << std::endl;
		}
	}*/
}

#endif /* RHECC_SEP_IS_FOUND */

#endif //BCP_RCSP_IS_FOUND
