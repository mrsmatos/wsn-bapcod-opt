/**
 *
 * This file dvcpHeuristics.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dvcpData.hpp"
#include "dvcpHeuristics.hpp"

#ifndef _MSC_VER
extern "C"
{
#include "ls_greedy.h"
}
#endif


void runSimpleGreedyHeuristic(DvcpData & data, const int printLevel, VCPsolution & solution)
{
  std::vector<int> vertexIsTaken(data.numV, false);
  int numTakenVertices = 0;
  int columnIndex = 0;
  Set::iterator setIt;

  if (printLevel >= 1)
      std::cout << "Simple greedy heuristic is run" << std::endl;

  while (numTakenVertices < data.numV)
  {
	  std::vector<int> hereVertexIsTaken(vertexIsTaken);
	  std::vector<int> verticesInCurrentStableSet;

	  int hereNumTakenVertices = 0;
	  int currentIndex = 0;
	  while (hereNumTakenVertices < data.numV - numTakenVertices)
	  {
		  while (hereVertexIsTaken[currentIndex])
		      currentIndex += 1;
		  verticesInCurrentStableSet.push_back(currentIndex);
		  hereNumTakenVertices += 1;
	      hereVertexIsTaken[currentIndex] = true;

		  for (setIt = data.N_v[currentIndex].begin(); setIt != data.N_v[currentIndex].end(); ++setIt)
		  {
		      hereVertexIsTaken[*setIt] = true;
			  hereNumTakenVertices += 1;
		  }
		  currentIndex += 1;
	  }

      solution.push_back(std::set<int>());
      if (printLevel >= 1)
        std::cout << "Color " << columnIndex + 1 << ":";
	  for (std::vector<int>::iterator vecIt = verticesInCurrentStableSet.begin();
           vecIt != verticesInCurrentStableSet.end(); ++vecIt)
	  {
          if (printLevel >= 1)
  	        std::cout << " " << *vecIt;
	      solution[columnIndex].insert(*vecIt);
	      vertexIsTaken[*vecIt] = true;
	      numTakenVertices += 1;
	  }
	  columnIndex += 1;
      if (printLevel >= 1)
  	    std::cout << std::endl;
  }
}

#ifndef _MSC_VER

void runDsaturGreedyHeuristic(DvcpData & data, const int printLevel, VCPsolution & solution)
{
    int ncount = data.numV;
    int ecount = data.numE;
    std::vector<int> elist(ecount * 2, 0);
    COLORset * colorclasses = (COLORset *)NULL;
    int ncolors = 0;

    if (printLevel >= 1)
        std::cout << "DSATUR  heuristic is run" << std::endl;

    for (int e = 0; e < ecount; e++)
    {
        elist[e*2] = data.T_e[e];
        elist[e*2+1] = data.H_e[e];
    }
    
    COLORdsatur(ncount, ecount, &elist[0], &ncolors, &colorclasses);
    
    for (int colorIndex = 0; colorIndex < ncolors; ++colorIndex)
    {
        if (printLevel >= 1)
            std::cout << "Color " << colorIndex + 1 << ":";
        solution.push_back(std::set<int>());
        for (int i = 0; i < colorclasses[colorIndex].count; i++)
        {
            if (printLevel >= 1)
                std::cout << " " << colorclasses[colorIndex].members[i];
            solution[colorIndex].insert(colorclasses[colorIndex].members[i]);
        }
        if (printLevel >= 1)
            std::cout << std::endl;
    }
    COLORfree_sets(&colorclasses, &ncolors);
}

#endif

