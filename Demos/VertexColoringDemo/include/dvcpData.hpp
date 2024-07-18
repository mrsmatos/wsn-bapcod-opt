/**
 *
 * This file dvcpData.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef DvcpData_h
#define DvcpData_h

#include "bcUsefulHeadFil.hpp"

typedef std::set<int> Set;
typedef std::map<int, Set> IndexedSet;
typedef std::map<int, int> IndexedValue;

typedef std::vector<std::set<int> > VCPsolution;

struct DvcpData
{
  //number of nodes
  int numV;

  //number of edges
  int numE;

  // //set of adjacent nodes for node v
  IndexedSet N_v;

  // //set of non-adjacent nodes for node v
  IndexedSet A_v;

  //tail node of edge e
  IndexedValue T_e;

  //head node of edge e
  IndexedValue H_e;

  /**
   * Read instances of DIMACS standard format.
   */
  void readData(const std::string & inputFileName);
};

#endif // DvcpData_h
