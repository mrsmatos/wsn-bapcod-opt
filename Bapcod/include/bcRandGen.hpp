/**
 *
 * This file bcRandGen.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef RANDGEN_H
#define RANDGEN_H

#include "bcUsefulHeadFil.hpp"

/**
 * Generate unique random integer 
 * numbers from lb to ub:
 */
class URandGen 
{
  std::set<int> used;
  int modulus;
  int ifloor;
 public:
  URandGen(const int& ub, const int& lb = 1,const long& seed = (long) clock);
  int operator()();
};

/**
 * RandGen returns a random integer number uniformly distributed in [lb,ub];
 * when constructing this functor the users specify lb, ub and the seed for the random generator.
 * For a fixed seed, the function will always generate the same sequence of random number.
 */
class RandGen
{
  int modulus;
  int ifloor;
 public:
  RandGen(const int& ub, const int& lb = 1, const long& seed = (long)  clock);
  int operator()();
};

#endif // RANDGEN_H
