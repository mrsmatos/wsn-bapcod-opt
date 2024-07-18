/**
 *
 * This file bcRandGen.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcUsefulHeadFil.hpp"
#include "bcRandGen.hpp"
#include "bcPrintC.hpp"

/**
 * Generate unique random integer numbers from lb to ub:
 */
URandGen::URandGen(const int& ub, const int& lb, const long& seed) :
  modulus(ub - lb + 1), ifloor(lb)
{
  srand(seed);

  return;
}

int URandGen::operator()()
{
  while(true)
    {
      int i = (int) rand() % modulus;
      if (used.find(i) == used.end())
        {
          used.insert(i);
          return(ifloor + i);
        }
    }
}

/**
 * Generate random integer numbers from lb to ub:
 */

/// RandGen() returns a random integer number uniformly distributed in [lb,ub]; 
/// when constructing this functor the users specify lb, ub and the seed for the random generator. 
/// For a fixed seed, the function will always generate the same sequence of random number.

RandGen::RandGen(const int& ub, const int& lb, const long& seed) :
  modulus(ub - lb + 1), ifloor(lb)
{
  srand(seed);
  
  return;
}

int RandGen::operator()()
{
  return ifloor + (int) rand() % modulus;
}
