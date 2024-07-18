/**
 *
 * This file bcBoundLevC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcBoundLevC.hpp"
#include "bcControlParameters.hpp"

const Double computeOptimalityGap(const Double & primalBd, const Double & dualBd)
{

  if (primalBd.isZero() && dualBd.isZero()) return Double::staticZero;

  Double product(primalBd);
  product *= dualBd;
  if (product.negative())  return Double::staticOne;

  if (primalBd > dualBd)
    {
      Double gap(primalBd - dualBd);
      gap/= Fabs(primalBd);

      if (printL(1))
	std::cout << "computeOptimalityGap() = " << gap  << std::endl;
      return gap;
    }
  Double gap(dualBd - primalBd);
  gap/= Fabs(dualBd);

  if (printL(1))
    std::cout << "computeOptimalityGap() = " << gap  << std::endl;
  return gap;

}

bool gapSmallerThanTol(const Bound & dualBound, const Bound & primalBound, const ControlParameters & param)
{
  if ( primalBound <= dualBound )
    return true;

  if ( Fabs(primalBound - dualBound) <= param.optimalityGapTolerance() )
    return true;

  return computeOptimalityGap(primalBound, dualBound).val() <= param.relOptimalityGapTolerance();
}
