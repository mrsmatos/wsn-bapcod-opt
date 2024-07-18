/**
 *
 * This file bcIndexC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcIndexC.hpp"
#include <boost/container_hash/hash.hpp>

std::size_t ihash::operator()(MultiIndex const & e) const
{
    std::size_t seed = 0;
	int f;
	for (int i=e.endPosition-1; i>=0; i--)
	{
		f=e.index(i);
		boost::hash_combine(seed,f);
	}
	return seed;
}
