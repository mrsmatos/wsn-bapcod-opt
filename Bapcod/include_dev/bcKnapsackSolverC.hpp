/**
 *
 * This file bcKnapsackSolverC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef KNAPSACKSOLVERCLASS_H
#define KNAPSACKSOLVERCLASS_H

#include "bcUsefulHeadFil.hpp"

class Variable;

struct KnpItem
{
  Variable * varPtr;
  int ref;
  double  p; /* profit                  */
  double  w; /* weight                  */
  int m; /* multiplicity */
  int  icl; /* reference to prod   */
  double ratio;
  int x;
  
  KnpItem(); 
  KnpItem(Variable * varPointer,
	  const int & lItemRef,
	  const double & lProfit,
	  const double & lWeight,
	  const int & lMultiplicity,
	  const int & lRefToProd = -1,
	  const double & lRatio = 0,
	  const int & lxVal = 0);

  KnpItem(const KnpItem & i);

  virtual std:: ostream& print(std::ostream& os =std::cout ) const;


  friend bool operator<(const KnpItem &it1, const KnpItem &it2);

  const KnpItem & operator-=(const KnpItem &it1);
  virtual ~KnpItem();
};

inline std::ostream& operator<<(std::ostream& os,  KnpItem & that)
{return that.print(os);}

extern bool binknap(const int & n, const int & ni,  const std::vector<double> & p, const std::vector<double> & w,  std::vector<int> & ci,
                    std::vector<KnpItem> &  vit, const double & minw, const double & cap,
                    int &countUBeval, int &countBBnode, double &inc);

extern bool MCbinknap(const int & n, const int & ni,  const std::vector<double> & p, const std::vector<double> & w, std::vector<int> & ci, std::vector<int> & maxci,
                      std::vector<KnpItem> &  vit, const double & minw,
                      const double &  cap, int &countUBeval, int &countBBnode, double &inc, const bool & PijEqualMijtimesPi = false);


extern bool BBfor01knap(const int & n, const  double * profit, const double * weight, const double & cap, double &inc, int * xsol, int &countUBeval, int &countBBnode, const double & EPS = 0.00001, const bool & firstBoundInc = false); // Assumes that the items are sorted in order of non-increasing ratio (profit/weight)


#endif // KNAPSACKSOLVERCLASS_H
