/**
 *
 * This file bcKnapsackSolverC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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
#include "bcPrintC.hpp"
#include "bcKnapsackSolverC.hpp"
#include "bcDoubleC.hpp"

/**
 * Generic code
 */
using namespace std;

#define P(i) #i << "=" << i


KnpItem::KnpItem(): varPtr(NULL), ref(0), p(0), w(0), m(0), icl(-1), ratio(0), x(0)
{
}

KnpItem::KnpItem(Variable * varPointer,
        const int & lItemRef,
        const double & lProfit,
        const double & lWeight,
        const int & lMultiplicity,
        const int & lRefToProd,
        const double & lRatio,
        const int & lxVal):
    varPtr(varPointer),
    ref(lItemRef),
    p(lProfit),
    w(lWeight),
    m(lMultiplicity),
    icl(lRefToProd),
    ratio(lRatio),
    x(lxVal)
{
}

KnpItem::KnpItem(const KnpItem & i)
{
  ref  = i.ref;
  varPtr = i.varPtr;
  p = i.p;
  w = i.w;
  m = i.m;
  icl  = i.icl;
  ratio = i.ratio;
  x = i.x;
}

KnpItem::~KnpItem()
{
}

std:: ostream& KnpItem::print(std::ostream& os) const
{
  os << "ref = " << ref <<  ", p = " << p << ", w = " << w << ", m = " <<m
     << ", icl = " << icl << ", ratio = " << ratio
     << ", x = " << x  << std::endl;
  return os;
}


const KnpItem & KnpItem::operator-=(const KnpItem &it1)
{
  w -= it1.w;
  p -= it1.p;
  ratio = p / w;
  return *this;
}

bool operator<(const KnpItem &it1, const KnpItem &it2)
{
  if (it1.ratio > it2.ratio) return true;
  if (it1.ratio < it2.ratio) return false;
  if (it1.w > it2.w) return true; // if same ratio, order by non-increasing weight
  if (it1.w < it2.w) return false;
  if (it1.m > it2.m) return true; // if same ratio and weight, order by non-increasing multiplicity
  return false;

}

bool BBfor01knap(const int & n, 
		 const  double * profit, 
		 const double * weight, 
		 const double & cap, 
		 double & inc, 
		 int * xsol,
		 int &countUBeval, 
		 int &countBBnode, 
		 const double & EPS, 
		 const bool & firstBoundInc)
{

  if (printL(6)) 
    std::cout << " BB01knap"<< std::endl;

  /** 
   * k is the level (depth) in the B&B tree 
   * z is a lower bound on the optimal value at the current node 
   */
  int  k, l, i;
  int* y = new int[n];
  for (l=0; l < n; l++) 
    xsol[l] = y[l] = 0;
  
  double z(0);
  double c(cap);
  double clp(0);
  double ub(0);
  //double inc(returnInc);
  k = 0;
  countBBnode = 0;
  countUBeval = 0;

  /// While problem not solved
  while (true) 
    {
      /// Sequence of Forward moves while don't need backtracking
      while (true)
        {
          /// Compute the LP upper bound for the current node
          countUBeval++;
          ub = z;
          clp = c;
          for(i= k; i < n; i++)
            {
	      // (weight[i] >  clp)
              if (weight[i] > clp + EPS)
	      //	if (weight[i] >= clp)
                {
                  ub += clp / weight[i] * profit[i];
                  break;
                }
              clp -= weight[i];
              ub += profit[i];
            }

          if (printL(6)) 
	    std::cout << "Upper Bound[" << countUBeval 
		      << "] = "  << ub 
		      << std::endl;

	  /// Prune node by bound if (ub <= inc): then go to backtracking
          if (ub < inc + EPS) 
	    //if (ub <= inc) 
	    break;

	  ///  While (next exists and fits: weight[k] <=  c) forward move
          while((k < n) && (weight[k] < c + EPS))
            {
              y[k] = 1;
              c -= weight[k];
              z += profit[k];
              if (printL(6)) 
		std::cout <<  "FORWARD MOVE take item " << k 
			  << "(w=" << weight[k] 
			  << ", p=" << profit[k] 
			  << ") ; cap =  " << c 
			  << " z = " << z 
			  << std::endl;
	      
              countBBnode++;
              k++;
            }
	  
	  /// Leaf node: record solution
          if ((k >= n) || ((double)c == 0.0))
            {
              ub = z;
              
	      /// Compare current int sol to incumbent
              if (z > inc + EPS) // (z > inc)
                {
                  inc = z;
                  if (printL(6)) 
		    std::cout << "new incumbent " << inc << std::endl;

                  for (l= 0; l < n; l++) 
		    xsol[l] = y[l];
		  
                  if (firstBoundInc)
                    {
                      delete[] y;
                      return(false);
                    }
                }
	      
	      /// Go to bactracking
              break;
            } 
	  // if (k >= n)
	  
	  /// Else branch on y[k] = 0
          y[k] = 0;
          
	  if (printL(6)) 
	    std::cout << "set item to zero" << k 
		      <<"(w="<< weight[k] 
		      << ", p=" << profit[k] 
		      << ") ; cap =  " << c 
		      << " z = " << z 
		      << std::endl;
          
	  countBBnode++;
          k++;
	  
	  /// Leaf node: record solution
          if (k >= n)
            {
              ub = z;
              
	      /// Compare current int sol to incumbent
              if (z > inc + EPS) //(z > inc) 
		{
                  inc = z;
                  if (printL(6)) 
		    std::cout << "new incumbent " << inc << std::endl;

                  for (l= 0; l < n; l++) 
		    xsol[l] = y[l];
		  
                  if (firstBoundInc)
                    {
                      delete[] y;
                      return(false);
                    }
                }
              if (y[n-1] != 0)
                {
                  std::cerr << "ERROR (y[k] != 0) " << std::endl;
                  exit(1);
                }
	      
	      /// Go to bactracking
              break; 
            } 
        }

      /// Backtracking
      l=k-1;

      /** 
       *  (y[k] == 1) the branch y[k] = 0 needs not be explored
       *  since if the current node has been pruned by bound,
       *  the  branch y[k] = 0 will also be pruned by bound
       *  (as the UB on that branch is at least as big as the
       *  ub on the branch y[k] = 1), and if at the node y[k] = 1
       *  an integer solution has been bound
       *  then no better int sol can be found at the node y[k] = 0.
       *  Therefore backtrack further 
       */
      if (y[l] == 1)
        {
          y[l] = 0;
          z -= profit[l];
          c += weight[l];
          if (printL(6)) 
	    std::cout << "untake last item " << l <<"(w="<< weight[l] 
		      << ", p=" << profit[l] 
		      << ") ; cap =  " << c 
		      << " z = " << z 
		      << std::endl;
          l--;
        }

      /// Go back up to the last level at which the branch  y[l] = 1 was chosen
      for (; l >= 0; l--)
        {
          if (y[l] == 0)
            {
              if (printL(6)) 
		std::cout <<  "unset from zero item " << l 
			  <<"(w="<< weight[l] 
			  << ", p=" << profit[l] 
			  << ") ; cap =  " << c 
			  << " z = " << z 
			  << std::endl;
            }
          else break;
        }
      
      /// No more bactracking possible, algorithm terminates
      if (l == -1) break;

      /// Now explore the branch y[l] = 0
      y[l] = 0;
      z -= profit[l];
      c += weight[l];
      if (printL(6)) 
	std::cout << "untake item and set it to zero  " << l 
		  <<"(w="<< weight[l] << ", p=" << profit[l] 
		  << ") ; cap =  " << c 
		  << " z = " << z 
		  << std::endl;
      
      /// Branch on y[l] = 0
      k = l + 1; 
      countBBnode++;
    }
  
  delete[] y;
  return(false);
}

bool binknap(const int & n, 
	         const int & ni,
	         const vector<double> & p,
	         const vector<double> & w,
	         vector<int> & ci,
             vector<KnpItem> &  vit, 
	         const double & minw,
	         const double & cap,
             int &countUBeval, 
	         int &countBBnode,
	         double &inc)
{
  /**
   * k is the level (depth) in the B&B tree
   * z is a lower bound on the optimal value at the current node
   */
  int  k, l, i;
  
  // Output("binknap");
  int* y = new int[n];
  for (l=0; l < n; l++) y[l] = 0;
  double z(0);
  double c(cap);
  double clp(0);
  double ub(0);

  std::string mes;

  k = 0;
  countBBnode = 0;
  countUBeval = 0;
  long start(clock());

  /// While problem not solved
  while (true)
    {
      /// Sequence of Forward moves while don't need backtracking
      while (true)
        {
          /// Compute the LP upper bound for the current node
          countUBeval++;
          if (printL(6))
            {
              mes = "bound eval ";
              std::cout << mes 
			<< P(k) 
			<< P(vit[k].icl) 
			<< P(ci[vit[k].icl])  
			<< P(c) 
			<< P(z) 
			<< P(countUBeval) 
			<< std::endl;
            }
	  
          ub = z;
          clp = c;
          for(i= vit[k].icl; i < ni; i++) 
	    if (ci[i] > 0)
	      {
		if (w[i] * ci[i] > clp)
		  {
		    ub += clp / w[i] * p[i];
		    if (printL(6))
		      {
			mes = "compute ub adding ";
			std::cout << mes 
				  << P(k) 
				  << P(i) 
				  << P(ci[i]) 
				  << P(clp) 
				  << P(ub) 
				  << std::endl;
		      }
		    break;
		  }
		clp -= w[i] * ci[i];
		ub += p[i] * ci[i];
		if (printL(6))
		  {
		    mes = "compute ub adding ";
		    std::cout << mes 
			      << P(k) 
			      << P(i) 
			      << P(ci[i]) 
			      << P(clp) 
			      << P(ub) 
			      << std::endl;
		  }
	      }
          if (printL(6))
	      std::cout << "BK Upper Bound ub = " << ub
		      << " nb " << countUBeval 
		      << "  clock = " << clock() - start 
		      << std::endl;

	  /// Prune node by bound: go to backtracking
          if (ub <= inc) 
	    break;

	  /// While (next exists and fits) FORWARD MOVE
          while((k < n) && (vit[k].w <= c)) 
            {
              y[k] = 1;
              c -= vit[k].w;
              ci[vit[k].icl] -= vit[k].m;
              z += vit[k].p;
              if (printL(6))
                {
                  mes = "FORWARD MOVE take item ";
                  std::cout << mes
			    << P(k) 
			    << P(y[k]) 
			    << P(vit[k].icl) 
			    << P(ci[vit[k].icl]) 
			    << P(c) 
			    << P(z) 
			    << P(countBBnode) 
			    << std::endl;
                }
              countBBnode++;
              k++;
            }

	  /// Leaf node: record solution
          if ((k >= n) || (zeroTest(c)))
            {
              ub = z;
              
	      /// Compare current int sol to incumbent
              if (z > inc)
                {
                  inc = z;
                  if (printL(6))
                    {
                      mes = "new incumbent ";
                      std::cout << mes << P(inc) << std::endl;
                    }
                  for (l= 0; l < n; l++)
                    {
                      vit[l].x = y[l];
                      if (printL(6))
                        {
                          mes = "contains ";
                          std::cout << mes 
				    << P(l) 
				    << P(y[l]) 
				    << std::endl;
                        }
                    }
                }
	      
	      /// Go to bactracking 
              break; 
            }
	  
          /// Else branch on y[k] = 0
          y[k] = 0;
          ci[vit[k].icl] -= vit[k].m;

          if (printL(6))
            {
              mes = "set item to zero";
              std::cout << mes 
			<< P(k)  
			<< P(y[k]) 
			<< P(vit[k].icl) 
			<< P(ci[vit[k].icl]) 
			<< P(countBBnode) 
			<< std::endl;
            }
          countBBnode++;
          k++;
	  
	  /// Leaf node: record solution
          if (k >= n)
            {
              ub = z;
          
	      /// Compare current int sol to incumbent 
              if (z > inc)
                {
                  inc = z;
                  if (printL(6))
                    {
                      mes = "new incumbent ";
                      std::cout << mes 
				<< P(inc) 
				<< std::endl;
                    }
                  for (l= 0; l < n; l++)
                    {
                      vit[l].x = y[l];
                      if (printL(6))
                        {
                          mes = "contains ";
                          std::cout << mes 
				    << P(l) 
				    << P(y[l]) 
				    << std::endl;
                        }
                    }
                }
              if (y[n-1] != 0)
                {
                  ///@todo: check if those lines are ok.
                  //std::cout << "ERROR (y[k] != 0) " << std::endl;
                  std::cout << "binknap(): ERROR (y[k] != 0)" << std::endl;
                  //check(1, "binknap(): ERROR (y[k] != 0)");
                }
	      
	      /// Go to bactracking 
              break;
            }
        } 
      /// End while forward sequence: exited when node pruned by optimality 

      /// Backtracking 
      l=k-1;

      /** (y[k] == 1) the branch y[k] = 0 needs not be explored
       *  since if the current node has been pruned by bound,
       *  the  branch y[k] = 0 will also be pruned by bound
       *  (as the UB on that branch is at least as big as the
       *  ub on the branch y[k] = 1), and if at the node y[k] = 1
       *  an integer solution has been bound
       *  then no better int sol can be found at the node y[k] = 0.
       *  Therefore backtrack further 
       */
      if (y[l] == 1)
        {


          y[l] = 0;
          z -= vit[l].p;
          c += vit[l].w;
          ci[vit[l].icl] += vit[l].m;

          if (printL(6))
            {
              mes = "untake last item ";
              std::cout << mes 
			<< P(l) 
			<< P(y[l]) 
			<< P(vit[l].icl) 
			<< P(ci[vit[l].icl])  
			<< P(c) 
			<< P(z) 
			<< P(countBBnode) 
			<< std::endl;
            }
          l--;
        }

      /// Go back up to the last level at which the branch  y[l] = 1 was chosen 
      for (; l >= 0; l--)
        {
          if (y[l] == 0)
            {
              ci[vit[l].icl] += vit[l].m;
              if (printL(6))
                {
                  mes = "unset from zero item ";
                  std::cout << mes 
			    <<P(l)  
			    << P(y[l]) 
			    << P(vit[l].icl) 
			    << P(ci[vit[l].icl]) 
			    << std::endl;
                }
            }
          else 
	    break;
        }
      /// No more bactracking possible, algorithm terminates
      if (l == -1) 
	break;

      /// Now explore the branch y[l] = 0 
      y[l] = 0;
      z -= vit[l].p;
      c += vit[l].w;
      /// Do not add multiplicty back into icl capacity because branch on y[l] = 0
      // ci[vit[l].icl] += vit[l].m;
      
      /// Branch on y[l] = 0
      k = l + 1;
      countBBnode++;
      if (printL(6))
        {
          mes = "untake item and set it to zero ";
          std::cout << mes 
		    <<P(l)  
		    << P(y[l]) 
		    << P(vit[l].icl) 
		    << P(ci[vit[l].icl])
		    << P(c) 
		    << P(z) 
		    << P(countBBnode) 
		    << std::endl;
        }
    } 
  
  /// End while problem not solved 
  delete[] y;
  return(false);
}

bool MCbinknap(const int & n, 
	           const int & ni,
	           const vector<double> & p,
	           const vector<double> & w,
               vector<int> & classBound, 
	           vector<int> & maxci,
               vector<KnpItem> &  vit, 
	           const double & minw,
               const double &  cap,
	           int &countUBeval,
	           int &countBBnode,
               double &inc,
	           const bool & PijEqualMijtimesPi)
{
  int  k, l, i;
  std::string mes;

  /** 
   * k is the level (depth) in the B&B tree 
   * z is a lower bound on the optimal value at the current node 
   */

  int* y = new int[n];
  
  /// Current solution being constructed
  for (l=0; l < n; l++) 
    y[l] = 0;

  /// Primal solution value associated to best solution x
  double z(0);
  
  /// Current residual capacity given partial solution y
  double c(cap);
  
  // Current residual capacity while computing lp dual bound
  double clp(0);
  
  // Current lp dual bound
  double ub(0);
  
  /**
   *i = product index
   * k = item index
   * Current item that is considered
   */
  k = 0;
  
  /**
   * ni = number of product
   * classBound = residual upper bound on the number of copies of product i given fixation to 1 in partial solution y
   * maxci = residual maximum number of copies of product i that could be achieved given fixation to 1 or 0 in partial solution y
   */
  countBBnode = 0;
  countUBeval = 0;

  /// While problem not solved 
  while (true) 
    {
      /// Sequence of forward moves while don't need backtracking 
      while (true)
        {
          /**
	   * Compute the LP upper bound for the current node 
	   * Begin 
	   */
          countUBeval++;
          if (printL(6))
            {
              mes = "bound eval ";
              std::cout << mes 
			<<  P(k) 
			<< P(vit[k].icl)  
			<< P(classBound[vit[k].icl]) 
			<< P(maxci[vit[k].icl]) 
			<<P(c) 
			<< P(z) 
			<< P(countUBeval) 
			<< std::endl;
            }
          if(PijEqualMijtimesPi)
            {
	      /**
	       * Minimum of the upper bounds on the number 
	       * of copies that can be taken for current class i
	       */
              int minci(0);
              ub = z;
              clp = c;
              l = k + 1;
              for(i = vit[k].icl; i < ni;)
                {
                  minci = intMin(classBound[i],maxci[i]);
                  if (minci > 0)
                    {

                      if (w[i] * minci > clp)
                        {
                          ub += clp / w[i]  * p[i];
                          if (printL(6))
                            {
                              mes = "compute ub adding ";
                              std::cout << mes 
					<< P(k) 
					<< P(i) 
					<< P(minci) 
					<< P(clp) 
					<< P(ub) 
					<< std::endl;
                            }
			  
                          break;
                        }
		      
                      clp -= w[i] * minci;
                      ub += p[i] *  minci;
                      if (printL(6))
                        {
                          mes = "compute ub adding ";
                          std::cout <<mes 
				    << P(k) 
				    << P(i) 
				    << P(minci) 
				    << P(clp) 
				    << P(ub) 
				    << std::endl;
                        }
		      
                    }
                  for (; l < n; l++) 
		    {
		      if (vit[l].icl != i) 
			break;
		    }
		  
                  if (l == n) 
		    /// Exit for (i ...
		    break;
		  
                  i = vit[l].icl;
                }
            }
          else
            {
              double ub2(z);
              int fm(0);
              int *cli = new int[ni];

              for(i= 0; i < ni; i++) 
		cli[i] = intMin(classBound[i],maxci[i]) ;

              for(clp = c, l = k; l < n ; l++) 
		if (cli[vit[l].icl] > 0)
		  {
		    if (vit[l].m <= cli[vit[l].icl])
		      {
			if (vit[l].w <= clp) 
			  // && (vit[l].m <= cli[vit[l].icl])
			  {
			    clp -= vit[l].w;
			    cli[vit[l].icl] -= vit[l].m;
			    ub2 += vit[l].p;
			    if (printL(6))
			      {
				mes = "compute ub2 adding ";
				std::cout << mes 
					  << P(l) 
					  << P(vit[l].icl) 
					  << P(cli[vit[l].icl]) 
					  << P(clp) 
					  << P(ub2) 
					  << std::endl;
			      }
			  }
			/// (vit[l].w > clp) but (vit[l].m <= cli[vit[l].icl])
			else
			  {
			    ub2 +=  clp / vit[l].w * vit[l].p;
			    if (printL(6))
			      {
				mes = "compute ub2 adding ";
				std::cout << mes 
					  << P(l) 
					  << P(vit[l].icl) 
					  << P(cli[vit[l].icl]) 
					  << P(clp) 
					  << P(ub2) 
					  << std::endl;
			      }
			    break;
			  }
		      }
		    /// (vit[l].m > cli[vit[l].icl])
		    else
		      {
			fm = cli[vit[l].icl];

			/// Frac fits
			if (fm * w[vit[l].icl] <= clp)
			  {
			    ub2 +=  vit[l].p * fm / vit[l].m; // fm * p[vit[l].icl];
			    clp -= fm * w[vit[l].icl];
			    cli[vit[l].icl] -= fm;
			    if (printL(6))
			      {
				mes = "compute ub2 adding ";
				std::cout << mes 
					  << P(l) 
					  << P(vit[l].icl) 
					  << P(cli[vit[l].icl]) 
					  << P(clp) 
					  << P(ub2) 
					  << std::endl;
			      }
			  }
			/// (vit[l].m > cli[vit[l].icl]) && (fm * vit[l].w > clp)
			else
			  {
			    ub2 += (double) clp / (double) vit[l].w * vit[l].p;
			    if (printL(6))
			      {
				mes = "compute ub2 adding ";
				std::cout << mes 
					  << P(l) 
					  << P(vit[l].icl) 
					  << P(cli[vit[l].icl]) 
					  << P(clp) 
					  << P(ub2) 
					  << std::endl;
			      }
			    
			    break;
			  }
		      }
		  }
	      
              delete[] cli;
              ub = ub2;
            }


	  /// Prune node by bound: go to backtracking 
          if (ub <= inc) break;

	  /// While (next exists and fits) FORWARD MOVE 
          while((k < n) && (vit[k].w <= c) && (vit[k].m <= classBound[vit[k].icl]))
            {
              y[k] = 1;
              c -= vit[k].w;
              // minci[vit[k].icl] -= vit[k].m;
              classBound[vit[k].icl] -= vit[k].m;
              maxci[vit[k].icl] -= vit[k].m;
              z += vit[k].p;
              if (printL(6))
                {
                  mes = "FORWARD MOVE take item ";
                  std::cout << mes
			    << P(k)  
			    << P(y[k]) 
			    << P(vit[k].icl) 
			    << P(classBound[vit[k].icl]) 
			    << P(c) 
			    << P(z) 
			    << P(countBBnode)
			    << std::endl;
                }
              countBBnode++;
              k++;
            }

	  /// Leaf node: record solution
          if ((k >= n) || (zeroTest(c)))
            {
              ub = z;
	      
              /// Compare current int sol to incumbent 
              if (z > inc)
                {
                  inc = z;
                  if (printL(6))
                    {
                      mes = "new incumbent ";
                      std::cout << mes << P(inc) << std::endl;
                    }
		  
                  for (l= 0; l < n; l++)
                    {
                      vit[l].x = y[l];
                      if (printL(6))
                        {
                          mes = "contains ";
                          std::cout << mes << P(l) << P(y[l]) << std::endl;
                        }
                    }
                }

	      /// Go to bactracking 
              break; 
            }

          /// else branch on y[k] = 0
          y[k] = 0;
          maxci[vit[k].icl] -= vit[k].m;
          if (printL(6))
            {
              mes = "set item to zero";
              std::cout << mes 
			<< P(k)  
			<< P(y[k]) 
			<< P(vit[k].icl) 
			<< P(classBound[vit[k].icl]) 
			<< P(countBBnode) 
			<< std::endl;
            }
	  
          countBBnode++;
          k++;
	  
	  /// Leaf node: record solution
          if (k >= n)
            {
              ub = z;
              /// Compare current int sol to incumbent 
              if (z > inc)
                {
                  inc = z;
                  if (printL(6))
                    {
                      mes = "new incumbent ";
                      std::cout << mes 
				<< P(inc) 
				<< std::endl;
                    }
                  for (l= 0; l < n; l++)
                    {
                      vit[l].x = y[l];
                      if (printL(6))
                        {
                          mes = "contains ";
                          std::cout << mes 
				    << P(l) 
				    << P(y[l]) 
				    << std::endl;
                        }
                    }
                }
              if (y[n-1] != 0)
                {
                  std::cout << "MCbinknap(): ERROR (y[k] != 0)" << std::endl;
                }
	      
	      /// Go to bactracking 
              break; 
            } // if (k >= n)
        } /// End while forward sequence: exited when node pruned by optimality 

      /// Backtracking 
      l=k-1;

      /// Dominance Test
      
      /**
       * (y[k] == 1) the branch y[k] = 0 needs not be explored
       * since if the current node has been pruned by bound,
       * the  branch y[k] = 0 will also be pruned by bound
       * (as the UB on that branch is at least as big as the
       * ub on the branch y[k] = 1), and if at the node y[k] = 1
       * an integer solution has been bound
       * then no better int sol can be found at the node y[k] = 0.
       * Therefore backtrack further 
       */
      if (y[l] == 1)
        {
          y[l] = 0;
          z -= vit[l].p;
          c += vit[l].w;
          classBound[vit[l].icl] += vit[l].m;
          maxci[vit[l].icl] += vit[l].m;
          if (printL(6))
            {
              mes = "untake last item ";
              std::cout << mes 
			<< P(l) 
			<< P(y[l]) 
			<< P(vit[l].icl) 
			<< P(classBound[vit[l].icl])  
			<< P(c) 
			<< P(z) 
			<< P(countBBnode) 
			<< std::endl;
            }
          l--;
        }

      /// Go back up to the last level at which the branch  y[l] = 1 was chosen 
      for (; l >= 0; l--)
        {
          if (y[l] == 0)
            {
              maxci[vit[l].icl] += vit[l].m;
              if (printL(6))
                {
                  mes = "unset from zero item ";
                  std::cout << mes 
			    << P(l) 
			    << P(y[l])  
			    << P(vit[l].icl) 
			    << P(classBound[vit[l].icl]) 
			    << std::endl;
                }
            }
          else 
	    break;
        }
      if (l == -1)
	/// No more bactracking possible, algorithm terminates
	break; 

      /// Now explore the branch y[l] = 0 
      y[l] = 0;
      z -= vit[l].p;
      c += vit[l].w;
      classBound[vit[l].icl] += vit[l].m;
      
      /// Do not add multiplicty back into icl capacity because branch on y[l] = 0
      // maxci[vit[l].icl] += vit[l].m;
      
      /// Branch on y[l] = 0
      k = l + 1;
      countBBnode++;
      if (printL(6))
        {
          mes = "untake item and set it to zero ";
          std::cout << mes 
		    << P(l) 
		    << P(y[l])  
		    << P(vit[l].icl) 
		    << P(classBound[vit[l].icl]) 
		    << P(c) 
		    << P(z) 
		    << P(countBBnode) 
		    << std::endl;
        }
    } 
  /// End while problem not solved 
  
  delete[] y;
  return(0);
}
