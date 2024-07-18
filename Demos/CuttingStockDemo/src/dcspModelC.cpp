/**
 *
 * This file dcspModelC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "dcspModelC.hpp"
#include "bcModelingLanguageC.hpp"
#include "dcspColGenSubproblem.hpp"
#include "dcspCutSeparationRoutines.hpp"

// a/**
// *  \file
// *  Problem:<br/>
// *  "Given set of items, \f$ I=\{1, \ldots, n\} \f$, with given width \f$ w_i \f$ and demand \f$ d_i \f$, the Cutting 
// *  Stock Problem (CSP)
// *  consists in assigning items to bins in such a way that total width of item assigned to a bin does not exceed 
// *  the bin capacity  \f$ C \f$, 
// *  the number of copy of each item equal its demand, and the number of bins used is minimized."
// *  \f{eqnarray*}{
// *     min \sum_{=1}^{n} y_{k} \\
// *   \sum_{k=1}^{n} x_{ik} = d_i \qquad i \in V \\
// *    \sum_{i=1}^{n} w_i x_{ik} \leq C  y_{k} \qquad k = 1, ..., n \\
// *   x_{ik} \in \{0, 1, \ldots , \min\{C / w_i , d_i\}\} \qquad i \in V, k = 1, ..., n \\
// *     y_{k} \in  \{0, 1\} \qquad k = 1, ..., n
// *  \f}
// * 
// *  A Cardinality Cut can be added of the form \f$ \sum_{=1}^{n} y_{k} \geq L \f$, 
// *  where \f$ L \f$  is computed as the rounded up value of  \f$ \sum_{=1}^{n} y_{k}  \f$ in the LP solution.
// */ 
// a/// ***************************************************************************************************************

/*!
 \file
 
 this file, we demonstrate how to formulate the cutting stock problem in Bapcod in order to use the power of this 
 framework to solve it using column generation. The following lecture will hopefully facilitate the understanding 
 of the code.
 
 First off, we will give the description of this problem and provide 
 the classical formulation of it. Then we will show how one could efficiently solve this problem using the column 
 generation method. Finally we will explain how by simply noticing that the Dantzig-Wolfe decomposition applies 
 to the classical formulation, We can let the generic framework Bapcod solve it using a hidden column generation 
 method that is pretty similar to the one we describe.
 
 Consider a mill of paper that has \f$K\f$ paper rolls of fixed width \f$W\f$. 
 The clients are requesting paper rolls of different widths. 
 The objective is to cut the paper rolls in order to minimize the left over. 
 Let us consider that every client \f$i \in \{1, ..., M\}\f$ is requesting 
 a demand \f$D_i\f$ of rolls having a width \f$W_i\f$ with \f$ W_i < W \f$. 
 Using two variables, we can write the classical IP formulation of this problem. 
 For an initial roll of paper \f$k \in \{1, ..., K\}\f$, let \f$y^k\f$ be the 
 binary variable that is set to \f$1\f$ if the initial paper roll is cut and \f$0\f$ otherwise. 
 Let the integer variable \f$x^k_i\f$ determines how many times item \f$i\f$ is cut on roll \f$k\f$. 
 The following IP1 describes formally our problem:
 
 \f{eqnarray*}{
 (IP1)\min \sum_{k=1}^K y^k \\
 s.t. \sum_{k=1}^K x_i^k \ge D_i  \qquad i = 1, ..., M  \\
 \sum_{i=1}^M W_i x_i^k \le W y^k  \qquad k = 1, ..., K  \\
 x_i^k \in Z^+, y^k \in \{0,1\} \qquad i = 1, ..., M \qquad k = 1, ..., K.
 \f} 
  
 This formulation is known to be computationally hard to solve because of its poor relaxation. 
 Hence it is preferable to use an alternative formulation with a column generation schema where each 
 column corresponds to a pattern for cutting the initial rolls. For a pattern \f$j\f$ let us consider 
 the variable \f$a_{ij}\f$ to be the number of times item \f$i\f$ is cut in pattern \f$j\f$. 
 Let \f$u_j\f$ be the number of initial rolls cut according to the pattern \f$j\f$. 
 If \f$\mathcal{J}\f$ is a set of patterns that is “sufficient” to reach the optimal, 
 the cutting stock problem can be described by the following IP2:
 
 \f{eqnarray*}{
 (IP2) \min \sum_{j \in \mathcal{J}} u_j \\
 s.t. \sum_{j \in \mathcal{J}} a_{ij} u_j \ge D_i  \qquad i = 1, ..., M  \\
 u_j \in Z^+ \qquad j \in \mathcal{J}.
 \f} 
 
 The issue now is that the total number of possible patterns can be exponentially large which makes it hard 
 to even solve the linear relaxation of it that we call LP2. However we know that the number of non-zero variables 
 (the basis variables) is bounded by the number of constraints, i.e., \f$M\f$. Hence, even if the number of possible 
 variables (columns) may be large, we only need a small subset of these columns in the optimal
 solution of LP2. The main idea of column generation is to start with a small 
 subset \f$\mathcal{P} \subset \mathcal{J}\f$ such that the following subproblem RLP2 is feasible 
 and generate ``best'' candidate patterns (eventually remove bad ones) as far as some which are not 
 in \f$\mathcal{P}\f$ have negative reduced costs.
 
 \f{eqnarray*}{
 (RLP2) \min \sum_{j \in \mathcal{J}} u_j \\
 s.t. \sum_{j \in \mathcal{J}} a_{ij} u_j \ge D_i  \qquad i = 1, ..., M  \\
 u_j \in R^+ \qquad j \in \mathcal{P}.
 \f}
 
 Given the optimal dual solution \f$q\f$ of (RLP2), the reduced cost of a pattern (column) \f$j \in \mathcal{J} 
 \setminus \mathcal{P}\f$ is \f$1 -  \sum_{i=1}^M a_{ij} q_i\f$, 
 and the best pattern to add is the one minimizing this value. 
 The components of such a pattern can be found by solving the following knapsack problem (IP3).
 
 \f{eqnarray*}{
 \max \sum_{i=1}^M  q_i a_i \\
 s.t. \sum_{i=1}^M W_i a_{i} \le W,  \\
 a_i  \in Z^+ \qquad i = 1, ..., M.
 \f}
 
 Now let us go back to IP1 and see how we need to implement it in Bapcod in order to solve it with the column 
 generation method (see the function ::buildModel). One can notice that after the \f$ M \f$ first constraints, 
 we have \f$ K \f$ constraints 
 that are pair-wise independant as they don't share any variables in common. 
 One can also notice that those constraints together with their variables are different instances
 of the same subroblem. Therefore we only define one subproblem "colGenSp[0]", and how each of its solutions
 is involved in the global objective value as well as in the \f$ M \f$ linking constraints. 
 
 Remark that those subproblems are similar to IP3. Hence, we provide an oracle to solve 
 IP3 (that belong to the well-known class of knapsack probklems). 
 However, this is not necessary as otherwise Bapcod uses the defined solver (ex: Cplex) 
 for solving those subproblems. This is the case of the similar demonstration that you can find under the name: 
 SimpleCuttingStock
 
 Finally, We added a Cardinality Cut (implemented in dcspCutSeparationRoutines.cpp) 
 of the form \f$ \sum_{=1}^{n} y_{k} \geq L \f$, 
 where \f$ L \f$  is computed as the rounded up value of  \f$ \sum_{=1}^{n} y_{k}  \f$ in the LP solution.
 ***************************************************************************************************************
 */ 

/*! 
 \mainpage 
 @copydoc dcspModelC.cpp 
*/

using namespace std;

void buildModel(DcspData& data, BcModel & dcspModel)
{
  int nbOfOrders = data.orderPtrVector().size();
  Double maxNbBin = 0;
  Double minNbBin = 0;
  for (int orderIndex = 0; orderIndex < nbOfOrders; ++orderIndex)
    {
      maxNbBin += data.orderPtrVector()[orderIndex]->demand();
      minNbBin += data.orderPtrVector()[orderIndex]->demand() * data.orderPtrVector()[orderIndex]->width();
    }
  minNbBin = Dceil (minNbBin / data.stockSheet().width());

 
  std::cout << " nbOfOrders = " << nbOfOrders << std::endl; 
  std::cout << " minNbBin = " << minNbBin << std::endl; 
  std::cout << " maxNbBin = " << maxNbBin << std::endl; 


  /// create the  model objective and assigned its order of magnitude
  /// **************************************************************
  BcObjective objective(dcspModel);
  objective.setStatus(BcObjStatus::minInt);
  objective <= maxNbBin;
  objective.setArtCostValue((maxNbBin + minNbBin) / 2);

  /// create a handle on model  master
  /// *********************************
  BcMaster master(dcspModel);

  /// create a master constraint object
  /// ********************************
  BcConstrArray covConstr(master, "COV");

  /// Assign a default sense and rhs to constraints of type  covConstr 
  /// ****************************************************************
  covConstr == 1;

  /// Define new object of column generation subproblem
  /// ************************************************
  BcColGenSpArray colGenSp(dcspModel);

  /// Define a pricing problem solver
  /// *******************************
  IntKnapsackSubProb * oraclePtr = new IntKnapsackSubProb(&data);

  /// Create a new column generation supbroblem of index '0' and assign to it an solver (the above oracle)
  /// *********************************************************************************
  colGenSp(0).attach(oraclePtr);

  /// Define bound on the number of subproblem '0' solutions used in master
  /// *********************************************************************
  colGenSp[0] <= maxNbBin;
  colGenSp[0] >= minNbBin;

  /// Create subproblem variable class atteached to subproblem '0'
  /// ***********************************************************
  BcVarArray xVar(colGenSp[0], "X");
  BcVarArray yVar(colGenSp[0], "Y");

  /// set the defalut type and index reference letter for Variables of type xVar
  /// (letter assigned to index improves readability of variable name)
  /// **************************************************************************
  xVar.type('I');
  xVar.defineIndexNames(MultiIndexNames('i'));
  xVar.priorityForMasterBranching(-1);
  xVar.priorityForSubproblemBranching(1);

  /// set the default type and lower bound for Variables of type yVar and set a negative branching priority
  /// *****************************************************************************************************
  yVar.type('B');
  yVar >= 1;
  yVar.priorityForMasterBranching(-1);
  yVar.priorityForSubproblemBranching(-1);

  /// Create a variable of class yVar with index 0 
  /// *********************************************
  yVar(0);

  /// create subproblem constraint class attached to a subproblem
  /// ************************************************************
  BcConstrArray knpConstr(colGenSp[0], "KNP");

  /// Create a constraint of class knpConstr with index 0 and set its rhs
  /// *******************************************************************
  knpConstr(0) <= 0;

  /// get constraint of class knpConstr with index 0 and set  yVar[0] coeffeicient in the constraint
  /// **********************************************************************************************
  knpConstr[0] -=  data.stockSheet().width() * yVar[0];

  for (int orderIndex = 0; orderIndex < nbOfOrders; ++orderIndex)
    {
      /// Create an instanciation of a variable of class xVar with index 'jobIndex' and set its upper bound
      /// ************************************************************************************************
      xVar(orderIndex) <= Dmin( data.orderPtrVector()[orderIndex]->demand(),
                                Dfloor(data.stockSheet().width()/data.orderPtrVector()[orderIndex]->width()));

      /// set branching priority on  xVar 
      /// *******************************
      xVar[orderIndex].branchingPriority(data.orderPtrVector()[orderIndex]->width());

      /// add coefficient of  variable xVar[orderIndex] in the constraint of type knpConstr with index 0
      /// ***********************************************************************************************
      knpConstr[0] += data.orderPtrVector()[orderIndex]->width() * xVar[orderIndex];

      /// create a constraint of type covConstr with index 'orderIndex' and assigned to it a default dual value
      /// *****************************************************************************************************
      covConstr(orderIndex).dualVal(data.orderPtrVector()[orderIndex]->width() / data.stockSheet().width());

      /// add coefficient of variable  xVar[orderIndex] in the constraint covConstr[orderIndex]
      /// *************************************************************************************
      covConstr[orderIndex] += xVar[orderIndex];
 
      /// set the rhs of constraint covConstr[orderIndex]
      /// ***********************************************
      covConstr[orderIndex] == data.orderPtrVector()[orderIndex]->demand();
    }

  /// Create an instanciation of a variable of class yVar with index '0' and set it in the objective with coefficient 1
  /// *****************************************************************************************************************
  objective += yVar[0];


  /// Define a cut template 
  /// *********************
  BcCutConstrArray cardinalityConstr(master, "CARD", 'F', 3.0, 1.0);

  /// Define a separation routine and associate it with the cut template
  /// ******************************************************************
  CardinalityCutSeparationRoutine * separationRoutinePtr = new CardinalityCutSeparationRoutine(&data);
  cardinalityConstr.attach(separationRoutinePtr);

}
