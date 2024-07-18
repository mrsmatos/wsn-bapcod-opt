/**
 *
 * This file bcNonPublicCglCut.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifdef _CGL_FOUND

#define CGLCutsPrintLevel 1

#include "bcModelNonPublicCuts.hpp"
#include <chrono>
#include <vector>
#include <list>

#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglSimpleRounding.hpp"
//#include "CglPreProcess.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglOddHole.hpp"
#include "CglKnapsackCover.hpp"
//#include "CglGomory.hpp"
//#include "CglProbing.hpp"
//#include "CglRedSplit.hpp"
//#include "CglRedSplit2.hpp"
//#include "CglGMI.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglTwomir.hpp"
//#include "CglDuplicateRow.hpp"
//#include "CglStored.hpp"
//#include "CglLandP.hpp"
#include "CglResidualCapacity.hpp"
#include "CglZeroHalf.hpp"
//#include "CglLiftAndProject.hpp"

#include "bcInstanciatedVarConstrC.hpp"
#include "bcMasterConfC.hpp"
#include "bcSpVarConstrC.hpp"
#include "bcColGenSpConfC.hpp"

CGLSeparationFunctor::CGLSeparationFunctor() :
   last_cut_index(0), enabledCuts(BcCutType::_BcCutTypeMax, false)
{
    enableCutAll();
}

struct noderange {
    int startnode;
    int endnode;
};
typedef struct noderange NODERANGE;

int filterCuts(OsiSolverInterface & osi, OsiCuts & cs, double effectivenessLb, std::vector<std::string> & cutsName)
{
    int nbcuts(0);

    /// Loop once for each row cut
    /// here we just have clear the cutsName for bad cuts
    for (int cutId = 0; cutId < cs.sizeRowCuts(); cutId++)
    {
        if (cs.rowCut(cutId).effectiveness() < effectivenessLb)
        {
            cutsName[cutId].clear();
            continue;
        }
        if (!cs.rowCut(cutId).consistent())
        {
            cutsName[cutId].clear();
            continue;
        }
        if (!cs.rowCut(cutId).consistent(osi))
        {
            cutsName[cutId].clear();
            continue;
        }
        if (cs.rowCut(cutId).infeasible(osi))
        {
            cutsName[cutId].clear();
            continue;
        }
        nbcuts++;
        //applyRowCut(cs.rowCut(cutId));
    }

    return nbcuts;
}

void printCuts(OsiCuts & cuts, int start_row_idx)
{
    for (int cutId = start_row_idx; cutId < cuts.sizeRowCuts(); cutId++)
    {
        const OsiRowCut * cut = cuts.rowCutPtr(cutId);
        cut->print();
    }
}

void CGLSeparationFunctor::enableCut(BcCutType::Cuts cutType)
{
    enabledCuts[cutType] = true;
}

void CGLSeparationFunctor::disableCut(BcCutType::Cuts cutType)
{
    enabledCuts[cutType] = false;
}

void CGLSeparationFunctor::enableCutAll()
{
    for (int cutType = 0; cutType < BcCutType::_BcCutTypeMax; cutType++)
        enabledCuts[cutType] = true;
}

void CGLSeparationFunctor::disableCutAll()
{
    for (int cutType = 0; cutType < BcCutType::_BcCutTypeMax; cutType++)
        enabledCuts[cutType] = false;
}

int tryCutGeneration(CglCutGenerator & cgl, OsiClpSolverInterface & osi, OsiCuts & cuts,
                     std::vector<std::string> & cutsName)
{
    int old_nb_cuts = cuts.sizeCuts();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    cgl.generateCuts(osi, cuts);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    int micro = std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count();

    double timeInSeconds = (double)micro * 1e-6;

    int new_nb_cuts = cuts.sizeCuts();
    int nb_cuts = (new_nb_cuts - old_nb_cuts);
    if (nb_cuts > 0)
    {
        std::cout << "bcCglCut: " << " generate " << nb_cuts << " cut(s)" << " in " << timeInSeconds << " from "
                  << typeid (cgl).name() << std::endl;

        for (int i = 0; i < nb_cuts; i++)
            cutsName.push_back(typeid(cgl).name());

        if (printL(2))
            printCuts(cuts, old_nb_cuts);

    } else {
        std::cout << "bcCglCut: " << " generate no cut" << " in " << timeInSeconds << " from " << typeid (cgl).name()
                  << std::endl;
    }

    return nb_cuts;
}

void printOsiProb(const OsiClpSolverInterface & osi) {
    int ncol = osi.getNumCols();
    int row_b_size = osi.getNumRows();
    const char *rowSense = osi.getRowSense();
    const double *rowLower = osi.getRowLower();
    const double *rowUpper = osi.getRowUpper();
    const double *colLower = osi.getColLower();
    const double *colUpper = osi.getColUpper();
    const double *coef = osi.getObjCoefficients();
    const double *colSol = osi.getColSolution();

    const char *colType = osi.getColType();
    char convType[] = {'C', 'B', 'I', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?'};

    std::cout << "osi row lsu \n";
    for (int i = 0; i < row_b_size; i++)
        std::cout << "(" << rowLower[i] << "," << rowSense[i] << "," << rowUpper[i] << ") ";
    std::cout << std::endl;

    std::cout << "osi col type \n";
    for (int i = 0; i < ncol; i++)
        std::cout << "(" << colLower[i] << "," << convType[colType[i]] << "," << colUpper[i] << ") ";
    std::cout << std::endl;

    std::cout << "osi obj \n";
    for (int i = 0; i < ncol; i++) {
        std::cout << coef[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "osi sol \n";
    for (int i = 0; i < ncol; i++) {
        std::cout << colSol[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "osi matrix \n";

    const CoinPackedMatrix *rmat = osi.getMatrixByCol();

    for (int y = 0; y < rmat->getNumRows(); y++) {
        for (int x = 0; x < rmat->getNumCols(); x++) {
            rmat->printMatrixElement(y, x);
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}

void nicePrintOsiProb(const OsiClpSolverInterface & osi, const std::vector<InstanciatedVar *> & instVarPtr)
{
    int numOsiVars = osi.getNumCols();
    const char* rowSense = osi.getRowSense();
    const double* rowLower = osi.getRowLower();
    const double* rowUpper = osi.getRowUpper();
    const double* colLower = osi.getColLower();
    const double* colUpper = osi.getColUpper();
    const double* coef = osi.getObjCoefficients();
    const double* colSol = osi.getColSolution();

    const char* colType = osi.getColType();
    char convType[] = {'C', 'B', 'I', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?'};

    std::cout << "Vars in Osi :" << std::endl;
    for (int varId = 0; varId < numOsiVars; varId++)
    {
        std::string varName = (varId < (int)instVarPtr.size()) ? instVarPtr[varId]->name() : "NN";
        std::cout << " " << varName << ": ref=" << varId << ", bounds=[" << colLower[varId] << ","
                  << colUpper[varId]  << "]" << ", type=" << convType[colType[varId]] << ", obj=" << coef[varId]
                  << ", sol=" << colSol[varId] << std::endl;
    }

    std::cout << "Constrs in Osi :" << std::endl;

    const CoinPackedMatrix * rmat = osi.getMatrixByRow();

    const double * elements = rmat->getElements();
    const int * indices = rmat->getIndices();
    int elemId = 0;
    for (int rowId = 0; rowId < rmat->getNumRows(); rowId++)
    {
        int rowSize = rmat->getVectorSize(rowId);
        if (rowSize > 0) {
            std::cout << "row " << rowId << " :";
            std::cout << rowLower[rowId] << " <= ";
            for (int rowElemId = 0; rowElemId < rowSize; ++rowElemId)
            {
                int varId = indices[elemId];
                double coeff = elements[elemId];
                std::string varName = (varId < (int)instVarPtr.size()) ? instVarPtr[varId]->name() : "NN";

                std::cout << ((coeff >= 0) ? "+" : "-") << abs(coeff) << "*" << varName;
                elemId += 1;
            }
            std::cout <<  " <= " << rowUpper[rowId] << " (t:" << rowSense[rowId] << ")";
            std::cout << std::endl;
        }
    }
}

void smallTest()
{
    // test coupe knapsack
    OsiClpSolverInterface kosi;

    const double zeroes[7] = {0, 0, 0, 0, 0, 0, 0};
    const double ones[7] = {1, 1, 1, 1, 1, 1, 1};
    const double kobj[7] = {100, 100, 20, 5, 2, 1, 1};
    const double row[7] = {11, 6, 6, 6, 5, 6, 1};
    //  const double   ksol[7]={1,1,0,0,0,0.333,0};
    //  const double   ksol[7]={1,0.333,1,0,0,0,0};
    //const double   ksol[7]={1,1,0.3333333,0,0,0,0};
    //  const double   ksol[7]={1,1,0.333,0.000333,0,0,0};
    const double ksol[7] = {1, 1, 0, 0, 0, 0.333333, 0};
    const double rhs[1] = {19};

    std::cout << "Create kmatrix..." << std::endl;

    CoinPackedMatrix kmatrix(false, 0, 0);
    //kmatrix->setDimensions(1, 7); // rows x columns   (indices max)
    //kmatrix->setDimensions(0, 6); // rows x columns   (indices max)

    std::cout << "Create vector..." << std::endl;

    CoinPackedVector kcol;
    for (int i = 0; i < 7; i++)
        kcol.insert(i, row[i]);

    std::cout << "Append Kcol..." << std::endl;
    kmatrix.appendRow(kcol);

    std::cout << "Loading kmatrix..." << std::endl;

    kosi.loadProblem(kmatrix,&zeroes[0], &ones[0],&kobj[0],NULL, &rhs[0]); // row l/u bound

    for (int i = 0; i < 7; i++)
        kosi.setInteger(i);

    kosi.setObjSense(-1);
    kosi.setColSolution(&ksol[0]);

    printOsiProb(kosi);


    std::cout << " Coucou Issam : " << std::setprecision(9) << std::endl << std::endl;
    std::cout << " Cgl Osi Tolerance " << kosi.getIntegerTolerance() << std::endl;

    CglKnapsackCover kcgk;
    kosi.setColSolution(&ksol[0]);

    OsiCuts kcuts;
    kcgk.generateCuts(kosi, kcuts);

    OsiSolverInterface::ApplyCutsReturnCode kacRc = kosi.applyCuts(kcuts, 0.0);

    // Print applyCuts return code
    std::cout << " Knap Test : " << std::endl << std::endl;
    std::cout << kcuts.sizeCuts() << " cuts were generated" << std::endl;
    std::cout << "  " << kacRc.getNumInconsistent() << " were inconsistent" << std::endl;
    std::cout << "  " << kacRc.getNumInconsistentWrtIntegerModel()
         << " were inconsistent for this problem" << std::endl;
    std::cout << "  " << kacRc.getNumInfeasible() << " were infeasible" << std::endl;
    std::cout << "  " << kacRc.getNumIneffective() << " were ineffective" << std::endl;
    std::cout << "  " << kacRc.getNumApplied() << " were applied" << std::endl;
    std::cout << std::endl << std::endl;

    kcuts.printCuts();

    printOsiProb(kosi);
}


void addVariables(ProbConfig * probConfPtr,
                  std::vector<double> & varLb,
                  std::vector<double> & varUb,
                  std::vector<double> & varObjCoeff,
                  std::vector<char> & varType,
                  std::vector<InstanciatedVar *> & instVarPtr,
                  std::vector<int> & varRefToVarId,
                  int & numOsiVars)
{
    VarIndexManager & varSet = probConfPtr->probPtr()->probVarSet();
    for (VarIndexManager::iterator varIt = varSet.begin(VcIndexStatus::Active, 's');
         varIt != varSet.end(VcIndexStatus::Active, 's'); varIt++)
    {
        if ((*varIt)->isTypeOf(VcId::InstanciatedVarMask))
        {
            instVarPtr.push_back(static_cast<InstanciatedVar *>(*varIt));
            varUb.push_back((*varIt)->globalUb());
            varLb.push_back((*varIt)->globalLb());
            varObjCoeff.push_back((*varIt)->costrhs());
            varType.push_back((*varIt)->type());
            varRefToVarId.resize((*varIt)->ref() + 1, -1);
            varRefToVarId[(*varIt)->ref()] = numOsiVars;
            numOsiVars += 1;
        }
    }
}

void addConstraint(Constraint * constrPtr,
                   const std::vector<int> & varRefToVarId,
                   std::vector<double> & constrLb,
                   std::vector<double> & constrUb,
                   CoinPackedMatrix & matrix,
                   int & numOsiConstrs)
{

}

void addConstraints(ProbConfig * probConfPtr,
                    const std::vector<int> & varRefToVarId,
                    std::vector<double> & constrLb,
                    std::vector<double> & constrUb,
                    CoinPackedMatrix & matrix,
                    int & numOsiConstrs)
{
    int maxVarRef = (int)varRefToVarId.size() - 1;
    ConstrIndexManager & constrSet = probConfPtr->probPtr()->probConstrSet();
    for (ConstrIndexManager::iterator constrPtrIt = constrSet.begin(VcIndexStatus::Active, 's');
         constrPtrIt != constrSet.end(VcIndexStatus::Active, 's'); constrPtrIt++)
    {
        if (!(*constrPtrIt)->isTypeOf(VcId::Base4NonLinearConstraintMask) && (*constrPtrIt)->inCurForm())
        {
            CoinPackedVector row;
            for (ConstVarConstrPtr2Double::const_iterator membIt = (*constrPtrIt)->member2coefMap().begin();
                 membIt != (*constrPtrIt)->member2coefMap().end(); membIt++)
            {
                int varFlag = membIt->first->flag();
                if (varFlag == 's')
                {
                    int varRef = membIt->first->ref();
                    if ((varRef <= maxVarRef) && (varRefToVarId[varRef] >= 0))
                        row.insert(varRefToVarId[varRef], (double) membIt->second);
                }
            }
            if ((*constrPtrIt)->isTypeOf(VcId::InstMasterConstrMask))
            {
                InstMasterConstr * iMastConstrPtr = static_cast<InstMasterConstr * >(*constrPtrIt);
                MapSubProbVariablePtr2Double::const_iterator membIt;
                for (membIt = iMastConstrPtr->subProbVarMember2coefMap().begin();
                     membIt != iMastConstrPtr->subProbVarMember2coefMap().end(); membIt++)
                {
                    int varFlag = membIt->first->flag();
                    if (varFlag == 's')
                    {
                        int varRef = membIt->first->ref();
                        if ((varRef <= maxVarRef) && (varRefToVarId[varRef] >= 0))
                            row.insert(varRefToVarId[varRef], (double) membIt->second);
                    }
                }
            }
            if (row.getNumElements() > 0)
            {
                switch ((*constrPtrIt)->sense()) {
                    case 'L':
                    {
                        constrLb.push_back((*constrPtrIt)->lb());
                        constrUb.push_back((*constrPtrIt)->rhs());
                        break;
                    }
                    case 'G':
                    {
                        constrLb.push_back((*constrPtrIt)->rhs());
                        constrUb.push_back((*constrPtrIt)->ub());
                        break;
                    }
                    case 'E':
                    {
                        constrLb.push_back((*constrPtrIt)->rhs());
                        constrUb.push_back((*constrPtrIt)->rhs());
                        break;
                    }
                    default:
                    {
                        constrLb.push_back((*constrPtrIt)->lb());
                        constrUb.push_back((*constrPtrIt)->ub());
                        if (printL(2))
                            std::cout << " unknown constraint sense ERROR : " << (*constrPtrIt)->sense() << std::endl;
                    }
                }
                matrix.appendRow(row);
                numOsiConstrs += 1;
            }
        }
    }
}

int CGLSeparationFunctor::operator()(BcFormulation formPtr,
                                     BcSolution & primalSol,
                                     double & maxViolation,
                                     std::list< BcConstr > & cutList)
{
    OsiClpSolverInterface osi;

    std::vector<double> varUb, varLb, varObjCoeff;
    std::vector<char> varType;
    std::vector<InstanciatedVar *> instVarPtr;
    std::vector<int> varRefToVarId;

    ProbConfig * probConfPtr = formPtr.probConfPtr();
    int numOsiVars = 0;
    addVariables(probConfPtr, varLb, varUb, varObjCoeff, varType, instVarPtr,varRefToVarId,
                 numOsiVars);

    for (std::vector< ColGenSpConf * >::const_iterator spConfPtrIt = probConfPtr->colGenSubProbConfPts().begin();
         spConfPtrIt != probConfPtr->colGenSubProbConfPts().end(); spConfPtrIt++)
        addVariables(*spConfPtrIt, varLb, varUb, varObjCoeff, varType, instVarPtr,varRefToVarId,
                     numOsiVars);

    int numOsiConstrs = 0;
    std::vector<double> constrLb, constrUb;
    CoinPackedMatrix matrix(false, 0, 0);

    addConstraints(probConfPtr, varRefToVarId, constrLb, constrUb, matrix, numOsiConstrs);

    for (std::vector< ColGenSpConf * >::const_iterator spConfPtrIt = probConfPtr->colGenSubProbConfPts().begin();
         spConfPtrIt != probConfPtr->colGenSubProbConfPts().end(); spConfPtrIt++)
    {
         if (*(*spConfPtrIt)->upperBoundPtr() <= 1)
             addConstraints(*spConfPtrIt, varRefToVarId, constrLb, constrUb, matrix, numOsiConstrs);
    }

    /// set solution variable vector to Zero
    std::vector<double> sol(numOsiVars, 0.0);

    /// copy non zero solution variables into solution vector
    std::set< BcVar > varList;
    primalSol.extractVar(varList);
    for (std::set< BcVar >::const_iterator varIt = varList.begin(); varIt != varList.end(); ++varIt)
    {
        int varRef = ((InstanciatedVar *)(*varIt))->ref();
        int varFlag = ((InstanciatedVar *)(*varIt))->flag();
        if ((varRef < (int)varRefToVarId.size()) && (varFlag == 's') && (varRefToVarId[varRef] >= 0))
        {
            sol[varRefToVarId[varRef]] = ((InstanciatedVar *)(*varIt))->solVal();
        }
    }

    if (printL(2))
        std::cout << " Loading matrix..." << std::endl;

    osi.loadProblem(matrix, &varLb[0], &varUb[0], &varObjCoeff[0],&constrLb[0],&constrUb[0]);
    osi.setColSolution(&sol[0]);
    osi.setObjSense(1); /// minimization objective

    /// set variable type
    for (int varId = 0; varId < numOsiVars; varId++) {
        switch (varType[varId]) {
            case 'C':
            {
                osi.setContinuous(varId);
                break;
            }
            case 'I':
            {
                osi.setInteger(varId);
                break;
            }
            case 'B':
            {
                osi.setInteger(varId);
                break;
            }
            default:
            {
                std::cout << " unknown variable type ERROR : " << varType[varId] << std::endl;
            }
        }
    }

    if (printL(CGLCutsPrintLevel))
        nicePrintOsiProb(osi, instVarPtr);

    OsiCuts cuts;

    osi.setWarmStart(0); // to prevent Gomory cut generator to segfault...

    //CglPreProcess cgpp;
    // // cgl.setTolerance(0.0000001);
    //  cgl.generateCuts(*osi,cuts);
    // stolen from ClpSOlver

//    CglProbing probingGen;
//    probingGen.setUsingObjective(1);
//    probingGen.setMaxPass(1);
//    probingGen.setMaxPassRoot(1);
//    // Number of unsatisfied variables to look at
//    probingGen.setMaxProbe(10);
//    probingGen.setMaxProbeRoot(50);
//    // How far to follow the consequences
//    probingGen.setMaxLook(10);
//    probingGen.setMaxLookRoot(50);
//    probingGen.setMaxLookRoot(10);
//    // Only look at rows with fewer than this number of elements
//    probingGen.setMaxElements(200);
//    probingGen.setMaxElementsRoot(300);
//    probingGen.setRowCuts(3);
//    // set default action (0=off,1=on,2=root)
//    int probingAction = 1;
    // probingGen.generateCuts(*osi,cuts);

    std::vector<std::string> cutsName;

    if (enabledCuts[BcCutType::cglClique])
    {
        CglClique cgc;
        cgc.setRowCliqueReport(false);
        cgc.setStarCliqueReport(false);
        tryCutGeneration(cgc, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglKnapsackCoverCut])
    {
        CglKnapsackCover cgk;
        cgk.switchOnExpensive();
        tryCutGeneration(cgk, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglOddHoleCut])
    {
        CglOddHole cgoh;
        tryCutGeneration(cgoh, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglZeroHalf])
    {
        CglZeroHalf cgzh;
        tryCutGeneration(cgzh, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglFlowCover])
    {
        CglFlowCover cgfc;
        tryCutGeneration(cgfc, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglMixedIntegerRoundingCut])
    {
        CglMixedIntegerRounding cgmir;
        tryCutGeneration(cgmir, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglMixedIntegerRounding2Cut])
    {
        CglMixedIntegerRounding2 cgmir2(1, true, 1);
        cgmir2.setDoPreproc(1);
        tryCutGeneration(cgmir2, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglTwoMir])
    {
        CglTwomir cgtm;
        cgtm.setCutTypes(true, true, false, true);
        tryCutGeneration(cgtm, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglResidualCapacity])
    {
        CglResidualCapacity cgrc;
        tryCutGeneration(cgrc, osi, cuts, cutsName);
    }

    if (enabledCuts[BcCutType::cglSimpleRoundingCut])
    {
        CglSimpleRounding cgsr;
        tryCutGeneration(cgsr, osi, cuts, cutsName);
    }




    if (filterCuts(osi, cuts, 0, cutsName))
    {
        std::string maxViolationCutName;

        BcCutConstrArray CGLCutArray(formPtr, "CGL");

        maxViolation = 0;
        /// Loop once for each row cut
        for (int cutId = 0; cutId < cuts.sizeRowCuts(); cutId++)
        {
            if (cutsName[cutId].length() > 0)
            {
                OsiRowCut rc = cuts.rowCut(cutId);

                int constr_low(0), constr_high(0);
                double lb = rc.lb();
                double ub = rc.ub();

                /// this constraint is ok
                if (printL(2))
                    std::cout << "Cgl: osi cut dump " << cutId << " " << cutsName[cutId] << " " << rc.sense() << " ";

                switch (rc.sense())
                {
                    case 'E':
                    {
                        constr_low = last_cut_index++;
                        CGLCutArray(constr_low) == lb;
                        break;
                    }

                    case 'G':
                    {
                        constr_low = last_cut_index++;
                        CGLCutArray(constr_low) >= lb;
                        break;
                    }

                    case 'L':
                    {
                        constr_low = last_cut_index++;
                        CGLCutArray(constr_low) <= ub;
                        break;
                    }

                    case 'R':
                    {
                        constr_low = last_cut_index++;
                        CGLCutArray(constr_low) >= lb;
                        constr_high = last_cut_index++;
                        CGLCutArray(constr_high) <= ub;
                        break;
                    }
                    case 'N':
                    {
                        if (printL(-1))
                            std::cout << "BaPCod WARNING : Cgl : unmanaged case - inequality cut generated... "
                                      << "look at the code ! TODO\n";
                        std::cerr << "BaPCod WARNING : Cgl : unmanaged case - inequality cut generated... "
                                  << "look at the code ! TODO\n";
                    }
                }

                double mhs = 0.0; /// evaluation of expression to mesure violation

                const double * coeffs = rc.row().getVectorElements();
                const int * varIds = rc.row().getVectorIndices();

                for (int elemId = 0; elemId < rc.row().getVectorNumElements(); elemId++)
                {
                    if (coeffs[elemId] != 0)
                    {
                        if (printL(2))
                            std::cout << " " << coeffs[elemId] << "*@" << varIds[elemId];

                        BcVar myVar(instVarPtr[varIds[elemId]]);

                        if (myVar.isDefined())
                        {
                            if (constr_low)
                                CGLCutArray(constr_low) += coeffs[elemId] * myVar; // TODO slow extraction ?
                            if (constr_high)
                                CGLCutArray(constr_high) += coeffs[elemId] * myVar; // TODO slow extraction ?

                            if (printL(2) && (coeffs[elemId] > 1e-5))
                                std::cout << " " << myVar.name() << " * " << coeffs[elemId] << " + ";

                            mhs += coeffs[elemId] * sol[varIds[elemId]];
                        }
                        else if (printL(2))
                        {
                            std::cout << "Cgl: Possible ERROR refToIVarPtrMap[elemId] for elemId=" << varIds[elemId];
                        }
                    }
                }

                if (printL(2))
                    std::cout << std::endl;

                if (mhs > ub)
                {
                    double violation = abs(mhs - ub);
                    if (maxViolation < violation)
                    {
                        maxViolation = violation;
                        if (printL(2))
                            std::cout << "Cgl: new maxViolation (ub) : " << maxViolation << " from "
                                      << (maxViolationCutName = cutsName[cutId]);
                    }
                }

                if (mhs < lb) {
                    double violation = abs(mhs - lb);
                    if (maxViolation < violation)
                    {
                        maxViolation = violation;
                        if (printL(2))
                            std::cout << "Cgl: new maxVIolation (lb) : " << maxViolation << " from "
                                      << (maxViolationCutName = cutsName[cutId]);
                    }
                }

                if (printL(2))
                    std::cout << std::endl;

                if (constr_low)
                {
                    cutList.push_back(CGLCutArray(constr_low));
                    if (printL(CGLCutsPrintLevel))
                        CGLCutArray(constr_low).nicePrint();
                }
                if (constr_high)
                {
                    cutList.push_back(CGLCutArray(constr_high));
                    if (printL(CGLCutsPrintLevel))
                        CGLCutArray(constr_low).nicePrint();
                }
            }
        }

        if (printL(2))
            std::cout << "Cgl: CutSeparationRoutine violation = " << maxViolation << " from " << maxViolationCutName
                      << std::endl;

        return cutList.size();
    }

    if (printL(2))
        cuts.printCuts();

    return 0;
}

#endif /// _CGL_FOUND