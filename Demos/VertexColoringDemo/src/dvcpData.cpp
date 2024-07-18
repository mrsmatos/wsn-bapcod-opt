/**
 *
 * This file dvcpData.cpp is a part of BaPCod - a generic Branch-And-Price Code.
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

void DvcpData::readData(const std::string & inputFileName)
{
    std::ifstream inFile(inputFileName.c_str(), std::ios::in);
    if (!inFile)
    {
        std::cout << "cannot find input file: " << inputFileName << std::endl;
    }
    if (inFile.eof())
    {
        std::cout << "empty input file" << std::endl;
    }

    char line[256];
    std::string str = "";
    int e=0;
    while (!inFile.eof())
    {
        inFile.getline(line, 256);
        std::string stLine(line);
        std::istringstream lineStream(stLine);
        int tail, head;
        switch (line[0])
        {
            case 'p' :
                lineStream >> str;
                lineStream >> str;
                lineStream >> numV;
                lineStream >> numE;
                //std::cout << numV << " vertices and " << numE << " edges" << std::endl;
                break;
            case 'e' :
                lineStream >> str;
                lineStream >> tail;
                lineStream >> head;
                //E.insert(order);
                T_e[e] = tail - 1;
                H_e[e++] = head - 1;
                N_v[tail-1].insert(head-1);
                N_v[head-1].insert(tail-1);
                //std::cout << " edge " << tail << "--" << head << std::endl;
                break;
            default :
                /// this is comment, do nothing
                break;
        }
    }

    /// construct the sets of non-adjacent nodes
    for (int i = 0; i < numV; i++)
    {
        for (int j = 0; j < numV; j++)
            if (i != j)
                A_v[i].insert(j);
        for (Set::iterator setIt = N_v[i].begin(); setIt != N_v[i].end(); ++setIt)
            A_v[i].erase(*setIt);
    }
}
