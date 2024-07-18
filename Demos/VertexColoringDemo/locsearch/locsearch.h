#ifndef __MWIS_H
#define __MWIS_H
/**
    This file is part of exactcolors.

    exactcolors is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exactcolors is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ls_graph.h"
#include "ls_color_defs.h"

typedef struct _MWISls_env   MWISls_env;


extern
int COLORstable_LS(MWISls_env** env,
                   COLORset** newsets, int* nnewsets, int ncount,
                   int ecount, const int elist[],
                   const COLORNWT nweights[],COLORNWT cutoff);

int COLORstable_init_LS(MWISls_env** env,
                        int ncount,
                        int ecount, const int elist[],
                        const COLORNWT nweights[], COLORNWT cutoff);

extern int COLORstable_free_ls_env(MWISls_env** env);

#endif
