#ifndef __COLOR_H
#define __COLOR_H
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

#include "ls_color_defs.h"


/* find a greedy coloring.*/
extern int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
                 COLORset **colorclasses);
/* find a greedy coloring using the DSATUR algorithm.*/
extern int COLORdsatur(int ncount, int ecount, int *elist, int *ncolors,
                COLORset **colorclasses);


#endif  /* __COLOR_H */
