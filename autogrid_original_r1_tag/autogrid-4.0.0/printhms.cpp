/*

 $Id: printhms.cpp,v 1.6 2007/05/03 20:46:06 garrett Exp $

 AutoGrid 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* printhms.cc */

    #include <stdio.h>
    #include "printhms.h"

extern FILE *logFile;

void printhms( Real t )

{
    int   h,
          m;
    Real T, s;
    Real min = 60.,
	  hrs = 3600.;

    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,       "%.2fs",       s );
        else
            fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }
}
/* EOF */
