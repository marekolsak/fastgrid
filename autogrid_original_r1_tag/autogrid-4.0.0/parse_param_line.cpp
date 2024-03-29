/*

 $Id: parse_param_line.cpp,v 1.2 2007/05/03 20:46:06 garrett Exp $

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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "partokens.h"
#include "parse_param_line.h"

extern int debug;
extern FILE *logFile;

int parse_param_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_param_line                                                */
/*  Function: Parse the docking parameter file line                           */
/* Copyright: (C) 2005, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 08/03/05                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 08/03/05 GMM     Entered code.                                             */
/******************************************************************************/

{
    int j, i, token = PAR_;               /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    // tokentablesize should be set to the length of the tokentable
    // 
    const int tokentablesize = 6;

    const struct {
       char *lexeme;
       int tokenvalue;
    } tokentable[] = {{"FE_coeff_vdW", PAR_VDW}, // 1
                      {"FE_coeff_hbond", PAR_HBOND}, // 2
                      {"FE_coeff_estat", PAR_ESTAT}, // 3
                      {"FE_coeff_desolv", PAR_DESOLV}, // 4
                      {"FE_coeff_tors", PAR_TORS}, // 5
                      {"atom_par", PAR_ATOM_PAR} // 6
              }; // 6 tokens  // remember to set tokentablesize earlier

    c[0] = '\0';
    for (j=0; ((line[j]!='\0')&&(line[j]!=' ')&&(line[j]!='\t')&&(line[j]!='\n')); j++) {
        /*  Ignore case */
        c[j] = (char)tolower((int)line[j]);
        if (debug > 0) {
            (void)fprintf(logFile,"%c",c[j]);
        }
    }
    if (debug > 0) {
        (void)fprintf(logFile,"\nj = %d\n",j);
    }

    /*  Recognize one character tokens  */

    if ((c[0]=='\n') || (c[0]=='\0')) {
        token = PAR_NULL;
    } else if (c[0]=='#') {
        token = PAR_COMMENT;
    }

    /*  Recognize token strings  */

    for (i=0;  (i < tokentablesize) && (token == PAR_);  i++) {
        if (debug > 0) {
            (void)fprintf(logFile,"i = %d, tokentable[i].lexeme = %s, tokentable[i].value = %d, c = %s\n",i,tokentable[i].lexeme,tokentable[i].tokenvalue,c);
        }
        if (strncasecmp(tokentable[i].lexeme, c, j) == 0) {
            token = tokentable[i].tokenvalue;
        }
    }
    return(token);
}
/* EOF */
