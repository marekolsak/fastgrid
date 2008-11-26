/*

 $Id: processProgramParameters.cpp,v 1.8 2007/05/03 20:46:06 garrett Exp $

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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "programparameters.h"
#include "utils.h"
#include "exceptions.h"

ProgramParameters::ProgramParameters(int argc, char **argv): debug(0)
{
    programName[0] = 0;
    gridParameterFilename[0] = 0;
    logFilename[0] = 0;

    parse(argc, argv);
}

const char *ProgramParameters::getProgramName() const
{
    return programName;
}

const char *ProgramParameters::getGridParameterFilename() const
{
    return gridParameterFilename;
}

const char *ProgramParameters::getLogFilename() const
{
    return logFilename;
}

int ProgramParameters::getDebugLevel() const
{
    return debug;
}

int ProgramParameters::parse(int argc, char **argv)
{
    strncpy(programName, argv[0], MAX_CHARS);

    // Loop over arguments
    int argindex = 1;
    while((argc > 1) && (argv[1][0] == '-'))
    {
        switch(argv[1][1])
        {
        case 'd':
            debug++;
            break;
        case 'u':
            fprintf(stderr, "usage: %s -p parameter_filename\n"
                            "                -l log_filename\n"
                            "                -d (increment debug level)\n"
                            "                -u (display this message)\n\n", argv[0]);
            throw ExitProgram(0);
            break;
        case 'l':
            strncpy(logFilename, argv[2], MAX_CHARS);
            argv++;
            argc--;
            argindex++;
            break;
        case 'p':
            strncpy(gridParameterFilename, argv[2], MAX_CHARS);
            argv++;
            argc--;
            argindex++;
            break;
        default:
            fprintf(stderr,"%s: unknown switch -%c\n",argv[0],argv[1][1]);
            throw ExitProgram(1);
            break;
        }
        argindex++;
        argc--;
        argv++;
    }
    return argindex;
}
