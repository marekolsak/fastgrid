/*
    AutoGrid

    Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
    All Rights Reserved.
    Copyright (C) 2008-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.

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

#include "parameterlibrary.h"
#include "utils.h"
#include "exceptions.h"
#include <cstring>
#include <cctype>

// Uh, this is ugly, but... we want to make the param_string variable (declared
// in default_parameters.h) const without changing the autodock source code
const
#include "../autodock-4.0.1/default_parameters.h"

ParameterLibrary::ParameterLibrary(LogFile *logFile, int debug, int outputLevel): logFile(logFile), debug(debug), outputLevel(outputLevel)
{
    memset(dictionary, 0, sizeof(dictionary));

    // Setup default parameters
    // These are set up in "default_parameters.h"
    // and stored in the param_string[MAX_LINES] array
    logFile->print("Setting up parameter library with factory defaults.\n\n\n");
    for (int i = 0; param_string[i]; i++)
        readLine(param_string[i]);
}

ParameterLibrary::~ParameterLibrary()
{
    for (int i = 0; i < MAXKEY; i++)
        if (dictionary[i])
            delete dictionary[i];
}

unsigned int ParameterLibrary::hash(const char *key)
{
    switch (strlen(key))
    {
    case 0:
        return 0;
    case 1:
        return key[0];
    default:
        return (unsigned int)key[0] + 256*(unsigned int)key[1];
    }
}

void ParameterLibrary::insertAtomParameter(const char *key, const ParameterEntry &value)
{
    unsigned int hashKey = hash(key);
    if (!dictionary[hashKey])
        dictionary[hashKey] = new ParameterEntry();
    *dictionary[hashKey] = value;
}

ParameterEntry *ParameterLibrary::findAtomParameter(const char *key)
{
    return dictionary[hash(key)];
}

const ParameterEntry *ParameterLibrary::findAtomParameter(const char *key) const
{
    return dictionary[hash(key)];
}

int ParameterLibrary::getAtomParameterRecIndex(const char *key) const
{
    const ParameterEntry *foundParam = findAtomParameter(key);
    if (foundParam)
        return foundParam->recIndex;
    return -1;
}

void ParameterLibrary::load(const char *filename)
{
    logFile->print("Using read_parameter_library\n");

    // Open and read the parameter library
    FILE *parameterLibraryFile;
    if ((parameterLibraryFile = boincOpenFile(filename, "r")) == 0)
    {
         fprintf(stderr,"Sorry, I can't find or open %s\n", filename);
         throw ExitProgram(-1);
    }

    char line[MAX_CHARS];
    while (fgets(line, sizeof(line), parameterLibraryFile) != 0)
        readLine(line);
}

void ParameterLibrary::readLine(const char *line)
{
    ParameterEntry thisParameter;
    int int_hbond_type = 0;
    int nfields;
    ParserTokens paramKeyword = parseParamLine(line);

    if (debug > 0)
        logFile->printFormatted("DEBUG: line = %sDEBUG: paramKeyword          = %d\n", line, paramKeyword);

    switch (paramKeyword)
    {
    case PAR_:
    case PAR_NULL:
    case PAR_COMMENT:
        break;

    case PAR_VDW:
        nfields = sscanf(line, "%*s %lf", &this->coeff_vdW);
        if (nfields < 1)
        {
            logFile->printTitled("WARNING:  Please supply a coefficient as a floating point number.\n\n");
            return; // skip any line without enough info
        }
        logFile->printFormatted("Free energy coefficient for the van der Waals term = \t%.4lf\n\n", this->coeff_vdW);
        break;

    case PAR_HBOND:
        nfields = sscanf(line, "%*s %lf", &this->coeff_hbond);
        if (nfields < 1)
        {
            logFile->printTitled("WARNING:  Please supply a coefficient as a floating point number.\n\n");
            return; // skip any line without enough info
        }
        logFile->printFormatted("Free energy coefficient for the H-bonding term     = \t%.4lf\n\n", this->coeff_hbond);
        break;

    case PAR_ESTAT:
        nfields = sscanf(line, "%*s %lf", &this->coeff_estat);
        if (nfields < 1)
        {
            logFile->printTitled("WARNING:  Please supply a coefficient as a floating point number.\n\n");
            return; // skip any line without enough info
        }
        logFile->printFormatted("Free energy coefficient for the electrostatic term = \t%.4lf\n\n", this->coeff_estat);
        break;

    case PAR_DESOLV:
        nfields = sscanf(line, "%*s %lf", &this->coeff_desolv);
        if (nfields < 1)
        {
            logFile->printTitled("WARNING:  Please supply a coefficient as a floating point number.\n\n");
            return; // skip any line without enough info
        }
        logFile->printFormatted("Free energy coefficient for the desolvation term   = \t%.4lf\n\n", this->coeff_desolv);
        break;

    case PAR_TORS:
        nfields = sscanf(line, "%*s %lf", &this->coeff_tors);
        if (nfields < 1)
        {
            logFile->printTitled("WARNING:  Please supply a coefficient as a floating point number.\n\n");
            return; // skip any line without enough info
        }
        logFile->printFormatted("Free energy coefficient for the torsional term     = \t%.4lf\n\n", this->coeff_tors);
        break;

    case PAR_ATOM_PAR:
        // Read in one line of atom parameters;
        // NB: scanf doesn't try to write missing fields
        nfields = sscanf(line, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                            thisParameter.autogridType,
                            &thisParameter.Rij,
                            &thisParameter.epsij,
                            &thisParameter.vol,
                            &thisParameter.solpar,
                            &thisParameter.RijHB,
                            &thisParameter.epsijHB,
                            &int_hbond_type,
                            &thisParameter.recIndex,
                            &thisParameter.mapIndex,
                            &thisParameter.bondIndex);
        if (nfields < 2)
            return; // skip any line without enough info

        if (int_hbond_type == 0)
            thisParameter.hbond = NON;
        else if (int_hbond_type == 1)
            thisParameter.hbond = DS;
        else if (int_hbond_type == 2)
            thisParameter.hbond = D1;
        else if (int_hbond_type == 3)
            thisParameter.hbond = AS;
        else if (int_hbond_type == 4)
            thisParameter.hbond = A1;
        else if (int_hbond_type == 5)
            thisParameter.hbond = A2;
        else
            thisParameter.hbond = NON;

        thisParameter.epsij   *= this->coeff_vdW;
        thisParameter.epsijHB *= this->coeff_hbond;

        insertAtomParameter(thisParameter.autogridType, thisParameter);
        logFile->printFormatted("Parameters for the atom type named \"%s\" were read in from the parameter library as follows:\n", thisParameter.autogridType);

        if (outputLevel > 2)
            logFile->printFormatted("\tR-eqm = %5.2f Angstrom\n\tweighted epsilon = %5.3f\n\tAtomic fragmental volume = %5.3f\n\tAtomic solvation parameter = %5.3f\n\tH-bonding R-eqm = %5.3f\n\tweighted H-bonding epsilon = %5.3f\n\tH-bonding type = %d,  bond index = %d\n\n",
                    thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                    thisParameter.RijHB, thisParameter.epsijHB, thisParameter.hbond, thisParameter.bondIndex);
        else
            logFile->printFormatted("\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,  At.frag.vol. = %.3f,  At.solv.par. = %.3f, \n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,  Hb type = %d,  bond index = %d\n\n",
                    thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                    thisParameter.RijHB, thisParameter.epsijHB, thisParameter.hbond, thisParameter.bondIndex);
        break;
    }
}

ParameterLibrary::ParserTokens ParameterLibrary::parseParamLine(const char *line)
{
    int j, i;
    ParserTokens token = PAR_; // return -1 if nothing is recognized.
    char c[LINE_LEN];

    // tokentablesize should be set to the length of the tokentable
    const int tokentablesize = 6;

    const struct
    {
       const char *lexeme;
       ParserTokens tokenvalue;
    } tokentable[] = {
        {"FE_coeff_vdW", PAR_VDW},          // 1
        {"FE_coeff_hbond", PAR_HBOND},      // 2
        {"FE_coeff_estat", PAR_ESTAT},      // 3
        {"FE_coeff_desolv", PAR_DESOLV},    // 4
        {"FE_coeff_tors", PAR_TORS},        // 5
        {"atom_par", PAR_ATOM_PAR}          // 6
    }; // 6 tokens  // remember to set tokentablesize earlier

    c[0] = '\0';
    for (j = 0; line[j]!='\0' && line[j]!=' ' && line[j]!='\t' && line[j]!='\n'; j++)
    {
        //  Ignore case
        c[j] = (char)tolower((int)line[j]);
        if (debug > 0)
            logFile->printFormatted("%c",c[j]);
    }
    if (debug > 0)
        logFile->printFormatted("\nj = %d\n",j);

    //  Recognize one character tokens
    if ((c[0]=='\n') || (c[0]=='\0'))
        token = PAR_NULL;
    else if (c[0]=='#')
        token = PAR_COMMENT;

    //  Recognize token strings
    for (i = 0; i < tokentablesize && token == PAR_; i++)
    {
        if (debug > 0)
            logFile->printFormatted("i = %d, tokentable[i].lexeme = %s, tokentable[i].value = %d, c = %s\n",i,tokentable[i].lexeme,tokentable[i].tokenvalue,c);
        if (strncasecmp(tokentable[i].lexeme, c, j) == 0)
            token = tokentable[i].tokenvalue;
    }
    return token;
}
