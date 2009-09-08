/*
    FastGrid (formerly AutoGrid)

    Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
    Copyright (C) 2009 Masaryk University. All rights reserved.

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

#include <cstring>
#include <cctype>
#include "ParameterLibrary.h"
#include "Utils.h"
#include "Exceptions.h"

#define char const char // ugly hack to calm down G++
#include "../autodock/default_parameters.h"
#undef char

ParameterLibrary::ParameterLibrary(LogFile *logFile, const char *modelText, Unbound_Model unboundModel, int debug, int outputLevel)
    : logFile(logFile), debug(debug), outputLevel(outputLevel)
{
    filename[0] = 0;
    memset(dictionary, 0, sizeof(dictionary));

    // Setup default parameters
    // These are set up in "default_parameters.h"
    // and stored in the param_string_VERSION_NUM[MAX_LINES] array
    // so far we have param_string_4_0 and param_string_4_1
    const char **paramString;
    if (unboundModel == Extended)
        paramString = param_string_4_0;
    else if (unboundModel == Unbound_Same_As_Bound)
        paramString=param_string_4_1;
    else
    {
        logFile->printFormatted("DEBUG: cannot determine %s parameter values \n", modelText);
        throw ExitProgram(-1);
    }

    logFile->printFormatted("Setting up parameter library with AutoDock %s values.\n\n\n", modelText);
    for (int i = 0; paramString[i]; i++)
        readLine(paramString[i]);
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
    logFile->printFormatted("Using read_parameter_library() to try to open and read \"%s\".\n\n", filename);
    strncpy(this->filename, filename, MAX_CHARS);

    // Open and read the parameter library
    FILE *parameterLibraryFile;
    if ((parameterLibraryFile = boincOpenFile(filename, "r")) == 0)
    {
        logFile->printFormatted("Sorry, I can't find or open %s\n", filename);
        throw ExitProgram(-1);
    }

    char line[LINE_LEN];
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

    case PAR_UNBOUND:
        logFile->printError(WARNING, "the unbound model cannot be specified in the parameter library file.\n\n");
        logFile->print("Use the DPF parameter 'unbound_model' instead.\n");
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

        thisParameter.hbond = int_hbond_type <= A2 ? HBondType(int_hbond_type) : NON;
        thisParameter.epsij *= this->coeff_vdW;
        thisParameter.epsijHB *= this->coeff_hbond;

        insertAtomParameter(thisParameter.autogridType, thisParameter);
        if (filename[0])
            logFile->printFormatted("Parameters for the atom type \"%s\" were read in from \"%s\" as follows:\n\n", thisParameter.autogridType, filename);
        else
            logFile->printFormatted("Parameters for the atom type \"%s\" were initialised with the following default values:\n\n", thisParameter.autogridType);

        if (outputLevel > 2)
            logFile->printFormatted("\tR-eqm = %5.2f Angstrom\n\tweighted epsilon = %5.3f\n\tAtomic fragmental volume = %5.3f\n\tAtomic solvation parameter = %5.3f\n\tH-bonding R-eqm = %5.3f\n\tweighted H-bonding epsilon = %5.3f\n\tH-bonding type = %d,  bond index = %d\n\n",
                    thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                    thisParameter.RijHB, thisParameter.epsijHB, thisParameter.hbond, thisParameter.bondIndex);
        else
            logFile->printFormatted("\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,\n\tAt.frag.vol. = %.3f,  At.solv.par. = %.3f,\n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,\n\tHb type = %d,  bond index = %d\n\n",
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
