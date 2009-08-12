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

#pragma once
#include "LogFile.h"

#define MAX_LEN_AUTOGRID_TYPE 7

// hbonding character:
enum HBondType
{
    NON,    // none
    DS,     // spherical donor
    D1,     // directional donor
    AS,     // spherical acceptor
    A1,     // acceptor of 1 directional hbond
    A2      // acceptor of 2 directional hbonds
};

struct ParameterEntry
{
    char autogridType[MAX_LEN_AUTOGRID_TYPE + 1];    // autogridType is based on babel_types assigned by PyBabel
    double Rij;         // Lennard-Jones equilibrium separation
    double epsij;       // Lennard-Jones energy well-depth
    double vol;         // solvation volume
    double solpar;      // solvation parameter
    HBondType hbond;    // hbonding character
    double RijHB;       // 12-10 Lennard-Jones equilibrium separation
    double epsijHB;     // 12-10 Lennard-Jones energy well-depth
    int recIndex;       // used to set up receptor atom_types
    int mapIndex;       // used to set up map atom_types
    int bondIndex;      // used to set up bonds; corresponds to the enum in mdist.h
};

// Free energy coefficients and atom parameters
class ParameterLibrary
{
public:
    // Linear free energy model
    double coeff_vdW;                 // Free energy coefficient for van der Waals term
    double coeff_hbond;               // Free energy coefficient for H-bonding term
    double coeff_estat;               // Free energy coefficient for electrostatics term
    double coeff_desolv;              // Free energy coefficient for desolvation term
    double coeff_tors;                // Free energy coefficient for torsional term

    double stderr_vdW;                // Free energy standard error for van der Waals term
    double stderr_hbond;              // Free energy standard error for H-bonding term
    double stderr_estat;              // Free energy standard error for electrostatics term
    double stderr_desolv;             // Free energy standard error for desolvation term
    double stderr_tors;               // Free energy standard error for torsional term

    // Atom parameters
    void insertAtomParameter(const char *key, const ParameterEntry &value);
    ParameterEntry *findAtomParameter(const char *key);
    const ParameterEntry *findAtomParameter(const char *key) const;
    int getAtomParameterRecIndex(const char *key) const;

    ParameterLibrary(LogFile *logFile, const char *modelText, Unbound_Model unboundModel, int debug, int outputLevel = -1);
    ~ParameterLibrary();
    void load(const char *filename);

private:
    enum { MAXKEY = 256*256 };

    // Define tokens for parsing AutoDock atomic parameter files
    enum ParserTokens
    {
        PAR_ = -1,
        PAR_NULL = 0,
        PAR_VDW,
        PAR_HBOND,
        PAR_ESTAT,
        PAR_DESOLV,
        PAR_TORS,
        PAR_ATOM_PAR,
        PAR_COMMENT,
        PAR_UNBOUND
    };

    char filename[MAX_CHARS];
    ParameterEntry *dictionary[MAXKEY];
    LogFile *logFile;
    int debug, outputLevel;

    void readLine(const char *line);
    ParserTokens parseParamLine(const char *line);
    static unsigned int hash(const char *key);
};
