#pragma once
#include "logfile.h"

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
    PAR_COMMENT
};

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

    ParameterLibrary(LogFile *logFile, int debug, int outputLevel = -1);
    ~ParameterLibrary();
    void load(const char *filename);

private:
    enum { MAXKEY = 256*256 };

    ParameterEntry *dictionary[MAXKEY];
    LogFile *logFile;
    int debug, outputLevel;

    void readLine(const char *line);
    int parseParamLine(const char *line);
    static unsigned int hash(const char *key);
};
