#pragma once

#define MAX_LEN_AUTOGRID_TYPE 7

// hbonding character:
enum HBondType
{
    NON,
    DS,
    D1,
    AS,
    A1,
    A2
};

// was "parm_info" in earlier AutoGrid 4 code
struct ParameterEntry
{
    char autogridType[MAX_LEN_AUTOGRID_TYPE + 1];    // autogridType is based on babel_types assigned by PyBabel
    double Rij;            /* Lennard-Jones equilibrium separation */
    double epsij;            /* Lennard-Jones energy well-depth */
    double vol;            /* solvation volume */
    double solpar;        /* solvation parameter */
    HBondType hbond;        /* hbonding character:
                   NON: none,
                   DS: spherical donor
                   D1: directional donor
                   AS: spherical acceptor
                   A1: acceptor of 1 directional hbond
                   A2: acceptor of 2 directional hbonds */
    double RijHB;        /* 12-10 Lennard-Jones equilibrium separation */
    double epsijHB;        /* 12-10 Lennard-Jones energy well-depth */
    int recIndex;        /* used to set up receptor atom_types */
    int mapIndex;        /* used to set up map atom_types */
    int bondIndex;        /* used to set up bonds; corresponds to the enum in mdist.h */
};

class AtomParameterManager
{
public:
    AtomParameterManager();
    ~AtomParameterManager();
    void insert(const char *key, const ParameterEntry &value);
    ParameterEntry *find(const char *key) const;
    int getRecIndex(const char *key) const;

private:
    enum { MAXKEY = 256*256 };
    ParameterEntry *dictionary[MAXKEY];

    static unsigned int hash(const char *key);
};
