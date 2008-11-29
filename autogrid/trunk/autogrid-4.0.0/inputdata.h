#pragma once
#include "gridmap.h"
#include "programparameters.h"

struct InputData
{
    // if the first char is equal to '\0', the filename is not specified
    char fldFilenameAVS[MAX_CHARS];
    char floatingGridFilename[MAX_CHARS];
    char receptorFilename[MAX_CHARS];
    char xyzFilename[MAX_CHARS];
    char parameterLibraryFilename[MAX_CHARS];   // the AD4 parameters .dat file name

    // variables for RECEPTOR:
    // each type is now at most two characters, eg 'NA\0'
    // NB: these are sparse arrays, some entries are not set
    char receptorTypes[NUM_RECEPTOR_TYPES][3];
    int numReceptorTypes; // number of different receptor atom types actually found in receptor PDBQT

    int numGridPointsPerMap;
    int numReceptorAtoms;

    double charge[AG_MAX_ATOMS];
    double vol[AG_MAX_ATOMS];
    double solpar[AG_MAX_ATOMS];
    int atomType[AG_MAX_ATOMS];
    HBondType hbond[AG_MAX_ATOMS];
    double coord[AG_MAX_ATOMS][XYZ];

    double cgridmin[XYZ];
    double center[XYZ];
    double covpos[XYZ];         // Cartesian-coordinate of covalent affinity well.
    int ne[XYZ];
    int n1[XYZ];
    int nelements[XYZ];

    double epsilon[MAX_DIST];

    // for NEW3 desolvation terms
    double solparQ;   // unweighted value restored 3:9:05
    double invDielCal;
    double rSmooth;
    double spacing;     // One quarter of a C-C bond length.
    double covHalfWidth;
    double covBarrier;

    bool distDepDiel, disorderH;
};

class InputDataLoader : public InputData
{
public:
    InputDataLoader(LogFile *logFile);
    void load(const char *gridParameterFilename, GridMapList &gridmaps, ParameterLibrary &parameterLibrary);

private:
    LogFile *logFile;

    int checkSize(int nelements, char axischar);
    static int parseGPFLine(const char *line);
    static double calculateDDDMehlerSolmajer(double distance, double approx_zero);
    static int parseTypes(char *line, char **words, int maxwords);
    static int strIndex(const char *s, const char *t);
};
