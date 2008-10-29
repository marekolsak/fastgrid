#pragma once
#include <cstdio>
#include "program_parameters.h"
#include "parameters.h"
#include "structs.h"

struct MapObject
{
    int atom_type;          // corresponds to receptor numbers????
    int map_index;
    int is_covalent;
    int is_hbonder;
    FILE *map_fileptr;
    char map_filename[MAX_CHARS];
    char type[3];           // eg HD or OA or NA or N
    double constant;        // this will become obsolete
    double energy_max;
    double energy_min;
    double energy;
    double vol_probe;
    double solpar_probe;
    // new 6/28
    double Rij;
    double epsij;
    hbond_type hbond;       // hbonding character:
    double Rij_hb;
    double epsij_hb;
    // per receptor type parameters, ordered as in receptor_types
    double nbp_r[NUM_RECEPTOR_TYPES];   // radius of energy-well minimum
    double nbp_eps[NUM_RECEPTOR_TYPES]; // depth of energy-well minimum
    int xA[NUM_RECEPTOR_TYPES]; // generally 12
    int xB[NUM_RECEPTOR_TYPES]; // 6 for non-hbonders 10 for h-bonders
    int hbonder[NUM_RECEPTOR_TYPES];
};

struct GridParameterInfo
{
    // variables for RECEPTOR:
    // each type is now at most two characters, eg 'NA\0'
    // NB: these are sparse arrays, some entries are not set
    char receptor_types[NUM_RECEPTOR_TYPES][3];

    // number of different receptor atom types actually found in receptor PDBQT
    int receptor_types_ct;

    // AG_MAX_ATOMS
    double charge[AG_MAX_ATOMS];
    double vol[AG_MAX_ATOMS];
    double solpar[AG_MAX_ATOMS];

    int atom_type[AG_MAX_ATOMS];
    hbond_type hbond[AG_MAX_ATOMS];

    int rexp[AG_MAX_ATOMS];
    double coord[AG_MAX_ATOMS][XYZ];

    // canned atom type number
    int hydrogen, carbon, arom_carbon, oxygen, nitrogen;
    int nonHB_hydrogen, nonHB_nitrogen, sulphur, nonHB_sulphur;

    double cgridmin[XYZ];
    double center[XYZ];
    double covpos[XYZ]; // Cartesian-coordinate of covalent affinity well.

    int ne[XYZ];
    int n1[XYZ];
    int nelements[XYZ];

    // MAX_CHARS
    char AVS_fld_filename[MAX_CHARS];
    char floating_grid_filename[MAX_CHARS];
    char receptor_filename[MAX_CHARS];

    double epsilon[MAX_DIST];

    // for NEW3 desolvation terms
    double solpar_q;
    double invdielcal;
    double percentdone;
    double r_smooth;
    double spacing;     // One quarter of a C-C bond length.
    double covhalfwidth;
    double covbarrier;

    double factor; // Used to convert between calories and SI units

    int num_maps;
    int num_atom_maps;
    int floating_grid, dddiel, disorder_h;
    int elecPE;
    int dsolvPE;

    int outlev;

    int num_grid_points_per_map;
    int i_smooth;

    int map_index;
    int num_receptor_atoms;

    GridParameterInfo();
};

void read_grid_parameter_file(ProgramParameters &programParams, FILE *logFile, Linear_FE_Model &AD4, MapObject *(&gridmap), GridParameterInfo &params);
