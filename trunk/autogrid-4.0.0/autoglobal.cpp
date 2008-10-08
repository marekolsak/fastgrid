#include "autoglobal.h"

char    *programname;
char    *AutoDockHelp = "           -p parameter_filename\n           -l log_filename\n           -o (Use old PDBQ format, charge q in columns 55-61)\n"
                        "           -k (Keep original residue numbers)\n           -i (Ignore header-checking)\n           -t (Parse the PDBQ file to check torsions, then stop.)\n"
                        "           -c < command_file (Command mode, by file)\n           -c | control_program (Command mode, by control_program)\n\n";

char    *AutoGridHelp = "-p parameter_filename\n-l log_filename\n-o (old PDBQ format)\n-d (increment debug level)\n-u (display this message)\n";

char    dock_param_fn[MAX_CHARS];
char    grid_param_fn[MAX_CHARS];

int     command_mode = FALSE;
int     debug = 0;
int	    ElecMap = 0;
int	    DesolvMap = 0;
int     ignore_errors = FALSE;
int     keepresnum = 1;
int     oldpdbq = FALSE;
int     parse_tors_mode = FALSE;
int	    true_ligand_atoms = 0;

FILE    *command_in_fp;
FILE    *command_out_fp;
FILE    *parFile;
FILE    *logFile;

#ifdef USE_DOUBLE
Real	idct = 1.0L;
#else
Real	idct = 1.0;
#endif

Linear_FE_Model AD4;

FILE    *stateFile;
int     write_stateFile = FALSE;
