#include "pairwiseinteractionenergies.h"
#include <cmath>

void PairwiseInteractionEnergies::calculate(const GridMapList &gridmaps, LogFile &logFile,
                                            int numReceptorTypes, const char (&receptorTypes)[NUM_RECEPTOR_TYPES][3], double rSmooth)
{
    double energySmooth[MAX_DIST];
    double dxA;
    double dxB;
    double rA;
    double rB;
    double Rij, epsij;
    double cA, cB, tmpconst;
    int xA, xB;

    // Angstrom is divided by A_DIVISOR in look-up table.
    // Typical value of rSmooth is 0.5 Angstroms
    // so iSmooth = 0.5 * 100 / 2 = 25
    int iSmooth = int(rSmooth * A_DIVISOR / 2);

    logFile.print("\n\nCalculating Pairwise Interaction Energies\n"
                  "=========================================\n\n");

    // do the map stuff here:
    // set up xA, xB, npb_r, npb_eps and hbonder
    // before this pt
    for (int ia = 0; ia < gridmaps.getNumAtomMaps(); ia++)
        if (!gridmaps[ia].isCovalent)
        {
            // i is the index of the receptor atom type, that the ia type ligand probe will interact with. *//* GPF_MAP
            for (int i = 0; i < numReceptorTypes; i++)
            {
                // for each receptor_type
                xA = gridmaps[ia].xA[i];
                xB = gridmaps[ia].xB[i];
                Rij = gridmaps[ia].nbpR[i];
                epsij = gridmaps[ia].nbpEps[i];

                // for each receptor_type get its parms and fill in tables
                cA = (tmpconst = epsij / (double)(xA - xB)) * pow(Rij, (double)xA) * (double)xB;
                cB = tmpconst * pow(Rij, (double)xB) * (double)xA;
                if (isnan(cA))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cA is not a number.  AutoGrid must exit.");
                if (isnan(cB))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cB is not a number.  AutoGrid must exit.");
                dxA = (double)xA;
                dxB = (double)xB;
                if (xA == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xA is 0.  AutoGrid must exit.");
                if (xB == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xB is 0.  AutoGrid must exit.");

                logFile.printFormatted("\n             %9.1lf       %9.1lf \n"
                                       "    E    =  -----------  -  -----------\n"
                                       "     %s, %s         %2d              %2d\n"
                                       "                r               r \n\n"
                                       "Calculating energies for %s-%s interactions.\n",
                                       cA, cB, gridmaps[ia].type, receptorTypes[i], xA, xB, gridmaps[ia].type, receptorTypes[i]);

                // loop over distance index, indx_r, from 0 to MAX_DIST
                for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                {
                    double r = angstrom(indx_r);
                    rA = pow(r, dxA);
                    rB = pow(r, dxB);
                    energyLookup[i][indx_r][ia] = min(EINTCLAMP, (cA / rA - cB / rB));
                }               // for each distance
                energyLookup[i][0][ia] = EINTCLAMP;
                energyLookup[i][MAX_DIST-1][ia] = 0;

                // smooth with min function *//* GPF_MAP
                if (iSmooth > 0)
                {
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                    {
                        energySmooth[indx_r] = 100000;
                        for (int j = max(0, indx_r - iSmooth); j < min(MAX_DIST, indx_r + iSmooth); j++)
                            energySmooth[indx_r] = min(energySmooth[indx_r], energyLookup[i][j][ia]);
                    }
                    for (int indx_r = 1; indx_r < MAX_DIST; indx_r++)
                        energyLookup[i][indx_r][ia] = energySmooth[indx_r];
                }               // endif smoothing
            }                   // for i in receptor types: build energy table for this map

            // Print out a table, of distance versus energy...
            // GPF_MAP
            logFile.printFormatted("\n\nFinding the lowest pairwise interaction energy within %.1f Angstrom (\"smoothing\").\n\n  r ", rSmooth);
            for (int iat = 0; iat < numReceptorTypes; iat++)
                logFile.printFormatted("    %s    ", receptorTypes[iat]);
            logFile.print("\n ___");
            for (int iat = 0; iat < numReceptorTypes; iat++)
                logFile.print(" ________");                   // iat
            logFile.print("\n");
            for (int j = 0; j <= 500; j += 10)
            {
                logFile.printFormatted("%4.1lf", angstrom(j));
                for (int iat = 0; iat < numReceptorTypes; iat++)
                    logFile.printFormatted((energyLookup[iat][j][ia] < 100000) ? "%9.2lf" : "%9.2lg", energyLookup[iat][j][ia]);               // iat
                logFile.print("\n");
            }                   // j
            logFile.print("\n");
        }
        else
            // parsing for intnbp not needed for covalent maps
            logFile.print("\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
}
