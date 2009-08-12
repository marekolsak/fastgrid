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

#include "PairwiseInteractionEnergies.h"

PairwiseInteractionEnergies::PairwiseInteractionEnergies()
{
    energyLookup = new LookupTable();
}

PairwiseInteractionEnergies::~PairwiseInteractionEnergies()
{
    delete energyLookup;
}

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
                cA = (tmpconst = epsij / (xA - xB)) * pow(Rij, xA) * xB;
                cB = tmpconst * pow(Rij, xB) * xA;
                if (isnan(cA))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cA is not a number.  " APPNAME " must exit.");
                if (isnan(cB))
                    logFile.printError(FATAL_ERROR, "Van der Waals coefficient cB is not a number.  " APPNAME " must exit.");
                dxA = xA;
                dxB = xB;
                if (xA == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xA is 0.  " APPNAME " must exit.");
                if (xB == 0)
                    logFile.printError(FATAL_ERROR, "Van der Waals exponent xB is 0.  " APPNAME " must exit.");

                logFile.printFormatted("\n             %9.1lf       %9.1lf \n"
                                       "    E    =  -----------  -  -----------\n"
                                       "     %s, %s         %2d              %2d\n"
                                       "                r               r \n\n"
                                       "Calculating energies for %s-%s interactions.\n",
                                       cA, cB, gridmaps[ia].type, receptorTypes[i], xA, xB, gridmaps[ia].type, receptorTypes[i]);

                // loop over distance index, indexR, from 0 to MAX_DIST
                for (int indexR = 1; indexR < MAX_DIST; indexR++)
                {
                    double r = indexToAngstrom<double>(indexR);
                    rA = pow(r, dxA);
                    rB = pow(r, dxB);
                    energyLookup->table[i][indexR][ia] = Mathd::Min(EINTCLAMP, (cA / rA - cB / rB));
                }               // for each distance
                energyLookup->table[i][0][ia] = EINTCLAMP;
                energyLookup->table[i][MAX_DIST-1][ia] = 0;

                // smooth with min function *//* GPF_MAP
                if (iSmooth > 0)
                {
                    for (int indexR = 1; indexR < MAX_DIST; indexR++)
                    {
                        energySmooth[indexR] = 100000;
                        for (int j = Mathi::Max(0, indexR - iSmooth); j < Mathi::Min(MAX_DIST, indexR + iSmooth); j++)
                            energySmooth[indexR] = Mathd::Min(energySmooth[indexR], energyLookup->table[i][j][ia]);
                    }
                    for (int indexR = 1; indexR < MAX_DIST; indexR++)
                        energyLookup->table[i][indexR][ia] = energySmooth[indexR];
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
                logFile.printFormatted("%4.1lf", indexToAngstrom<double>(j));
                for (int iat = 0; iat < numReceptorTypes; iat++)
                    logFile.printFormatted((energyLookup->table[iat][j][ia] < 100000) ? "%9.2lf" : "%9.2lg", energyLookup->table[iat][j][ia]);               // iat
                logFile.print("\n");
            }                   // j
            logFile.print("\n");
        }
        else
            // parsing for intnbp not needed for covalent maps
            logFile.print("\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
}
