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
#include "BondVectors.h"

BondVectors::BondVectors(int numReceptorAtoms, LogFile *logFile): logFile(logFile)
{
    disorder = new bool[numReceptorAtoms];
    rexp = new int[numReceptorAtoms];
    rvector = new Vec3d[numReceptorAtoms];
    rvector2 = new Vec3d[numReceptorAtoms];
    
    memset(rexp, 0, sizeof(sizeof(int) * numReceptorAtoms));
}

BondVectors::~BondVectors()
{
    delete [] disorder;
    delete [] rexp;
    delete [] rvector;
    delete [] rvector2;
}

void BondVectors::calculate(const InputData *input, const ParameterLibrary &parameterLibrary)
{
    Vec3d d;
    Vec3d dc;
    double rdot;
    int from, to;
    int nbond;
    int i1 = 0, i2 = 0, i3 = 0;
    double inv_rd, rd2;

    // Loop over all RECEPTOR atoms to
    // calculate bond vectors for directional H-bonds
    // setup the canned atom types here....
    // at this point set up hydrogen, carbon, oxygen and nitrogen
    int hydrogen = parameterLibrary.getAtomParameterRecIndex("HD");
    int nonHB_hydrogen = parameterLibrary.getAtomParameterRecIndex("H");
    int carbon = parameterLibrary.getAtomParameterRecIndex("C");
    int arom_carbon = parameterLibrary.getAtomParameterRecIndex("A");
    int oxygen = parameterLibrary.getAtomParameterRecIndex("OA");
    int sulphur = parameterLibrary.getAtomParameterRecIndex("SA");

    // These are not used
    /*int nitrogen = parameterLibrary.getAtomParameterRecIndex("NA");
    int nonHB_nitrogen = parameterLibrary.getAtomParameterRecIndex("N");
    int nonHB_sulphur = parameterLibrary.getAtomParameterRecIndex("S");*/

    // 7:CHANGE HERE: scan the 'mapIndex' from the input
    for (int ia = 0; ia < input->numReceptorAtoms; ia++)
    {
        //** ia = i_receptor_atom_a **
        disorder[ia] = false;   // initialize disorder flag.
        bool warned = false;

        // Set scan limits looking for bonded atoms
        from = Mathi::Max(ia - 20, 0);
        to = Mathi::Min(ia + 20, input->numReceptorAtoms - 1);

        // If 'ia' is a hydrogen atom, it could be a
        // RECEPTOR hydrogen-BOND DONOR,
        // TODO: 8:CHANGE HERE: fix the input->atomType vs atom_types problem in following
        if (input->hbond[ia] == D1) // D1 hydrogen bond donor
        {
            for (int ib = from; ib <= to; ib++)
                if (ib != ia) // ib = i_receptor_atom_b
                {
                    // =>  NH-> or OH->
                    // if ((input->atomType[ib] == nitrogen) || (input->atomType[ib]==nonHB_nitrogen) ||(input->atomType[ib] == oxygen)||(input->atomType[ib] == sulphur)||(input->atomType[ib]==nonHB_sulphur)) {

                    // Calculate the square of the N-H or O-H bond distance, rd2,
                    //                            ib-ia  ib-ia
                    d = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[ib]);
                    rd2 = d.MagnitudeSqr();
                    // If ia & ib are less than 1.3 A apart -- they are covalently bonded,
                    if (rd2 < 1.90)
                    {
                        // INCREASED for H-S bonds
                        if (rd2 < APPROX_ZERO)
                        {
                            if (rd2 == 0)
                                logFile->printErrorFormatted(WARNING,
                                    "While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n",
                                    ia + 1, ib + 1);
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = Mathd::Rsqrt(rd2);

                        // N-H: Set exponent rexp to 2 for m/m H-atom,
                        // if (input->atomType[ib] == nitrogen) rexp[ia] = 2;
                        if ((input->atomType[ib] != oxygen) && (input->atomType[ib] != sulphur))
                            rexp[ia] = 2;

                        // O-H: Set exponent rexp to 4 for m/m H-atom,
                        // and flag disordered hydroxyls
                        if ((input->atomType[ib] == oxygen) || (input->atomType[ib] == sulphur))
                        {
                            rexp[ia] = 4;
                            if (input->disorderH)
                                disorder[ia] = true;
                        }

                        // Normalize the vector from ib to ia, N->H or O->H...
                        rvector[ia] = d * inv_rd;

                        // First O-H/N-H H-bond-donor found; Go on to next atom,
                        break;
                    }           // Found covalent bond.
                    // } Found NH or OH in receptor.
                }
            // Finished scanning for the NH or OH in receptor.
            // If 'ia' is an Oxygen atom, it could be a
            // RECEPTOR H_BOND ACCEPTOR,
        }
        else if (input->hbond[ia] == A2)
        {
            // A2
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia + 5)th m/m-atom
            // determine number of atoms bonded to the oxygen
            nbond = 0;
            int ib = from;
            for (; ib <= to; ib++)
                if (ib != ia)
                {
                    dc = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[ib]);
                    rd2 = dc.MagnitudeSqr();

                    if (((rd2 < 3.61) && ((input->atomType[ib] != hydrogen) && (input->atomType[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((input->atomType[ib] == hydrogen) || (input->atomType[ib] == nonHB_hydrogen))))
                    {
                        if (nbond == 2)
                            logFile->printErrorFormatted(WARNING, "Found an H-bonding atom with three bonded atoms, atom serial number %d\n", ia + 1);
                        if (nbond == 1)
                        {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0)
                        {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }               // (ib != ia)

            // if no bonds, something is wrong
            if (nbond == 0)
                logFile->printErrorFormatted(WARNING, "Oxygen atom found with no bonded atoms, atom serial number %d, input->atomType %d\n", ia + 1, input->atomType[ia]);

            // one bond: Carbonyl Oxygen O=C-X
            if (nbond == 1)
            {
                // calculate normalized carbonyl bond vector rvector[ia][]
                rvector[ia] = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[i1]);
                rd2 = rvector[ia].MagnitudeSqr();
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && !warned)
                    {
                        logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, i1 + 1);
                        warned = true;
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = Mathd::Rsqrt(rd2);
                rvector[ia] *= inv_rd;

                // find a second atom (i2) bonded to carbonyl carbon (i1)
                for (int i2 = from; i2 <= to; i2++)
                    if ((i2 != i1) && (i2 != ia))
                    {
                        dc = Vec3d(input->receptorAtom[i1]) - Vec3d(input->receptorAtom[i2]);
                        rd2 = dc.MagnitudeSqr();
                        if (((rd2 < 2.89) && (input->atomType[i2] != hydrogen)) || ((rd2 < 1.69) && (input->atomType[i2] == hydrogen)))
                        {
                            // found one
                            // d[i] vector from carbon to second atom
                            d = Vec3d(input->receptorAtom[i2]) - Vec3d(input->receptorAtom[i1]);
                            rd2 = d.MagnitudeSqr();
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0) && !warned)
                                {
                                    logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                                    warned = true;
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = Mathd::Rsqrt(rd2);
                            d *= inv_rd;

                            // C=O cross C-X gives the lone pair plane normal
                            rvector2[ia] = Vec3d::Cross(rvector[ia], d);
                            rd2 = rvector2[ia].MagnitudeSqr();
                            if (rd2 < APPROX_ZERO)
                            {
                                if ((rd2 == 0) && !warned)
                                {
                                    logFile->printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                                    warned = true;
                                }
                                rd2 = APPROX_ZERO;
                            }
                            inv_rd = Mathd::Rsqrt(rd2);
                            rvector2[ia] *= inv_rd;
                        }
                    }
            }                   // endif nbond==1

            // two bonds: Hydroxyl or Ether Oxygen X1-O-X2
            if (nbond == 2)
                // disordered hydroxyl
                if ((input->atomType[i1] == hydrogen || input->atomType[i2] == hydrogen) && input->atomType[i1] != input->atomType[i2] && input->disorderH)
                {
                    if ((input->atomType[i1] == carbon) || (input->atomType[i1] == arom_carbon))
                        ib = i1;
                    if ((input->atomType[i2] == carbon) || (input->atomType[i1] == arom_carbon))
                        ib = i2;
                    disorder[ia] = true;
                    rvector[ia] = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[ib]);
                    rd2 = rvector[ia].MagnitudeSqr();
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && !warned)
                        {
                            logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                            warned = true;
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = Mathd::Rsqrt(rd2);
                    rvector[ia] *= inv_rd;
                }
                else
                {
                    // not a disordered hydroxyl
                    // normalized X1 to X2 vector, defines lone pair plane
                    rvector2[ia] = Vec3d(input->receptorAtom[i2]) - Vec3d(input->receptorAtom[i1]);
                    rd2 = rvector2[ia].MagnitudeSqr();
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && !warned)
                        {
                            logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                            warned = true;
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = Mathd::Rsqrt(rd2);
                    rvector2[ia] *= inv_rd;

                    // vector pointing between the lone pairs:
                    // front of the vector is the oxygen atom,
                    // X1->O vector dotted with normalized X1->X2 vector plus
                    // coords of X1 gives the point on the X1-X2 line for the
                    // back of the vector.
                    rdot = Vec3d::Dot(Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[i1]), rvector2[ia]);
                    rvector[ia] = Vec3d(input->receptorAtom[ia]) - (rvector2[ia]*rdot + Vec3d(input->receptorAtom[i1]));
                    rd2 = rvector[ia].MagnitudeSqr();
                    if (rd2 < APPROX_ZERO)
                    {
                        if ((rd2 == 0) && !warned)
                        {
                            logFile->printError(WARNING, "Attempt to divide by zero was just prevented.\n\n");
                            warned = true;
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = Mathd::Rsqrt(rd2);
                    rvector[ia] *= inv_rd;
                }               // end disordered hydroxyl
        }
        else if (input->hbond[ia] == A1)
        {                       // A1
            // Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
            //        to (ia+5)th m/m-atom
            // determine number of atoms bonded to the oxygen
            nbond = 0;
            int ib = from;
            for (; ib <= to; ib++)
                if (ib != ia)
                {
                    dc = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[ib]);
                    rd2 = dc.MagnitudeSqr();

                    if (((rd2 < 2.89) && ((input->atomType[ib] != hydrogen) && (input->atomType[ib] != nonHB_hydrogen))) ||
                        ((rd2 < 1.69) && ((input->atomType[ib] == hydrogen) || (input->atomType[ib] == nonHB_hydrogen))))
                    {
                        if (nbond == 2)
                        {
                            nbond = 3;
                            i3 = ib;
                        }
                        if (nbond == 1)
                        {
                            nbond = 2;
                            i2 = ib;
                        }
                        if (nbond == 0)
                        {
                            nbond = 1;
                            i1 = ib;
                        }
                    }
                }               // (ib != ia)

            // if no bonds, something is wrong
            if (nbond == 0)
                logFile->printErrorFormatted(WARNING, "Nitrogen atom found with no bonded atoms, atom serial number %d\n", ia);

            // one bond: Azide Nitrogen :N=C-X
            if (nbond == 1)
            {
                // calculate normalized N=C bond vector rvector[ia][]
                rvector[ia] = Vec3d(input->receptorAtom[ia]) - Vec3d(input->receptorAtom[i1]);
                rd2 = rvector[ia].MagnitudeSqr();
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && !warned)
                    {
                        logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = true;
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = Mathd::Rsqrt(rd2);
                rvector[ia] *= inv_rd;
            }                   // endif nbond==1

            // two bonds: X1-N=X2
            if (nbond == 2)
            {
                // normalized vector from Nitrogen to midpoint between X1 and X2
                rvector[ia] = Vec3d(input->receptorAtom[ia]) - (Vec3d(input->receptorAtom[i2]) + Vec3d(input->receptorAtom[i1])) / 2;
                rd2 = rvector[ia].MagnitudeSqr();
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && !warned)
                    {
                        logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = true;
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = Mathd::Rsqrt(rd2);
                rvector[ia] *= inv_rd;
            }                   // end two bonds for nitrogen

            // three bonds: X1,X2,X3
            if (nbond == 3)
            {
                // normalized vector from Nitrogen to midpoint between X1, X2, and X3
                rvector[ia] = Vec3d(input->receptorAtom[ia]) - (Vec3d(input->receptorAtom[i1]) + Vec3d(input->receptorAtom[i2]) + Vec3d(input->receptorAtom[i3])) / 3;
                rd2 = rvector[ia].MagnitudeSqr();
                if (rd2 < APPROX_ZERO)
                {
                    if ((rd2 == 0) && !warned)
                    {
                        logFile->printErrorFormatted(WARNING, "Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        warned = true;
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = Mathd::Rsqrt(rd2);
                rvector[ia] *= inv_rd;
            }                   // end three bonds for Nitrogen
            // endNEW directional N Acceptor
        }                       // end test for atom type
    }                           // Do Next receptor atom...
}
