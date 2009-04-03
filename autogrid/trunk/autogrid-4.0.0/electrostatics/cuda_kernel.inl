/*
    AutoGrid

    Copyright (C) 1989-2007, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
    All Rights Reserved.
    Copyright (C) 2008-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.

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

// The generic kernel source code
{
    int x = blockIdx.x;
    int y = blockIdx.y;

    float gridPosX = (x - numGridPointsDiv2.x) * gridSpacing;
    float gridPosY = (y - numGridPointsDiv2.y) * gridSpacing;

    float energy = 0;

    //  Do all Receptor (protein, DNA, etc.) atoms...
    for (int ia = 0; ia < numAtoms; ia++)
    {
        // Get the distance from current grid point to this receptor atom (|receptorAtomCoord - gridPos|)
        float dx = atoms[ia].x - gridPosX;
        float dy = atoms[ia].y - gridPosY;
        float rSq = dx*dx + dy*dy + atoms[ia].z;
        float r = sqrt(rSq);
        float invR = 1 / r;

        // The estat forcefield coefficient/weight is premultiplied
        energy += atoms[ia].w * min(invR, 2.f)
#if defined(DDD)
            * epsilon[min(int(r * A_DIVISOR), int(MAX_DIST-1))]
#endif
            ;
    }

    // Round to 3 decimal places
    int outputIndex = outputIndexZBase + y * numGridPointsX + x;
    outEnergies[outputIndex] += abs(energy) < 0.0005f ? 0.f : rint(energy * 1000.f) * 0.001f;
}
