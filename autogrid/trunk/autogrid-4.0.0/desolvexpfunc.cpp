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

#include "desolvexpfunc.h"
#include "autogrid.h"

DesolvExpFunc::DesolvExpFunc(double coeffDesolv)
{
    func = new double[MAX_DIST];
    double minusInvTwoSigmaSquared;
    double sigma;

    // exponential function for receptor and ligand desolvation
    // note: the solvation term will not be smoothed
    sigma = 3.6;
    minusInvTwoSigmaSquared = -1 / (2 * sigma * sigma);
    for (int indexR = 1; indexR < MAX_DIST; indexR++)
    {
        double r = angstrom(indexR);
        func[indexR] = exp(sq(r) * minusInvTwoSigmaSquared);
        func[indexR] *= coeffDesolv;
    }
}

DesolvExpFunc::~DesolvExpFunc()
{
    delete [] func;
}
