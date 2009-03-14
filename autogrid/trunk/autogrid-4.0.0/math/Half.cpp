/*
    Auxiliary Math library

    Copyright (C) 2003-2009, Marek Olsak (maraeo@gmail.com), All Rights Reserved.
    Copyright (C) 2003-2005, Tomas Pastorek (tomas@tomaspastorek.cz), All Rights Reserved.

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

#include "All.h"

namespace Rune
{
    float half::HalfToFloat(uint16 val)
    {
        ieee_half h;
        h.bits = val;

        if (!h.ieee.e)
        {
            if (!h.ieee.m) return 0;
            return ((h.ieee.s)? -1 : 1) * h.ieee.m * 0.000000059604644775390625f;
        }

        ieee_single f;
        f.ieee.s = h.ieee.s;

        if (h.ieee.e == 31)
        {
            f.ieee.e = 0xFF;
            f.ieee.m = !!h.ieee.m;
        }
        else
        {
            f.ieee.e = h.ieee.e + 112;
            f.ieee.m = h.ieee.m << 13;
        }

        return f.f;
    }

    uint16 half::FloatToHalf(float val)
    {
        ieee_single f;
        ieee_half h;
        f.f = val;

        h.ieee.s = f.ieee.s;

        if ((f.ieee.e==0) && (f.ieee.m==0))
        {
            h.ieee.m = 0;
            h.ieee.e = 0;
        }
        else if ((f.ieee.e==0) && (f.ieee.m!=0))
        {
            h.ieee.m = 0;
            h.ieee.e = 0;
        }
        else if ((f.ieee.e==0xff) && (f.ieee.m==0))
        {
            h.ieee.m = 0;
            h.ieee.e = 31;
        }
        else if ((f.ieee.e==0xff) && (f.ieee.m!=0))
        {
            h.ieee.m = 1;
            h.ieee.e = 31;
        }
        else
        {
            int new_exp = f.ieee.e-127;
            if (new_exp<-24)
            {
                h.ieee.m = 0;
                h.ieee.e = 0;
            }

            if (new_exp<-14)
            {
                uint32 exp_val = uint32(-14 - new_exp); // 2^-exp_val
                h.ieee.e = 0;
                switch (exp_val)
                {
                    case 0: h.ieee.m = 0; break;
                    case 1: h.ieee.m = 512 + (f.ieee.m>>14); break;
                    case 2: h.ieee.m = 256 + (f.ieee.m>>15); break;
                    case 3: h.ieee.m = 128 + (f.ieee.m>>16); break;
                    case 4: h.ieee.m = 64 + (f.ieee.m>>17); break;
                    case 5: h.ieee.m = 32 + (f.ieee.m>>18); break;
                    case 6: h.ieee.m = 16 + (f.ieee.m>>19); break;
                    case 7: h.ieee.m = 8 + (f.ieee.m>>20); break;
                    case 8: h.ieee.m = 4 + (f.ieee.m>>21); break;
                    case 9: h.ieee.m = 2 + (f.ieee.m>>22); break;
                    case 10: h.ieee.m = 1; break;
                }
            }
            else if (new_exp>15)
            {
                h.ieee.m = 0;
                h.ieee.e = 31;
            }
            else
            {
                h.ieee.e = new_exp+15;
                h.ieee.m = (f.ieee.m >> 13);
            }
        }
        return h.bits;
    }
}
