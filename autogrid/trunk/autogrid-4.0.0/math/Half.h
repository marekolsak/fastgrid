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

#pragma once

namespace Rune
{
    /**
        Datovy typ s plovouci desetinnou carkou o vel. 2 bajty (polovicni float).

        Vyhodou tohoto typu je mala velikost, a proto se hodi jako nahrada za float pri praci se
        streamy (prace se soubory, sitovy prenos).
        POZOR: Tento typ neni vhodny pro vypocty, protoze je emulovan a tudiz je extremne pomaly.
        Rychle s nim umi pracovat jen nektere graficke karty.
    **************************************************************************************************/
    class half
    {
    public:
        half() {}
        half(uint16 s): m_var(s) {}
        half(float f): m_var(FloatToHalf(f)) {}

        void operator =(half f)         { m_var = f.m_var; }
        void operator +=(half f)        { m_var = FloatToHalf(HalfToFloat(m_var) + HalfToFloat(f.m_var)); }
        void operator -=(half f)        { m_var = FloatToHalf(HalfToFloat(m_var) - HalfToFloat(f.m_var)); }
        void operator *=(half f)        { m_var = FloatToHalf(HalfToFloat(m_var) * HalfToFloat(f.m_var)); }
        void operator /=(half f)        { m_var = FloatToHalf(HalfToFloat(m_var) / HalfToFloat(f.m_var)); }

        half operator +(half h) const   { return half(HalfToFloat(m_var) + HalfToFloat(h.m_var)); }
        half operator -(half h) const   { return operator +(-h); }
        half operator *(half h) const   { return half(HalfToFloat(m_var) * HalfToFloat(h.m_var)); }
        half operator /(half h) const   { return half(HalfToFloat(m_var) / HalfToFloat(h.m_var)); }
        half operator -() const         { return half(uint16((m_var & 0x8000)? m_var & ~0x8000 : m_var | 0x8000)); }

        bool operator ==(half h) const  { return m_var == h.m_var; }
        bool operator !=(half h) const  { return !operator ==(h); }
        bool operator <(half h) const   { return HalfToFloat(m_var) < HalfToFloat(h.m_var); }
        bool operator >(half h) const   { return HalfToFloat(m_var) > HalfToFloat(h.m_var); }
        bool operator <=(half h) const  { return !operator >(h); }
        bool operator >=(half h) const  { return !operator <(h); }

        operator float() const          { return HalfToFloat(m_var); }
        uint16 Get() const              { return m_var; }

    private:
        uint16 m_var;

        RUNEMATH_API static float HalfToFloat(uint16 val);
        RUNEMATH_API static uint16 FloatToHalf(float val);

        union ieee_half
        {
            uint16 bits;
            struct ieee_bits
            {
                uint32 m : 10;   // mantissa ... rozsah 0-1023, presnost prakticky na 3 desetinna mista (dost malo)
                uint32 e : 5;    // exponent ... rozsah 0-31, pro vypocty kolem 3D grafiky je to dostacujici
                uint32 s : 1;    // sign
            } ieee;
        };

        union ieee_single // ieee-754 single floating point type
        {
            float f;
            struct ieee_bits
            {
                uint32 m : 23;   // mantissa ... rozsah 0-8388608, presnost nejcasteji na 7 desetinnych mist
                uint32 e : 8;    // exponent ... rozsah 0-255, bohate staci
                uint32 s : 1;    // sign
            } ieee;
        };
    };
}
