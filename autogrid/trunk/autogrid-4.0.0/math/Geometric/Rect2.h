/*
    Linear Algebra / Math library

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
    template<typename T>
    class Rect2
    {
    public:
        Vec2<T> leftTop;
        Vec2<T> rightBottom;

        Rect2() {}
        Rect2(T x1, T y1, T x2, T y2): leftTop(x1, y1), rightBottom(x2, y2) {}
        Rect2(const Vec2<T> &LeftTop, const Vec2<T> &RightBottom): leftTop(LeftTop), rightBottom(RightBottom) {}
        T GetWidth() const              { return rightBottom.x - leftTop.x; }
        T GetHeight() const             { return rightBottom.y - leftTop.y; }
        Vec2<T> GetCenter() const       { return Vec2<T>::GetCenter(leftTop, rightBottom); }
        Vec2<T> GetSize() const         { return rightBottom - leftTop; }
        void SetPos(T x, T y)           { Set(x, y, x + GetWidth(), y + GetHeight()); }
        void SetPos(const Vec2<T> &p)   { SetPos(p.x, p.y); }
        void SetSize(T x, T y)          { rightBottom.Set(leftTop.x + x, leftTop.y + y); }
        void SetSize(const Vec2<T> &s)  { SetSize(s.x, s.y); }
        void Set(T x1, T y1, T x2, T y2) { leftTop.Set(x1, y1); rightBottom.Set(x2, y2); }

        bool operator ==(const Rect2 &r) const      { return leftTop == r.leftTop && rightBottom == r.rightBottom; }
        bool operator !=(const Rect2 &r) const      { return leftTop != r.leftTop && rightBottom != r.rightBottom; }
        Rect2 operator +(const Vec2<T> &p) const    { return Rect2(leftTop + p, rightBottom + p); }
    };

    template<typename T>
    std::ostream &operator <<(std::ostream &s, const Rect2<T> &r)
    {
        s << Vec4<T>(r.leftTop, r.rightBottom);
        return s;
    }

    template<typename T>
    std::istream &operator >>(std::istream &s, Rect2<T> &r)
    {
        Vec4<T> v;
        s >> v;
        r.Set(v.x, v.y, v.z, v.w);
        return s;
    }

    typedef Rect2<float> Rect2f;
    typedef Rect2<double> Rect2d;
    typedef Rect2<int32> Rect2i;
}