#pragma once

namespace Rune
{
    /**
        Trida kruhu ve 2D
    **************************************************************************************************/
    template<typename T>
    class Circle2
    {
    public:
        Vec2<T> pos;
        T radius;

        Circle2() {}
        Circle2(const Vec2<T> &position, const T circleradius): pos(position), radius(circleradius) {}
    };

    typedef Circle2<float> Circle2f;
    typedef Circle2<double> Circle2d;
}
