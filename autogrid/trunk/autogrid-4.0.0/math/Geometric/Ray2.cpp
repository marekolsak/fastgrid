#include "../All.h"

namespace Rune
{
    /**
        Pro bod - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz lezi na primce
    **************************************************************************************************/
    template<typename T>
    int Ray2<T>::GetSide(const Vec2<T> &point) const
    {
        T f = GetDistance(point);
        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Ray2<float>::GetSide(const Vec2<float> &point) const;
    template int Ray2<double>::GetSide(const Vec2<double> &point) const;

    /**
        Pro kruh - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz se dotyka primky
    **************************************************************************************************/
    template<typename T>
    int Ray2<T>::GetSide(const Circle2<T> &circle) const
    {
        T f = GetDistance(circle.pos);
        if (f < circle.radius) return 0;

        if (f > 0) f -= circle.radius;
        else f += circle.radius;

        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Ray2<float>::GetSide(const Circle2<float> &circle) const;
    template int Ray2<double>::GetSide(const Circle2<double> &circle) const;
}
