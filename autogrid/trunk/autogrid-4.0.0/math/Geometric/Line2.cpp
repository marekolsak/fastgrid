#include "../All.h"

namespace Rune
{
    /**
        Nastavi body usecky
    **************************************************************************************************/
    template<typename T>
    void Line2<T>::SetPoints(const Vec2<T> &p1, const Vec2<T> &p2)
    {
        origin = Vec2<T>::GetCenter(p1, p2);
        ray.Set((p2-p1).GetNormalized().GetNormal(), origin);
        extent = Vec2<T>::GetDistance(origin, p1);
    }

    template void Line2<float>::SetPoints(const Vec2<float> &p1, const Vec2<float> &p2);
    template void Line2<double>::SetPoints(const Vec2<double> &p1, const Vec2<double> &p2);

    /**
        Pro bod - vraci -1, kdyz je za primkou, 1 kdyz pred primkou a 0, kdyz lezi na primce
    **************************************************************************************************/
    template<typename T>
    int Line2<T>::GetSide(const Vec2<T> &_point) const
    {
        T f = ray.GetDistance(_point);
        if (f > Math<T>::Epsilon()) return 1;
        if (f < -Math<T>::Epsilon()) return -1;
        return 0;
    }

    template int Line2<float>::GetSide(const Vec2<float> &_point) const;
    template int Line2<double>::GetSide(const Vec2<double> &_point) const;
}
