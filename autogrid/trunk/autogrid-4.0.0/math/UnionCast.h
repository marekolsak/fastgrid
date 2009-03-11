#pragma once

namespace Rune
{
    // Sablona pro pretypovani pres union, vhodne pro bezpecnou reinterpretaci vsech datovych typu
    // krome pointeru a referenci (pro tyto se pouziva reinterpret_cast).
    template<typename T1, typename T2>
    inline T1 union_cast(const T2 &a)
    {
        union
        {
            T2 a;
            T1 b;
        } myUnion;
        myUnion.a = a;
        return myUnion.b;
    }
}
