#pragma once

/*#if defined(_MSC_VER)
    #if !defined(RUNE_STATIC)
        #if defined(RUNEMATH_EXPORTS)
            #define RUNEMATH_API __declspec(dllexport)
        #else
            #define RUNEMATH_API __declspec(dllimport)
        #endif
    #else
        #define RUNEMATH_API
    #endif
#else*/
    #define RUNEMATH_API
//#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include "Int.h"
#include "UnionCast.h"

// Zaklad
#include "Half.h"
#include "Math.h"

// Vektory
#include "Geometric/Vec2.h"
#include "Geometric/Vec3.h"
#include "Geometric/Vec4.h"

// 2D
#include "Geometric/Circle2.h"
#include "Geometric/Ray2.h"
#include "Geometric/Line2.h"
#include "Geometric/Rect2.h"
#include "Geometric/Intersection2.h"

// 3D
#include "Geometric/Quaternion.h"
#include "Geometric/Matrix3.h"
#include "Geometric/Matrix4.h"
#include "Geometric/Vertex3.h"
#include "Geometric/Ray3.h"
#include "Geometric/Sphere3.h"
#include "Geometric/AxisAlignedBox3.h"
#include "Geometric/Plane3.h"
#include "Geometric/Frustum3.h"
#include "Geometric/Intersection3.h"
