#pragma once
#include "math.h"

template<typename T>    // T should be an integer (short, int, ...)
class Octree
{
public:
    Octree(): elements(0), spheres(0), nodes(0) {}
    ~Octree()
    {
        if (elements)
            delete [] elements;
        if (spheres)
            delete [] spheres;
        if (nodes)
            delete [] nodes;
    }

    void create(const Vec3d &gridSize, const Vec3d &centerPos, int maxNestedLevel, int maxElements, int maxElementsInLastLevel)
    {
        this->gridSize = gridSize;
        this->centerPos = centerPos;
        this->maxNestedLevel = maxNestedLevel;
        this->maxElements = maxElements;
        this->maxElementsInLastLevel = maxElementsInLastLevel;

        box.SetCenterAndExtents(centerPos, gridSize * 0.5);

        numElements = maxNestedLevel > 0 ? maxElements : maxElementsInLastLevel;
        numAdded = 0;
        elements = new T[numElements];

        // no need to store spheres in the last leaf
        spheres = maxNestedLevel > 0 ? new Sphere[numElements] : 0;
        nodes = 0;
    }

    const Octree *getNodeByPos(const Vec3d &pos) const
    {
        const Octree *p = this;
        while (p->nodes)
        {
            int index;
            if (pos.x < p->centerPos.x)
                if (pos.y < p->centerPos.y)
                    if (pos.z < p->centerPos.z)
                        index = 0;
                    else
                        index = 1;
                else
                    if (pos.z < p->centerPos.z)
                        index = 2;
                    else
                        index = 3;
            else
                if (pos.y < p->centerPos.y)
                    if (pos.z < p->centerPos.z)
                        index = 4;
                    else
                        index = 5;
                else
                    if (pos.z < p->centerPos.z)
                        index = 6;
                    else
                        index = 7;

            p = p->nodes + index;
        }
        return p;
    }

    int getNumElements() const
    {
        return numAdded;
    }

    T getElement(int index) const
    {
        return elements[index];
    }

    void insertSphere(const Sphere3d &s, T id, bool skipTest = false)
    {
        if (!skipTest)
        {
            // if the sphere is outside the box, return
            bool in = Intersect(box, SOLID, s, SOLID);
            if (!in)
                return;

            // if the box is inside the sphere while not intersecting its surface, we can assume
            // the test does not need to be carried out in the following nodes
            skipTest = !Intersect(box, SOLID, s, HOLLOW);
        }

        // If this is a leaf
        if (elements)
        {
            // Insert a new element if possible
            if (numAdded < numElements)
            {
                elements[numAdded] = id;
                if (spheres)
                {
                    spheres[numAdded].s = s;
                    spheres[numAdded].skipTest = skipTest;
                }
                ++numAdded;
            }
            else if (maxNestedLevel > 0)
            {
                Vec3d halfGridSize = gridSize * 0.5;
                Vec3d quarterGridSize = gridSize * 0.25;

                const Vec3d shift[] = {
                    Vec3d(-1, -1, -1),
                    Vec3d(-1, -1,  1),
                    Vec3d(-1,  1, -1),
                    Vec3d(-1,  1,  1),
                    Vec3d( 1, -1, -1),
                    Vec3d( 1, -1,  1),
                    Vec3d( 1,  1, -1),
                    Vec3d( 1,  1,  1)
                };

                // Otherwise, try to create 8 new leaves
                nodes = new Octree[8];
                for (int i = 0; i < 8; i++)
                {
                    nodes[i].create(halfGridSize, centerPos + quarterGridSize * shift[i], maxNestedLevel-1, maxElements, maxElementsInLastLevel);

                    // Move all elements to each leaf
                    for (int j = 0; j < numElements; j++)
                        nodes[i].insertSphere(spheres[j].s, elements[j], spheres[j].skipTest);

                    // And insert the new one
                    nodes[i].insertSphere(s, id, skipTest);
                }

                // Make this leaf a node by removing the elements
                delete [] elements;
                if (spheres)
                    delete [] spheres;
                elements = 0;
                spheres = 0;
            }
        }
        else
            // This is a node
            for (int i = 0; i < 8; i++)
                nodes[i].insertSphere(s, id, skipTest);
    }

private:
    // parameters of the create functions
    Vec3d gridSize, centerPos;
    int maxNestedLevel, maxElements, maxElementsInLastLevel;

    // Node properties
    AxisAlignedBox3d box;

    // Elements
    int numElements, numAdded;
    T *elements;

    struct Sphere
    {
        Sphere3d s;
        bool skipTest;
    } *spheres;

    // Nodes
    Octree *nodes;
};
