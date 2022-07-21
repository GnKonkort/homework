#include "triangle.h"

Point Triangle::center_of_mass()
{
    Point result;
    result.x = (a.x + b.x + c.x) / 3;
    result.y = (a.y + b.y + c.y) / 3;
    result.z = (a.z + b.z + c.z) / 3;
    return result;
}
