#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "point.h"

class Triangle
{
public:
    Point a,b,c;
    Point center_of_mass();
    Triangle() = default;
};

#endif // TRIANGLE_H
