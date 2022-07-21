#ifndef GRID_H
#define GRID_H
#include "point.h"
#include "triangle.h"
class Grid
{
public:
    int nx, ny;
    Point a,b;
    Triangle* triangles;
    Point* points;
    void rebuild_domain(int nx, int ny, double x0, double y0, double x1, double y1);
    Grid(int nx, int ny, double x0, double y0, double x1, double y1);
};

#endif // GRID_H
