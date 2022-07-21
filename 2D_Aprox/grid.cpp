#include "grid.h"
//Grid::Grid()
//{

//}

void Grid::rebuild_domain(int nx, int ny, double x0, double y0, double x1, double y1)
{
    delete [] this->triangles;
    delete [] this->points;
    this->nx = nx;
    this->ny = ny;
    this->triangles = new Triangle[nx * ny * 2];
    this->points = new Point[(nx + 1) * (ny + 1)];
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    for(int i = 0; i < nx + 1; i++){
        for(int j = 0; j < ny + 1; j++){
            points[i * (nx + 1) + j].x = x0 + j * hx;
            points[i * (nx + 1) + j].y = y0 + i * hy;
            points[i * (nx + 1) + j].z = 0;
        }
    }
    //fill traingles with filled points
    for(int i = 0; i < 2 * ny; i++){
        for(int j = 0; j < nx; j++){
            triangles[ny * i + j].a = points[j + i / 2 * (ny + 1)];
            triangles[ny * i + j].b = points[j + 1 + i / 2 * (ny + 1)];
            triangles[ny * i + j].c = points[j + (i / 2 + 1) * (ny + 1)];

            triangles[ny * (i + 1) + j].a = points[j + 1 + i / 2 * (ny + 1)];
            triangles[ny * (i + 1) + j].b = points[j + (i / 2 + 1) * (ny + 1)];
            triangles[ny * (i + 1) + j].c = points[j + 1 + (i / 2 + 1) * (ny + 1)];
        }
    }
}

Grid::Grid(int nx, int ny, double x0, double y0, double x1, double y1)
{
    this->nx = nx;
    this->ny = ny;
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    this->triangles = new Triangle[nx * ny * 2];
    this->points = new Point[(nx + 1) * (ny + 1)];
    //fill points
    for(int i = 0; i < nx + 1; i++){
        for(int j = 0; j < ny + 1; j++){
            points[i * (nx + 1) + j].x = x0 + j * hx;
            points[i * (nx + 1) + j].y = y0 + i * hy;
            points[i * (nx + 1) + j].z = 0;
        }
    }
    //fill traingles with filled points
    for(int i = 0; i < 2 * ny; i+=2){
        for(int j = 0; j < nx; j++){
            triangles[ny * i + j].a = points[j + i / 2 * (ny + 1)];
            triangles[ny * i + j].b = points[j + 1 + i / 2 * (ny + 1)];
            triangles[ny * i + j].c = points[j + (i / 2 + 1) * (ny + 1)];

            triangles[ny * (i + 1) + j].a = points[j + 1 + i / 2 * (ny + 1)];
            triangles[ny * (i + 1) + j].b = points[j + (i / 2 + 1) * (ny + 1)];
            triangles[ny * (i + 1) + j].c = points[j + 1 + (i / 2 + 1) * (ny + 1)];
        }
    }
}
