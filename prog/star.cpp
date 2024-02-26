#include "star.h"

Star::~Star()
{
    lights.clear();
}

void Star::calc_lights()
{
    int size = GetEdgesNum();
    lights = std::vector<Point> (size);
    Edge* edges = GetEdges();
    Point* points = GetPoints();
    Point a, b, c, mid;
    for (int i = 0; i < size; i++)
    {
        a = points[edges[i].GetP1()];
        b = points[edges[i].GetP2()];
        c = points[edges[i].GetP3()];

        mid.SetX((a.GetX() + b.GetX() + c.GetX()) / 3);
        mid.SetY((a.GetY() + b.GetY() + c.GetY()) / 3);
        mid.SetZ((a.GetZ() + b.GetZ() + c.GetZ()) / 3);

        lights[i] = mid;
    }
}

Star Star::operator=(Star p)
{
    this->lights = p.lights;
    this->m = p.m;
    return *this;
}

void Star::AddOrbitPoint(Point p)
{
    if (orbits_v.size() > 1000)
        orbits_v.clear();
    orbits_v.push_back(Point(p.GetX(), p.GetY(), p.GetZ()));
}

std::vector<Point> Star::GetOrbitNum(void)
{
    std::vector<Point> tmp = orbits_v;
    return tmp;
}

void Star::MashOrbit(double k)
{
    this->orbits_v.clear();
}
