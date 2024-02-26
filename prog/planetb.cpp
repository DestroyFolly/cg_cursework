#include "planetb.h"

PlanetB::PlanetB(int distance, int a_xy, int a_yz, Point star_center,
               int radius_a, int radius_b, int net_nodes_num,
               Point *lights,
               QColor surface_color) : Sphere(CalcCenter(distance, a_xy, a_yz, star_center), radius_a, radius_b, net_nodes_num)
{
    this->surface_color = surface_color;
    CalcOrbit(distance, star_center, a_xy, a_yz);
    this->lights = lights;
    CalcIntensities(lights);
}

Point PlanetB::CalcCenter(int dist, int a_xy, int a_yz, Point star_center)
{
    Point planet_c = star_center;
    planet_c.SetX(star_center.GetX() + dist);
    planet_c.RotateZ(a_xy, star_center);
    planet_c.RotateY(a_yz, star_center);
    return planet_c;
}

void PlanetB::CalcOrbit(double dist, Point sun_center, double a_xy, double a_yz)
{
    int n = GetNetNodesParam();
    orbit = new Point[4 + n * 4];
    orbit[0] = Point(sun_center.GetX() - dist, sun_center.GetY(), sun_center.GetZ());
    double step = 90 / (n + 1), angle = step;
    for (int i = 1; i <= n; i++)
    {
        double z = sin(angle * M_PI / 180) * dist;
        double g = dist * sqrt(2 * (1 - cos(angle * M_PI / 180)));
        double x = sqrt(g * g - z * z);
        orbit[i] = Point(sun_center.GetX() - dist + x, sun_center.GetY(), sun_center.GetZ() - z);
        angle += step;
    }
    orbit[1 + n] = Point(sun_center.GetX(), sun_center.GetY(), sun_center.GetZ() - dist);
    for (int i = 1; i <= n; i++)
    {
        orbit[i + n + 1] = orbit[n + 1 - i];
        orbit[i + n + 1].SetX(2 * sun_center.GetX() - orbit[n + 1 - i].GetX());

        orbit[3 + n * 3 - i] = orbit[n + 1 + i];
        orbit[3 + n * 3 - i].SetZ(sun_center.GetZ() + (sun_center.GetZ() - orbit[n + 1 + i].GetZ()));

        orbit[3 + n * 3 + i] = orbit[n + 1 - i];
        orbit[3 + n * 3 + i].SetZ(sun_center.GetZ() + (sun_center.GetZ() - orbit[n + 1 - i].GetZ()));
    }
    orbit[2 + 2 * n] = Point(sun_center.GetX() + dist, sun_center.GetY(), sun_center.GetZ());
    orbit[3 + 3 * n] = Point(sun_center.GetX(), sun_center.GetY(), sun_center.GetZ() + dist);

    for (int i = 0; i < 4 + 4 * n; i++)
    {
        orbit[i].RotateZ(a_xy, sun_center);
    }
}

void PlanetB::CalcIntensities(Point *lights)
{
    Edge *edges = GetEdges();
    Point *points = GetPoints();
    int n = GetNetNodesParam();

    Vector *normals = new Vector[GetEdgesNum()];
    Point a, b, c;
    Vector ab, ac, ra;
    int I_0 = 1, I;
    double cos_lamb, cos_r;
    double len_n, len_a, len_r;
    double sn, sr;

    Point center = GetCenter();

    for (int i = 0; i < GetEdgesNum(); i++)
    {
        intensities[i] = 0;
        for (int j = 0; j < GetEdgesNum(); j++)
        {
            a = points[edges[i].GetP1()];
            b = points[edges[i].GetP2()];
            c = points[edges[i].GetP3()];

            ab.Set(a, b);
            ac.Set(a, c);
            normals[i] = ab.VectProd(ac);

            ra.Set(a, center);
            len_r = ra.GetLen();
            sr = normals[i].ScalarProd(ra);
            len_n = normals[i].GetLen();
            len_a = a.GetDistance(lights[j]);

            cos_r = sr / len_r / len_n;
            if (cos_r > 0)
                normals[i].Neg();

            sn = normals[i].ScalarProd(Vector(a, lights[j]));
            cos_lamb = sn / len_n / len_a;

            intensities[i] += I_0 * cos_lamb;
        }
    }
    delete [] normals;
}

void PlanetB::RotateX(int degree, Point c)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateX(degree, c);
    }
    for (int i = 0; i < GetOrbitNodesNum(); i++)
    {
        orbit[i].RotateX(degree, c);
    }
    SetCenter(center.RotateX(degree, c));
}

void PlanetB::RotateY(int degree, Point c)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateY(degree, c);
    }
    for (int i = 0; i < GetOrbitNodesNum(); i++)
    {
        orbit[i].RotateY(degree, c);
    }
    SetCenter(center.RotateY(degree, c));
}

void PlanetB::RotateZ(int degree, Point c)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateZ(degree, c);
    }
    for (int i = 0; i < GetOrbitNodesNum(); i++)
    {
        orbit[i].RotateZ(degree, c);
    }
    SetCenter(center.RotateZ(degree, c));
}

void PlanetB::Scale(double kx, double ky, double kz, Point c)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Scale(kx, ky, kz, c);
    }
    for (int i = 0; i < GetOrbitNodesNum(); i++)
    {
        orbit[i].Scale(kx, ky, kz, c);
    }
    SetCenter(center.Scale(kx, ky, kz, c));
}

void PlanetB::Scale(double k, Point c)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Scale(k, c);
    }
    for (int i = 0; i < GetOrbitNodesNum(); i++)
    {
        orbit[i].Scale(k, c);
    }
    SetCenter(center.Scale(k, c));
}

void PlanetB::Move(double dx, double dy, double dz)
{
    Point *points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Move(dx, dy, dz);
    }
    Point tmp = orbit[0];
    for (int i = 0; i < GetOrbitNodesNum() - 1; i++)
    {
        orbit[i] = orbit[i + 1];
    }
    orbit[GetOrbitNodesNum() - 1] = tmp;
    SetCenter(center.Move(dx, dy, dz));
    CalcIntensities(this->lights);
}

int PlanetB::GetOrbitNodesNum()
{
    int n = GetNetNodesParam();
    return 4 + 4 * n;
}

double *PlanetB::GetIntensities()
{
    return intensities;
}

Point *PlanetB::GetOrbit()
{
    return orbit;
}
