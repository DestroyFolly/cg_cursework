#include "planeta.h"
#include <iostream>
#include <cmath>

PlanetA::PlanetA(int radius_a, int radius_b, int net_nodes_num, std::vector<Point> lights,
               QColor surface_color, int type) : Sphere(Point(0.0, 0.0, 0.0), radius_a, radius_b, net_nodes_num, type)
{
    this->surface_color = surface_color;
    this->lights = lights;
    CalcIntensities(lights);
}

void PlanetA::CalcIntensities(std::vector<Point> lights)
{
    int initialIntensity = 1, intensity;
    double cosLambda, cosR, lengthN, lengthA, lengthR, scalarN, scalarR;
    Edge* edgeArray = GetEdges();
    Point* pointArray = GetPoints();
    int netNodes = GetNetNodesParam();

    Vector *normalVectors = new Vector[GetEdgesNum()];
    Point pointA, pointB, pointC;
    Vector vectorAB, vectorAC, vectorRA;

    Point centerPoint = GetCenter();

    for (int edgeIndex = 0; edgeIndex < GetEdgesNum(); edgeIndex++)
    {
        intensities[edgeIndex] = 0;
        for (int lightIndex = 0; lightIndex < GetEdgesNum(); lightIndex++)
        {
            pointA = pointArray[edgeArray[edgeIndex].GetP1()];
            pointB = pointArray[edgeArray[edgeIndex].GetP2()];
            pointC = pointArray[edgeArray[edgeIndex].GetP3()];

            vectorAB.Set(pointA, pointB);
            vectorAC.Set(pointA, pointC);
            normalVectors[edgeIndex] = vectorAB.VectProd(vectorAC);

            vectorRA.Set(pointA, centerPoint);
            lengthR = vectorRA.GetLen();
            lengthN = normalVectors[edgeIndex].GetLen();
            lengthA = pointA.GetDistance(lights[lightIndex]);
            scalarR = normalVectors[edgeIndex].ScalarProd(vectorRA);

            cosR = scalarR / lengthR / lengthN;
            if (cosR > 0)
                normalVectors[edgeIndex].Neg();

            scalarN = normalVectors[edgeIndex].ScalarProd(Vector(pointA, lights[lightIndex]));
            cosLambda = scalarN / lengthN / lengthA;

            intensities[edgeIndex] += initialIntensity * cosLambda;
        }
    }
    delete [] normalVectors;
}

void PlanetA::RotateX(int degree, Point c)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateX(degree, c);
    }
    SetCenter(center.RotateX(degree, c));
}

void PlanetA::RotateY(int degree, Point c)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateY(degree, c);
    }
    SetCenter(center.RotateY(degree, c));
}

void PlanetA::RotateZ(int degree, Point c)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].RotateZ(degree, c);
    }
    SetCenter(center.RotateZ(degree, c));
}

void PlanetA::Scale(double kx, double ky, double kz, Point c)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Scale(kx, ky, kz, c);
    }

    SetCenter(center.Scale(kx, ky, kz, c));
}

void PlanetA::Scale(double k, Point c)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Scale(k, c);
    }

    SetCenter(center.Scale(k, c));
}

void PlanetA::Move(double dx, double dy, double dz)
{
    Point* points = GetPoints();
    int n = GetNetNodesParam();
    Point center = GetCenter();
    for (int i = 0; i < GetPointsNum(); i++)
    {
        points[i].Move(dx, dy, dz);
    }

    SetCenter(center.Move(dx, dy, dz));
    CalcIntensities(this->lights);
}

double *PlanetA::GetIntensities()
{
    return intensities;
}

void PlanetA::AddOrbitPoint(Point p)
{
    if (orbits_v.size() > 1000)
        orbits_v.clear();
    orbits_v.push_back(Point(p.GetX(), p.GetY(), p.GetZ()));
}

std::vector<Point> PlanetA::GetOrbitNum(void)
{
    std::vector<Point> tmp = orbits_v;
    return tmp;
}

void PlanetA::MashOrbit(double k)
{
    this->orbits_v.clear();
}

