#include "vector.h"

Vector::Vector(Point p1, Point p2)
{
    p.SetX(p2.GetX() - p1.GetX());
    p.SetY(p2.GetY() - p1.GetY());
    p.SetZ(p2.GetZ() - p1.GetZ());
}

void Vector::Set(Point p1, Point p2)
{
    p.SetX(p2.GetX() - p1.GetX());
    p.SetY(p2.GetY() - p1.GetY());
    p.SetZ(p2.GetZ() - p1.GetZ());
}

Point Vector::Get()
{
    return p;
}

double Vector::GetLen()
{
    double x = p.GetX(), y = p.GetY(), z = p.GetZ();
    int len = sqrt(x * x + y * y + z * z);
    return len;
}

void Vector::Neg()
{
    p.SetX(-p.GetX());
    p.SetY(-p.GetY());
    p.SetZ(-p.GetZ());
}

Vector Vector::VectProd(Vector a)
{
    Vector res;
    res.p.SetX(p.GetY() * a.p.GetZ() - p.GetZ() * a.p.GetY());
    res.p.SetY(p.GetZ() * a.p.GetX() - p.GetX() * a.p.GetZ());
    res.p.SetZ(p.GetX() * a.p.GetY() - p.GetY() * a.p.GetX());
    return res;
}

double Vector::ScalarProd(Vector v)
{
    return p.GetX() * v.p.GetX() + p.GetY() * v.p.GetY() + p.GetZ() * v.p.GetZ();
}
