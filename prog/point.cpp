#include "point.h"

double Point::GetDistance(Point p)
{
    return sqrt((x - p.GetX()) * (x - p.GetX()) + (y - p.GetY()) * (y - p.GetY()) + (z - p.GetZ()) * (z - p.GetZ()));
}

bool Point::operator ==(Point p) const
{
    return p.GetX() == x && p.GetY() == y && p.GetZ() == z;
}

Point Point::operator-(Point p)
{
    return Point(this->GetX() - p.GetX(), this->GetY() - p.GetY(), this->GetZ() - p.GetZ());
}

Point Point::Scale(double kx, double ky, double kz,  double xc, double yc, double zc)
{
    x = (x - xc) * kx + xc;
    y = (y - yc) * ky + yc;
    z = (z - zc) * kz + zc;
    return *this;
}

Point Point::Scale(double kx, double ky, double kz, Point center)
{
    x = (x - center.GetX()) * kx + center.GetX();
    y = (y - center.GetY()) * ky + center.GetY();
    z = (z - center.GetZ()) * kz + center.GetZ();
    return *this;
}

Point Point::Scale(double k, Point center)
{
    x = (x - center.GetX()) * k + center.GetX();
    y = (y - center.GetY()) * k + center.GetY();
    z = (z - center.GetZ()) * k + center.GetZ();
    return *this;
}

Point Point::Move(double dx, double dy, double dz)
{
    this->x += dx;
    this->y += dy;
    this->z += dz;
    return *this;
}

Point Point::RotateX(int degree, Point center)
{
    double radians = degree * M_PI / 180;
    Move(-center.x, -center.y, -center.z);
    double temp = y;
    y = cos(radians) * y - sin(radians) * z;
    z = cos(radians) * z + sin(radians) * temp;
    Move(center.x, center.y, center.z);
    return *this;
}

Point Point::RotateY(int degree, Point center)
{
    double radians = degree * M_PI / 180;
    Move(-center.x, -center.y, -center.z);
    double temp = x;
    x = cos(radians) * x - sin(radians) * z;
    z = cos(radians) * z + sin(radians) * temp;
    Move(center.x, center.y, center.z);
    return *this;
}

Point Point::RotateZ(int degree, Point center)
{
    double radians = degree * M_PI / 180;
    Move(-center.x, -center.y, -center.z);
    double temp = x;
    x = cos(radians) * x - sin(radians) * y;
    y = cos(radians) * y + sin(radians) * temp;
    Move(center.x, center.y, center.z);
    return *this;
}
