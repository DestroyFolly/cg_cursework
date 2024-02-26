#ifndef POINT_H
#define POINT_H

#include <cmath>

class Point
{
public:
    Point(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z) {};

    double GetX(void) {return this->x;};
    double GetY(void) {return this->y;};
    double GetZ(void) {return this->z;};

    void SetX(double x) {this->x = x;};
    void SetY(double y) {this->y = y;};
    void SetZ(double z) {this->z = z;};

    Point Get(void) {return *this;};
    double GetDistance(Point p);

    bool operator == (Point p) const;
    Point operator - (Point p);

    Point Scale(double kx, double ky, double kz, double xc, double yc, double zc);
    Point Scale(double kx, double ky, double kz, Point center);
    Point Scale(double k, Point center);

    Point Move(double dx, double dy, double dz);

    Point RotateX(int degree, Point center);
    Point RotateY(int degree, Point center);
    Point RotateZ(int degree, Point center);

private:
    double x;
    double y;
    double z;
};

#endif // POINT_H
