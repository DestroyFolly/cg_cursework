#ifndef VECTOR_H
#define VECTOR_H

#include "point.h"

class Vector: public Point
{
public:
    Vector(): p(0,0,0) {};
    Vector(Point p1, Point p2);

    void Set(Point p1, Point p2);
    Point Get();
    double GetLen();

    void Neg();

    Vector VectProd(Vector a);
    double ScalarProd(Vector v);
private:
    Point p;
};

#endif // VECTOR_H
