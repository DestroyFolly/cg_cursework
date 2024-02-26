#ifndef SPHERE_H
#define SPHERE_H

#include "point.h"
#include "edge.h"
#include <vector>

class Sphere
{
public:
    Sphere(Point center, int radius_a, int radius_b, int net_nodes_param = 0, int type = 0);
    ~Sphere();

    Point GetCenter(void) {return this->center;};
    int GetRadius(void) {return this->radius_a;};
    int GetNetNodesParam(void) {return this->net_nodes_param;};

    int GetPointsNum();
    int GetEdgesNum();
    Sphere operator = (Sphere p);

    Point* GetPoints() {return points;};
    Edge* GetEdges() {return edges;};

    void SetCenter(Point center) {this->center = center;};

    void RotateX(int degree, Point c);
    void RotateY(int degree, Point c);
    void RotateZ(int degree, Point c);

    void Scale(double kx, double ky, double kz, Point c);

    void Move(double dx, double dy, double dz);

    void ChangeR(double k);
    void ChangeR(int r_1, int r_2);

private:
    Point center;
    int radius_a;
    int radius_b;
    int m_radius_a;
    int m_radius_b;

    int net_nodes_param;
    int type;

    Edge* edges;
    Point* points;
    void CreateNet(void);
    void CreatePoints(void);
    void CreateEdges(void);
};

#endif // SPHERE_H
