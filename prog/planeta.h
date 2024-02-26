#ifndef PLANETA_H
#define PLANETA_H

#include "sphere.h"
#include <QColor>
#include "vector.h"
#include <vector>


class PlanetA: public Sphere
{
public:
    PlanetA(int radius_a, int radius_b, int net_nodes_num, std::vector<Point> lights, QColor surface_color, int type);

    QColor GetSurfaceColor() {return surface_color;};

    void RotateX(int degree, Point center);
    void RotateY(int degree, Point center);
    void RotateZ(int degree, Point center);
    void Scale(double kx, double ky, double kz, Point c);
    void Scale(double k, Point c);
    void Move(double dx, double dy, double dz);

    double *GetIntensities();

    void SetPhisicalNumb(double w, double vx_n, double vy_n, double vz_n){weight = w; vx = vx_n; vy = vy_n; vz = vz_n;};
    void SetNewPos(Point n) {Point new_p = n - this->GetCenter(); this->Move(new_p.GetX(), new_p.GetY(), new_p.GetZ());};
    double GetWeight(void) {return this->weight;}
    double GetDist(void) {return this->dist;}
    double GetVx(void) {return vx;}
    double GetVy(void) {return vy;}
    double GetVz(void) {return vz;}
    Point GetPos(void) {return p_center;}
    void SetPos(Point p) {p_center = p;};

    void AddOrbitPoint(Point p);
    std::vector<Point> GetOrbitNum(void);
    void MashOrbit(double k);

private:
    Point CalcCenter(int dist, int a_xy, int a_yz, Point star_center);
    void CalcIntensities(std::vector<Point> lights);

    QColor surface_color;
    double *intensities = new double[GetEdgesNum()];
    std::vector<Point> orbits_v;
    std::vector<Point> lights;

    double weight = 0;
    double dist = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    Point p_center;
};

#endif // PLANETA_H
