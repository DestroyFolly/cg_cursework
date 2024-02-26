#ifndef STAR_H
#define STAR_H

#include "sphere.h"
#include "vector.h"

class Star: public Sphere
{
public:
    Star(Point center, int radius, int net_nodes_num, double w): Sphere(center, radius, radius, net_nodes_num), m(w) {calc_lights();};
    ~Star();
    std::vector<Point> GetLights(void) {return lights;};
    void SetWeight(double w) {m = w;};
    void SetV(double vx, double vy, double vz) {m_vx = vx; m_vy = vy; m_vz = vz;}
    double GetVx(void) {return m_vx;}
    double GetVy(void) {return m_vy;}
    double GetVz(void) {return m_vz;}

    Point GetPos(void) {return p_center;}
    void SetPos(Point p) {p_center = p;}
    void SetNewPos(Point n) {Point new_p = n - this->GetCenter(); this->Move(new_p.GetX(), new_p.GetY(), new_p.GetZ());};

    double GetWeight(void) {return m;};
    Star operator=(Star p);

    void AddOrbitPoint(Point p);
    std::vector<Point> GetOrbitNum(void);
    void MashOrbit(double k);
private:
    std::vector<Point> lights;
    double m;
    double m_vx = 0;
    double m_vy = 0;
    double m_vz = 0;
    std::vector<Point> orbits_v;
    Point p_center = Point(0,0,0);
    void calc_lights();
};

#endif // STAR_H
