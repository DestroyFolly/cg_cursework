#ifndef MYGRAPHICSVIEW_H
#define MYGRAPHICSVIEW_H

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QWidget>
#include <QGraphicsItem>
#include <QTimer>
#include <QHBoxLayout>
#include <QDebug>
#include "point.h"
#include "planeta.h"
#include "star.h"

#define MAX_OBJECTS_SIZE 100
#define Z_MATRIX_ROW_SIZE 2

//Класс окна графического представления QGraphicsView
class myGraphicsView: public QGraphicsView
{
    Q_OBJECT
public:
    myGraphicsView(QWidget *parent = 0);
    void Start(void);
    void StartRun(void);
    void Stop(void);
    void Clean(void);
    void UpdateStar(double m, int R, bool flag);
    void AddPlanet(Point pos, double Vx, double Vy, double Vz, double m, int R, int R_1, int type, bool flag, QPen p);
    void MinusMash(void);
    void PlusMash(void);
    void DownloadSolarSystem(void);
    void DownloadBistarSystem(void);
    void GenImage(int n);
    void CalcNewPos();
    std::vector<PlanetA> GetPlanets(void) {return this->objects;};
private:
    void CalcStars(void);
    //void CalcNewPos();
    void calculateForce(double xs, double ys, double zs, double m,int i, int i_s, double &gx, double &gy, double &gz);
    void DrawSystem();
    void DrawStars();
    void DrawStar();
    void DrawPlanet(PlanetA &p);
    void DrawPoint(Point p);
    void DrawLine(Point p1, Point p2);
    void DrawPolygon(Point p1, Point p2, Point p3, double k, Point c);
    QGraphicsScene *scene;
    QTimer *myRectTimer;

    QPen pen;
    QBrush brush;

    int n = 7;
    Point star_c = Point(600,400, 0);
    Star Sun = Star(star_c, 0, n, 0);
    std::vector<Point> lights;

    std::vector<PlanetA> objects;
    std::vector<Star> star;
    int angle = 200;
    double k_num = 1.0;
    std::vector<Point> stars;

    double dt = 60*60*24.0;

private slots:
    void myTimerRect();
};

#endif // MYGRAPHICSVIEW_H
