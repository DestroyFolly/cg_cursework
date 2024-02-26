#include "testg.h"
#include <QtTest>
#include <vector>
#include <iostream>

TestG::TestG(QObject* parent, myGraphicsView* g)
    : QObject { parent }
{
    this->m_g = g;
}

void TestG::test_g_01()
{
    this->m_g->UpdateStar(2e30, 20, false);
    this->m_g->AddPlanet(Point(1.5e11,0,0), 0,0,30000, 5e24, 10, 10, 0, false, QPen(Qt::yellow));

    this->m_g->CalcNewPos();

    std::vector<PlanetA> planets = this->m_g->GetPlanets();

    Point new_pos = planets[0].GetPos();

    QCOMPARE(qFuzzyCompare(new_pos.GetZ(), 2.592e+09), true);
}

void TestG::test_g_02()
{
    this->m_g->UpdateStar(2e30, 20, false);
    this->m_g->AddPlanet(Point(1.5e11,0,0), 0,0,30, 10e24, 10, 10, 0, false, QPen(Qt::yellow));

    this->m_g->CalcNewPos();

    std::vector<PlanetA> planets = this->m_g->GetPlanets();

    Point new_pos = planets[0].GetPos();

    QCOMPARE(qFuzzyCompare(new_pos.GetY(), 0.0), true);
}

void TestG::test_g_03()
{
    Star s (Point(600,400,5), 10, 0, 2e30) ;
    Point * test = s.GetPoints();
    std::vector<Point> test_tmp = {Point(590, 400, 5),
                                   Point(600, 410, 5),
                                   Point(610, 400, 5),
                                   Point(600, 400, 15),
                                   Point(600, 390, 5),
                                   Point(600, 400, -5),
                                   Point(600, 410, 5),
                                   Point(0, 0, 0),
                                   Point(0, 0, 0),
                                   Point(0, 0, 0),
                                   Point(0, 0, 0)};
    for(int i = 0; i < s.GetPointsNum(); i++)
    {
        QCOMPARE(qFuzzyCompare(test[i].GetX(), test_tmp[i].GetX()), true);
        QCOMPARE(qFuzzyCompare(test[i].GetY(), test_tmp[i].GetY()), true);
        QCOMPARE(qFuzzyCompare(test[i].GetZ(), test_tmp[i].GetZ()), true);
    }
}

void TestG::test_g_04()
{
    Star s (Point(600,400,5), 10, 0, 2e30) ;
    std::vector<Point> l = s.GetLights();
    PlanetA p (10, 10, 0, l, QColor(Qt::blue), 0);
    p.SetNewPos(Point(500,500, 20));
    Point *test = p.GetPoints();

    std::vector<Point> test_tmp = {Point(490, 500, 20),
                                   Point(500, 510, 20),
                                   Point(510, 500, 20),
                                   Point(500, 500, 30),
                                   Point(500, 490, 20),
                                   Point(500, 500, 10),
                                   Point(500, 510, 20),
                                   Point(500, 500, 20),
                                   Point(500, 500, 20),
                                   Point(500, 500, 20),
                                   Point(500, 500, 20)};
    for(int i = 0; i < s.GetPointsNum(); i++)
    {
        QCOMPARE(qFuzzyCompare(test[i].GetX(), test_tmp[i].GetX()), true);
        QCOMPARE(qFuzzyCompare(test[i].GetY(), test_tmp[i].GetY()), true);
        QCOMPARE(qFuzzyCompare(test[i].GetZ(), test_tmp[i].GetZ()), true);
    }
}

