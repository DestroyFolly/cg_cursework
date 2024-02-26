#include "mygraphicsview.h"
#include <iostream>
#include <QDir>
#include <QMessageBox>

myGraphicsView::myGraphicsView(QWidget *parent) : QGraphicsView(parent)
{
    scene = new QGraphicsScene(0,0, 1000, 800, this);
    scene->addRect(scene->sceneRect());
    setScene(scene);
    pen = QPen(QColor(255,140,0), 2);
    this->setBackgroundBrush(QBrush(Qt::black));
    brush = QBrush(Qt::yellow);

    CalcStars();
    myRectTimer = new QTimer(this);
    connect(myRectTimer, SIGNAL(timeout()), this, SLOT(myTimerRect()));
}

void myGraphicsView::Start(void)
{
    this->dt = 60*60*24.0;
    myRectTimer->start(1000/50);
}

void myGraphicsView::StartRun(void)
{
    this->dt *= 5.0;
    myRectTimer->start(1000/50);
}

void myGraphicsView::Stop()
{
    myRectTimer->stop();
}

void myGraphicsView::Clean()
{
    myRectTimer->stop();
    this->objects.clear();
    this->star.clear();
    this->scene->clear();
    this->DrawStars();
}

void myGraphicsView::UpdateStar(double m, int R, bool flag)
{
    Point c = Point(600, 400, 0);
    if (star.size() == 0)
        this->star.push_back(Star(c, R, 7, m));
    else
    {
        this->star[0].ChangeR(R, R);
        this->star[0].SetWeight(m);
        this->star[0].ChangeR(this->k_num);
        this->star[0].SetPos(Point(0,0,0));
    }
    this->star[0].SetWeight(m);
    this->lights = star[0].GetLights();
    star[0].RotateX(angle, star[0].GetCenter());
    star[0].ChangeR(this->k_num);
    if (flag)
        DrawSystem();
}

void myGraphicsView::AddPlanet(Point pos, double Vx, double Vy, double Vz, double m, int R, int R_1, int type, bool flag, QPen p)
{
    this->objects.push_back(PlanetA(R, R_1, 7, lights, p.color(), type));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(m, Vx, Vy, Vz);
    this->objects[this->objects.size() - 1].SetPos(pos);
    if (flag)
    {
        this->CalcNewPos();
        DrawSystem();
        std::cout << this->objects[this->objects.size() - 1].GetCenter().GetX() << std::endl;
        std::cout << this->objects[this->objects.size() - 1].GetCenter().GetY() << std::endl;
        std::cout << this->objects[this->objects.size() - 1].GetCenter().GetZ() << std::endl;
        if (this->objects[this->objects.size() - 1].GetCenter().GetX() < 0 || this->objects[this->objects.size() - 1].GetCenter().GetX() > 1000)
        {
            QMessageBox::critical(NULL, "Упс", "Система с такими данными не существует");
        }
    }
}

void myGraphicsView::MinusMash()
{
    this->k_num -=0.05;
    star[0].ChangeR(this->k_num);
    for(int i = 0; i < objects.size(); i++)
    {
        objects[i].ChangeR(this->k_num);
        objects[i].MashOrbit(this->k_num);
    }

    CalcNewPos();
    DrawSystem();
}

void myGraphicsView::PlusMash()
{
    this->k_num +=0.05;
    star[0].ChangeR(this->k_num);
    for(int i = 0; i < objects.size(); i++)
    {
        objects[i].ChangeR(this->k_num);
        objects[i].MashOrbit(this->k_num);
    }
    CalcNewPos();
    DrawSystem();
}

void myGraphicsView::DownloadSolarSystem()
{
    this->objects.clear();
    this->k_num = 1.0;
    this->star.clear();
    // Sun
    Point c = Point(600, 400, 0);
    this->star.push_back(Star(c, 40, 7, 2e30));
    this->lights = star[0].GetLights();
    star[0].RotateX(angle, star[0].GetCenter());
    star[0].ChangeR(this->k_num);
    // Mercury
    this->objects.push_back(PlanetA(10, 10, 5, lights, QColor(220, 92, 0), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(3.33e23, 0, 0, 38700);
    this->objects[this->objects.size() - 1].SetPos(Point(6.9817e10, 0, 0));
    // Venus
    this->objects.push_back(PlanetA(12, 12, 5, lights, QColor(202, 194, 0), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(4.8675e24, 0, 0, 35020);
    this->objects[this->objects.size() - 1].SetPos(Point(1.0894e11, 0, 0));
    //Earth
    this->objects.push_back(PlanetA(14, 14, 5, lights, QColor(8, 169, 227), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(5.972e24, 0, 0, 29790);
    this->objects[this->objects.size() - 1].SetPos(Point(1.52505e11, 0, 0));
    //Mars
    this->objects.push_back(PlanetA(12, 12, 5, lights, QColor(243, 49, 5), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(6.4171e23, 0, 0, 24007);
    this->objects[this->objects.size() - 1].SetPos(Point(2.49232e11, 0, 0));
    //Jupiter
    this->objects.push_back(PlanetA(30, 30, 5, lights, QColor(179, 143, 123), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(1.8986e27, 0, 0, 16070);
    this->objects[this->objects.size() - 1].SetPos(Point(3.785e11, 0, 0));
    //Saturn
    this->objects.push_back(PlanetA(20, 20, 5, lights, QColor(255, 228, 116), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(0.5986e26, 0, 0, 13070);
    this->objects[this->objects.size() - 1].SetPos(Point(4.785e11, 0, 0));
    //Uranus
    this->objects.push_back(PlanetA(20, 20, 5, lights, QColor(173, 213, 230), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(5.8936e25, 0, 0, 10070);
    this->objects[this->objects.size() - 1].SetPos(Point(5.785e11, 0, 0));
    //Neptune
    this->objects.push_back(PlanetA(20, 20, 5, lights, QColor(26, 73, 111), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(2.8986e25, 0, 0, 9070);
    this->objects[this->objects.size() - 1].SetPos(Point(6.785e11, 0, 0));
    this->CalcNewPos();
    DrawSystem();
}



void myGraphicsView::GenImage(int n)
{
    this->DownloadSolarSystem();

    for (int i = 0; i < n; i++)
    {
        this->CalcNewPos();
        this->DrawSystem();
        const QPixmap pixmap = this->grab();
        const QString fileName(QString::fromStdString(std::to_string(i)) +".png");
        qDebug() << pixmap.save(fileName);
    }
}

void myGraphicsView::DownloadBistarSystem()
{
    this->objects.clear();
    this->k_num = 1.0;
    this->star.clear();
    // Sun
    Point c = Point(600, 400, 0);
    this->star.push_back(Star(c, 40, 7, 2e30));
    this->lights = star[0].GetLights();
    star[0].RotateX(angle, star[0].GetCenter());
    star[0].ChangeR(this->k_num);
    //Earth
    this->objects.push_back(PlanetA(14, 14, 7, lights, QColor(8, 169, 227), 0));
    this->objects[this->objects.size() - 1].SetPhisicalNumb(5.972e24, 0, 0, 29790);
    this->objects[this->objects.size() - 1].SetPos(Point(1.52505e11, 0, 0));
    this->CalcNewPos();
    DrawSystem();
}


void myGraphicsView::myTimerRect()
{
    this->CalcNewPos();
    DrawSystem();
}

void myGraphicsView::DrawPlanet(PlanetA &p)
{
    pen = QPen(Qt::yellow);
    brush = QBrush(p.GetSurfaceColor());
    Edge *net = p.GetEdges();
    Point* points = p.GetPoints();
    double *intenses = p.GetIntensities();

    int n = p.GetNetNodesParam();
    int cz = p.GetCenter().GetZ();
    Point c = p.GetCenter();

    Point p1, p2, p3;
    Point sc = star[0].GetCenter();
    double sr = star[0].GetRadius();

    std::vector<Point> orbit = p.GetOrbitNum();
    int no = orbit.size();
    int j;
    for (int i = 0; i < no - 1; i++)
    {
        DrawLine(orbit[i], orbit[i+1]);
    }

    int r, g, b;
    double t = p.GetCenter().GetZ() / 450;
    for (int i = 0; i < p.GetEdgesNum(); i++)
    {
        p.GetSurfaceColor().getRgb(&r, &g, &b);
        double k = intenses[i] / 1 / ((n * 2 + 1) * 2 * (n + 1) * 4 + 1);
        r += k * 255;
        if (r > 255) r = 255;
        if (r < 0) r = 0;
        g += k * 255;
        if (g > 255) g = 255;
        if (g < 0) g = 0;
        brush = QBrush(QColor(r, g, b));
        pen = QPen(QColor(r, g, b), 1);

        p1 = points[net[i].GetP1()];
        p2 = points[net[i].GetP2()];
        p3 = points[net[i].GetP3()];

        if (p1.GetZ() >= cz || p2.GetZ() >= cz || p3.GetZ() >= cz)
            DrawPolygon(p1, p2, p3, 1 + t, c);
    }
}


void myGraphicsView::DrawPoint(Point p)
{
    this->scene->addLine(p.GetX(), 800 - p.GetY(), p.GetX(), 800 - p.GetY(), this->pen);
}

void myGraphicsView::DrawLine(Point p1, Point p2)
{
    this->scene->addLine(p1.GetX(), 800 - p1.GetY(), p2.GetX(), 800 - p2.GetY(),this->pen);
}


void myGraphicsView::CalcStars()
{
    for (int i = 0; i < 1200; i += 1 + std::rand()/((RAND_MAX + 1u)/60))
    {
        int k = 1 + std::rand()/((RAND_MAX + 1u)/15);
        for (int j = 0; j < 800; j += 1 + std::rand()/((RAND_MAX + 1u)/60))
        {
            k *= -1;
            stars.push_back(Point(i+k, j ,0.0));
            scene->addLine(i + k, j, i + k, j, QPen(Qt::white, 1));
        }
    }
}

void myGraphicsView::DrawStars()
{
    for (int i = 0; i < stars.size(); i ++)
    {
        scene->addLine(stars[i].GetX(), stars[i].GetY(), stars[i].GetX(), stars[i].GetY(), QPen(Qt::white, 1));
    }
}

void myGraphicsView::DrawStar()
{
    pen = QPen(Qt::yellow);
    brush = QBrush(Qt::yellow);
    Edge* net = star[0].GetEdges();
    Point* points = star[0].GetPoints();
    int cz = star[0].GetCenter().GetZ();
    Point p1, p2, p3;
    double t = star[0].GetCenter().GetZ() / 450;

    std::vector<Point> orbit = star[0].GetOrbitNum();
    int no = orbit.size();
    for (int i = 0; i < no - 1; i++)
    {
        DrawLine(orbit[i], orbit[i+1]);
    }
    for (int i = 0; i < star[0].GetEdgesNum(); i++)
    {
        p1 = points[net[i].GetP1()];
        p2 = points[net[i].GetP2()];
        p3 = points[net[i].GetP3()];
        if (p1.GetZ() >= cz || p2.GetZ() >= cz || p3.GetZ() >= cz)
            DrawPolygon(p1, p2, p3, 1 + t, star[0].GetCenter());
    }
}

void myGraphicsView::DrawPolygon(Point p1, Point p2, Point p3, double k, Point c)
{
    QPolygonF polygon;

    p1.Scale(k, c);
    p2.Scale(k, c);
    p3.Scale(k, c);

    polygon << QPointF(p1.GetX(), 800 - p1.GetY());
    polygon << QPointF(p2.GetX(), 800 - p2.GetY());
    polygon << QPointF(p3.GetX(), 800 - p3.GetY());

    this->scene->addPolygon(polygon, pen, brush);
}

void myGraphicsView::DrawSystem()
{
    scene->clear();
    DrawStars();

    int z[MAX_OBJECTS_SIZE][Z_MATRIX_ROW_SIZE];
    for (int i = 0; i < 8; i++)
    {
        z[i][0] = (objects[i]).GetCenter().GetZ();
        z[i][1] = i;
    }

    if (objects.size() > 0)
    {
        for (int i = 0; i < 8 - 1; i++)
            for (int j = 0; j < objects.size() - i - 1; j++)
                if (z[j][0] > z[j + 1][0])
                {
                    std::swap(z[j][0], z[j + 1][0]);
                    std::swap(z[j][1], z[j + 1][1]);
                }
    }
    int i = 0;

    while (i < objects.size() && z[i][0] < star[0].GetCenter().GetZ())
        DrawPlanet(objects[z[i++][1]]);
    DrawStar();
    while (i < objects.size())
        DrawPlanet(objects[z[i++][1]]);
}

void myGraphicsView::calculateForce(double xCoord, double yCoord, double zCoord, double mass, int objIndex, int starIndex, double &forceX, double &forceY, double &forceZ)
{
    double gravityConst = 6.67e-11;
    forceX = 0;
    forceY = 0;
    forceZ = 0;

    for(int index = 0; index < objects.size(); index++)
    {
        if (index == objIndex)
            continue;
        double distance = pow(pow(xCoord -  objects[index].GetPos().GetX(), 2) + pow(yCoord -  objects[index].GetPos().GetY(), 2) + pow(zCoord -  objects[index].GetPos().GetZ(), 2), 1.5);
        double dgv = -(gravityConst*objects[index].GetWeight() * mass);
        double deltaX = ((xCoord - objects[index].GetPos().GetX()) * (dgv)) / distance;
        double deltaY = ((yCoord - objects[index].GetPos().GetY()) * (dgv))/distance;
        double deltaZ = ((zCoord - objects[index].GetPos().GetZ()) * (dgv))/distance;

        forceX += deltaX;
        forceY += deltaY;
        forceZ += deltaZ;
    }
    for(int index = 0; index < star.size(); index++)
    {
        if (index == starIndex)
            continue;
        double distance = pow(pow(xCoord - star[index].GetCenter().GetX(), 2) + pow(yCoord - star[index].GetCenter().GetY(), 2) + pow(zCoord - star[index].GetCenter().GetZ(), 2), 1.5);
        double dgv = -(gravityConst*star[index].GetWeight() * mass);
        double deltaX = ((xCoord - star[index].GetPos().GetX()) * (dgv)) / distance;
        double deltaY = ((yCoord - star[index].GetPos().GetY()) * (dgv))/distance;
        double deltaZ = ((zCoord - star[index].GetPos().GetZ()) * (dgv))/distance;

        forceX += deltaX;
        forceY += deltaY;
        forceZ += deltaZ;
    }

}

void myGraphicsView::CalcNewPos(void)
{
    std::vector<Point> buf_pos;
    std::vector<Point> buf_speed;
    std::vector<Point> buf_pos_star;
    std::vector<Point> buf_speed_star;
    for(int i = 0; i < objects.size(); i++)
    {
        double xs = objects[i].GetPos().GetX();
        double ys = objects[i].GetPos().GetY();
        double zs = objects[i].GetPos().GetZ();
        double m = objects[i].GetWeight();
        double vx = objects[i].GetVx();
        double vy = objects[i].GetVy();
        double vz = objects[i].GetVz();
        double gx, gy, gz;
        calculateForce(xs, ys, zs, m, i, -1, gx, gy, gz);
        vx += (gx*dt)/m;
        vy += (gy*dt)/m;
        vz += (gz*dt)/m;

        xs += vx*dt;
        ys += vy*dt;
        zs += vz*dt;
        buf_pos.push_back(Point(xs, ys, zs));
        buf_speed.push_back(Point(vx, vy, vz));
    }
    for(int i = 0; i < star.size(); i++)
    {
        double xs = star[i].GetPos().GetX();
        double ys = star[i].GetPos().GetY();
        double zs = star[i].GetPos().GetZ();
        double m = star[i].GetWeight();
        double vx = star[i].GetVx();
        double vy = star[i].GetVy();
        double vz = star[i].GetVz();
        double gx, gy, gz;
        calculateForce(xs, ys, zs, m, -1, i, gx, gy, gz);
        vx += (gx*dt)/m;
        vy += (gy*dt)/m;
        vz += (gz*dt)/m;

        xs += vx*dt;
        ys += vy*dt;
        zs += vz*dt;
        buf_pos_star.push_back(Point(xs, ys, zs));
        buf_speed_star.push_back(Point(vx, vy, vz));
    }
    for (int k = 0; k < objects.size(); k++)
    {
        Point cur_c ((buf_pos[k].GetX()/328348571.4285) * k_num + star_c.GetX(),
                     (buf_pos[k].GetY()/ 328348571.4285) * k_num + star_c.GetY(),
                     (buf_pos[k].GetZ()/328348571.4285) * k_num+ star_c.GetZ());
        cur_c.RotateX(angle, star_c);
        objects[k].SetNewPos(cur_c);
        objects[k].SetPos(Point(buf_pos[k].GetX(), buf_pos[k].GetY(), buf_pos[k].GetZ()));
        objects[k].AddOrbitPoint(cur_c);
        objects[k].SetPhisicalNumb(objects[k].GetWeight(), buf_speed[k].GetX(), buf_speed[k].GetY(), buf_speed[k].GetZ());
    }
    for (int k = 0; k < star.size(); k++)
    {
        Point cur_c ((buf_pos_star[k].GetX()/328348571.4285) * k_num + star_c.GetX(),
                     (buf_pos_star[k].GetY()/ 328348571.4285) * k_num + star_c.GetY(),
                     (buf_pos_star[k].GetZ()/328348571.4285) * k_num + star_c.GetZ());
        cur_c.RotateX(angle, star_c);
        star[k].SetNewPos(cur_c);
        star[k].SetPos(Point(buf_pos_star[k].GetX(), buf_pos_star[k].GetY(), buf_pos_star[k].GetZ()));
        star[k].AddOrbitPoint(cur_c);
        star[k].SetV(buf_speed_star[k].GetX(), buf_speed_star[k].GetY(), buf_speed_star[k].GetZ());
    }
}
