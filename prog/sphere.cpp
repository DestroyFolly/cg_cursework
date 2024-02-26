#include "sphere.h"
#include <iostream>

Sphere::Sphere(Point center, int radius_a, int radius_b, int net_nodes_num, int type)
{
    this->center = center;
    this->radius_a = radius_a;
    this->radius_b = radius_b;
    this->m_radius_a = radius_a;
    this->m_radius_b = radius_b;
    this->net_nodes_param = net_nodes_num;
    this->type = type;
    CreateNet();
}

Sphere Sphere::operator=(Sphere p)
{
    this->center = p.center;
    this->radius_a = p.radius_a;
    this->radius_b = p.radius_b;
    this->net_nodes_param = p.net_nodes_param;
    this->type = p.type;
    this->edges = p.edges;
    this->points = p.points;
    return *this;
}

Sphere::~Sphere()
{
}

void Sphere::CreateNet()
{
    CreatePoints();
    CreateEdges();
}

void Sphere::CreatePoints()
{
    points = new Point [GetPointsNum()];
    int k = 0;
    int r = radius_a;
    int cx = center.GetX(); int cy = center.GetY(); int cz = center.GetZ();

    points[k++] = Point(cx - r, cy, cz);

    int step_degree = 90 / (net_nodes_param + 1);
    double radians, a, b, c;
    for (int i = 1; i <= net_nodes_param; i++)
    {
        radians = (90 - step_degree * i) * M_PI / 180;
        if (type == 0)
        {
            a = sqrt(2 * r * r * (1 - cos(radians)));
            b = sin(radians) * r;
            c = sqrt(a * a - b * b);
            points[k++] = Point(cx - b, cy + r - c, cz);
        }
        else if (type == 1)
        {
            int r_tmp = round(sqrt((radius_a*radius_a*radius_b*radius_b) / (cos(radians)*cos(radians)*radius_a*radius_a+sin(radians)*sin(radians)*radius_b*radius_b)));
            a = sqrt(2 * r_tmp * r_tmp * (1 - cos(radians)));
            b = sin(radians) * r_tmp;
            c = sqrt(a * a - b * b);
            points[k++] = Point(cx - b, cy + radius_a - c, cz);
        }
        //points[k++] = Point(cx - b, cy + r - c, cz);
    }
    points[k++] = Point(cx, cy + r, cz);

    for (int i = 0; i <= net_nodes_param; i++)
    {
        points[k++] = Point(2 * cx - points[net_nodes_param - i].GetX(),
                          points[net_nodes_param - i].GetY(),
                          points[net_nodes_param - i].GetZ());
    }
}

int Sphere::GetPointsNum()
{
    int n = net_nodes_param;
    return (3 + n * 2) + (n * 2 + 1) * 2 * (360 / 90 * (n + 1));
}

void Sphere::CreateEdges()
{
    int n = 3 + net_nodes_param * 2;
    int m = 2 + net_nodes_param * 2 * 2;

    Edge tempEdges[m];
    edges = new Edge[GetEdgesNum()];

    int edgeIndex = 0;
    int topIndices[n];
    int bottomIndices[n];
    double degreeStep = 90 / (net_nodes_param + 1);
    int k = n;
    int currentIndex, rightIndex, leftIndex;

    for (int i = 0; i < n; i++)
    {
        topIndices[i] = i;
    }

    int s = n;

    for (double degree = degreeStep; degree <= 360; degree += degreeStep)
    {
        for (int i = 1; i < n - 1; i++)
        {
            Point tempPoint = points[i];
            points[s++] = tempPoint.RotateX(degree, center);
        }
        for (int i = 0; i < n - 2; i++)
        {
            bottomIndices[i] = k + i;
        }
        for (int i = 0; i < m; i++)
        {
            tempEdges[i] = Edge();
        }

        tempEdges[0].Append(topIndices[0]);
        tempEdges[n - 2].Append(topIndices[n - 1]);

        for (int i = 1; i < n - 1; i++)
        {
            tempEdges[i - 1].Append(topIndices[i]);
            tempEdges[i].Append(topIndices[i]);
        }

        for (int i = 2; i < n - 1; i++)
        {
            tempEdges[n - 1 + i - 2].Append(topIndices[i]);
        }
        for (int i = 0; i < n - 2; i++)
        {
            currentIndex = (n + i) % (n - 1);
            rightIndex = (n - 1) + i;
            leftIndex = rightIndex - 1;
            if (currentIndex == 1)
            {
                leftIndex = 0;
            }
            tempEdges[currentIndex].Append(bottomIndices[i]);
            tempEdges[leftIndex].Append(bottomIndices[i]);
            if (currentIndex != n - 2)
            {
                tempEdges[rightIndex].Append(bottomIndices[i]);
            }
        }
        for (int i = edgeIndex; i < edgeIndex + m; i++)
        {
            edges[i] = tempEdges[i % m];
        }
        edgeIndex += m;
        topIndices[0] = 0;
        topIndices[n - 1] = n - 1;
        for (int i = 0; i < n - 2; i++)
        {
            topIndices[i + 1] = bottomIndices[i];
        }
        k += n - 2;
    }

    for (int i = 0; i < n - 2; i++)
    {
        bottomIndices[i] = topIndices[i + 1];
    }
    for (int i = 0; i < n; i++)
    {
        topIndices[i] = i;
    }
    for (int i = 0; i < m; i++)
    {
        tempEdges[i] = Edge();
    }

    tempEdges[0].Append(topIndices[0]);
    tempEdges[n - 2].Append(topIndices[n - 1]);

    for (int i = 1; i < n - 1; i++)
    {
        tempEdges[i - 1].Append(topIndices[i]);
        tempEdges[i].Append(topIndices[i]);
    }

    for (int i = 2; i < n - 1; i++)
    {
        tempEdges[n - 1 + i - 2].Append(topIndices[i]);
    }
    for (int i = 0; i < n - 2; i++)
    {
        currentIndex = (n + i) % (n - 1);
        rightIndex = (n - 1) + i;
        leftIndex = rightIndex - 1;
        if (currentIndex == 1)
        {
            leftIndex = 0;
        }
        tempEdges[currentIndex].Append(bottomIndices[i]);
        tempEdges[leftIndex].Append(bottomIndices[i]);
        if (currentIndex != n - 2)
        {
            tempEdges[rightIndex].Append(bottomIndices[i]);
        }
    }
    for (int i = edgeIndex; i < edgeIndex + m; i++)
    {
        edges[i] = tempEdges[i % m];
    }

}

int Sphere::GetEdgesNum()
{
    int n = net_nodes_param;
    int num = (n * 2 * 2 + 2) * (360 / 90 * (n + 1) + 2);
    return num;
}

void Sphere::RotateX(int degree, Point c)
{
    int n = GetPointsNum();
    for (int i = 0; i < n; i++)
    {
        points[i].RotateX(degree, c);
    }
    center.RotateX(degree, c);
}

void Sphere::RotateY(int degree, Point c)
{
    int n = GetPointsNum();
    for (int i = 0; i < n; i++)
    {
        points[i].RotateY(degree, c);
    }
    center.RotateY(degree, c);
}

void Sphere::RotateZ(int degree, Point c)
{
    int n = GetPointsNum();
    for (int i = 0; i < n; i++)
    {
        points[i].RotateZ(degree, c);
    }
    center.RotateZ(degree, c);
}

void Sphere::Scale(double kx, double ky, double kz, Point c)
{
    int n = GetPointsNum();
    for (int i = 0; i < n; i++)
    {
        points[i].Scale(kx, ky, kz, c.GetX(), c.GetY(), c.GetZ());
    }
    center.Scale(kx, ky, kz, c.GetX(), c.GetY(), c.GetZ());
}

void Sphere::Move(double dx, double dy, double dz)
{
    int n = GetPointsNum();
    for (int i = 0; i < n; i++)
    {
        points[i].Move(dx, dy, dz);
    }
    center.Move(dx, dy, dz);
}

void Sphere::ChangeR(double k)
{
    this->radius_a = this->m_radius_a * k;
    this->radius_b = this->m_radius_b * k;
    delete []edges;
    delete []points;
    CreateNet();
}

void Sphere::ChangeR(int r_1, int r_2)
{
    this->radius_a = r_1;
    this->radius_b = r_2;
    delete []edges;
    delete []points;
    CreateNet();
}
