#ifndef EDGE_H
#define EDGE_H


class Edge
{
public:
    Edge(int P_1 = 0, int P_2 = 0, int P_3 = 0): p_1(P_1), p_2(P_2), p_3(P_3) {};

    void Append(int P);

    int GetP1(void) {return p_1;};
    int GetP2(void) {return p_2;};
    int GetP3(void) {return p_3;};

    void SetP1(int p) {this->p_1 = p;};
    void SetP2(int p) {this->p_2 = p;};
    void SetP3(int p) {this->p_3 = p;};

private:
    int p_1;
    int p_2;
    int p_3;

    int i = 0;
};

#endif // EDGE_H
