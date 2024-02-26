#include "edge.h"

void Edge::Append(int P)
{
    if (i == 0)
        this->p_1 = P;
    else if (i == 1)
        this->p_2 = P;
    else if (i == 2)
        this->p_3 = P;
    i++;
}
