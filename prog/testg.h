#ifndef TESTG_H
#define TESTG_H


#include <QObject>

#include "mygraphicsview.h"

class TestG : public QObject {
    Q_OBJECT
public:
    explicit TestG(QObject* parent = nullptr, myGraphicsView* g = nullptr);

private:
    myGraphicsView* m_g;
private slots:
    void test_g_01();
    void test_g_02();
    void test_g_03();
    void test_g_04();
};



#endif // TESTG_H
