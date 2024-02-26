#include "mainwindow.h"

#include <QApplication>
#include <getopt.h>
#include <iostream>
#include <unistd.h>
#include <fstream>


int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        QApplication a(argc, argv);
        MainWindow w;
        w.show();
        return a.exec();
    }

    return 0;
}
