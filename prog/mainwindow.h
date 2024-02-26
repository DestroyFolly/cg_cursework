#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "mygraphicsview.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE


//Класс окна QWidget
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void genImage(int n);
    myGraphicsView *GetView(void) {return this->myGraphics;}

private slots:
    void on_stop_clicked();

    void on_start_clicked();

    void on_updateStar_clicked();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_action_2_triggered();

    void on_action_3_triggered();

    void on_full_end_clicked();

    void on_action_triggered();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_color_p_clicked();

    void color_dialog(QColor& color);

private:
    Ui::MainWindow *ui;
    myGraphicsView *myGraphics;
    QHBoxLayout *myLayoutScene;

    QPen pen;
    QBrush brush;
    bool star_f = false;
};
#endif // MAINWINDOW_H
