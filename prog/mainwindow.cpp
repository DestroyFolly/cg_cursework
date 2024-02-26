#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <QPixmap>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <QMessageBox>
#include <QColorDialog>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    myLayoutScene = new QHBoxLayout(this);
    myGraphics = new myGraphicsView;
    myLayoutScene->addWidget(myGraphics);
    QWidget *placeholderWidget = new QWidget;
    placeholderWidget->setLayout(myLayoutScene);
    ui->widget->setLayout(myLayoutScene);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::genImage(int n)
{
    this->myGraphics->GenImage(n);
}

void MainWindow::on_stop_clicked()
{
    if (!star_f)
        QMessageBox::critical(NULL, "Ошибка", "Добавьте объекты на сцену");
    else
        myGraphics->Stop();
}


void MainWindow::on_start_clicked()
{
    if (!star_f)
        QMessageBox::critical(NULL, "Ошибка", "Добавьте объекты на сцену");
    else
        myGraphics->Start();
}







void MainWindow::on_pushButton_2_clicked()
{
    myGraphics->MinusMash();
}

void MainWindow::on_pushButton_3_clicked()
{
    myGraphics->PlusMash();
}


void MainWindow::on_action_2_triggered()
{
    myGraphics->DownloadSolarSystem();
    this->star_f = true;

}

void MainWindow::on_action_3_triggered()
{
    QMessageBox::information(NULL, "Информация", "Программа для визуализации солнечной системы \n Автор: Студент ИУ7-52Б Маслюков Павел");
}





void MainWindow::on_full_end_clicked()
{
    myGraphics->Clean();
    this->star_f = false;
}


void MainWindow::on_action_triggered()
{
    myGraphics->DownloadBistarSystem();
    this->star_f = true;
}


void MainWindow::on_pushButton_4_clicked()
{
    this->close();
}


void MainWindow::on_pushButton_5_clicked()
{
    if (!star_f)
        QMessageBox::critical(NULL, "Ошибка", "Добавьте объекты на сцену");
    else
        myGraphics->StartRun();
}




void MainWindow::on_updateStar_clicked()
{
    //    double m_tmp = ui->star_m->value();
    //    m_tmp *= pow(10.0, ui->star_m_e->value());
    //    int r_tmp = ui->star_r->value();
    //    std::cout << m_tmp;
    //    myGraphics->UpdateStar(m_tmp, r_tmp, true);
    //    this->star_f = true;
}


void MainWindow::on_pushButton_clicked()
{
    //    double m_tmp = ui->planet_m->value();
    //    m_tmp *= pow(10.0, ui->planet_m_e->value());
    //    int r_tmp = ui->planet_r->value();
    //    int r_tmp_2 = ui->planet_r_2->value();
    //    double X = ui->X_m->value();
    //    X *= pow(10.0, ui->X_m_E->value());
    //    double Y = ui->y_m->value();
    //    Y *= pow(10.0, ui->y_m_e->value());
    //    double Z = ui->z_m->value();
    //    Z *= pow(10.0, ui->z_m_e->value());

    //    double Vx = ui->Vx->value();
    //    double Vy = ui->Vy->value();
    //    double Vz = ui->Vz->value();

    //    int type = ui->comboBox->currentIndex();

    //    Point tmp (X, Y, Z);
    //    if (!star_f)
    //        QMessageBox::critical(NULL, "Ошибка", "Добавьте звезду");
//    else
    //        myGraphics->AddPlanet(tmp, Vx, Vy, Vz, m_tmp, r_tmp, r_tmp_2, type, true, this->pen);
}


void MainWindow::on_color_p_clicked()
{
    //    QColor tmpC = QColor(Qt::yellow);
    //    color_dialog(tmpC);
    //    QImage img = QImage(ui->color_l->geometry().width(), ui->color_l->geometry().height(), QImage::Format_RGB32);
    //    QPainter p(&img);
    //    p.setBrush(QBrush(tmpC));
    //    p.setPen(Qt::black);
    //    QRect rect = QRect(0, 0, ui->color_l->geometry().width(), ui->color_l->geometry().height());
    //    p.drawRect(rect);

    //    QPixmap pixmap = QPixmap::fromImage(img);
    //    ui->color_l->clear();
    //    ui->color_l->setPixmap(pixmap);
    //    pen.setColor(tmpC);
}

void MainWindow::color_dialog(QColor& color)
{
    //    QColorDialog dialog;
    //    dialog.setCurrentColor(color);
    //    dialog.show();
    //    dialog.exec();
    //    QColor tmp = dialog.selectedColor();
    //    if (!tmp.isValid())
    //        QMessageBox::critical(NULL, "Ошибка", "Выберите цвет");
    //    else
    //        color = tmp;
}




