#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>
#include <QDebug>
#include <QColor>
#include<omp.h>
#include <fstream>
#include<iomanip>
#include<iostream>
#include <QtCharts>
#include <QChartView>
#include <QtCharts/QLineSeries>

#define METHOD_ORDER 4
#define MINSTEP 1e-8

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow){
    ui->setupUi(this);
}

MainWindow::~MainWindow(){
    delete ui;
}

void MainWindow::Solve_Test_Eq(){
    unsigned int meth_char = pow(2, METHOD_ORDER);
    size_t C1= 0; size_t C2 = 0;
    double tstep = step;
    double local_error = 0;
    double local_error_abs = 0;
    double xcurr = X_Left;
    double ycurr = start_position_y;
    double yn = 0;
    double yt = 0;

    series1->append(xcurr, ycurr);
    series2->append(xcurr,Eq_F_Test_TrueSol(xcurr, start_position_y));

    array_of_xi.push_back(xcurr);
    array_of_vi.push_back(ycurr);
    array_of_v2i.push_back(ycurr);
    sub_vi_v2i.push_back(0);
    OLP.push_back(0);
    hi.push_back(tstep);
    C1_Vec.push_back(0);
    C2_Vec.push_back(0);
    Ui.push_back(Eq_F_Test_TrueSol(xcurr, start_position_y));
    Ui_Vi.push_back(abs(Eq_F_Test_TrueSol(xcurr, start_position_y)-ycurr));

    while(xcurr + tstep < X_Right)
    {
        yn = RK4_Test(xcurr, ycurr, tstep);
        yt = RK4_Test(xcurr + tstep/2, RK4_Test(xcurr, ycurr, tstep/2), tstep/2);
        local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
        local_error_abs = abs(local_error);

        if((local_error_abs >= (RTOL / (meth_char*2)))&&(local_error_abs <= RTOL)){
            ycurr = yn;
            xcurr += tstep;
        } else if(local_error_abs < (RTOL / (meth_char*2))){
            ycurr = yn;
            xcurr += tstep;
            tstep = tstep * 2;
            C2++;
        } else if(local_error_abs > RTOL){
            while((local_error_abs > RTOL) && (tstep > RTOL)){
                tstep = fmax(MINSTEP, tstep / 2);
                yn = RK4_Test(xcurr, ycurr, tstep);
                yt = RK4_Test(xcurr + tstep/2, RK4_Test(xcurr, ycurr, tstep/2), tstep/2);
                local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
                local_error_abs = abs(local_error);
                C1++;
            }
            ycurr = yn;
            xcurr += tstep;
        }

        C1_Vec.push_back(C1);
        C2_Vec.push_back(C2);

        series1->append(xcurr,ycurr);
        series2->append(xcurr,Eq_F_Test_TrueSol(xcurr, start_position_y));

        array_of_xi.push_back(xcurr);
        array_of_vi.push_back(ycurr);
        array_of_v2i.push_back(yt);
        sub_vi_v2i.push_back(ycurr-yt);
        OLP.push_back(local_error_abs);
        Ui.push_back(Eq_F_Test_TrueSol(xcurr, start_position_y));
        Ui_Vi.push_back(abs(Eq_F_Test_TrueSol(xcurr, start_position_y)-ycurr));
        hi.push_back(tstep);
        count_of_iter++;

        //std::cout<<tstep<<' '<<xcurr<<' '<<ycurr<<' '<<local_error_abs<<' '<<RTOL<<' '<<(RTOL / (meth_char*2))<<' '<<(local_error_abs < (RTOL / (meth_char*2)))<<' '<<(local_error_abs > RTOL)<<std::endl;
    }
    if (xcurr + tstep >= X_Right){
        tstep = X_Right - xcurr;
        yn = RK4_Test(xcurr, ycurr, tstep);
        yt = RK4_Test(xcurr + tstep/2, RK4_Test(xcurr, ycurr, tstep/2), tstep/2);
        local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
        local_error_abs = abs(local_error);
        hi.push_back(tstep);
        ycurr = yn;
        xcurr = X_Right;

        C1_Vec.push_back(C1);
        C2_Vec.push_back(C2);

        series1->append(xcurr,ycurr);
        series2->append(xcurr,Eq_F_Test_TrueSol(xcurr, start_position_y));

        array_of_xi.push_back(xcurr);
        array_of_vi.push_back(ycurr);
        array_of_v2i.push_back(yt);
        sub_vi_v2i.push_back(ycurr-yt);
        OLP.push_back(local_error_abs);
        Ui.push_back(Eq_F_Test_TrueSol(xcurr, start_position_y));
        Ui_Vi.push_back(abs(Eq_F_Test_TrueSol(xcurr, start_position_y)-ycurr));
        count_of_iter++;
    }

    Vec_Of_Vecs.push_back(array_of_xi);
    Vec_Of_Vecs.push_back(array_of_vi);
    Vec_Of_Vecs.push_back(array_of_v2i);
    Vec_Of_Vecs.push_back(sub_vi_v2i);
    Vec_Of_Vecs.push_back(OLP);
    Vec_Of_Vecs.push_back(hi);
    Vec_Of_Vecs.push_back(C1_Vec);
    Vec_Of_Vecs.push_back(C2_Vec);
    Vec_Of_Vecs.push_back(Ui);
    Vec_Of_Vecs.push_back(Ui_Vi);

    ui->tableWidget->setColumnCount(10);
    ui->tableWidget->setRowCount(count_of_iter+1);
    ui->tableWidget->setHorizontalHeaderLabels(QStringList()<<"X_i"<<"V_i"<<"V2i"<<"V_i - V2i"<<"OLP"<<"h_i"<<"C1"<<"C2"<<"Ui"<<"Ui - Vi");
    for(size_t i = 0 ; i < ui->tableWidget->columnCount();i++)
    {
        for(size_t j = 0 ; j < ui->tableWidget->rowCount();j++)
        {
            QTableWidgetItem *tbl = new QTableWidgetItem(QString::number(Vec_Of_Vecs[i][j]));

            ui->tableWidget->setItem(j,i,tbl);
        }
    }

    ui->aux_info->append(QString("Итераций: ") + QString::number(count_of_iter));
}

void MainWindow::Solve_Eq(){
    unsigned int meth_char = pow(2, METHOD_ORDER);
    size_t C1= 0; size_t C2 = 0;
    double tstep = step;
    double local_error = 0;
    double local_error_abs = 0;
    double xcurr = X_Left;
    double ycurr = start_position_y;
    double yn = 0;
    double yt = 0;

    series1->append(xcurr, ycurr);

    array_of_xi.push_back(xcurr);
    array_of_vi.push_back(ycurr);
    array_of_v2i.push_back(ycurr);
    sub_vi_v2i.push_back(0);
    OLP.push_back(0);
    hi.push_back(tstep);
    C1_Vec.push_back(0);
    C2_Vec.push_back(0);
    while(xcurr + tstep < X_Right)
    {
        yn = RK4(xcurr, ycurr, tstep);
        yt = RK4(xcurr + tstep/2, RK4(xcurr, ycurr, tstep/2), tstep/2);
        local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
        local_error_abs = abs(local_error);

        if((local_error_abs >= (RTOL / (meth_char*2)))&&(local_error_abs <= RTOL)){
            ycurr = yn;
            xcurr += tstep;
        } else if(local_error_abs < (RTOL / (meth_char*2))){
            ycurr = yn;
            xcurr += tstep;
            tstep = tstep * 2;
            C2++;
        } else if(local_error_abs > RTOL){
            while((local_error_abs > RTOL) && (tstep > RTOL)){
                tstep = fmax(MINSTEP, tstep / 2);
                yn = RK4(xcurr, ycurr, tstep);
                yt = RK4(xcurr + tstep/2, RK4(xcurr, ycurr, tstep/2), tstep/2);
                local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
                local_error_abs = abs(local_error);
                C1++;
            }
            ycurr = yn;
            xcurr += tstep;
        }

        C1_Vec.push_back(C1);
        C2_Vec.push_back(C2);

        series1->append(xcurr,ycurr);

        array_of_xi.push_back(xcurr);
        array_of_vi.push_back(ycurr);
        array_of_v2i.push_back(yt);
        sub_vi_v2i.push_back(ycurr-yt);
        OLP.push_back(local_error_abs);
        hi.push_back(tstep);
        count_of_iter++;

        //std::cout<<tstep<<' '<<xcurr<<' '<<ycurr<<' '<<local_error_abs<<' '<<RTOL<<' '<<(RTOL / (meth_char*2))<<' '<<(local_error_abs < (RTOL / (meth_char*2)))<<' '<<(local_error_abs > RTOL)<<std::endl;
    }
    if (xcurr + tstep >= X_Right){
        tstep = X_Right - xcurr;
        yn = RK4(xcurr, ycurr, tstep);
        yt = RK4(xcurr + tstep/2, RK4(xcurr, ycurr, tstep/2), tstep/2);
        local_error= ((yt - yn) / (meth_char - 1)) * (meth_char);
        local_error_abs = abs(local_error);
        hi.push_back(tstep);
        ycurr = yn;
        xcurr = X_Right;

        C1_Vec.push_back(C1);
        C2_Vec.push_back(C2);

        series1->append(xcurr,ycurr);

        array_of_xi.push_back(xcurr);
        array_of_vi.push_back(ycurr);
        array_of_v2i.push_back(yt);
        sub_vi_v2i.push_back(ycurr-yt);
        OLP.push_back(local_error_abs);
        count_of_iter++;
    }

    Vec_Of_Vecs.push_back(array_of_xi);
    Vec_Of_Vecs.push_back(array_of_vi);
    Vec_Of_Vecs.push_back(array_of_v2i);
    Vec_Of_Vecs.push_back(sub_vi_v2i);
    Vec_Of_Vecs.push_back(OLP);
    Vec_Of_Vecs.push_back(hi);
    Vec_Of_Vecs.push_back(C1_Vec);
    Vec_Of_Vecs.push_back(C2_Vec);

    ui->tableWidget->setColumnCount(8);
    ui->tableWidget->setRowCount(count_of_iter+1);
    ui->tableWidget->setHorizontalHeaderLabels(QStringList()<<"X_i"<<"V_i"<<"V2i"<<"V_i - V2i"<<"OLP"<<"h_i"<<"C1"<<"C2");
    for(size_t i = 0 ; i < ui->tableWidget->columnCount();i++)
    {
        for(size_t j = 0 ; j < ui->tableWidget->rowCount();j++)
        {
            QTableWidgetItem *tbl = new QTableWidgetItem(QString::number(Vec_Of_Vecs[i][j]));

            ui->tableWidget->setItem(j,i,tbl);
        }
    }

    ui->aux_info->append(QString("Итераций: ") + QString::number(count_of_iter));
}

void MainWindow::Solve_Sys(){
    unsigned int meth_char = pow(2, METHOD_ORDER);
    double tstep = step;
    double local_error = 0;
    double xcurr = X_Left;
    double y1curr = start_position_y;
    double y2curr = start_position_y_sys;
    double* ynn = new double[2];
    double* ynt1 = new double[2];
    double* ynt2 = new double[2];
    double local_error_abs = 0;
    std::cout<<a_param<<' '<<b_param<<std::endl;

    while((xcurr + tstep) < X_Right)
    {
        ynt1 = RK4_Sys(xcurr, y1curr, y2curr, tstep/2);
        ynt2 = RK4_Sys(xcurr + tstep/2, ynt1[0], ynt1[1], tstep/2);
        ynn = RK4_Sys(xcurr,  y1curr,  y2curr,  tstep);
        local_error= sqrt(pow((ynt2[0] - ynn[0]) / (meth_char - 1), 2) + pow((ynt2[1] - ynn[1]) / (meth_char - 1), 2));
        local_error_abs = fabs(local_error);
        if((local_error_abs >= (RTOL / (meth_char*2)))&&(local_error_abs <= RTOL)){
            y1curr = ynn[0];
            y2curr = ynn[1];
            xcurr += tstep;
        } else if(local_error_abs<(RTOL / (meth_char*2))){
            y1curr = ynn[0];
            y2curr = ynn[1];
            xcurr += tstep;
            tstep = tstep * 2;
        } else if(local_error_abs>RTOL){
            while((local_error_abs > RTOL) && (tstep > RTOL)){
                tstep = fmax(MINSTEP, tstep / 2);
                ynt1 = RK4_Sys(xcurr, y1curr, y2curr, tstep/2);
                ynt2 = RK4_Sys(xcurr + tstep/2, ynt1[0], ynt1[1], tstep/2);
                ynn = RK4_Sys(xcurr,  y1curr,  y2curr,  tstep);
                local_error= sqrt(pow((ynt2[0] - ynn[0]) / (meth_char - 1), 2) + pow((ynt2[1] - ynn[1]) / (meth_char - 1), 2));
                local_error_abs = fabs(local_error);
            }
            y1curr = ynn[0];
            y2curr = ynn[1];
            xcurr += tstep;
        }
        if(fabs(y1curr) > 1e+3 || fabs(y2curr) > 1e+3) break;

        switch (choice_gr)
        {
        case 0:
            series1->append(xcurr, y1curr);
            break;
        case 1:
            series1->append(xcurr, y2curr);
            break;
        case 2:
            series1->append(y1curr, y2curr);
            break;
        }
        //std::cout<<tstep<<' '<<xcurr<<' '<<y1curr<<' '<<y2curr<<' '<<local_error_abs<<' '<<RTOL<<' '<<(RTOL / (meth_char*2))<<' '<<(local_error_abs < (RTOL / (meth_char*2)))<<' '<<(local_error_abs > RTOL)<<std::endl;

        count_of_iter++;
    }

    delete[] ynn;
    delete[] ynt1;
    delete[] ynt2;

    ui->aux_info->append(QString("Итераций: ") + QString::number(count_of_iter));
}

void MainWindow::ClearScreen(){
    series1->detachAxis(axisX);
    series1->detachAxis(axisY);
    series2->detachAxis(axisX);
    series2->detachAxis(axisY);
    chart->removeAxis(axisX);
    chart->removeAxis(axisY);
    ui->horizontalLayout_10->removeWidget(chartView);
    chart->removeSeries(series1);
    chart->removeSeries(series2);
    disconnect(chartView);
    disconnect(chart);
    delete chart;
    delete chartView;
    delete axisX;
    delete axisY;

    if(series1->count() != 0) series1->clear();
    if(series2->count() != 0) series2->clear();

    chartView = new QChartView();
    chart = new QChart();
    axisX = new QValueAxis();
    axisY = new QValueAxis();

    ui->tableWidget->setRowCount(0);
    ui->tableWidget->setColumnCount(0);

    Vec_Of_Vecs.clear();

    array_of_xi.clear();
    array_of_vi.clear();
    array_of_v2i.clear();
    sub_vi_v2i.clear();
    OLP.clear();
    hi.clear();
    C1_Vec.clear();
    C2_Vec.clear();
    Ui.clear();
    Ui_Vi.clear();
    count_of_iter = 0;

    ui->aux_info->clear();
}

void MainWindow::GetData(){
    X_Left = ui->X_EDIT_MIN->text().toDouble();
    X_Right = ui->X_EDIT_MAX->text().toDouble();

    RTOL = ui->Acc_Edit_3->text().toDouble();
    step = ui->StepEdit_3->text().toDouble();
    start_position_y = ui->Start_Posit->text().toDouble();
    start_position_y_sys = ui->Start_y_for_sys->text().toDouble();

    a_param = ui->a_param_input->text().toDouble();
    b_param = ui->b_param_input->text().toDouble();

    if(ui->Test_task_button->isChecked())
    {
        choice = 0;
    }
    else if(ui->main_task->isChecked())
    {
        choice = 1;
    }
    else if(ui->main_task_syst->isChecked())
    {
        choice = 2;
    }
    else
        this->close();


    if(ui->xy1->isChecked())
    {
        choice_gr = 0;
    }
    else if(ui->xy2->isChecked())
    {
        choice_gr = 1;
    }
    else if(ui->y1y2->isChecked())
    {
        choice_gr = 2;
    }
    else
        if (choice == 2) this->close();
}

void MainWindow::on_Draw_but_3_clicked(){
    ClearScreen();
    GetData();
    std::cout<<"a"<<a_param<<std::endl;

    switch (choice)
    {
    case 0:
        {
            ui->horizontalLayout_10->addWidget(chartView);
            chart->legend()->hide();

            Solve_Test_Eq();

            QPen pen(QRgb(0xfdb157));
            pen.setWidth(1);
            series1->setPen(pen);
            chart->addSeries(series1);

            QPen pen1(QRgb(0x0000FF));
            pen.setWidth(1);
            series2->setPen(pen1);
            chart->addSeries(series2);
        }
        break;

    case 1:
        {
            ui->horizontalLayout_10->addWidget(chartView);
            chart->legend()->hide();

            Solve_Eq();
            QPen pen(QRgb(0xfdb157));
            pen.setWidth(1);
            series1->setPen(pen);
            chart->addSeries(series1);
        }
        break;

    case 2:
        {
            ui->horizontalLayout_10->addWidget(chartView);
            chart->legend()->hide();
            Solve_Sys();

            QPen pen(QRgb(0xfdb157));
            pen.setWidth(1);
            series1->setPen(pen);
            chart->addSeries(series1);
        }
        break;
    }
    chart->setTitle("Graphic");
    //QValueAxis *axisX = new QValueAxis();
    axisX->setTitleText("x");
    axisX->setLabelFormat("%f");
    axisX->setTickCount(1);
    chart->addAxis(axisX, Qt::AlignBottom);
    series1->attachAxis(axisX);


    // QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("y");
    axisY->setLabelFormat("%f");
    axisY->setTickCount(10);
    chart->addAxis(axisY, Qt::AlignLeft);
    series1->attachAxis(axisY);
    chartView->setChart(chart);
}

void MainWindow::on_Clear_But_3_clicked(){
    ClearScreen();
}


double MainWindow::RK4(double x,double y, double h){
    double k1, k2, k3, k4, k;
    double yn = 0;
    double x0 = x;
    double y0 = y;

    k1 = Eq_F(x0, y0);
    k2 = Eq_F((x0 + h/2), (y0 + (h/2) * k1));
    k3 = Eq_F((x0 + h/2), (y0 + (h/2) * k2));
    k4 = Eq_F((x0 + h), (y0 + h * k3));
    k = h*(k1 + 2 * k2 + 2 * k3 + k4)/6;
    yn = y0 + k;

    return yn;
}

double MainWindow::RK4_Test(double x,double y, double h){
    double k1, k2, k3, k4, k;
    double yn = 0;
    double x0 = x;
    double y0 = y;

    k1 = Eq_F_Test(x0, y0);
    k2 = Eq_F_Test((x0 + h/2), (y0 + (h/2) * k1));
    k3 = Eq_F_Test((x0 + h/2), (y0 + (h/2) * k2));
    k4 = Eq_F_Test((x0 + h), (y0 + h * k3));
    k = h*(k1 + 2 * k2 + 2 * k3 + k4)/6;
    yn = y0 + k;

    return yn;
}

double* MainWindow::RK4_Sys(double x, double y1, double y2, double h){
    double K[4][2];
    double* array_of_y = new double[2];

    K[0][0] = Sys_F1(x, y1, y2);
    K[0][1] = Sys_F2(x, y1, y2);
    K[1][0] = h * Sys_F1(x + h/2, y1+K[0][0]/2, y2 + K[0][1]/2);
    K[1][1] = h * Sys_F2(x + h/2, y1+K[0][0]/2, y2 + K[0][1]/2);
    K[2][0] = h * Sys_F1(x + h/2, y1+K[1][0]/2, y2 + K[1][1]/2);
    K[2][1] = h * Sys_F2(x + h/2, y1+K[1][0]/2, y2 + K[1][1]/2);
    K[3][0] = h * Sys_F1(x + h, y1+K[2][0], y2 + K[2][1]);
    K[3][1] = h * Sys_F2(x + h, y1+K[2][0], y2 + K[2][1]);
    array_of_y[0] =y1 + ((K[0][0] + 2 * K[1][0] + 2 * K[2][0] + K[3][0])/6);
    array_of_y[1] =y2 + ((K[0][1] + 2 * K[1][1] + 2 * K[2][1] + K[3][1])/6);

    return array_of_y;
}

double MainWindow::Sys_F1(double x, double y1, double y2){
    return y2;
}

double MainWindow::Sys_F2(double x, double y1, double y2){
    return (-a_param) * y2 * y2 - b_param * sin(y1);
}

double MainWindow::Eq_F(double x, double y){
    return ((log(x + 1) / (x * x + 1)) * y * y + y - y * y * y * sin(10 * x));
}

double MainWindow::Eq_F_Test(double x, double y){
    return ((-5)/2 * y);
}

double MainWindow::Eq_F_Test_TrueSol(double x, double y0){
    return (y0 * exp(((-5)/2) * x));
}
