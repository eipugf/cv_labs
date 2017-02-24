#include "mainwindow.h"
#include <QApplication>
#include <QDesktopWidget>
#include <memory>
#include "matrix.h"
#include "kernel.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    QRect rWindow = w.geometry();
    rWindow.moveCenter(a.desktop()->availableGeometry().center());
    w.setGeometry(rWindow);
    w.show();
    return a.exec();
}
