#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QImage>
#include <string>
#include <memory>
#include "matrix.h"

using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

    unique_ptr<Matrix> picture;

    QGraphicsScene * scene = new QGraphicsScene(this);

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_close_triggered();

    void on_openPicture_triggered();

    void on_savePicture_triggered();

   // void on_grayScale_triggered();

    void on_sobelAction_triggered();

    void on_gaussAction_triggered();

    void on_scaleSpace_triggered();

private:
    Ui::MainWindow *ui;

    void showPicture(Matrix & image);
    void save(const Matrix & level, const string & file) const;
    void showImage(QImage & image);
};

#endif // MAINWINDOW_H
