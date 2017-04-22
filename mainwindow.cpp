#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include "kernel.h"
#include "utils.h"
#include "scalespace.h"
#include <memory>
#include <QPainter>
#include <time.h>
#include <stdio.h>
#include "descriptor.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->image->setScene(scene);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_close_triggered()
{
    this->close();
}

void MainWindow::on_openPicture_triggered()
{
    QString path = QFileDialog::getOpenFileName(
                this,tr("Open file"),NULL,"PNG (*.png)");
    if(path.length() > 0){
        image = make_unique<QImage>(path);
        picture = make_unique<Matrix>(image->width(),image->height());
        for(int i = 0; i < image->width(); i++){
            for(int j = 0; j < image->height(); j++){
                picture->set(i,j,Utils::gray(qRed(image->pixel(i,j)),
                                       qGreen(image->pixel(i,j)),
                                       qBlue(image->pixel(i,j))));
            }
        }
        showImage(*image);
    }
}



void MainWindow::showImage(QImage &image)
{
    scene->clear();
    scene->addPixmap(QPixmap::fromImage(image));
    ui->image->invalidateScene();
}

void MainWindow::showPoints(vector<Point> &points)
{
    QImage sImage = QImage(*image);

    QPainter painter(&sImage);
    QPen pen;
    pen.setWidth(3);
    pen.setColor(Qt::red);
    painter.setPen(pen);

    for(vector<Point>::iterator p=points.begin(); p!=points.end(); p++){
        painter.drawPoint(p->x, p->y);
    }

    painter.end();
    showImage(sImage);
}

void MainWindow::showPicture(Matrix & img)
{
    auto mult = Utils::multiple(255);
    auto grayMatrix = std::move(img.compute(mult));
    auto image = make_unique<QImage>(img.width(),
                                     img.height(), QImage::Format_ARGB32);
    for(int i = 0; i < image->width(); i++){
        for(int j = 0; j < image->height(); j++){
            byte gray = grayMatrix.get(i,j);
            image->setPixel(i,j,qRgb(gray,gray,gray));
        }
    }
    showImage(*image);
}


void MainWindow::showPictureWithPoints(QImage & img,
                                       vector<pair<Point,Point>> & pairs)
{
    QPainter painter(&img);
    QPen pen;
    pen.setWidth(2);
    pen.setColor(Qt::red);
    painter.setPen(pen);

    for(auto &each:pairs){
        Point p1 = each.first;
        Point p2 = each.second;
        painter.drawLine(p1.x, p1.y,p2.x+img.width()/2, p2.y);
    }

    painter.end();

    showImage(img);
}



Matrix MainWindow::imageToMatrix(QImage &image)
{
    Matrix m = Matrix(image.width(),image.height());

    for(int i = 0; i<image.width(); i++){
        for(int j = 0; j < image.height(); j++){
            m.set(i,j,Utils::gray(qRed(image.pixel(i,j)),
                                   qGreen(image.pixel(i,j)),
                                   qBlue(image.pixel(i,j))));
        }
    }

    return m;
}

void MainWindow::on_savePicture_triggered()
{

}

void MainWindow::on_sobelAction_triggered()
{
    if(picture != nullptr) {
        auto resX(std::move(picture->
           convolution(KernelFactory::sobelX(),Matrix::Border::SIMPLE)));
        auto resY(std::move(picture->
           convolution(KernelFactory::sobelY(),Matrix::Border::SIMPLE)));
        auto hyp = Utils::hypotenuse();
        auto resS(std::move(resX.compute(resY,hyp).normalize()));
        showPicture(resS);
    }
}

void MainWindow::on_gaussAction_triggered()
{
    if(picture != nullptr) {
        auto resP = std::move(
           picture->convolution(KernelFactory::createGaussX(0.5),Matrix::Border::COPIED).
           convolution(KernelFactory::createGaussY(0.5),Matrix::Border::COPIED).
           normalize());
        showPicture(resP);
    }
}

void MainWindow::on_scaleSpace_triggered()
{
    ScaleSpace space = ScaleSpace(*picture,5);
    auto & octav = space.octavs();
    for(int i = 0; i < octav.size(); i++){
        for(int j = 0; j < space.octaveSize(); j++){
            auto & scale = octav[i][j];
            string fileName = "/home/eugene/qtprojects/scale_space/" +
                 to_string(i) + "_octave_g_"+to_string(scale.efectSigma)+"_"+
                 to_string(j)+"_level_"+to_string(scale.sigma)+".jpg";
            save(scale.matrix,fileName);
        }
    }

}

void MainWindow::save(const Matrix & level, const string & file) const
{
    QImage image = QImage(level.width(),level.height(),
                          QImage::Format_ARGB32);
    auto mul255 = Utils::multiple(255);
    auto pixImage = level.compute(mul255);
    for(int i = 0; i < pixImage.width(); i++){
        for(int j = 0; j < pixImage.height(); j++){
            double gray = pixImage.get(i,j);
            image.setPixel(i,j,qRgb(gray, gray, gray));
        }
    }
    image.save(QString(file.c_str()),"JPG");
}


void MainWindow::on_action_triggered()
{
    if(picture!=NULL){
        auto points = CornerDetectors().detect(*picture,Algorithm::MORAVEC);
        points = PointFileter(picture->width()*picture->height(),150).filter(points);
        showPoints(points);
    }
}

void MainWindow::on_haris_triggered()
{
    if(picture!=NULL){
        auto points = CornerDetectors().detect(*picture,Algorithm::HARIS);
        points = PointFileter(picture->width()*picture->height(),150).filter(points);
        showPoints(points);
    }
}

void MainWindow::on_simpleCompareAction_triggered()
{
    QImage image0 = QImage("/home/eugene/Lenna.png");
    QImage image1 = QImage("/home/eugene/Lenna1.png");

    Matrix m0 = imageToMatrix(image0);
    Matrix m1 = imageToMatrix(image1);
    auto pairs = PointMatcher(0.17).match(m0,m1);


    QImage pictures = QImage(m0.width()*2, m0.height(),QImage::Format_ARGB32);
    for(int i = 0; i<m0.width(); i++){
        for(int j = 0; j<m0.height(); j++){
            pictures.setPixel(i,j,image0.pixel(i,j));
            pictures.setPixel(i+image1.width(),j,image1.pixel(i,j));
        }
    }

    showPictureWithPoints(pictures,pairs);

}

void MainWindow::on_findBlobs_triggered()
{
    if(picture!=NULL){

        auto space = ScaleSpace(*picture,5);
        auto blobs = space.computeDiffs().searchBlobs();

        QImage img = QImage(*image);
        QPainter painter(&img);
        QPen rpen;
        rpen.setWidth(2);
        rpen.setColor(Qt::red);
        painter.setPen(rpen);
        for(Blob each:blobs){
            double k = pow(2,each.octav);
            double radius = each.sigma*k*M_SQRT2;
            painter.drawEllipse(QPointF(each.x*k,each.y*k), radius, radius);
        }
        painter.end();
        showImage(img);
    }
}

void MainWindow::on_scaleCompareAction_triggered()
{
    QImage image0 = QImage("/home/eugene/Lenna.png");
    QImage image1 = QImage("/home/eugene/Lenna1.png");

    Matrix m0 = imageToMatrix(image0);
    Matrix m1 = imageToMatrix(image1);
    auto pairs = PointMatcher(0.2).match(m0,m1,true);

    QImage pictures = QImage(m0.width()*2, m0.height(),QImage::Format_ARGB32);
    for(int i = 0; i<m0.width(); i++){
        for(int j = 0; j<m0.height(); j++){
            pictures.setPixel(i,j,image0.pixel(i,j));
            pictures.setPixel(i+image1.width(),j,image1.pixel(i,j));
        }
    }

    showPictureWithPoints(pictures,pairs);
}
