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
#include "haugh.h"
#include "ransac.h"

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
                this,tr("Open file"),NULL,"PNG (*.png | *.jpg)");

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
    pen.setWidth(1);
    pen.setColor(Qt::red);
    painter.setPen(pen);

    for(auto &each:pairs){
        Point p1 = each.first;
        Point p2 = each.second;
        double p2x = p2.x+img.width()/2;
        double p2y = p2.y;

        painter.drawLine(p1.x, p1.y,p2x, p2y);
        painter.drawEllipse(QPointF(p1.x,p1.y),p1.rad, p1.rad);

        painter.drawEllipse(QPointF(p2x,p2y),p2.rad, p2.rad);
        painter.drawLine(p1.x,p1.y,p1.x+p1.rad*sin(p1.angle),p1.y+p1.rad*cos(p1.angle));

        painter.drawLine(p2x,p2y,p2x+p2.rad*sin(p2.angle),p2y+p2.rad*cos(p2.angle));

    }
    painter.end();
    showImage(img);
}

void MainWindow::showPictureWithDescr(QImage &img, vector<Descriptor> &descr)
{
    QPainter painter(&img);
    QPen pen;
    pen.setWidth(1);
    pen.setColor(Qt::red);
    painter.setPen(pen);

    for(auto &each:descr){
        painter.drawEllipse(QPointF(each.x,each.y), each.rad, each.rad);
        double x2 = each.x + each.rad*sin(each.angle);
        double y2 = each.y + each.rad*cos(each.angle);
        painter.drawLine(each.x,each.y,x2,y2);
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
    auto pairs = PointMatcher(0.2).match(m0,m1);


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
        rpen.setWidth(1);
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
    QImage image0 = QImage("/home/eugene/rans.jpg");
    QImage image1 = QImage("/home/eugene/rans1.jpg");

    Matrix m0 = imageToMatrix(image0);
    Matrix m1 = imageToMatrix(image1);
    auto pairs = PointMatcher(0.05).match(m0,m1,true);

    QImage pictures = QImage(m0.width()*2, m0.height(),QImage::Format_ARGB32);
    for(int i = 0; i<m0.width(); i++){
        for(int j = 0; j<m0.height(); j++){
            pictures.setPixel(i,j,image0.pixel(i,j));
            pictures.setPixel(i+image1.width(),j,image1.pixel(i,j));
        }
    }

    showPictureWithPoints(pictures,pairs);
}

void MainWindow::on_seeDescriptors_triggered()
{
    if(picture!=NULL){
        auto descr = SIDiscrBuilder::build(*picture);
        showPictureWithDescr(*image,descr);
    }
}

void MainWindow::on_computeRansac_triggered()
{
    QImage image0 = QImage("/home/eugene/rans.jpg");
    QImage image1 = QImage("/home/eugene/rans1.jpg");
    Matrix m0 = imageToMatrix(image0);
    Matrix m1 = imageToMatrix(image1);
    auto pairs = PointMatcher(0.1).match(m0,m1,true);

    vector<pair<Point, Point>> swaped;
    for(auto & pair:pairs){
        Point p1(pair.first.x,pair.first.y,0);
        Point p2(pair.second.x,pair.second.y,0);

        swaped.emplace_back(std::pair<Point,Point>(p1,p2));
    }

    auto h = Ransac().searchTransform(swaped);

    QImage result((m0.width() + m1.width()) ,
                  (m0.height() + m1.height()) / 1.5,
                  QImage::Format_RGB32);



    QPainter painter(&result);
    painter.translate(300, 0);
    painter.drawImage(0, 0, image1);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::SmoothPixmapTransform);
    QTransform transform(h[0], h[3], h[6],
                         h[1], h[4], h[7],
                         h[2], h[5], h[8]
                         );
    painter.setTransform(transform, true);
    painter.drawImage(0, 0, image0);

    showImage(result);
}

void MainWindow::on_actionHough_triggered()
{
    QImage image0 = QImage("/home/eugene/hough.jpg");
    QImage image1 = QImage("/home/eugene/hough.jpg");
    Matrix m0 = imageToMatrix(image0);
    Matrix m1 = imageToMatrix(image1);
    auto pairs = PointMatcher(0.1).match(m0,m1,true);

    if(pairs.size() > 2){
        auto transf = Hough(m0.width(),m0.height(),m1.width(),
                            m1.height(),pairs).computeHaugh();
        QImage result = QImage(image1);
        QPainter painter(&result);
        int size1 = m1.width()*transf.scale;
        int size2 = m1.width()*transf.scale;
        painter.drawEllipse(QPointF(transf.x, transf.y), 10, 10);
        QRect rect = QRect(transf.x - m1.width()*transf.scale / 2,
                           transf.y - m1.height()*transf.scale / 2,
                           size1, size2);
        painter.drawRect(rect);
        showImage(result);
    } else {
        showImage(image1);
    }

}
