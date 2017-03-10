#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include "kernel.h"
#include "utils.h"
#include "scalespace.h"
#include <memory>

#include <time.h>

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
                this,tr("Open file"),NULL,"JPG/JPEG (*.jpg)");
    if(path.length() > 0){
        auto image = make_unique<QImage>(path);
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

void MainWindow::showPicture(Matrix & img)
{
    auto mult = Utils::multiple(255);
    auto grayMatrix = std::move(img.compute(mult));
    auto image = make_unique<QImage>(img.width(),
                                     img.hight(), QImage::Format_ARGB32);
    for(int i = 0; i < image->width(); i++){
        for(int j = 0; j < image->height(); j++){
            byte gray = grayMatrix.get(i,j);
            image->setPixel(i,j,qRgb(gray,gray,gray));
        }
    }
    showImage(*image);
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
        auto resP = std::move(picture->compress()
           .convolution(KernelFactory::createGaussX(0.5),Matrix::Border::CILINDER).
           convolution(KernelFactory::createGaussY(0.5),Matrix::Border::CILINDER).
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
    QImage image = QImage(level.width(),level.hight(),
                          QImage::Format_ARGB32);
    auto mul255 = Utils::multiple(255);
    auto pixImage = level.compute(mul255);
    for(int i = 0; i < pixImage.width(); i++){
        for(int j = 0; j < pixImage.hight(); j++){
            float gray = pixImage.get(i,j);
            image.setPixel(i,j,qRgb(gray, gray, gray));
        }
    }
    image.save(QString(file.c_str()),"JPG");
}
