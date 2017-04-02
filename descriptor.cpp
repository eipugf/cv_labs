#include "descriptor.h"
#include <stdio.h>

Descriptor::Descriptor(const int x, const int y):
    x(x),y(y){}

Descriptor::Descriptor(const int x, const int y, const vector<double> &data):
    Descriptor(x,y)
{
    for(double each:data){
        this->data.push_back(each);
    }

}

vector<Descriptor> DescrBuilder::build(const Matrix &m) const
{
    auto points = CornerDetectors().detect(m,Algorithm::HARIS);
    points = PointFileter(m.width()*m.height(),30).filter(points);
    auto sobleX = m.convolution(KernelFactory::sobelX(),Matrix::Border::CILINDER);
    auto sobleY = m.convolution(KernelFactory::sobelY(),Matrix::Border::CILINDER);
    auto descriptors = vector<Descriptor>();
    for(auto &each : points){
        auto descr = Descriptor(each.x,each.y,move(computeData(each,sobleX,sobleY)));
        descr = normilize(descr);
        descr = filterTrash(descr);
        descr = normilize(descr);
        descriptors.emplace_back(move(descr));
    }
    return descriptors;
}

vector<double> DescrBuilder::computeData(const Point &p,
                             const Matrix & sobelX,const Matrix & sobelY) const
{
    int size = sizeArea/sizeHist;
    auto data = vector<double>(size*size * numBins, 0.0);
    for(int i = 0; i < sizeArea; i++){
        for(int j = 0; j< sizeArea; j++){
            double dx = sobelX.get(p.x + i - sizeArea/2,
                          p.y + j - sizeArea/2, Matrix::Border::CILINDER);
            double dy = sobelY.get(p.x + i - sizeArea/2,
                          p.y + j - sizeArea/2, Matrix::Border::CILINDER);

            double magnitud = sqrt(dx*dx+dy*dy)*
                    Utils::gauss(i - sizeArea/2, j - sizeArea/2, sigma);
           // printf("%lg ",Utils::gauss(i - sizeArea/2, j - sizeArea/2, sigma));
            double phi = atan2(dx,dy)+M_PI;
            phi = fmod(phi+2*M_PI-M_PI/numBins,2*M_PI);

            double num = phi/(2*M_PI/numBins);

            int binIdx1 = ((int)floor(num))%numBins;
            int binIdx2 = ((int)(floor(num)+1))%numBins;

            int hIdx = ((i/sizeHist)*size+j/sizeHist)*numBins;
            data[hIdx + binIdx1] += magnitud*(ceil(num) - num);
            data[hIdx + binIdx2] += magnitud*(num - floor(num));
        }
       // printf("\n");
    }
    return data;
}

Descriptor DescrBuilder::filterTrash(const Descriptor &descr) const
{
    Descriptor result(descr.x,descr.y);
    for(double bin:descr.data){
        result.data.push_back(min(treshold, bin));
    }
    return result;
}

Descriptor DescrBuilder::normilize(const Descriptor &descr) const
{
    Descriptor result(descr.x,descr.y);
    double sum = 0;
    for(double each:descr.data){
        sum += each*each;
    }
    sum = sqrt(sum);
    for(double each:descr.data){
        result.data.push_back(each/sum);
    }
    return result;
}

vector<pair<Point, Point> > PointSearcher::findSamePoints(
        const Matrix &m1, const Matrix &m2) const
{
    auto descrM1 = DescrBuilder().build(m1);
    auto descrM2 = DescrBuilder().build(m2);

    vector<pair<Point, Point>> samePoints;
    for(Descriptor & each:descrM1){
        for(Descriptor & each1:descrM2){
             double sum = 0;
             for(int i = 0; i<each.data.size(); i++){
                 sum += (each.data[i]-each1.data[i])*(each.data[i]-each1.data[i]);
             }
             sum = sqrt(sum);

             if(sum < eps){
                 samePoints.emplace_back(
                      pair<Point,Point>(
                              Point(each.x,each.y,0),Point(each1.x,each1.y,0)));
                 break;
             }
        }
    }
    return samePoints;
}
