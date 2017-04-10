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

DescrBuilder::DescrBuilder(const Matrix &m)
{
    sobelX = m.convolution(KernelFactory::sobelX(),Matrix::Border::CILINDER);
    sobelY = m.convolution(KernelFactory::sobelY(),Matrix::Border::CILINDER);
    points = CornerDetectors().detect(m,Algorithm::HARIS);
    points = PointFileter(m.width()*m.height(),maxNumPoints).filter(points);
}

vector<Descriptor> DescrBuilder::build() const
{
    auto descriptors = vector<Descriptor>();
    for(auto &each : points){
        auto data = computeData(each,0,sizeArea,numBins*k);

        auto maxh = *max_element(data.begin(),data.end());

        double angles[3];
        int nangles = 0;
        for(int i = 0; i < data.size(); i++){
            auto h0 = data[i];
            auto hm = data[(i - 1 + data.size())%data.size()];
            auto hp = data[(i + 1 + data.size())%data.size()];

            if(h0 > 0.8*maxh && h0 > hm && h0 > hp){
                auto di = -0.5 * (hp - hm) / (hp + hm - 2 * h0);
                auto th = 2 * M_PI * (i + di + 0.5) / data.size();
                angles[nangles++] = th;
                if(nangles == 3)
                    break;
            }
        }

        //нашли слишком много пиков выбрасываем точку
        if(nangles > 1)
            continue;

        for(int i = 0; i<nangles; i++){
            descriptors.emplace_back(normilizeAll(
                Descriptor(each.x,each.y,computeData(each,angles[i],sizeHist,numBins))));
        }
    }
    return descriptors;
}

Descriptor DescrBuilder::normilizeAll(const Descriptor &descr) const
{
    return normilize(filterTrash(normilize(descr)));
}

vector<double> DescrBuilder::computeData(const Point &p,
    const double rotateAngle, const int sizeHist, const int numBins) const
{
    int size = sizeArea/sizeHist;
    auto data = vector<double>(size*size * numBins, 0.0);

    auto aSin = sin(rotateAngle);
    auto aCos = cos(rotateAngle);

    int radArea = sizeArea/2;
    for(int i = 0; i < sizeArea; i++){
        for(int j = 0; j< sizeArea; j++){

            int cx = i - radArea;
            int cy = j - radArea;

            int x = (cx * aCos - cy * aSin) + radArea;
            int y = (cx * aSin + cy * aCos) + radArea;

            //выпали из старой окрестности, выбрасываем эту точку
            if(x >= sizeArea|| y >= sizeArea || x < 0 || y < 0){
                continue;
            }

            double dx = sobelX.get(p.x + cx, p.y + cy, Matrix::Border::CILINDER);
            double dy = sobelY.get(p.x + cx, p.y + cy, Matrix::Border::CILINDER);

            double magnitud = sqrt(dx*dx+dy*dy)*Utils::gauss(cx, cy, sigma);

            double phi = atan2(dx,dy)+M_PI - rotateAngle;

            phi = fmod(phi+4*M_PI-M_PI/numBins,2*M_PI);

            double num = phi/(2*M_PI/numBins);

            int binIdx1 = ((int)floor(num))%numBins;
            int binIdx2 = ((int)(floor(num)+1))%numBins;

            int hIdx = ((x/sizeHist)*size+y/sizeHist)*numBins;

            data[hIdx + binIdx1] += magnitud*(ceil(num) - num);
            data[hIdx + binIdx2] += magnitud*(num - floor(num));
        }
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


PointMatcher::PointMatcher(const double eps):eps(eps){}

vector<pair<Point, Point> > PointMatcher::match(
        const Matrix &m1, const Matrix &m2) const
{
    auto descrM1 = DescrBuilder(m1).build();
    auto descrM2 = DescrBuilder(m2).build();
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

