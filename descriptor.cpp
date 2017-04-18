#include "descriptor.h"

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
    sobelX = m.convolution(KernelFactory::sobelX(),Matrix::Border::COPIED);
    sobelY = m.convolution(KernelFactory::sobelY(),Matrix::Border::COPIED);
    points = CornerDetectors().detect(m,Algorithm::HARIS);
    points = PointFileter(m.width()*m.height(),maxNumPoints).filter(points);
    sigma0 = 2.0;
    sigma = 2.0;
}

DescrBuilder::DescrBuilder(const Matrix &m,const double sigma,
                                        const vector<Point> & points)
{
    sobelX = m.convolution(KernelFactory::sobelX(),Matrix::Border::COPIED);
    sobelY = m.convolution(KernelFactory::sobelY(),Matrix::Border::COPIED);
    this->points = points;
    this->sigma0 = 1.6;
    this->sigma = sigma;
}


vector<Descriptor> DescrBuilder::build() const
{
    auto descriptors = vector<Descriptor>();
    for(auto &each : points){
        auto data = computeData(each,0,numRotateHist,sizeRotateHist,numRotateBins);

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

        if(nangles > 2)
            continue;

        for(int i = 0; i<nangles; i++){
            descriptors.emplace_back(normilizeAll(
               Descriptor(each.x,each.y,
                   computeData(each,angles[i],
                               numDescrHist, sizeDescrHist, numDescrBins))));
        }
    }
    return descriptors;
}

Descriptor DescrBuilder::normilizeAll(const Descriptor &descr) const
{
    return normilize(filterTrash(normilize(descr)));
}

vector<double> DescrBuilder::computeData(const Point &p,
    const double rotateAngle, const int numHist,
                                const int sizeHist, const int numBins) const
{
    int sizeArea = numHist*sizeHist*(sigma/sigma0);

    auto data = vector<double>(numHist * numHist * numBins, 0.0);

    auto aSin = sin(rotateAngle);
    auto aCos = cos(rotateAngle);

    int radArea = sizeArea/2;
    for(int i = 0; i < sizeArea; i++){
        for(int j = 0; j< sizeArea; j++){

            int cx = i - radArea;
            int cy = j - radArea;

            int x = (cx * aCos + cy * aSin) + radArea;
            int y = (-cx * aSin + cy * aCos) + radArea;

            if(x >= sizeArea|| y >= sizeArea || x < 0 || y < 0){
                continue;
            }

            double dx = sobelX.get(p.x + cx, p.y + cy, Matrix::Border::COPIED);
            double dy = sobelY.get(p.x + cx, p.y + cy, Matrix::Border::COPIED);

            double magnitud = sqrt(dx*dx+dy*dy)*Utils::gauss(cx, cy, sigma);

            double phi = atan2(dx,dy)+M_PI - rotateAngle;

            phi = fmod(phi+4*M_PI-M_PI/numBins,2*M_PI);

            double num = phi/(2*M_PI/numBins);

            int binIdx1 = ((int)floor(num))%numBins;
            int binIdx2 = ((int)(floor(num)+1))%numBins;

            int hIdx = ((x/sizeArea)*numHist+y/sizeArea)*numBins;

            assert(hIdx + binIdx1 < data.size() && hIdx + binIdx2 < data.size());

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


