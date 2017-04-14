#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H
#include "corner_detectors.h"

using namespace std;

struct Descriptor{
    int x;
    int y;
    vector<double> data;
    Descriptor(const int x,const int y);
    Descriptor(const int x, const int y, const vector<double> & data);
};


class DescrBuilder{
    const double treshold = 0.2;
    double sigma;

    const int maxNumPoints = 500;

    const int sizeArea = 16;
    const int sizeHist = 4;
    const int numBins = 8;

    const double k = 4.5;

    Matrix sobelX;
    Matrix sobelY;
    vector<Point> points;

public:
    DescrBuilder(const Matrix & m);
    DescrBuilder(const Matrix & m,const double sigma,const vector<Point> & points);
    vector<Descriptor> build() const;

private:


    vector<double> computeData(const Point &p,
            const double rotateAngle, const int sizeHist, const int numBins) const;

    Descriptor filterTrash(const Descriptor &descr) const;
    Descriptor normilize(const Descriptor & descr) const;
    Descriptor normilizeAll(const Descriptor & descr) const;
};

#endif // DESCRIPTOR_H
