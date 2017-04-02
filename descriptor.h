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
    const double sigma = 2.0;

    const int sizeArea = 16;
    const int sizeHist = 4;
    const int numBins = 8;

public:
    vector<Descriptor> build(const Matrix & m) const;

private:
    vector<double> computeData(const Point &p,
                   const Matrix & sobelX, const Matrix & sobelY) const;

    Descriptor filterTrash(const Descriptor &descr) const;
    Descriptor normilize(const Descriptor & descr) const;
};

class PointSearcher{

    const double eps = 0.05;

public:
    vector<pair<Point,Point>> findSamePoints( const Matrix &m1,const Matrix &m2) const;
};



#endif // DESCRIPTOR_H
