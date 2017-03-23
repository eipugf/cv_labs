#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H
#include "corner_detectors.h"

using namespace std;

struct Descriptor{
    int x;
    int y;
    vector<float> data;
    Descriptor(const int x,const int y);
    Descriptor(const int x, const int y, const vector<float> & data);
};


class DescrBuilder{
    const float treshold = 0.2;
    const float sigma = 2;

    const int sizeArea = 16;
    const int sizeHist = 4;
    const int numBeans = 8;

public:
    vector<Descriptor> build(const Matrix & m) const;

private:
    vector<float> computeData(const Point &p,
                   const Matrix & sobelX, const Matrix & sobelY) const;

    Descriptor filterTrash(const Descriptor &descr) const;
    Descriptor normilize(const Descriptor & descr) const;
};

#endif // DESCRIPTOR_H
