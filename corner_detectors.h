#ifndef CORNER_DETECTORS_H
#define CORNER_DETECTORS_H
#include "matrix.h"
#include "kernel.h"
#include "utils.h"
#include <cmath>
#include <limits>
#include <functional>
#include "kernel.h"
#include <vector>

using namespace std;
using namespace std::placeholders;

struct Point{
    int x;
    int y;
    float value;
    Point(const int x,const int y,const float v):x(x),y(y),value(v){};
};

enum Algorithm {MORAVEC, HARIS};

class CornerDetectors
{
    const float E = 1e-5;
    const float sigma = 1;

    const float trasholdMor = 0.05;
    const float trasholdHar = 0.25;
    const int winSize = 3;
    const int pSize = 3;

public:
    CornerDetectors();
    vector<Point> detect(const Matrix & m,const Algorithm alg = MORAVEC) const;

private:
    Matrix detectMoravec(const Matrix & m) const;
    Matrix detectHaris(const Matrix & m) const;
    vector<Point> localMinimums(const Matrix & m, const float tr) const;
    bool isMinimum(const Matrix & m,const int i,const int j, const float tr) const;
};

class PointFileter{
    const float factor = 0.9;
    int maxR;
    int maxPoints;

public:
    PointFileter(const int maxR, const int maxPoints);
    vector<Point> filter(const vector<Point> & points) const;
};


#endif // CORNER_DETECTORS_H
