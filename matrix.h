#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include "kernel.h"
#include <math.h>
#include <memory>
#include <cassert>

using namespace std;

class Matrix
{
    int _width;
    int _high;
    unique_ptr<float[]> matrix;

public:

    enum Border{DEFAULT,SIMPLE,CILINDER,COPIED,MIRROR};

    Matrix();
    Matrix(const int width,const int high);
    Matrix(const int width,const int high, const unique_ptr<float[]> matrix);
    Matrix(Matrix &&matrix) = default;

    int width() const;
    int hight() const;

    float get(const int i,const int j,const Border border = DEFAULT) const;
    void set(const int i,const int j,const float gray);

    Matrix normalize() const;

    Matrix convolution(const Kernel & kernel, const Border border = DEFAULT) const;
    Matrix compress() const;

    Matrix compute(const Matrix & other,const function<float(float,float)> & funk) const;
    Matrix compute(function<float(float)> & funk) const;

    Matrix& operator=(Matrix&& other) = default;

private:
    float convoluite(const int i,const int j,
                     const Kernel & kernel, const Border border = DEFAULT) const;
};

#endif // MATRIX_H
