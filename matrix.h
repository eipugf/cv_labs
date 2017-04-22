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
    int _height;
    unique_ptr<double[]> matrix;

public:

    enum Border{DEFAULT,SIMPLE,CILINDER,COPIED,MIRROR};

    Matrix();
    Matrix(const int width,const int height);
    Matrix(const int width,const int height, const unique_ptr<double[]> matrix);
    Matrix(Matrix &&matrix) = default;

    int width() const;
    int height() const;

    double get(const int i,const int j,const Border border = DEFAULT) const;
    void set(const int i,const int j,const double gray);

    Matrix normalize() const;

    Matrix convolution(const Kernel & kernel, const Border border = DEFAULT) const;
    Matrix compress() const;

    Matrix compute(const Matrix & other,const function<double(double,double)> & funk) const;
    Matrix compute(function<double(double)> & funk) const;

    Matrix& operator=(Matrix&& other) = default;

private:
    double convoluite(const int i,const int j,
                     const Kernel & kernel, const Border border = DEFAULT) const;
};

#endif // MATRIX_H
