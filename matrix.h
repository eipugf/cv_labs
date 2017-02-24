#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <kernel.h>
#include <math.h>
#include <memory>

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
    Matrix(Matrix &&matrix);

    int width() const;
    int hight() const;

    float get(const int i,const int j,const Border border = DEFAULT) const;
    void set(const int i,const int j,const float gray);

    Matrix normalize() const;
    Matrix denormalize() const;

    Matrix canvolution(const Kernel & kernel, const Border border = DEFAULT) const;

    float brightness(const int i,const int j,
                     const Kernel & kernel, const Border border = DEFAULT) const;

    Matrix compute(Matrix &, function<float(float,float)>);

    Matrix& operator=(Matrix&& other);
};

#endif // MATRIX_H
