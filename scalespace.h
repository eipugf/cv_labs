#ifndef SCALESPACE_H
#define SCALESPACE_H
#include "matrix.h"
#include <math.h>
#include <memory>
#include <vector>
#include <string>

using namespace std;

class ScaleSpace
{
public:
    struct ScaleLevel{
        Matrix matrix;
        float sigma;
        float efectSigma;
        ScaleLevel():sigma(0),efectSigma(0),matrix(Matrix(0,0)){}
        ScaleLevel(ScaleLevel &&) = default;
        ScaleLevel(Matrix & m,float sigma,float efectSigma):
            sigma(sigma),efectSigma(efectSigma),matrix(move(m)){}
    };
private:
    const int _layerSize = 7;
    const float startSigma = 0.5;
    const float sigmaA = 1.6;
    const float k = pow(2,1/(float)_layerSize);

    vector<vector<ScaleLevel>> _octavs;

public:
    ScaleSpace(const Matrix &initMatrix,const int numOctavs);

    int octaveSize() const;

    const vector<vector<ScaleLevel>> & octavs() const;

private:
    Matrix gauss(const Matrix & matrix,const float sigma) const;
    void calculate(const Matrix & m);
};

#endif // SCALESPACE_H
