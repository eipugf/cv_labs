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
        double sigma;
        double efectSigma;
        ScaleLevel():sigma(0),efectSigma(0),matrix(Matrix(0,0)){}
        ScaleLevel(ScaleLevel &&) = default;
        ScaleLevel(Matrix & m,double sigma,double efectSigma):
            sigma(sigma),efectSigma(efectSigma),matrix(move(m)){}
    };
private:
    const int _layerSize = 7;
    const double startSigma = 0.5;
    const double sigmaA = 1.6;
    const double k = pow(2,1/(double)_layerSize);

    vector<vector<ScaleLevel>> _octavs;

public:
    ScaleSpace(const Matrix &initMatrix,const int numOctavs);

    int octaveSize() const;

    const vector<vector<ScaleLevel>> & octavs() const;

private:
    Matrix gauss(const Matrix & matrix,const double sigma) const;
    void calculate(const Matrix & m);
};

#endif // SCALESPACE_H
