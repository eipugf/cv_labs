#ifndef SCALESPACE_H
#define SCALESPACE_H
#include "matrix.h"
#include <math.h>
#include <memory>
#include <vector>
#include <string>

using namespace std;

struct ScaleLevel{
    Matrix matrix;
    double sigma;
    double efectSigma;
    ScaleLevel():sigma(0),efectSigma(0),matrix(Matrix(0,0)){}
    ScaleLevel(ScaleLevel &&) = default;
    ScaleLevel(Matrix & m,double sigma,double efectSigma):
        sigma(sigma),efectSigma(efectSigma),matrix(move(m)){}
};

struct Blob{
    int level;
    int x;
    int y;
    bool min;
    double sigma;
    double efectSigma;
    Blob():min(false),level(0),x(0),y(0),sigma(0),efectSigma(0){}
    Blob(const int x, const int y, const int level,
            const double sigma, const double efectSigma, const bool min):
        level(level),x(x),y(y),sigma(0),efectSigma(efectSigma),min(min){}
};

class ScaleSpace
{
public:

private:
    const int _layerSize = 9;
    const double startSigma = 0.5;
    const double sigmaA = 1.6;
    const double k = pow(2,1/(double)_layerSize);

    vector<vector<ScaleLevel>> _octavs;

public:
    ScaleSpace(const int numOctavs);
    ScaleSpace(const Matrix &initMatrix,const int numOctavs);

    int octaveSize() const;

    const vector<vector<ScaleLevel>> & octavs() const;

    ScaleSpace computeDiffs() const;

    vector<Blob> searchBlobs() const;

private:
    int checkExtremum(const int octav, const int level, const int i, const int j) const;
    Matrix gauss(const Matrix & matrix,const double sigma) const;
    void calculate(const Matrix & m);
};

#endif // SCALESPACE_H
