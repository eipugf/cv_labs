#include "scalespace.h"
#include <QImage>

ScaleSpace::ScaleSpace(const Matrix &initMatrix,const int numOctavs):
    _numOctavs(numOctavs),_octavs(vector<vector<ScaleLevel>>(numOctavs))
{
    calculate(initMatrix);
}

vector<vector<ScaleSpace::ScaleLevel> > &ScaleSpace::octavs()
{
    return _octavs;
}

int ScaleSpace::numOctavs() const
{
    return _numOctavs;
}

int ScaleSpace::layerSize() const
{
    return _layerSize;
}

Matrix ScaleSpace::gauss(const Matrix &matrix,const float sigma) const
{
    return matrix.convolution(KernelFactory::createGaussX(sigma),Matrix::Border::CILINDER)
           .convolution(KernelFactory::createGaussY(sigma),Matrix::Border::CILINDER)
           .normalize();
}

void ScaleSpace::calculate(const Matrix & m)
{
    ScaleLevel initLevel;
    initLevel.sigma = sigmaA;
    initLevel.matrix = move(gauss(m, startSigma));
    _octavs[0].emplace_back(move(initLevel));
    for(int i = 0; i < _numOctavs; i++){
        for(int j = 0; j < _layerSize; j++){
            ScaleLevel level;
            if(j == 0){
                if(i == 0)
                    continue;
                level.sigma = _octavs[i - 1][_layerSize - 1].sigma * 2;
                level.matrix = move(_octavs[i-1][_layerSize -1].matrix.compress());
            } else {
                auto &prevL = _octavs[i][j - 1];
                level.sigma = k * prevL.sigma;
                float delta = sqrt(pow(level.sigma, 2) - pow(prevL.sigma, 2));
                level.matrix = move(gauss(prevL.matrix, delta));
            }
            _octavs[i].emplace_back(move(level));
        }
    }
}

