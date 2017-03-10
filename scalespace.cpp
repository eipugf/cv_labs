#include "scalespace.h"
#include <QImage>

ScaleSpace::ScaleSpace(const Matrix &initMatrix,const int octavNum):
    _octavs(vector<vector<ScaleLevel>>(octavNum))
{
    calculate(initMatrix);
}

int ScaleSpace::octaveSize() const
{
    return _layerSize;
}

const vector<vector<ScaleSpace::ScaleLevel> > &ScaleSpace::octavs() const
{
    return _octavs;
}


Matrix ScaleSpace::gauss(const Matrix &matrix,const float sigma) const
{
    return matrix
           .convolution(KernelFactory::createGaussX(sigma),Matrix::Border::CILINDER)
           .convolution(KernelFactory::createGaussY(sigma),Matrix::Border::CILINDER);
}

void ScaleSpace::calculate(const Matrix & m)
{
    ScaleLevel initLevel;
    initLevel.sigma = sigmaA;
    initLevel.matrix = move(gauss(m, sqrt(pow(sigmaA,2) - pow(startSigma,2))));
    initLevel.efectSigma = sigmaA;
    _octavs[0].emplace_back(move(initLevel));
    for(int i = 0; i < _octavs.size(); i++){
        for(int j = 0; j < _layerSize; j++){
            ScaleLevel level;
            if(j == 0){
                if(i == 0)
                    continue;
                level.sigma = sigmaA;
                level.efectSigma =k* _octavs[i-1][_layerSize -1].efectSigma;
                level.matrix = move(_octavs[i-1][_layerSize -1].matrix.compress());
            } else {
                auto &prevL = _octavs[i][j - 1];
                level.sigma = k * prevL.sigma;
                level.efectSigma = k * prevL.efectSigma;
                float delta = sqrt(pow(level.sigma, 2) - pow(prevL.sigma, 2));
                level.matrix = move(gauss(prevL.matrix, delta));
            }
            _octavs[i].emplace_back(move(level));
        }
    }
}

