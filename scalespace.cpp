#include "scalespace.h"
#include <QImage>
#include <corner_detectors.h>

ScaleSpace::ScaleSpace(const int numOctavs):
    _octavs(vector<vector<ScaleLevel>>(numOctavs)){}

ScaleSpace::ScaleSpace(const Matrix &initMatrix,const int octavNum):
    ScaleSpace(octavNum)
{
    calculate(initMatrix);
}

int ScaleSpace::octaveSize() const
{
    return _layerSize;
}

const vector<vector<ScaleLevel> > &ScaleSpace::octavs() const
{
    return _octavs;
}

ScaleSpace ScaleSpace::computeDiffs() const
{
    ScaleSpace diffs(_octavs.size());
    for(int octav = 0; octav < _octavs.size(); octav++){
        for(int layer = 0; layer < _layerSize - 1; layer++){
            ScaleLevel diff;
            auto & m1 = _octavs[octav][layer].matrix;
            auto & m2 = _octavs[octav][layer+1].matrix;
            diff.matrix = m1.compute(m2,std::minus<double>());
            diff.sigma = _octavs[octav][layer].sigma;
            diff.efectSigma = _octavs[octav][layer].efectSigma;
            diffs._octavs[octav].emplace_back(move(diff));
        }
    }
    return diffs;
}

vector<Blob> ScaleSpace::searchBlobs() const
{
    vector<Blob> result;
    for(int octav = 0; octav < _octavs.size(); octav++){
        for(int layer = 1; layer < _octavs[octav].size()-1; layer++){
            auto & oLayer = _octavs[octav][layer];
            for(int i = 0; i<_octavs[octav][layer].matrix.width(); i++){
                for(int j = 0; j < oLayer.matrix.height(); j++){
                    int type = checkExtremum(octav,layer,i,j);
                    bool min = type == -1;
                    bool max = type == 1;
                    if(min || max){
                        result.emplace_back(
                            Blob(i, j,
                                 layer, oLayer.sigma,oLayer.efectSigma, min));
                    }
                }
            }
        }
    }
    return result;
}

int ScaleSpace::checkExtremum(const int octav,
                              const int level, const int i, const int j) const
{
    bool min = true;
    bool max = true;
    double center = _octavs[octav][level].matrix.get(i,j);
    for(int x = -1; x < 2; x++){
        for(int y = -1; y < 2; y++){
            for(int z = -1; z < 2; z++){
                if(x == 0 && y == 0 && z == 0)
                    continue;
                double value = _octavs[octav][level+z].matrix.
                                get(i+x,j+y,Matrix::Border::COPIED);
                if(value < 1e-6)
                    continue;//считаем что тут мусор
                if(value > center) max = false;
                if(value < center) min = false;
                if(!min && !max)
                    return 0;
            }
        }
    }
    if(min && max)
        return 0;
    if(min)
        return -1;
    return 1;
}





Matrix ScaleSpace::gauss(const Matrix &matrix,const double sigma) const
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
                double delta = sqrt(pow(level.sigma, 2) - pow(prevL.sigma, 2));
                level.matrix = move(gauss(prevL.matrix, delta));
            }
            _octavs[i].emplace_back(move(level));
        }
    }
}


