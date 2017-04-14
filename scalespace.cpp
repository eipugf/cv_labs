#include "scalespace.h"

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

const vector<vector<ScaleLevel> > & ScaleSpace::octavs() const
{
    return _octavs;
}

ScaleSpace ScaleSpace::computeDiffs() const
{
    ScaleSpace diffs(_octavs.size());
    for(int octav = 0; octav < _octavs.size(); octav++){
        for(int layer = 0; layer < _layerSize + 2; layer++){
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
                    if(checkExtremum(octav,layer,i,j)){
                        result.emplace_back(
                            Blob(i, j, octav, layer, oLayer.sigma));
                    }
                }
            }
        }
    }
    return result;
}

bool ScaleSpace::checkExtremum(const int octav,
                              const int level, const int i, const int j) const
{
    bool min = true;
    bool max = true;
    double center = _octavs[octav][level].matrix.get(i,j);
    if(std::abs(center) < EPS)
        return false;
    for(int x = -1; x < 2; x++){
        for(int y = -1; y < 2; y++){
            for(int z = -1; z < 2; z++){
                if(!(x || y || z))
                    continue;
                double value = _octavs[octav][level+z].matrix.
                                get(i+x,j+y,Matrix::Border::COPIED);
                if(value > center) max = false;
                if(value < center) min = false;
                if(!(min || max))
                    return false;
            }
        }
    }
    return  min ^ max;
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
        for(int j = 0; j < _layerSize+3; j++){
            ScaleLevel level;
            if(j == 0){
                if(i == 0)
                    continue;
                level.sigma = sigmaA;
                level.efectSigma =k* _octavs[i-1][_layerSize -1].efectSigma;
                level.matrix = _octavs[i-1][_layerSize - 1].matrix.compress();
            } else {
                auto &prevL = _octavs[i][j - 1];
                level.sigma = k * prevL.sigma;
                level.efectSigma = k * prevL.efectSigma;
                double delta = sqrt(pow(level.sigma, 2) - pow(prevL.sigma, 2));
                level.matrix = gauss(prevL.matrix, delta);
            }
            _octavs[i].emplace_back(move(level));
        }
    }
}

vector<Blob> BlobFilter::filter(const vector<Blob> & blobs, ScaleSpace & space) const
{
    vector<Blob> result;
    const auto & oct = space.octavs();
    for(const Blob & each:blobs){
        double k = pow(2,each.octav);
        auto haris = CornerDetectors().computeHaris(
             oct[each.octav][each.layer].matrix,each.x,each.y,each.sigma*M_SQRT2);
        if(haris > treshold){
            result.emplace_back(Blob(each));
        }
    }
    return result;
}


vector<Descriptor> SIDiscrBuilder::build(const Matrix &m)
{
    auto space = ScaleSpace(m,5);
    auto blobs = BlobFilter().filter(space.computeDiffs().searchBlobs(),space);
    int curOct = -1;
    int curLayer = -1;
    vector<Point> points;
    vector<Descriptor> descriptors;
    for(const Blob & each:blobs){
        if(curLayer < 0 && curOct < 0){
            curLayer = each.layer;
            curOct = each.octav;
        }
        if(curOct != each.octav || curLayer != each.layer){
            auto descr = DescrBuilder(space.octavs()[curOct][curLayer].matrix,
                               space.octavs()[curOct][curLayer].sigma,points).build();
            descriptors.insert(descriptors.end(),descr.begin(),descr.end());
            points.clear();
            curOct = each.octav;
            curLayer = each.layer;
        }
        points.emplace_back(Point(each.x, each.y, 0));
    }
    return descriptors;
}

PointMatcher::PointMatcher(const double eps):eps(eps){}

vector<pair<Point, Point> > PointMatcher::match(
        const Matrix &m1, const Matrix &m2, const bool withScale) const
{
    auto descrM1 = withScale?SIDiscrBuilder::build(m1):DescrBuilder(m1).build();
    auto descrM2 = withScale?SIDiscrBuilder::build(m2):DescrBuilder(m2).build();
    vector<pair<Point, Point>> samePoints;
    for(Descriptor & each:descrM1){
        int idx = 0;
        for(Descriptor & each1:descrM2){
             double sum = 0;
             for(int i = 0; i<each.data.size(); i++){
                 sum += (each.data[i]-each1.data[i])*(each.data[i]-each1.data[i]);
             }
             sum = sqrt(sum);
             if(sum < eps){
                 samePoints.emplace_back(
                      pair<Point,Point>(
                              Point(each.x,each.y,0),Point(each1.x,each1.y,0)));
                 descrM2.erase(descrM2.begin()+idx);
                 break;
             }
             idx++;
        }
    }
    return samePoints;
}
