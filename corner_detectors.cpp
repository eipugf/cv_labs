#include "corner_detectors.h"

CornerDetectors::CornerDetectors(){}

vector<Point> CornerDetectors::detect(const Matrix &m,
                                      const Algorithm alg) const
{
    Matrix matrix;
    switch (alg) {
    case MORAVEC:
        matrix = std::move(detectMoravec(m));
        break;
    case HARIS:
        matrix = std::move(detectHaris(m));
    default:
        break;
    }
    return localMinimums(matrix);
}

Matrix CornerDetectors::detectMoravec(const Matrix &m) const
{
    auto errors = Matrix(m.width(),m.hight());
    for(int i = 0; i<m.width(); i++){
        for(int j = 0; j<m.hight(); j++){
            float min = std::numeric_limits<int>::max();
            for(int u = -1; u <= 1; u++){
                for(int v = -1; v <= 1; v++){
                    float err = 0;
                    for(int dx = -winSize; dx < winSize; dx++){
                        for(int dy = -winSize; dy < winSize; dy++){
                            err += pow(m.get(i+u,j+v,Matrix::Border::COPIED) -
                                       m.get(i+u+dx,j+v+dy,Matrix::Border::COPIED),2);
                        }
                    }
                    min = std::min(min, err);
                }
            }
            errors.set(i, j, min);
        }
    }
    return errors;
}

Matrix CornerDetectors::detectHaris(const Matrix &m) const
{
    auto lambdas = Matrix(m.width(),m.hight());
    auto derX = m.convolution(KernelFactory::sobelX(),Matrix::Border::SIMPLE);
    auto derY = m.convolution(KernelFactory::sobelY(),Matrix::Border::SIMPLE);
    auto w = KernelFactory::createGauss(sigma);
    int wSize = w.width/2;
    for(int i = 0; i < lambdas.width(); i++){
        for(int j = 0; j < lambdas.hight(); j++){
            float A = 0, B = 0, C = 0;
            for(int v = 0; v < w.high; v++){
                for(int u = 0; u < w.width; u++){
                    float Ix = derX.get(i+u-wSize, j+v-wSize, Matrix::Border::SIMPLE);
                    float Iy = derY.get(i+u-wSize, j+v-wSize, Matrix::Border::SIMPLE);
                    float k = w.matrix[u*w.high+v];
                    A += k*Ix*Ix;
                    B += k*Ix*Iy;
                    C += k*Iy*Iy;
                }
            }
            float descr = std::sqrt(std::pow(A-C,2) + 4*B*B);
            lambdas.set(i, j, std::min(abs((A+C-descr)/2),abs((A+C+descr)/2)));
        }
    }
    return lambdas;
}

vector<Point> CornerDetectors::localMinimums(const Matrix &m) const
{
    vector<Point> points;
    for(int i = 0; i < m.width(); i++){
        for(int j = 0; j < m.hight(); j++){
            if(m.get(i, j) < trashold)
                continue;
            bool foundMore = false;
            for(int k = -pSize; k < pSize && !foundMore; k++){
                for(int t = -pSize; t < pSize && !foundMore; t++){
                    if(k == 0 || t == 0)
                        continue;
                    if(m.get(i+k,j+t,Matrix::Border::SIMPLE) - m.get(i,j) > E){
                        foundMore = true;
                    }
                }
            }
            if(foundMore)
                continue;
            points.emplace_back(move(Point(i, j, m.get(i,j))));
        }
    }
    return points;
}

PointFileter::PointFileter(const int maxR, const int maxPoints):
    maxR(maxR),maxPoints(maxPoints){}

vector<Point> PointFileter::filter(const vector<Point> &points) const
{
    vector<Point> result(points);
    for(int rad = 0; rad < maxR && result.size() > maxPoints; rad++){
        for(int i = 0; i < result.size(); i++){
            for(int j = i+1; j < result.size(); j++){
                if(sqrt(pow(result[i].x - result[j].x,2) +
                       pow(result[i].y - result[j].y,2)) <= rad &&
                        factor*result[i].value > result[j].value){
                    result.erase(result.begin() + j);
                    break;
                }
            }
        }
    }
    return result;
}
