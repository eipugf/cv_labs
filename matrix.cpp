#include "matrix.h"

Matrix::Matrix():_width(0),_high(0),matrix(nullptr){}

Matrix::Matrix(const int width,const int high):
    _width(width),_high(high),matrix(make_unique<float[]>(_width*_high)){}

int Matrix::width() const{
    return _width;
}

int Matrix::hight() const{
    return _high;
}

float Matrix::get(const int i,const int j, const Border border) const
{
    switch (border) {
    case SIMPLE:
        return i < 0 || j < 0 || i >= _width || j >= _high?0:matrix[i*_high + j];
    case COPIED:
        return matrix[min(max(i, 0), _width - 1)*_high + min(max(j,0),_high-1)];
    case MIRROR:
        return matrix[(i >= _width ? 2*_width - i - 2 : abs(j))*_high +
                       (j >= _high ? 2*_high - j - 2 : abs(j))];
    case CILINDER:
        return matrix[((i + _width) % _width) * _high + (j + _high) % _high];
    default:
        return matrix[i*_high + j];
    }
}

void Matrix::set(const int i,const int j,const float gray)
{
    matrix[i*_high + j] = gray;
}

Matrix Matrix::normalize() const
{
    auto minMax = minmax_element(&matrix[0],&matrix[_width*_high]);
    float range = (*minMax.second - *minMax.first);
    auto funk = Utils::normalize(*minMax.first,range);
    return compute(funk);
}

Matrix Matrix::convolution(const Kernel &kernel, const Border border) const
{
    Matrix processed(_width,_high);
    for(int i = 0; i < _width; i++){
        for(int j = 0; j < _high; j++){
            processed.set(i,j,convoluite(i,j,kernel,border));
        }
    }
    return processed;
}

Matrix Matrix::compress() const
{
    Matrix compresed(_width/2,_high/2);
    for(int i = 0; i < compresed._width; i++){
        for(int j = 0; j < compresed._high; j++){
            float average = (get(i*2,j*2,CILINDER)+
                            get(i*2+1,j*2,CILINDER)+
                            get(i*2,j*2+1,CILINDER)+
                            get(i*2+1,j*2+1,CILINDER))/4.0;
            compresed.set(i, j, average);
        }
    }
    return compresed;
}

float Matrix::convoluite(const int elI, const int elJ, const Kernel &k,
                         const Border border) const
{
    float brightness = 0;
    for(int i = 0; i < k.width; i++){
        for(int j = 0; j < k.high; j++){
            int posX = elI + (i - (k.width / 2));
            int posY = elJ + (j - (k.high / 2));
            brightness += get(posX, posY,border)*k.matrix[i * k.high + j];
        }
    }
    return brightness;
}

Matrix Matrix::compute(function<float (float)> &funk) const
{
    Matrix computed(_width,_high);
    transform(matrix.get(),matrix.get()+_width*_high,computed.matrix.get(),funk);
    return computed;
}

Matrix Matrix::compute(const Matrix &other,
                       const function<float (float, float)> & funk) const
{
    assert(_width == other._width && _high == other._high);
    Matrix computed(_width,_high);
    transform(matrix.get(),matrix.get()+_width*_high,
                   other.matrix.get(),computed.matrix.get(),funk);
    return computed;
}
