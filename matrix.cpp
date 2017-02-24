#include "matrix.h"

Matrix::Matrix():_width(0),_high(0),matrix(nullptr){}

Matrix::Matrix(const int width,const int high):
    _width(width),_high(high),matrix(make_unique<float[]>(_width*_high)){}

Matrix::Matrix(Matrix && picture):_width(picture._width),
    _high(picture._high),matrix(move(picture.matrix))
{
    picture._high = picture._width = 0;
}

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
    Matrix normalized(_width,_high);
    auto minMax = minmax_element(&matrix[0],&matrix[_width*_high]);
    float range = (*minMax.second - *minMax.first);
    for(int i = 0; i < _width*_high; i++){
        normalized.matrix[i] = (matrix[i] - *minMax.first)/range;
    }
    return normalized;
}

Matrix Matrix::denormalize() const
{
    Matrix denormalized(_width,_high);
    for(int i = 0; i<_width*_high; i++){
        denormalized.matrix[i] = matrix[i]*255.0;
    }
    return denormalized;
}

Matrix Matrix::canvolution(const Kernel &kernel, const Border border) const
{
    Matrix processed(_width,_high);
    for(int i = 0; i < _width; i++){
        for(int j = 0; j < _high; j++){
            processed.set(i,j,brightness(i,j,kernel,border));
        }
    }
    return processed;
}

float Matrix::brightness(const int elI, const int elJ, const Kernel &k,
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

Matrix Matrix::compute(Matrix &other, function<float (float, float)> funk)
{
    Matrix computed(_width,_high);
    for(int i = 0; i<_width*_high; i++){
        computed.matrix[i] = funk(matrix[i],other.matrix[i]);
    }
    return computed;
}

Matrix & Matrix::operator=(Matrix &&other)
{
    _width = other._width;
    _high = other._high;
    matrix = move(other.matrix);
    other._high = other._width = 0;
    return *this;
}
