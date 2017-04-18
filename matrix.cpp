#include "matrix.h"

Matrix::Matrix():_width(0),_height(0),matrix(nullptr){}

Matrix::Matrix(const int width,const int height):
    _width(width),_height(height),matrix(make_unique<double[]>(_width*_height)){}

int Matrix::width() const{
    return _width;
}

int Matrix::height() const{
    return _height;
}

double Matrix::get(const int i,const int j, const Border border) const
{
    switch (border) {
    case SIMPLE:
        return i < 0 || j < 0 || i >= _width || j >= _height?0:matrix[i*_height + j];
    case COPIED:
        return matrix[min(max(i, 0), _width - 1)*_height + min(max(j,0),_height-1)];
    case MIRROR:
        return matrix[(i >= _width ? 2*_width - i - 2 : abs(j))*_height +
                       (j >= _height ? 2*_height - j - 2 : abs(j))];
    case CILINDER:
        return matrix[((i + _width) % _width) * _height + (j + _height) % _height];
    default:
        return matrix[i*_height + j];
    }
}

void Matrix::set(const int i,const int j,const double gray)
{
    matrix[i*_height + j] = gray;
}

Matrix Matrix::normalize() const
{
    auto minMax = minmax_element(&matrix[0],&matrix[_width*_height]);
    double range = (*minMax.second - *minMax.first);
    auto funk = Utils::normalize(*minMax.first,range);
    return compute(funk);
}

Matrix Matrix::convolution(const Kernel &kernel, const Border border) const
{
    Matrix processed(_width,_height);
    for(int i = 0; i < _width; i++){
        for(int j = 0; j < _height; j++){
            processed.set(i,j,convoluite(i,j,kernel,border));
        }
    }
    return processed;
}

Matrix Matrix::compress() const
{
    Matrix compresed(_width/2,_height/2);
    for(int i = 0; i < compresed._width; i++){
        for(int j = 0; j < compresed._height; j++){
            double average = (get(i*2,j*2,COPIED)+
                            get(i*2+1,j*2,COPIED)+
                            get(i*2,j*2+1,COPIED)+
                            get(i*2+1,j*2+1,COPIED))/4.0;
            compresed.set(i, j, average);
        }
    }
    return compresed;
}

double Matrix::convoluite(const int elI, const int elJ, const Kernel &k,
                         const Border border) const
{
    double brightness = 0;
    for(int i = 0; i < k.width; i++){
        for(int j = 0; j < k.height; j++){
            int posX = elI + (i - (k.width / 2));
            int posY = elJ + (j - (k.height / 2));
            brightness += get(posX, posY,border)*k.matrix[i * k.height + j];
        }
    }
    return brightness;
}

Matrix Matrix::compute(function<double (double)> &funk) const
{
    Matrix computed(_width,_height);
    transform(matrix.get(),matrix.get()+_width*_height,computed.matrix.get(),funk);
    return computed;
}

Matrix Matrix::compute(const Matrix &other,
                       const function<double (double, double)> & funk) const
{
    assert(_width == other._width && _height == other._height);
    Matrix computed(_width,_height);
    transform(matrix.get(),matrix.get()+_width*_height,
                   other.matrix.get(),computed.matrix.get(),funk);
    return computed;
}
