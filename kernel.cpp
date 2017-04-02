#include "kernel.h"
#include <stdio.h>


Kernel::Kernel():width(0),height(0),matrix(nullptr){}

Kernel::Kernel(int width,int height):width(width),height(height),
    matrix(make_unique<double []>(width*height)){}

Kernel::Kernel(int width, int height, unique_ptr<double[]> matrix):
    width(width),height(height),matrix(move(matrix)){}

Kernel KernelFactory::createGauss(double sigma)
{
    int kSize = (int)round(sigma * 3);
    Kernel kernel(kSize * 2 + 1,kSize * 2 + 1);
    int idx = 0;
    for (int x = -kSize; x <= kSize; x++){
        for (int y = -kSize; y <= kSize; y++){
            kernel.matrix[idx++] = Utils::gauss(x,y,sigma);
        }
    }
    return kernel;
}

Kernel KernelFactory::createGaussX(double sigma)
{
    int kSize = (int) round(sigma * 3);
    Kernel kernel(kSize*2+1,1);
    int idx = 0;
    double sum = 0;
    for(int i = -kSize; i <= kSize; i++){
        kernel.matrix[idx++] = Utils::gauss(i,sigma);
        sum += kernel.matrix[idx - 1];
    }
    for(int i = 0; i<kernel.width; i++){
        kernel.matrix[i] /= sum;
    }
    return kernel;
}

Kernel KernelFactory::createGaussY(double sigma)
{
    Kernel kernel = createGaussX(sigma);
    swap(kernel.height,kernel.width);
    return kernel;
}

Kernel KernelFactory::sobelX()
{
    return Kernel(3,3,unique_ptr<double[]>(new double[9]{-1,0,1,-2,0,2,-1,0,1}));
}

Kernel KernelFactory::sobelY()
{
    return Kernel(3,3,unique_ptr<double[]>(new double[9]{-1,-2,-1,0,0,0,1,2,1}));
}
