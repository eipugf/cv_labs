#ifndef KERNEL_H
#define KERNEL_H

#include <kernel.h>
#include <algorithm>
#include <memory>
#include "utils.h"

using namespace std;

struct Kernel
{
public:
    int width;
    int high;
    unique_ptr<float []> matrix;

public:
    Kernel();
    Kernel(int width,int high);
    Kernel(int width,int high,unique_ptr<float []>);
    Kernel(Kernel &&) = default;
};

//---------------------------------------------------------------
//ядерная фабрика
//---------------------------------------------------------------
class KernelFactory
{
public:
    static Kernel createGauss(float sigma);
    static Kernel createGaussX(float sigma);
    static Kernel createGaussY(float sigma);
    static Kernel sobelX();
    static Kernel sobelY();

    static Kernel testKernell();
};
#endif // KERNEL_H
