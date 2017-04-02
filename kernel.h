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
    int height;
    unique_ptr<double []> matrix;

public:
    Kernel();
    Kernel(int width,int height);
    Kernel(int width,int height,unique_ptr<double []>);
    Kernel(Kernel &&) = default;
};

//---------------------------------------------------------------
//ядерная фабрика
//---------------------------------------------------------------
class KernelFactory
{
public:
    static Kernel createGauss(double sigma);
    static Kernel createGaussX(double sigma);
    static Kernel createGaussY(double sigma);
    static Kernel sobelX();
    static Kernel sobelY();
};
#endif // KERNEL_H
