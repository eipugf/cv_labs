#include "utils.h"

function<double (double, double)> Utils::hypotenuse()
{
    return [](double x, double y){return sqrt(x*x + y*y);};
}

function<double (double)> Utils::multiple(double num)
{
    return std::bind1st(std::multiplies<double>(),num);
}

function<double (double)> Utils::div(double num)
{
    return std::bind2nd(std::divides<double>(),num);
}

function<double (double)> Utils::square()
{
    return [](double x){return x*x;};
}

function<double (double)> Utils::normalize(double min, double range)
{
    return std::bind([](double a,double b, double c)
    {return (a-b)/c;},_1,min,range);
}

double Utils::gray(byte r, byte g, byte b)
{
    return (r*0.299 + g*0.587 + b*0.114)/255.0;
}

double Utils::gauss(double x, double y, double sigma)
{
    return exp(-(x*x+y*y)/(2*sigma*sigma)) / 2 * M_PI * sigma * sigma;
}

double Utils::gauss(double x, double sigma)
{
    return exp(-x*x/(2*sigma*sigma)) / sqrt(2 * M_PI) * sigma;
}
