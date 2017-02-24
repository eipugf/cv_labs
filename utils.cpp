#include "utils.h"

function<float (float, float)> Utils::sqldiff()
{
    return [](float x, float y){return sqrt(x*x + y*y);};
}

float Utils::gray(byte r, byte g, byte b)
{
    return (r*0.299 + g*0.587 + b*0.114) / 255.0;
}

float Utils::gauss(float x, float y, float sigma)
{
    return exp(-(x*x+y*y)/(2*sigma*sigma)) / 2 * M_PI * sigma * sigma;
}

float Utils::gauss(float x, float sigma)
{
    return exp(-x*x/(2*sigma*sigma)) / sqrt(2 * M_PI) * sigma;
}
