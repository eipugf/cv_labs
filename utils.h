#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <memory>

using namespace std;

typedef unsigned char byte;

class Utils
{
public:

    static function<float(float,float)> sqldiff();

    static float gray(byte r, byte g, byte b);

    static float gauss(float x, float y, float sigma);
    static float gauss(float x, float sigma);
};

#endif // UTILS_H
