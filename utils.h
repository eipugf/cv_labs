#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <memory>
#include <functional>

using namespace std;
using namespace std::placeholders;

typedef unsigned char byte;

class Utils
{
public:

    static function<float(float,float)> hypotenuse();

    static function<float(float)> multiple(float num);

    static function<float(float)> normalize(float min, float range);

    static float gray(byte r, byte g, byte b);

    static float gauss(float x, float y, float sigma);
    static float gauss(float x, float sigma);


};


#endif // UTILS_H
