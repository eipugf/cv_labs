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

    static function<double(double,double)> hypotenuse();

    static function<double(double)> multiple(double num);
    static function<double (double)> div(double num);

    static function<double(double)> square();

    static function<double(double)> normalize(double min, double range);

    static double gray(byte r, byte g, byte b);

    static double gauss(double x, double y, double sigma);
    static double gauss(double x, double sigma);

};

#endif // UTILS_H
