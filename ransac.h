#ifndef RANSAC_H
#define RANSAC_H
#include <vector>
#include "descriptor.h"


using namespace std;

class Ransac
{

    const int M_WIDTH = 9;
    const int COUNT = 4;
    const int LIMIT_ITER = 2000;
    const double EPS = 5;

    vector<double> computeTrans(
        const vector<pair<Point,Point>> & inliers,
        const vector<int> & idx, const int size) const;
    vector<double> normalize(const vector<double> & v, const double r) const;
    double computeError(const vector <double> & h,
            const pair<Point,Point> & pairs, const int idx) const;
public:
    Ransac();

    vector<double> searchTransform(
            const vector<pair<Point,Point>> & pairs) const;
};

#endif // RANSAC_H
