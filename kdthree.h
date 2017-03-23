#ifndef KDTHREE_H
#define KDTHREE_H

#include "descriptor.h"
#include <vector>

using namespace std;

class KDThree
{
    class Item{
        Descriptor value;
        unique_ptr<Item> left;
        unique_ptr<Item> right;
    public:
        Item(const Descriptor & v);
        void add(const Descriptor & v,const unsigned depth);
        unique_ptr<Descriptor> find(const Descriptor & d,
                                    const unsigned depth,const float eps) const;
        bool compare(const Descriptor &d1,const Descriptor &d2,const float eps) const;
    };

    unique_ptr<Item> root;

public:
    KDThree();
    void add (const Descriptor & value);
    unique_ptr<Descriptor> findNear(const Descriptor & d,const float eps) const;
};

class CommonPoints{
    const float eps = 0.01;
public:
    vector<pair<Point,Point>> findSamePoints(
                                    const Matrix &m1,const Matrix &m2);
};

#endif // KDTHREE_H
