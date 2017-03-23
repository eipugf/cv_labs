#include "kdthree.h"

KDThree::KDThree():root(nullptr){}

KDThree::Item::Item(const Descriptor &v):value(v),left(nullptr),right(nullptr){}

void KDThree::Item::add(const Descriptor & v,const unsigned depth)
{
    unsigned cd = depth % v.data.size();
    if(v.data[cd] < value.data[cd]){
        if(left != nullptr)
            left->add(v, depth+1);
        else left = make_unique<Item>(v);
    } else {
        if(right != nullptr)
            right->add(v, depth + 1);
        else right = make_unique<Item>(v);
    }
}

bool KDThree::Item::compare(const Descriptor &d1,
                            const Descriptor &d2,const float eps) const
{
    float sum = 0;
    for(int i = 0; i< d1.data.size(); i++){
        sum += pow((d1.data[i] - d2.data[i]),2);
    }
    return sqrt(sum) < eps;
}

unique_ptr<Descriptor> KDThree::Item::find(const Descriptor &d,
                                           const unsigned depth,const float eps) const
{
    if(compare(value,d, eps)){
        return make_unique<Descriptor>(value);
    }

    unsigned cd = depth % value.data.size();

    if (d.data[cd] < value.data[cd]){
        if(left!=nullptr)
            return left->find(d, depth + 1, eps);
        else return nullptr;
    } else {
        if(right != nullptr)
            return right->find(d, depth + 1, eps);
        else return nullptr;
    }
}

void KDThree::add(const Descriptor &v)
{
    if(root == nullptr)
        root = make_unique<Item>(v);
    else
        root->add(v,0);
}

unique_ptr<Descriptor> KDThree::findNear(const Descriptor &d,const float eps) const
{
    if(root != nullptr)
        return root->find(d,0,eps);
    return nullptr;
}

vector<pair<Point,Point> > CommonPoints::findSamePoints(
        const Matrix &m1,const Matrix &m2)
{
    auto descrM1 = DescrBuilder().build(m1);
    auto descrM2 = DescrBuilder().build(m2);

    KDThree three;
    for(Descriptor & each:descrM1){
        three.add(each);
    }

    vector<pair<Point, Point>> samePoints;
    for(Descriptor & each:descrM2){
        auto same = move(three.findNear(each, eps));
        if(same == nullptr)
            continue;
        samePoints.emplace_back(pair<Point,Point>(
                                    Point(same->x,same->y,0),Point(each.x,each.y,0)));
    }
    return samePoints;
}
