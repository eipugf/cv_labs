#include "descriptor.h"

Descriptor::Descriptor(const int x, const int y):
    x(x),y(y){}

Descriptor::Descriptor(const int x, const int y, const vector<float> &data):
    Descriptor(x,y)
{
    for(float each:data){
        this->data.push_back(each);
    }

}

vector<Descriptor> DescrBuilder::build(const Matrix &m) const
{
    auto points = CornerDetectors().detect(m,Algorithm::HARIS);
    points = PointFileter(m.width()*m.height(),10).filter(points);
    auto sobleX = m.convolution(KernelFactory::sobelX(),Matrix::Border::CILINDER);
    auto sobleY = m.convolution(KernelFactory::sobelY(),Matrix::Border::CILINDER);
    auto descriptors = vector<Descriptor>();
    for(auto &each : points){
        auto descr = Descriptor(each.x,each.y,move(computeData(each,sobleX,sobleY)));
        descr = move(normilize(descr));
        descr = move(filterTrash(descr));
        descr = move(normilize(descr));
        descriptors.emplace_back(move(descr));
    }
    return descriptors;
}

vector<float> DescrBuilder::computeData(const Point &p,
                             const Matrix & sobelX,const Matrix & sobelY) const
{
    int size = sizeArea/sizeHist;
    auto data = vector<float>(size*size * numBeans, 0.0);
    for(int i = 0; i < sizeArea; i++){
        for(int j = 0; j< sizeArea; j++){
            float dx = sobelX.get(p.x + i - sizeArea/2,
                          p.y + j - sizeArea/2, Matrix::Border::CILINDER);
            float dy = sobelY.get(p.x + i - sizeArea/2,
                          p.y + j - sizeArea/2, Matrix::Border::CILINDER);

            float magnitud = sqrt(dx*dx+dy*dy)*
                        Utils::gauss(i - sizeArea/2, j - sizeArea/2, sigma);

            double phi = atan2(dx,dy);

            float beanIdx = (phi/M_PI + 1)*numBeans/2.0;

            int hIdx = ((i/sizeHist)*size+j/sizeHist)*numBeans+(int)beanIdx;

            data[hIdx] = magnitud*(beanIdx - ((int)beanIdx));
            data[(hIdx+1)%numBeans] = magnitud*(((int)beanIdx+1) - beanIdx);
        }
    }
    return data;
}

Descriptor DescrBuilder::filterTrash(const Descriptor &descr) const
{
    Descriptor result(descr.x,descr.y);
    for(float bin:descr.data){
        result.data.push_back(min(treshold, bin));
    }
    return result;
}

Descriptor DescrBuilder::normilize(const Descriptor &descr) const
{
    Descriptor result(descr.x,descr.y);
    float sum = 0;
    for(float each:descr.data){
        sum += each*each;
    }
    sum = sqrt(sum);
    for(float each:descr.data){
        result.data.push_back(each/sum);
    }
    return result;
}
