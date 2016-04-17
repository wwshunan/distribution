#ifndef GUARD_distribution_h
#define GUARD_distribution_h

#include <vector>
#include <TApplication.h>

class Distr_base_2D {
public:
    Distr_base_2D(): num(0), emit(0), R(0), dRdZ(0) {}
    Distr_base_2D(int n, double e, double r, double div): num(n), emit(e), R(r), dRdZ(div) {}
    virtual void create()=0;
    TH2F* display(TApplication*);
    //virtual void test() const=0;
protected:
    std::vector<double> array1;
    std::vector<double> array2;
    size_t num;
    double emit, R, dRdZ;
};

class Guass_2D: public Distr_base_2D {
public:
    Guass_2D() {}
    Guass_2D(int n, double e, double r, double div): Distr_base_2D(n, e, r, div) { create(); }
    void create();
};

#endif



    
