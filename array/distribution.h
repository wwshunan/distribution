#ifndef GUARD_distribution_h
#define GUARD_distribution_h

#include <vector>
#include <TApplication.h>

class Distr_base_2D {
public:
    Distr_base_2D(): array1(0), array2(0), num(0), emit(0), R(0), dRdZ(0) {}
    Distr_base_2D(int n, double e, double r, double div): array1(0), array2(0), num(n), emit(e), R(r), dRdZ(div) {}
    TH2F* display(TApplication*);
    //virtual void test() const=0;
protected:
    double* array1;
    double* array2;
    size_t num;
    double emit, R, dRdZ;
};

class Guass_2D: public Distr_base_2D {
public:
    Guass_2D() {}
    Guass_2D(int n, double e, double r, double div): Distr_base_2D(n, e, r, div) { create(); }
    ~Guass_2D();

private:
    void create();
};

#endif



    
