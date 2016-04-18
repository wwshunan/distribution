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

class Distr_base_4D {
public:
    Distr_base_4D(): array1(0), array2(0), array3(0), array4(0), num(0), Ex(0), Rx(0), dRxdZ(0), Ey(0), Ry(0), dRydZ(0) {}
    Distr_base_4D(int n, double emitX, double envX, double divX, double emitY, double envY, double divY): array1(0), array2(0), array3(0), array4(0), num(n), Ex(emitX), Rx(envX), dRxdZ(divX), Ey(emitX), Ry(envY), dRydZ(divY) {}
    TH2F* display_x_px(TApplication*);
    TH2F* display_y_py(TApplication*);
    TH2F* display_x_y(TApplication*);

protected:
    double* array1;
    double* array2;
    double* array3;
    double* array4;
    size_t num;
    double Ex, Rx, dRxdZ, Ey, Ry, dRydZ;
};

class Gauss_4D: public Distr_base_4D {
public:
    Gauss_4D() {}
    Gauss_4D(int n, double emitX, double envX, double divX, double emitY, double envY, double divY): Distr_base_4D(n, emitX, envX, divX, emitY, envY, divY) { create(); }
    ~Gauss_4D();

private:
    void create();
};

class KV_4D: public Distr_base_4D {
public:
    KV_4D() {}
    KV_4D(int n, double emitX, double envX, double divX, double emitY, double envY, double divY): Distr_base_4D(n, emitX, envX, divX, emitY, envY, divY) { create(); }
    ~KV_4D();

private:
    void create();
};

class Parabolic_4D: public Distr_base_4D {
public:
    Parabolic_4D() {}
    Parabolic_4D(int n, double emitX, double envX, double divX, double emitY, double envY, double divY): Distr_base_4D(n, emitX, envX, divX, emitY, envY, divY) { create(); }
    ~Parabolic_4D();

private:
    void create();
};

class Waterbag_4D: public Distr_base_4D {
public:
    Waterbag_4D() {}
    Waterbag_4D(int n, double emitX, double envX, double divX, double emitY, double envY, double divY): Distr_base_4D(n, emitX, envX, divX, emitY, envY, divY) { create(); }
    ~Waterbag_4D();

private:
    void create();
};

#endif



    
