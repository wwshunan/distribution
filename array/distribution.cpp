#include <TRandom1.h>
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TStatistic.h>
#include "distribution.h"

using namespace std;

Guass_2D::~Guass_2D()
{
    delete[] array1;
    delete[] array2;
}

TH2F* Distr_base_2D::display(TApplication* app)
{
    TH2F* ellipse = new TH2F("x-x'", "scatter plot", 10000, -1, 1, 10000, -10, 10);
    for (int i = 0; i != num; ++i) {
        ellipse->Fill(array1[i], array2[i]);
    }
    ellipse->Draw();
    return ellipse;
}

void Guass_2D::create()
{
    double sigma = R / sqrt(emit);
    double dSigma = dRdZ / sqrt(emit);
    double esqrt = 0.5  * sqrt(emit);

    array1 = new double[num];
    array2 = new double[num];

    TRandom1* myrand = new TRandom1();
    double x, dx;

    for (int i = 0; i != num; ++i) {
        x = myrand->Gaus(0, esqrt);
        dx = myrand->Gaus(0, esqrt);
        array1[i] = sigma * x;
        array2[i] = dx / sigma + x * dSigma;
    }
}

TH2F* Distr_base_4D::display_x_px(TApplication* app)
{
    TH2F* ellipse = new TH2F("x-x'", "scatter plot", 10000, -0.5, 0.5, 10000, -2, 2);
    for (int i = 0; i != num; ++i) {
        ellipse->Fill(array1[i], array3[i]);
    }
    ellipse->Draw();
    return ellipse;
}

TH2F* Distr_base_4D::display_y_py(TApplication* app)
{
    TH2F* ellipse = new TH2F("y-y'", "scatter plot", 10000, -0.5, 0.5, 10000, -2, 2);
    for (int i = 0; i != num; ++i) {
        ellipse->Fill(array2[i], array4[i]);
    }
    ellipse->Draw();
    return ellipse;
}

TH2F* Distr_base_4D::display_x_y(TApplication* app)
{
    TH2F* ellipse = new TH2F("x-y", "scatter plot", 10000, -1, 1, 10000, -10, 10);
    for (int i = 0; i != num; ++i) {
        ellipse->Fill(array1[i], array2[i]);
    }
    ellipse->Draw();
    return ellipse;
}

Gauss_4D::~Gauss_4D()
{
    delete[] array1;
    delete[] array2;
    delete[] array3;
    delete[] array4;
}

void Gauss_4D::create()
{
    double sigma_x = Rx / sqrt(Ex);
    double sigma_y = Ry / sqrt(Ey);
    double dSigma_x = dRxdZ / sqrt(Ex);
    double dSigma_y = dRydZ / sqrt(Ey);

    double ratio = Ex / Ey;
    double Ex2 = 2 * Ex;

    double rsqrt = 1 / sqrt(ratio);
    double esqrt = 0.5  * sqrt(Ex);

    array1 = new double[num];
    array2 = new double[num];
    array3 = new double[num];
    array4 = new double[num];

    TRandom1* myrand = new TRandom1();
    double x, dx, y, dy;

    for (int i = 0; i != num; ++i) {
        x = myrand->Gaus(0, esqrt);
        dx = myrand->Gaus(0, esqrt);
        y = myrand->Gaus(0, esqrt);
        dy = myrand->Gaus(0, esqrt);

        array1[i] = sigma_x * x;
        array2[i] = dx / sigma_x + x * dSigma_x;
        array3[i] = sigma_y * y * rsqrt;
        array4[i] = (dy / sigma_y + y * dSigma_y) * rsqrt;
    }
}

KV_4D::~KV_4D()
{
    delete[] array1;
    delete[] array2;
    delete[] array3;
    delete[] array4;
}

void KV_4D::create()
{
    double sigma_x = Rx / sqrt(Ex);
    double sigma_y = Ry / sqrt(Ey);
    double dSigma_x = dRxdZ / sqrt(Ex);
    double dSigma_y = dRydZ / sqrt(Ey);

    double ratio = Ex / Ey;
    double Ex2 = 2 * Ex;

    double rsqrt = 1 / sqrt(ratio);
    double esqrt = 0.5  * sqrt(Ex);

    double f = Ex;

    array1 = new double[num];
    array2 = new double[num];
    array3 = new double[num];
    array4 = new double[num];

    TRandom1* myrand = new TRandom1();
    double AyAy, AxAx, Ax, Ay;
    double phase_x, phase_y;
    double cosx, sinx, cosy, siny;

    for (int i = 0; i != num; ++i) {
        if (i % 2 == 0) {
            AyAy = myrand->Uniform(0, 1) * f / ratio;
            AxAx = f - AyAy * ratio;
        } else {
            AxAx = myrand->Uniform(0, 1) * f;
            AyAy = (f - AxAx) / ratio;
        }

        Ax = sqrt(AxAx);
        Ay = sqrt(AyAy);

        phase_x = 2 * TMath::Pi() * myrand->Uniform(0, 1);
        phase_y = 2 * TMath::Pi() * myrand->Uniform(0, 1);

		cosx = TMath::Cos(phase_x);
        cosy = TMath::Cos(phase_y);
        sinx = TMath::Sin(phase_x);
        siny = TMath::Sin(phase_y);

        array1[i] = Ax * sigma_x * cosx;
        array2[i] = Ay * sigma_y * cosy;
        array3[i] = Ax * (dSigma_x * cosx - sinx / sigma_x);
        array4[i] = Ay * (dSigma_y * cosy - siny / sigma_y);
    }
}

Parabolic_4D::~Parabolic_4D()
{
    delete[] array1;
    delete[] array2;
    delete[] array3;
    delete[] array4;
}

void Parabolic_4D::create()
{
    double sigma_x = Rx / sqrt(Ex);
    double sigma_y = Ry / sqrt(Ey);
    double dSigma_x = dRxdZ / sqrt(Ex);
    double dSigma_y = dRydZ / sqrt(Ey);

    double ratio = Ex / Ey;
    double Ex2 = 2 * Ex;

    double rsqrt = 1 / sqrt(ratio);
    double esqrt = 0.5  * sqrt(Ex);


    array1 = new double[num];
    array2 = new double[num];
    array3 = new double[num];
    array4 = new double[num];

    TRandom1* myrand = new TRandom1();
    double AyAy, AxAx, Ax, Ay;
    double phase_x, phase_y;
    double cosx, sinx, cosy, siny;
    double r, f, f1, f2, f3, alpha3;

    for (int i = 0; i != num; ++i) {
        r = myrand->Uniform(0, 1);
        alpha3 = (TMath::ACos(1 - 2 * r)) / 3;
        f1 = Ex * (1 + 2 * TMath::Cos(alpha3));
        f2 = Ex * (1 - 2 * TMath::Cos(alpha3 + TMath::Pi() * 2 / 3));
        f3 = Ex * (1 - 2 * TMath::Cos(alpha3 - TMath::Pi() * 2 / 3));

        if (f1 > 0 && f1 < Ex2) 
            f = f1;
        if (f2 > 0 && f2 < Ex2) 
            f = f2;
        if (f3 > 0 && f3 < Ex2) 
            f = f3;
        if (r == 0) 
            f = 0;
        if (r == 1) 
            f = Ex2;

        if (i % 2 == 0) {
            AyAy = myrand->Uniform(0, 1) * f / ratio;
            AxAx = f - AyAy * ratio;
        } else {
            AxAx = myrand->Uniform(0, 1) * f;
            AyAy = (f - AxAx) / ratio;
        }

        Ax = sqrt(AxAx);
        Ay = sqrt(AyAy);

        phase_x = 2 * TMath::Pi() * myrand->Uniform(0, 1);
        phase_y = 2 * TMath::Pi() * myrand->Uniform(0, 1);

		cosx = TMath::Cos(phase_x);
        cosy = TMath::Cos(phase_y);
        sinx = TMath::Sin(phase_x);
        siny = TMath::Sin(phase_y);

        array1[i] = Ax * sigma_x * cosx;
        array2[i] = Ay * sigma_y * cosy;
        array3[i] = Ax * (dSigma_x * cosx - sinx / sigma_x);
        array4[i] = Ay * (dSigma_y * cosy - siny / sigma_y);
    }
}

Waterbag_4D::~Waterbag_4D()
{
    delete[] array1;
    delete[] array2;
    delete[] array3;
    delete[] array4;
}

void Waterbag_4D::create()
{
    double sigma_x = Rx / sqrt(Ex);
    double sigma_y = Ry / sqrt(Ey);
    double dSigma_x = dRxdZ / sqrt(Ex);
    double dSigma_y = dRydZ / sqrt(Ey);

    double ratio = Ex / Ey;
    double Ex2 = 2 * Ex;

    double rsqrt = 1 / sqrt(ratio);
    double esqrt = 0.5  * sqrt(Ex);


    array1 = new double[num];
    array2 = new double[num];
    array3 = new double[num];
    array4 = new double[num];

    TRandom1* myrand = new TRandom1();
    double AyAy, AxAx, Ax, Ay;
    double phase_x, phase_y;
    double cosx, sinx, cosy, siny;
    double r, f;

    for (int i = 0; i != num; ++i) {
        f = Ex * 1.5 * sqrt(myrand->Uniform(0, 1));

        if (i % 2 == 0) {
            AyAy = myrand->Uniform(0, 1) * f / ratio;
            AxAx = f - AyAy * ratio;
        } else {
            AxAx = myrand->Uniform(0, 1) * f;
            AyAy = (f - AxAx) / ratio;
        }

        Ax = sqrt(AxAx);
        Ay = sqrt(AyAy);

        phase_x = 2 * TMath::Pi() * myrand->Uniform(0, 1);
        phase_y = 2 * TMath::Pi() * myrand->Uniform(0, 1);

		cosx = TMath::Cos(phase_x);
        cosy = TMath::Cos(phase_y);
        sinx = TMath::Sin(phase_x);
        siny = TMath::Sin(phase_y);

        array1[i] = Ax * sigma_x * cosx;
        array2[i] = Ay * sigma_y * cosy;
        array3[i] = Ax * (dSigma_x * cosx - sinx / sigma_x);
        array4[i] = Ay * (dSigma_y * cosy - siny / sigma_y);
    }
}

int main(int argc, char* argv[])
{
    TApplication* app = new TApplication("ellipse", &argc, argv);
    TH2F* ellipse;
    //Guass_2D gaus(100000, 0.2, 0.2, 1);
    Waterbag_4D pb(80000, 0.2, 0.2, 1, 0.25, 0.2, 1);
    //KV_4D pb(20000, 0.2, 0.2, 1, 0.25, 0.2, 1);
    //ellipse = gaus.display(app);
    ellipse = pb.display_x_px(app);
    app->Run();
    delete app;
    delete ellipse;
    return 0;
}


                
