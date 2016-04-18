#include <TRandom1.h>
#include <cmath>
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

int main(int argc, char* argv[])
{
    TApplication* app = new TApplication("ellipse", &argc, argv);
    TH2F* ellipse;
    Guass_2D gaus(100000, 0.2, 0.2, 1);
    ellipse = gaus.display(app);
    app->Run();
    delete app;
    delete ellipse;
    return 0;
}


                
