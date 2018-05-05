/// Soliton Initial Condition.
/** @file
*/

#include "Soliton.h"

std::vector<double>     soliton(

                        int                   nx,
                        int                   ny,
                        double                alpha,
                        double                beta,
                        double                a,
                        double&               delx,
                        double&               dely,
                        double                xmin,
                        double                ymin,
                        std::vector<double>&  eta){

    int n = nx * ny;
    eta.resize(n, 0.0);

    double theta = 0.0;
    double k = 0.5 * sqrt((a / 3.0) * (alpha / beta));
    double l = tan(theta) * k; // k_y = 0 for the soliton wave to be parallel to the x-axis

    //==================================================
    // SOLITON GENERATION
    //==================================================

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double x = xmin + delx * i;
            double y = ymin + dely * j;
            eta[i + j * nx] = a * pow((1.0 / cosh(k * x - l * y)), 2.0);
        }
    }
    return eta;
}
