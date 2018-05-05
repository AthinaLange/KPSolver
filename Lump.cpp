/// Lump Initial Condition.
/** @file
*/

#include "Lump.h"

void    lump(

        int                   nx,
        int                   ny,
        double                a,
        double                h,
        double&               delx,
        double&               dely,
        double                xmin,
        double                ymin,
        std::vector<double>&  eta) {

    int n = nx * ny;
    eta.resize(n, 0.0);

    //=================================================
    // LUMP GENERATION
    //=================================================
    std::cout << delx << " " << dely << std::endl;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double x = xmin + delx * i;
            double y = ymin + dely * j;
            double num = (3/pow(h, 2)) + pow(h, 2) * pow(y, 2) - pow((x + a * y), 2);
            double den = pow((3/pow(h, 2)) + pow(h, 2) * pow(y, 2) + pow((x + a * y), 2), 2);
            eta[i + j * nx] = 4.0 * num/den;
        }
    }
}
