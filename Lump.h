/// Lump Initial Condition.
/** @file
*/

#ifndef SEA_STATES_IC_H
#define SEA_STATES_IC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

/// The initial free surface is generated via formula \f$\eta(x,y) = 4 * [-(x-h_0^2y)^2+h_0^2y^2 + 3/h_0^2]/[(x-h_0^2y)^2+h_0^2y^2 + 3/h_0^2]^2\f$
/**
Parameters for KPI are defined as:
\f$c_0 = 0\f$ - linear nondispersive wave velocity
\f$\gamma = -1.0\f$ - coeff in front of second derivative in y
Distance of grid: \f$60 x 100\f$

With simulation parameter defaults (in order to satisfy stability condition)
Time interval: 0.0005
Run time : 1,000,000 for 500s

@param[in]	nx	256 - meshgrid dimensions.
@param[in]	ny	256 - meshgrid dimensions.
@param[in]	a	\f$a = 1.0\f$ - lump amplitude.
@param[in]	h	\f$h = 5.0\f$ - water depth 
@param[in]	delx	\f$\delta x\f$ - spacing in x. 
@param[in]	dely	\f$\delta y\f$ - spacing in y.
@param[in]	xmin	Minimum x grid value.
@param[in]	ymin	Minimum y grid value.
@param[out]	eta	Surface height \f$\eta\f$.
*/
void    lump(

        int                   nx,
        int                   ny,
        double                a,
        double                h,
        double&               delx,
        double&               dely,
        double                xmin,
        double                ymin,
        std::vector<double>&  eta);

#endif //SEA_STATES_IC_H
