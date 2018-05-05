/// Soliton Initial Condition.
/** @file
*/

#ifndef SEA_STATES_SOLITON_H
#define SEA_STATES_SOLITON_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

/// The initial free surface is generated via \f$\eta = a*cosh^{-2}(k_x*x - k_y * y)\f$.
/**
Parameters for KPII are defined as:
\f$c_0 =\sqrt{g*h_0}\f$ - linear nondispersive wave velocity
\f$\gamma = c_0/2\f$ - coeff in front of second derivative in y
Distance of grid: \f$400 x 200\f$

With simulation parameter defaults (in order to satisfy stability condition)
Time interval: 0.001
Run time : 500,000 for 500s

@param[in]	nx	128 - dimensions of meshgrid.
@param[in]	ny	64 - dimensions of meshgrid.
@param[in]	alpha	\f$\alpha = 3c_0/2h_0\f$ - nonlinear coefficient.
@param[in]	beta	\f$\beta = h_0^2c_0/6\f$ - dispersive coefficient.
@param[in]	a	\f$a = 0.1\f$ - soliton amplitude.
@param[in]	delx	\f$\delta x\f$ - spacing in x.
@param[in]	dely	\f$\delta y\f$ - spacing in y.
@param[in]	xmin	Minimum x grid value.
@param[in]	ymin	Minimum y grid value.
@param[out]	eta	Surface height \f$\eta\f$.
*/

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
                        std::vector<double>&  eta);

#endif //SEA_STATES_SOLITON_H
