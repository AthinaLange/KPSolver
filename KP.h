/// KP Solver.
/** @file
The code is pseudo-spectral and it assumes periodic boundary conditions, both in x and y.
It is a natural extension of the KdV code using the method developed in Fornberg, Bengt,
and G. B. Whitham, "A numerical and theoretical study of certain nonlinear wave phenomena".
Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and
Engineering Sciences 289.1361 (1978): 373-404
*/

#ifndef SEA_STATES_KP_H
#define SEA_STATES_KP_H

#include <vector>
#include <complex>

/// Assigns wavenumber values on a grid
/**
@param[in]	ni	Meshgrid dimensions in x or y - assumed to be even.
@param[in]	dki	\f$\Delta k_x\f$ or \f$\Delta k_y\f$ - spacing between grid points.
@param[out]	deltaki	\f$\delta k_x\f$ or \f$\delta k_y\f$.
*/
void  wavenumbers(

      int                   ni,
      double                dki,
      std::vector<double>&  deltaki);

/// Calculates the dispersion factor
/**
@param[in]	nx	Dimension of meshgrid: x.
@param[in]	ny	Dimension of meshgrid: y.
@param[in]	deltat  Time step.

@param[in]	beta	Dispersive coeff.
@param[in]	gama	Coeff in front of second derivative in y.

@param[in]	deltakx	Spacing in x.
@param[in]	deltaky	Spacing in y.

@param[out]	adel	\f$sin(\omega \Delta t)\f$ - inverse small number approximation, assuming \f$\Delta t \rightarrow 0\f$, corresponding to 3rd derivative of nonlinear term.
*/
// !!! \f$sin(\omega \Deltat)\f$ - corresponds to \f$\mu^3 \Deltat + \mathfrak{O}(\Deltat^3)\f$ in the approximation \f$\Deltat \rightarrow 0\f$ comes from 3rd derivative of nonlinear term.

void  disprelation(

      int nx,
      int ny,
      double deltat,
      double beta,
      double gama,
      const std::vector<double>& deltakx,
      const std::vector<double>& deltaky,
      std::vector<double>& adel);

/// Completes a FT of the input data, to generate an ODE from the original PDE
/**
@param[in]	isign	Indicates wether FT or inverse FT.
@param[in]	ndim	Indicates how many dimensions are present (1 or 2d).
@param[in]	nn	Contains nx and ny values.
@param[out]	data	Array of values that will be FT.
*/
void  fourier(

      int isign,
      int ndim,
      const std::vector<int>& nn,
      std::vector<std::complex<double>>& data);


/// Iterates the Fornberg, Bengt, and G. B. Whitham method through time, outputing a csv file with the free surface height \f$\eta\f$ at every nvis
/**
@param[in]	nx	Number of points in the x-direction (must be even).
@param[in]	ny	Number of points in the y-direction (must be even).
@param[in]	h	Water depth.
@param[in]	alpha	Nonlinear coeff.
@param[in]	beta	Dispersive coeff.
@param[in]	gama	Coeff in front of second derivative in y.
@param[in]	delx	Spacing in x.
@param[in]	dely	Spacing in y.
@param[in]	deltat  Time step.
@param[in]	nstep  	Number of steps.
@param[in]	nvis  	Visualize output after nvis steps.
@param[out]	eta	Surface height \f$\eta\f$.
*/
std::vector<double>  KP(

                     int nx, 
                     int ny, 
                     double h, 
                     double alpha,
                     double beta, 
                     double gama,
                     double delx,
                     double dely, 
                     double deltat, 
                     int nstep,
                     int nvis, 
                     std::vector<double>& eta);

#endif //SEA_STATES_KP_H
