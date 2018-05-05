/// KP Solver.
/** @file
*/

#include <iostream>
#include <iomanip>
#include <fstream>

#include "global.h"
#include "csv.h"
#include "KP.h"

void wavenumbers(int ni,
                 double dki,
                 std::vector<double>& deltaki) { //calculates wavenumber at meshgrid points

    deltaki.resize(ni);

    for (int i = 0; i < ni; ++i) {
        if (i <= (ni/2)) { //ni assumed to be even
            deltaki[i] = dki * i;
        } else {
            deltaki[i] = dki * (i-ni);
        }
    }
}

void disprelation(int nx,
                  int ny,
                  double deltat,
                  double beta,
                  double gama,
                  const std::vector<double>& deltakx,
                  const std::vector<double>& deltaky,
                  std::vector<double>& adel) {
//const so that no changes occur within vector deltakx, and & to pass as reference and save memory bc not copying completely
    int n = nx * ny;
    adel.resize(n);

    for (int i = 0; i < n; ++i) {
        int j = i % nx;
        int l = i / nx;

        if (deltakx[j] == 0.0) {
            adel[i] = 0.0; //omega = 0.0;
        } else {
            double x = deltakx[j];
            double y = deltaky[l];
            double omega = (double)(beta * std::pow((long double)(x), (long double)(3.0))
                                    - (gama * std::pow((long double)(y), (long double)(2.0))) / x); //if statement, gives x =/0
            adel[i] = sin(omega * deltat); //from nonlinear 3rd derivative
        }
    }
}

void fourier(int isign,
             int ndim,
             const std::vector<int>& nn,
             std::vector<std::complex<double>>& data) {

    int ntot = 1;
    for (int idim = 0; idim < ndim; ++idim) {
        ntot *= nn[idim];
    }

    int nprev = 1;
    for (int idim = 0; idim < ndim; ++idim) {
        int s = nn[idim];
        int nrem = ntot / (s * nprev);

        int ip1 = nprev;
        // changed indexing of data to complex numbers > took off: 2*
        int ip2 = ip1 * s;
        int ip3 = ip2 * nrem;
        int i2rev = 0;

        for (int i2 = 0; i2 < ip2; i2 += ip1) {
            if (i2 < i2rev) {
                for (int i1 = i2; i1 < i2 + ip1; ++i1) {
                    // changed indexing of data to complex numbers > changed: -2 to -1 and +2 to +1
                    for (int i3 = i1; i3 <= ip3; i3 += ip2) {;
                        int i3rev = i2rev + i3 - i2;
                        std::swap(data[i3], data[i3rev]);
                        // changed indexing of data to complex numbers
                    }
                }
            }
            int ibit = ip2 / 2;
            while ((ibit >= ip1) && (i2rev >= ibit)) {
                i2rev -= ibit;
                ibit /= 2;
            }
            i2rev += ibit;
        }

        int ifp1 = ip1;

        while (ifp1 < ip2) {
            int ifp2 = 2 * ifp1;
            double theta = (isign * 2.0 * PI) / (ifp2 / ip1);
            std::complex<double> wp = {-2.0 * pow(sin(0.5 * theta), 2.0), sin(theta)};
            std::complex<double> w = {1.0, 0.0};
            for (int i3 = 0; i3 < ifp1; i3 += ip1) {
                for (int i1 = i3; i1 < i3 + ip1; ++i1) {
                    // changed indexing of data to complex numbers > changed: -2 to -2 and +2 to +1
                    for (int i2 = i1; i2 < ip3; i2 += ifp2) {
                        int k1 = i2;
                        int k2 = k1 + ifp1;
                        std::complex<double> temp = {(double)((long double)(w.real()) * (long double)(data[k2].real())
                                                              - (long double)(w.imag()) * (long double)(data[k2].imag())),
                                                     (double)((long double)(w.real()) * (long double)(data[k2].imag())
                                                              + (long double)(w.imag()) * (long double)(data[k2].real()))};
                        data[k2] = data[k1] - temp;
                        data[k1] = data[k1] + temp;
                    }
                }
                w = {w.real() * wp.real() - w.imag() * wp.imag() + w.real(),
                     w.imag() * wp.real() + w.real() * wp.imag() + w.imag()};
            }
            ifp1 = ifp2;
        }
        nprev = s * nprev;
    }
}


std::vector<double> KP(int nx, // number of points in the x-direction (must be even)
                       int ny, // number of points in the y-direction (must be even)
                       double h, // water depth
                       double alpha, // nonlinear coeff
                       double beta, // dispersive coeff
                       double gama, // coeff in front of second derivative in y
                       double delx, // spacing in x
                       double dely, // spacing in y
                       double deltat, // time step
                       int nstep, // number of steps
                       int nvis, // visualize output after nvis steps
                       std::vector<double>& eta) {

    std::cout << std::setprecision(7) << std::endl;
    std::cout << std::endl;
    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    std::cout << "h: " << h << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "gamma: " << gama << std::endl;
    std::cout << "deltat: " << deltat << std::endl;
    std::cout << "nstep: " << nstep << std::endl;
    std::cout << "nvis: " << nvis << std::endl;
    std::cout << std::endl;

    //=================================================
    // VARIABLE DECLARATION
    //=================================================
    int n; // total number of grid points (nx*ny)
    std::vector<int> nn; // (nx, ny)
    nn.resize(2);
    int ndim; // dimension of grid

    double dkx;
    double dky;

    // Value assigned in wavenumber
    std::vector<double> deltakx;
    std::vector<double> deltaky;

    // Value assigned in disrelation
    std::vector<double> adel;

    int isign; // value of 1 or -1, depending on direction of Fourier Transform
    double time;

    n = nx * ny;
    nn[0] = nx;
    nn[1] = ny;
    ndim = 2;

    //multiple storage arrays for complex eta for leap-frog scheme
    std::vector<std::complex<double>> eta0(n, {0.0, 0.0});
    std::vector<std::complex<double>> eta1(n, {0.0, 0.0});
    std::vector<std::complex<double>> eta2(n, {0.0, 0.0});

    std::vector<std::complex<double>> nonlin(n, {0.0, 0.0});
    std::vector<std::complex<double>> dis(n, {0.0, 0.0});
    std::vector<std::complex<double>> dam(n, {0.0, 0.0});

    //=================================================
    // ASSIGNING VARIABLE VALUES
    //=================================================

    dkx = 2.0 * PI / (delx * nx); // spacing between wave numbers
    dky = 2.0 * PI / (dely * ny); //

    wavenumbers(nx, dkx, deltakx);
    wavenumbers(ny, dky, deltaky);

    disprelation(nx, ny, deltat, beta, gama, deltakx, deltaky, adel);

    //=================================================
    // INTEGRATE KP
    //=================================================

    //initialize values of eta0 and eta1 for later use
    for (int i = 0; i < n; ++i) {
        eta0[i] = {eta[i], 0.0};
        eta1[i] = eta0[i];
    }

    bool calc_blow_up = false; //setting an exit parameter to terminate program

    //Time iterations
    for (int k = 0; k < nstep; ++k) {
        time = k * deltat; //to output total run time (s)

        for (int i = 0; i < n; ++i) {
            dam[i] = {eta1[i].real() / n,
                      eta1[i].imag() / n}; // put surface height in a complex vector for FT use
        }

        //=================================================
        // EVALUATE DERIVATIVES - FT
        //=================================================

        isign = 1;
        fourier(isign, ndim, nn, dam); //calculates eta into fourier space


        for (int i = 0; i < n; ++i) {
            int j = i % nx;
            nonlin[i] = dam[i] * deltakx[j]; // allows calculation of derivative of nonlinear term
            dis[i] = dam[i] * adel[i]; // advances the linear part in Fourier space
        }

        isign = -1;
        fourier(isign, ndim, nn, nonlin); //calculating the nonlinear part into normal space
        fourier(isign, ndim, nn, dis); //calculating the


        for (int i = 0; i < n; ++i) {
            eta2[i] = eta0[i] //original value
                      + 2.0 * deltat * std::complex<double>(0.0, 1.0) * alpha * (nonlin[i] * eta[i]) //nonlinear term
                      - 2.0 * std::complex<double>(0.0, 1.0) * dis[i]; //dispersive term

            eta0[i] = eta1[i];
            eta1[i] = eta2[i];
            eta[i]  = eta1[i].real();

            if (eta2[i].real() > 100.0) {
                calc_blow_up = true;
                break; // check to avoid calculation blowing up
            }
        }

        std::cout << "iteration = " << k << std::endl;

        if (calc_blow_up) {
            std::cout << "Calculation blow up at iteration: " << k << "!" << std::endl;
            break;
        }

        if ((k % nvis) == 0) {
        //write value of eta to csv file
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(5) << k/nvis + 1; // +1 to conform to fortran indexing
            write_wave_csv("Wave_" + ss.str() + ".csv", nx, ny, delx, dely, eta);

        }


    //std::ofstream output_file(fn);
    /*if (output_file.is_open()) {

        for (int j = 0; j < nx * ny; ++j) {
            output_file << eta[j] << std::endl;
        }
        output_file.close();
    }*/
} // end loop in time

    std::cout << "Program ended at time: " << time << std::endl;
    return eta;
}
