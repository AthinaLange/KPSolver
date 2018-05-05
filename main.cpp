/// Program KP
/** @file

\mainpage 
This code solves the KP equation in dimensional units (SI) for the water wave
problem in a frame of reference moving with velocity c_0 = sqrt(g*h).

The code is pseudo-spectral and it assumes periodic boundary conditions, both in x and y.
It is a natural extension of the KdV code using the method developed in Fornberg, Bengt,
and G. B. Whitham, "A numerical and theoretical study of certain nonlinear wave phenomena".
Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and
Engineering Sciences 289.1361 (1978): 373-404.

The integral term is done in Fourier space as 1/k F(k).

The stability condition required is \f$|\frac{\Delta t}{\Delta x}| |\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}| < \frac{1}{\pi ^3}\f$

C++ code, based off a Fortran code developed by Dr. Miguel Onorato.

Initial condition and simulation parameters set: Choice of Soliton, Lump or input for KPII via eta.csv file.
Output: free surface height \f$\eta\f$ via eta.csv for every nvis time step

Generating initial conditions done in MATLAB directly when running with MEX-functions

\author Athina Lange, based off code from Dr. Miguel Onorato
\date April 8th, 2018 
*/

#include "global.h"
#include "Soliton.h"
#include "Lump.h"
#include "KP.h"
#include "csv.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[]){

    try {

        //=================================================
        // VARIABLE DECLERATION
        //=================================================

	std::string initial_condition; // initial condition: "Soliton", "Lump" or "SeaStates" - default: "Soliton"

        int nx; // dimension of meshgrid: x - default: 128
        int ny; // dimension of meshgrid: y - default: 64
        double c0; // linear nondispersive wave velocity c0 = sqrt(g*h)
        double a; // soliton amplitude
        double h; // water depth  - default: 5.0
        double alpha; // nonlinear coeff
        double beta; // dispersive coeff
        double gama; // coeff in front of second derivative in y

        double delx; // spacing in x
        double dely; // spacing in y

        double deltat; // time step - default: 0.001
        int nstep; // number of steps - default: 500000 (500s)
        int nvis; // visualize output after nvis steps - default: 1000

        // Read eta from csv file (no header, single value per line)
        std::vector<double> eta;

        //=================================================
        // VARIABLE INITIALIZING
        //=================================================

        po::options_description desc("Options");
        desc.add_options()
                ("help,h", "help message")
                ("initial-condition", po::value<std::string>(&initial_condition)->default_value("soliton"), "initial condition: soliton, lump or seastates")
                ("nx", po::value<int>(&nx)->default_value(128), "number of points in x direction (even number)")
                ("ny", po::value<int>(&ny)->default_value(64), "number of points in y direction (even number)")
                ("h", po::value<double>(&h)->default_value(5.0), "water depth")
                ("deltat", po::value<double>(&deltat)->default_value(0.001), "time step")
                ("nstep", po::value<int>(&nstep)->default_value(500000), "number of steps")
                ("nvis", po::value<int>(&nvis)->default_value(1000), "visualize output after nvis steps");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }

        if (initial_condition == "soliton") {

            a = 0.1;
            c0 = sqrt(G * h);
            alpha = (3.0 / 2.0) * c0 / h; //assumed to be >=0
            beta = h * h * c0 / 6.0; //assumed to be > 0
            gama = 0.5 * c0;

            double xmin = 0.0;
            double xmax = 400.0;
            double ymin = 0.0;
            double ymax = 200.0;

            delx = (xmax - xmin) / nx;
            dely = (ymax - ymin) / ny;

            soliton(nx, ny, alpha, beta, a, delx, dely, xmin, ymin, eta); // Run soliton initial condition

        } else if (initial_condition == "lump") {

            a = 1.0;
            c0 = 0.0;
            alpha = 6.0;
            beta = 1.0;
            gama = -1.0;

            //Change parameters for meshgrid (256 x 256), time interval (0.0005) and run time (1000000) to satisfy stability condition
            nx = 256;
            ny = 256;
            deltat = 0.0005;
            nstep = 1000000;

            double xmin = -30.0;
            double xmax = 30.0;
            double ymin = -50.0;
            double ymax = 50.0;

            delx = (xmax - xmin) / nx;
            dely = (ymax - ymin) / ny;

            lump(nx, ny, a, h, delx, dely, xmin, ymin, eta); // Run lump initial condition

        } else if (initial_condition == "seastates") {

            a = 0.1;
            c0 = sqrt(G * h);
            alpha = (3.0 / 2.0) * c0 / h; //assumed to be >=0
            beta = h * h * c0 / 6.0; //assumed to be > 0
            gama = 0.5 * c0;

            double xmin = 0.0;
            double xmax = 400.0;
            double ymin = 0.0;
            double ymax = 200.0;

            delx = (xmax - xmin) / nx;
            dely = (ymax - ymin) / ny;

            read_wave_csv("eta.csv", eta); // Read value of free surface eta from file

        } else {
            std::cout << "initial-condition needs to be soliton, lump or seastates, not " << initial_condition << std::endl;
            return 1;
        }

        write_soliton_csv("Soliton.csv", nx, ny, eta); // Write initial condition data to file

        //=================================================
        // RUN PROPAGATOR
        //=================================================

        KP(nx, ny, h, alpha, beta, gama, delx, dely, deltat, nstep, nvis, eta); // Run KP propagator

        write_wave_csv("Wave_1000.csv", nx, ny, delx, dely, eta); // Write final free surface eta to file
    }

    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!" << std::endl;
        return 1;
    }

    return 0;
}
